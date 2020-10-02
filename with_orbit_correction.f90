program tracking

  use orbit_correction_mod
  use get_pseudo_inverse_mod
  use runge_kutta_mod
  use bmad
  use random_mod
  use beam_mod
  implicit none
  type (lat_struct), target :: lat
  type (coord_struct), allocatable :: orb0(:), orb0_corrected(:),track_turn(:), CO_orbit(:)
  type (coord_struct) :: orbit
  integer, parameter :: i_dim = 6
  integer ( kind = 4 ) mx, my, nx, ny
  integer i_part, i_turn, i_turn_acc
  integer ( kind = 4 ) loc_steer, loc_steer_copy, count
  real ( kind = 8 ), allocatable, dimension ( : , : ) :: orm
  real ( kind = 8 ), allocatable, dimension ( : ) :: deltaOrbit, deltaOrbit_corrected
  real ( kind = 8 ), allocatable, dimension ( : ) :: kicks
  real(rp) :: s_body
  type (ele_struct), pointer :: ele, slave
  real(rp) :: loc_vec(3)
! Beam
  type (beam_init_struct) beam_init
  type (beam_struct) beam
  type (bunch_params_struct) bunch_params
  logical err_beam, err_co
! for time information
  real :: start, finish

! ------------------- changes are allowed only here ------------------------  
! Parameters
! general menu for misalignments, kickers, orbit correction and EDM
  logical, parameter :: orbit_correction = .false. ! false -> NO correction (N,C)
  logical, parameter :: kickers = .false.  ! false -> standard kicker (S,D) 
  logical, parameter :: edm = .true. ! false -> EDM=0 (N,E)
  logical, parameter :: qmisalignments = .false. ! false -> NO quad misalignments (N,Q)
  logical, parameter :: dmisalignments = .false. ! false -> NO dipole misalignments (N,D)
  real,    parameter :: edmvalue=1. ! EDM value
  integer, parameter :: N_turn = 1000 ! number of turns
  integer, parameter :: N_seeds = 1 ! number of seeds
! ------------------- end of allowed changes -------------------------------  
  
  integer, parameter :: & ! for quadrupoles and dipoles misalignments
       qox = 1     , & ! Q offset x
       qoy = 1     , & ! Q offset y
       qoz = 0     , & ! Q offset z
       qtx = 1     , & ! Q tilt x
       qty = 1     , & ! Q tilt y
       qtz = 0     , & ! Q tilt z
       dox = 0     , & ! D offset x
       doy = 0     , & ! D offset y
       doz = 0     , & ! D offset z
       dtx = 0     , & ! D tilt x
       dty = 1     , & ! D tilt y
       dtz = 1     , & ! D tilt z
       dro = 0         ! D roll z
  
  integer ( kind = 4 ) icont
!
  integer, parameter :: &
       N_parts = 1 , &
       N_bunch = 1 , &
       outfile = 20  
  real(rp), parameter :: &
       EPSx = 0._rp , &
       EPSy = 0._rp , &
       sigmaZ = 0.*10._rp/5  , &
       sigmadp = 0.*2.05e-3/5 

! misalignments
  real ( kind = 8 ) muq, mud
  real ( kind = 8 ) sigmaq, sigmad
  real ( kind = 8 ) r8_normal_ab
  real ( kind = 8 ) r_x
  real ( kind = 8 ) r_y
  real ( kind = 8 ) r_z
  real ( kind = 8 ) r_tilt_x
  real ( kind = 8 ) r_tilt_y
  real ( kind = 8 ) r_tilt_z
  real(rp) ::  rev_freq_ref
  integer ( kind = 4 ) n_ele
  integer ( kind = 8 ) seed_xq
  integer ( kind = 8 ) seed_yq
  integer ( kind = 8 ) seed_zq
  integer ( kind = 8 ) seed_tilt_xq
  integer ( kind = 8 ) seed_tilt_yq
  integer ( kind = 8 ) seed_tilt_zq
  integer ( kind = 8 ) seed_xd
  integer ( kind = 8 ) seed_yd
  integer ( kind = 8 ) seed_zd
  integer ( kind = 8 ) seed_tilt_xd
  integer ( kind = 8 ) seed_tilt_yd
  integer ( kind = 8 ) seed_tilt_zd
  integer ( kind = 8 ) in_seed
  integer :: i, j, fid
  character(100) :: corr="N", kickern="S", edmy="N", quad="N", dip="N", test, fname
  real(dp) :: vec0(6)
  icont=0
  if(edm) icont=10
  if(dmisalignments .or. qmisalignments) icont=icont+1

! Set mu and sigma for quadrupole misalignments
  muq = 0.0
  sigmaq = 0.0002
! Set mu and sigma for dipole misalignments
  mud = 0.0
  sigmad = 0.0002
! Set seed
  seed_xq=67357
  seed_yq=856031219
  seed_zq=317456131
  seed_tilt_xq=9815709
  seed_tilt_yq=945161
  seed_tilt_zq=1398145
  seed_xd=5673451
  seed_yd=34853941
  seed_zd=9458671
  seed_tilt_xd=2314834
  seed_tilt_yd=2981375
  seed_tilt_zd=856418
  loc_steer = 0
  
! Parsing and initialization
  bmad_com%spin_tracking_on = .true.
  bmad_com%auto_bookkeeper = .false.
  call bmad_parser ("/home/anjali/bmad/Anjali/COSY_steerer_simulations/lattice/prof_lattice/lattice_after_quad_dipole_overlap/COSY_default_ideal.bmad", lat)   ! Read in a lattice.
  call set_on_off (rfcavity$, lat, on$) !!!!!!!!!!!!! changed !!!!!!!!!!!!!!!!!!
  ! Set EDM value
  n_ele = lat%n_ele_max
  bmad_com%electric_dipole_moment=edmvalue*bmad_com%electric_dipole_moment
  if( .not. edm ) bmad_com%electric_dipole_moment=0.
  call twiss_at_start (lat)
  call twiss_propagate_all (lat)
  call closed_orbit_calc (lat, orb0, i_dim) !!!!!!!!!!!!! changed !!!!!!!!!!!!!!!!!!
  call cpu_time(start)
  do i = 1, n_ele
    ele => lat%ele(i)
    write(*,*) key_name(ele%key) 
  enddo
  do i = 1,N_seeds ! Here you define how many random sets you want to take
! "test" defined below will be part of each file name, so you can easily get the seeds. Be aware: "test" is not directly the seed. See below how they are calculated. If you really need the misalignments, better write them to a file each time.
     write(*,*) i
10 continue

! reset all kickers     
     do j = 0, lat%n_ele_max
        ele => lat%ele(j)
        if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"H",.false.) == 1 .and. index(ele%name, "TOROID", .false.) /= 1 .and. index(ele%name, "MBU", .false.) /= 1 ) then
           ele%value(kick$) = 0.
           call set_flags_for_changed_attribute (ele, ele%value(kick$))
        endif
        if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"V",.false.) == 1 ) then
           ele%value(kick$) = 0.
           call set_flags_for_changed_attribute (ele, ele%value(kick$))
        endif
        call lattice_bookkeeper (lat)
     enddo

! quadrupoles misalignments
     if( qmisalignments ) then
        do j = 0, lat%n_ele_max
           r_x = r8_normal_ab ( muq, sigmaq, seed_xq )
           r_y = r8_normal_ab ( muq, sigmaq, seed_yq )
           r_z = r8_normal_ab ( muq, sigmaq, seed_zq )
           r_tilt_x = r8_normal_ab ( muq, sigmaq, seed_tilt_xq )
           r_tilt_y = r8_normal_ab ( muq, sigmaq, seed_tilt_yq )
           r_tilt_z = r8_normal_ab ( muq, sigmaq, seed_tilt_zq )
           ele => lat%ele(j)
           if (index(ele%name,"QU",.false.) == 1 .or. index(ele%name,"QT",.false.) == 1 ) then
!           if (index(ele%name,"QT3 ",.false.) == 1 ) then
              if(qox==1) ele%value(x_offset$) = r_x
              if(qoy==1) ele%value(y_offset$) = r_y
              if(qoz==1) ele%value(z_offset$) = r_z
              if(qtx==1) ele%value(x_pitch$) = r_tilt_x
              if(qty==1) ele%value(y_pitch$) = r_tilt_y
              if(qtz==1) ele%value(tilt$) = r_tilt_z
              if(qox==1) call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
              if(qoy==1) call set_flags_for_changed_attribute (ele, ele%value(y_offset$))
              if(qoz==1) call set_flags_for_changed_attribute (ele, ele%value(z_offset$))
              if(qtx==1) call set_flags_for_changed_attribute (ele, ele%value(x_pitch$))
              if(qty==1) call set_flags_for_changed_attribute (ele, ele%value(y_pitch$))
              if(qtz==1) call set_flags_for_changed_attribute (ele, ele%value(tilt$))
              call lattice_bookkeeper (lat) !!!!!!!!!!!!! changed !!!!!!!!!!!!!!!!!!
           endif
        enddo
     endif
     ! dipoles misalignments
     if( dmisalignments ) then
        do j = 0, lat%n_ele_max
           r_x = r8_normal_ab ( muq, sigmaq, seed_xd )
           r_y = r8_normal_ab ( muq, sigmaq, seed_yd )
           r_z = r8_normal_ab ( muq, sigmaq, seed_zd )
           r_tilt_x = r8_normal_ab ( muq, sigmaq, seed_tilt_xd )
           r_tilt_y = r8_normal_ab ( muq, sigmaq, seed_tilt_yd )
           r_tilt_z = r8_normal_ab ( muq, sigmaq, seed_tilt_zd )
           ele => lat%ele(j)
           if (index(ele%name,"BE",.false.) == 1 .and. (index(ele%name,"BEG",.false.) == 0)) then
              if(dox==1) ele%value(x_offset$) = r_x
              if(doy==1) ele%value(y_offset$) = r_y
              if(doz==1) ele%value(z_offset$) = r_z
              if(dtx==1) ele%value(x_pitch$) = r_tilt_x
              if(dty==1) ele%value(y_pitch$) = r_tilt_y
              if(dtz==1) ele%value(ref_tilt$) = r_tilt_z
              if(dro==1) ele%value(roll$) = r_tilt_z
              if(dox==1) call set_flags_for_changed_attribute (ele, ele%value(x_offset$))
              if(doy==1) call set_flags_for_changed_attribute (ele, ele%value(y_offset$))
              if(doz==1) call set_flags_for_changed_attribute (ele, ele%value(z_offset$))
              if(dtx==1) call set_flags_for_changed_attribute (ele, ele%value(x_pitch$))
              if(dty==1) call set_flags_for_changed_attribute (ele, ele%value(y_pitch$))
              if(dtz==1) call set_flags_for_changed_attribute (ele, ele%value(ref_tilt$))
              if(dro==1) call set_flags_for_changed_attribute (ele, ele%value(roll$))
              call lattice_bookkeeper (lat)
           endif
        enddo
     endif
     call closed_orbit_calc (lat, CO_orbit, i_dim, 1, 0, err_co, .false.)
     if( err_co ) then
        write(*,*) 'bad orbit 1'
        go to 10
     endif
     call lat_make_mat6 (lat)
     call track_all (lat, track_turn)
     
     ! orbit correction
     if( orbit_correction ) then
        !vec0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] !Initial phase spac vector for the particle
        count = 1
        vec0 = CO_orbit(0)%vec !Initial phase spac vector for the particle
        track_turn(0)%vec=CO_orbit(0)%vec !!!!!!!!!!!!! changed !!!!!!!!!!!!!!!!!!
        call calc_ORM(lat, CO_orbit, loc_steer, mx, my, nx, ny ,orm)
        call get_inverse_orbit(lat, CO_orbit, mx, my, deltaOrbit)
        call get_steerer_kicks(orm, deltaOrbit, mx, my, nx, ny, kicks)
!        print *, "Round", count, " finished."
        if (maxval(abs(kicks)) > 0.01) then
           count = count + 1
           if(count>5) then
              write(*,*) 'count ',count
              go to 10 ! allowed to exclude 5 steerers
           endif
!           print *, "Maximum steerer kick exceeds 0.01 radian."
!           print *, "Steerer at position", maxloc(abs(kicks)), "will be excluded."
!           print *, "Recalculation of ORM starts."
           
           loc_steer = maxloc(abs(kicks), 1)
           loc_steer_copy = loc_steer
           
           call calc_ORM(lat, CO_orbit, loc_steer, mx, my, nx, ny ,orm)
           call get_inverse_orbit(lat, CO_orbit, mx, my, deltaOrbit)
           call get_steerer_kicks(orm, deltaOrbit, mx, my, nx, ny, kicks)
!           print *, "Round", count, " finished."
        end if
!        print *, "Steerer values okay."
        
        call apply_kicks(lat, loc_steer_copy, kicks, nx, ny)
        
        call closed_orbit_calc (lat, CO_orbit, i_dim, 1, 0, err_co, .false.)
        if( err_co ) then
           write(*,*) 'bad orbit 2'
           go to 10
        endif
     
        call get_inverse_orbit(lat, CO_orbit, mx, my, deltaOrbit_corrected) 
        call lat_make_mat6 (lat)
     endif
     
     !horizontal tune
     !print *, 'horizontal tune = ',lat%ele(lat%n_ele_track)%a%phi/twopi
     !vertical tune
     !print *, 'vertical tune = ',lat%ele(lat%n_ele_track)%b%phi/twopi
     !rev_freq_ref = 1./(CO_orbit(lat%n_ele_track)%t-CO_orbit(0)%t)
     !print *, 'revolution frequency = ',rev_freq_ref
     !print *,'horizontal = ',rev_freq_ref*lat%ele(lat%n_ele_track)%a%phi/twopi
     !print *,'vertical = ',rev_freq_ref*lat%ele(lat%n_ele_track)%b%phi/twopi
     
     beam_init%n_particle = N_parts
     beam_init%n_bunch = N_bunch
     beam_init%a_emit = EPSx 
     beam_init%b_emit = EPSy
     
     beam_init%random_sigma_cutoff = 4.0
     beam_init%init_spin = .true.
     beam_init%spin = [0.,0.,1.]
     
     if(i<10) write(test,'(i1)') i 
     if(9<i.and.i<100) write(test,'(i2)') i 
     if(99<i.and.i<1000) write(test,'(i3)') i 
     if(999<i.and.i<10000) write(test,'(i4)') i 
     if(9999<i.and.i<100000) write(test,'(i5)') i 
     fid = outfile + i
     if(qmisalignments) quad="Q"
     if(dmisalignments) dip="D"
     if(edm) edmy="E"
     if(kickers) kickern="D"
     if(orbit_correction) corr="C"
     ! Here the value of "test" is added to the filename in order to identify the seeds
     write(fname,"(5A1,13I1,'_',A)") corr,kickern,edmy,quad,dip,qox,qoy,qoz,qtx,qty,qtz,dox,doy,doz,dtx,dty,dtz,dro,trim(test)
     open(fid, file=fname,status="replace")
     write(fid,"(i2,1x,i5,1x,i3,1x,13i1,1x,f6.0,1x,4(f10.4,1x),es8.1e2)") icont,N_turn,N_parts,qox,qoy,qoz,qtx,qty,qtz,dox,doy,doz,dtx,dty,dtz,dro,edmvalue,muq,sigmaq,mud,sigmad,-bmad_com%electric_dipole_moment/((m_deuteron*2)/(h_bar_planck *c_light)*1.E-2)
     !write(fid,"(i2,1x,i5,1x,i3,1x,13i1,1x,f6.0,1x,es7.1e2)") icont,N_turn,N_parts,qox,qoy,qoz,qtx,qty,qtz,dox,doy,doz,dtx,dty,dtz,dro,edmvalue,-bmad_com%electric_dipole_moment/((m_deuteron*2)/(h_bar_planck *c_light)*1.E-2)
     call init_beam_distribution (lat%ele(0), lat%param, beam_init, beam)
     call calc_bunch_params(beam%bunch(1), bunch_params, err_beam, .true. )
     ! Start tracking  
     do i_part=1, bunch_params%n_particle_tot
        call reallocate_coord (track_turn, lat%n_ele_max)
        !        track_turn(0) = beam%bunch(1)%particle(i_part)
        track_turn(0)%vec=CO_orbit(0)%vec
        track_turn(0)%spin=[0.,0.,1.]
        i_turn_acc=0
        do i_turn=1, N_turn
           call track_all (lat, track_turn)
           track_turn(0) = track_turn(lat%n_ele_track)
           loc_vec = track_turn(lat%n_ele_track)%spin
           do j=1,lat%n_ele_track
              ele => lat%ele(j)
              !if (index(ele%name,"EDDATGT",.false.) == 1 .or. index(ele%name,"WASATGT",.false.) == 1) 
              write(fid,*) track_turn(j)%t,track_turn(j)%spin(2) ,lat%ele(j)%name  
           enddo
           i_turn_acc=i_turn          
        enddo
        if(i_turn_acc/=N_turn) then
           write(*,*) 'liczba obrotow = ',i_turn_acc
           go to 10
        endif
     enddo
     close(fid)
  enddo
  call cpu_time(finish)
!  write(*,'(a7,f10.6,a8)') "Time = ",finish-start," seconds"
end program tracking

