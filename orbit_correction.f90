module orbit_correction_mod

contains
    subroutine calc_ORM(lat, CO_orbit, loc_steer, mx, my, nx, ny, orm)

    !*****************************************************************************
    !
    !! CALC_ORM gets the orbit response matrix.
    !
    !  Parameters:
    !
    !    Input, lat_struct lat, contains all TWISS information about the underlying lattice.
    !
    !    Input, coord_struct CO_orbit, contains closed orbit values. 
    !
    !    Input, integer (kind = 4), loc_steer, location of steerer, that will be excluded in next round ORM calculation. Default is 0 (all steerer used, no second round needed).
    !
    !    Output integer (kind = 4) mx, my, nx, ny, the number of BPMS (m) and correctors (n) in the given direction (x or y).
    ! 
    !    Output, real ( kind = 8 ) orm(M,N), the orbit response matrix for the given direction (in m/rad). 
    !
        use bmad 

        implicit none

        type (lat_struct) lat
        type (coord_struct) :: CO_orbit(:)

        integer ( kind = 4 ) mx
        integer ( kind = 4 ) my
        integer ( kind = 4 ) nx
        integer ( kind = 4 ) ny
        integer ( kind = 4 ) n_ele
        integer ( kind = 4 ) loc_steer, i,j, count1, count2, count3, count4

        real ( kind = 8 ) length
        real ( kind = 8 ) gamma
        real ( kind = 8 ) tune_x
        real ( kind = 8 ) tune_y
        real ( kind = 8 ) mom_comp_fac
        real ( kind = 8 ) gammaTR

        type (normal_modes_struct) :: modes
        type (rad_int_all_ele_struct) ::  ele_rad_int

        real ( kind = 8 ), allocatable, dimension ( : , : ) :: orm
        real ( kind = 8 ), allocatable, dimension ( : ) :: beta_BPM
        real ( kind = 8 ), allocatable, dimension ( : ) :: beta_Corr
        real ( kind = 8 ), allocatable, dimension ( : ) :: phi_BPM
        real ( kind = 8 ), allocatable, dimension ( : ) :: phi_Corr
        real ( kind = 8 ), allocatable, dimension ( : ) :: disp_BPM
        real ( kind = 8 ), allocatable, dimension ( : ) :: disp_Corr


        type (ele_struct), pointer :: ele,slave

        !****************************************************************************
        count1 = 0
        count2 = 0
        count3 = 0
        count4 = 0
        mx = 0
        my = 0
        nx = 0
        ny = 0
        n_ele = lat%n_ele_max

        !****************************************************************************
        
        length = lat%param%total_length
        gamma = sqrt((lat%ele(3)%value(p0c$)/(mass_of(lat%param%particle)))**2+1)
       
        tune_x = lat%ele(lat%n_ele_track)%a%phi/twopi
        tune_y = lat%ele(lat%n_ele_track)%b%phi/twopi
        call radiation_integrals (lat, CO_orbit, modes, rad_int_by_ele = ele_rad_int)
        gammaTR = 1/sqrt(modes%synch_int(1)/ lat%param%total_length)

        !****************************************************************************

        do i = 1, n_ele
            ele => lat%ele(i)
            if (key_name(ele%key) .EQ. 'Monitor' .and. index(ele%name,"H",.false.) > 1) then

                mx = mx+1  

            else if (key_name(ele%key) .EQ. 'Monitor' .and. index(ele%name,"V",.false.) > 1) then

                my = my+1 

            else if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"H",.false.) == 1 .and. index(ele%name, "TOROID", .false.) /= 1 .and. index(ele%name, "MBU", .false.) /= 1 ) then

                nx = nx+1

            else if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"V",.false.) == 1 ) then

                ny = ny+1

            end if
        enddo

        if (loc_steer <= nx .and. loc_steer /= 0) then
            nx = nx -1
        else if (loc_steer > nx .and. loc_steer /= 0) then
            ny = ny-1
        end if

        if (allocated(orm)) then
            deallocate(orm)
        end if
        allocate ( orm(1:(mx+my),1:(nx+ny)) )
        allocate ( beta_BPM(1:(mx+my)) )
        allocate ( beta_Corr(1:(nx+ny)) )
        allocate ( phi_BPM(1:(mx+my)) )
        allocate ( phi_Corr(1:(nx+ny)) )
        allocate ( disp_BPM(1:(mx+my)) )
        allocate ( disp_Corr(1:(nx+ny)) )

        do i = 1, n_ele
            ele => lat%ele(i)
            if (key_name(ele%key) .EQ. 'Monitor' .and. index(ele%name,"H",.false.) > 1) then

                count1 = count1 + 1
                beta_BPM(count1) = ele%a%beta
                phi_BPM(count1) = ele%a%phi
                disp_BPM(count1) = ele%a%eta

            else if (key_name(ele%key) .EQ. 'Monitor' .and. index(ele%name,"V",.false.) > 1) then

                count2 = count2 + 1
                beta_BPM(count2+mx) = ele%b%beta
                phi_BPM(count2+mx) = ele%b%phi
                disp_BPM(count2+mx) = 0.        

            else if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"H",.false.) == 1 .and. index(ele%name, "TOROID", .false.) /= 1 .and. index(ele%name, "MBU", .false.) /= 1) then

                count3 = count3 + 1
                if (count3 .EQ. loc_steer .and. loc_steer /= 0) then
                    count3 = count3 - 1
                    loc_steer = 0
                else
                    beta_Corr(count3) = ele%a%beta
                    phi_Corr(count3) = ele%a%phi
                    disp_Corr(count3) = ele%a%eta

                end if

            else if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"V",.false.) == 1 ) then

                count4 = count4 + 1
                if (count4+nx .EQ. loc_steer .and. loc_steer /= 0) then
                    count4 = count4 - 1
                    loc_steer = 0
                else
                    beta_Corr(count4+nx) = ele%b%beta
                    phi_Corr(count4+nx) = ele%b%phi
                    disp_Corr(count4+nx) = 0.   
                end if     

            end if

        enddo

        !****************************************************************************

        do i = 1, (mx+my)
            do j = 1, (nx+ny)
                if (i <= mx .and. j <= nx) then
                    orm(i,j)=sqrt(beta_BPM(i)*beta_Corr(j))/(2*sin(pi*tune_x))*cos(abs(phi_BPM(i)-phi_Corr(j))-pi*tune_x)-disp_BPM(i)*disp_Corr(j)/(length*((1/gammaTR)**2-(1/gamma)**2))
                else if (i > mx .and. j > nx) then
                    orm(i,j)=sqrt(beta_BPM(i)*beta_Corr(j))/(2*sin(pi*tune_y))*cos(abs(phi_BPM(i)-phi_Corr(j))-pi*tune_y)-disp_BPM(i)*disp_Corr(j)/(length*((1/gammaTR)**2-(1/gamma)**2))
                else
                    orm(i,j) = 0.
                end if
            end do
        end do

        return

    end subroutine calc_ORM

    subroutine get_inverse_orbit(lat, CO_orbit, mx, my, deltaX)

    !*****************************************************************************
    !
    !! GET_INVERSE_ORBIT gives the inverse closed orbit in a vector at the positions of the BPMs.
    !
    !  Parameters:
    !
    !    Input, lat_struct lat, contains all TWISS information about the underlying lattice.
    !
    !    Input, coord_struct CO_orbit, contains closed orbit values. 
    !
    !    Input integer (kind = 4) mx, my, number of BPMs in horizontal (x) and vertical (y) direction. 
    !
    !    Output, real ( kind = 8 ) deltaX(mx+my), the inverse orbit vector (in m). 
    !
    
        use bmad 

        implicit none

        type (lat_struct) lat

        type (coord_struct) :: CO_orbit(:)

        integer ( kind = 4 ) mx, my, count1, count2, i

        real ( kind = 8 ), allocatable, dimension ( : ) :: deltaX

        count1 = 0
        count2 = 0

        if (allocated(deltaX)) then
            deallocate(deltaX)
        end if
        allocate ( deltaX(1:(mx+my)) )

        do i = 1, lat%n_ele_max
            if (key_name(lat%ele(i)%key) .EQ. 'Monitor' .and. index(lat%ele(i)%name,"H",.false.) > 1) then

                count1 = count1 + 1
                deltaX(count1) = -1.*CO_orbit(i)%vec(1)

            else if (key_name(lat%ele(i)%key) .EQ. 'Monitor' .and. index(lat%ele(i)%name,"V",.false.) > 1) then
                   
                count2 = count2 + 1
                deltaX(count2 + mx) = -1.*CO_orbit(i)%vec(3)

            end if
        end do

        return

    end subroutine get_inverse_orbit

    subroutine get_steerer_kicks(orm, deltaX, mx, my, nx, ny, kicks)

    !*****************************************************************************
    !
    !! GET_STEERER_KICKS gives the inverse closed orbit in a vector at the positions of the BPMs.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8) orm, orbit response matrix.
    !
    !    Input, real (kind = 8) deltaX, inverse orbit vector. 
    !
    !    Input integer (kind = 4) mx, my, number of BPMs in horizontal (x) and vertical (y) direction.
    !
    !    Input integer (kind = 4) nx, ny, number of steerers in horizontal (x) and vertical (y) direction.
    !
    !    Output, real ( kind = 8 ) kicks, vector containing the necessary steerer kicks (in rad). 
    !

        use bmad
        use get_pseudo_inverse_mod

        implicit none

        real ( kind = 8 ), allocatable, dimension ( : , : ) :: orm
        real ( kind = 8 ), allocatable, dimension ( : , : ) :: inv_orm
        real ( kind = 8 ), allocatable, dimension ( : ) :: deltaX
        real ( kind = 8 ), allocatable, dimension ( : ) :: kicks
        integer ( kind = 4 ) mx, my, nx, ny, m, n
        real ( kind = 8 ), allocatable, dimension ( :, : ) :: s
        real ( kind = 8 ), allocatable, dimension ( :, : ) :: u
        real ( kind = 8 ), allocatable, dimension ( :, : ) :: v


        m = mx+my
        n = nx+ny
        if (allocated(kicks)) then
            deallocate(kicks)
        end if
        allocate ( inv_orm(1:n,1:m) )
        allocate ( kicks(1:(nx+ny)) )
        allocate ( u(1:m,1:m) )
        allocate ( s(1:m,1:n) )
        allocate ( v(1:n,1:n) )

        call get_SVD(m, n, orm, u, s, v)
        call svd_product_test ( m, n, orm, u, s, v )
        call pseudo_inverse(m, n, u, s, v, inv_orm)


        kicks = matmul(inv_orm,deltaX)

        return

    end subroutine get_steerer_kicks

    subroutine apply_kicks(lat, loc_steer, kicks, nx, ny)

    !*****************************************************************************
    !
    ! APPPLY_KICKS applys the previously calculated steerer kicks to the steerer.
    !
    !  Parameters:
    !
    !    Input, lat_struct lat, contains all TWISS information about the underlying lattice.
    !
    !    Input, real ( kind = 8 ) kicks, vector containing the necessary steerer kicks (in rad).
    !
    !    Input integer (kind = 4) nx, ny, number of steerers in horizontal (x) and vertical (y) direction.
    !

        use bmad

        implicit none

        type (lat_struct) lat
        real ( kind = 8 ), allocatable, dimension ( : ) :: kicks
        integer ( kind = 4 ) nx
        integer ( kind = 4 ) ny
        integer ( kind = 4 ) n_ele
        integer ( kind = 4 ) i, count1, count2, loc_steer
        real ( kind = 8 ) a
        type (ele_struct), pointer :: ele

        n_ele = lat%n_ele_max
        count1 = 0
        count2 = 0
        a = 0.01

        do i = 1, n_ele
            ele => lat%ele(i)
    
            if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"H",.false.) == 1 .and. index(ele%name, "TOROID", .false.) /= 1 .and. index(ele%name, "MBU", .false.) /= 1 ) then

                count1 = count1 + 1
                if (count1 .EQ. loc_steer .and. loc_steer /= 0) then
                    count1 = count1 - 1
                    loc_steer = 0
                else if (abs(kicks(count1)) >= 0.000000001) then
                    ele%value(kick$) = ele%value(kick$)+kicks(count1)
                    call set_flags_for_changed_attribute (ele, ele%value(kick$))
                end if
            
            else if (index(key_name(ele%key),"kicker",.false.) == 2 .and. index(key_name(ele%key),"V",.false.) == 1 ) then

                count2 = count2 + 1
                if (count2+nx .EQ. loc_steer .and. loc_steer /= 0) then
                    count2 = count2 - 1
                    loc_steer = 0
                else if (abs(kicks(count2+nx)) >= 0.000000001) then
                    ele%value(kick$) = ele%value(kick$)+kicks(count2+nx)
                    call set_flags_for_changed_attribute (ele, ele%value(kick$))
                end if
        
            end if

        enddo

    end subroutine apply_kicks

    subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************
    !
    !! R8MAT_PRINT_SOME prints some of an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is a matrix of real ( kind = 8 ) values.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
    !
    !    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
    !
    !    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
        implicit none

        integer ( kind = 4 ), parameter :: incx = 5
        integer ( kind = 4 ) m
        integer ( kind = 4 ) n

        real ( kind = 8 ) a(m,n)
        character ( len = 14 ) ctemp(incx)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) i2hi
        integer ( kind = 4 ) i2lo
        integer ( kind = 4 ) ihi
        integer ( kind = 4 ) ilo
        integer ( kind = 4 ) inc
        integer ( kind = 4 ) j
        integer ( kind = 4 ) j2
        integer ( kind = 4 ) j2hi
        integer ( kind = 4 ) j2lo
        integer ( kind = 4 ) jhi
        integer ( kind = 4 ) jlo
        character ( len = * ) title

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) trim ( title )

        do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

            j2hi = j2lo + incx - 1
            j2hi = min ( j2hi, n )
            j2hi = min ( j2hi, jhi )

            inc = j2hi + 1 - j2lo

!            write ( *, '(a)' ) ' '

            do j = j2lo, j2hi
                j2 = j + 1 - j2lo
                write ( ctemp(j2), '(i8,6x)' ) j
            end do

!            write ( *, '(''  Col   '',5a14)' ) ctemp(1:inc)
!            write ( *, '(a)' ) '  Row'
!            write ( *, '(a)' ) ' '

            i2lo = max ( ilo, 1 )
            i2hi = min ( ihi, m )

            do i = i2lo, i2hi

                do j2 = 1, inc

                    j = j2lo - 1 + j2

                    if ( a(i,j) == real ( int ( a(i,j) ), kind = 8 ) ) then
                        write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
                    else
                        write ( ctemp(j2), '(g14.6)' ) a(i,j)
                    end if

                end do

!                write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

            end do

        end do

        return

    end subroutine r8mat_print_some

end module orbit_correction_mod
