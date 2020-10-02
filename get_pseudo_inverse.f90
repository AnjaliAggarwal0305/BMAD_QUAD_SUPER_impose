module get_pseudo_inverse_mod

contains

    subroutine get_SVD(m, n, a, u, s, v)

    !*****************************************************************************
    !
    !! R8MAT_SVD_LAPACK gets the SVD of a matrix using a call to LAPACK.
    !
    !  Discussion:
    !
    !    The singular value decomposition of a real MxN matrix A has the form:
    !
    !      A = U * S * V'
    !
    !    where
    !
    !      U is MxM orthogonal,
    !      S is MxN, and entirely zero except for the diagonal;
    !      V is NxN orthogonal.
    !
    !    Moreover, the nonzero entries of S are positive, and appear
    !    in order, from largest magnitude to smallest.
    !
    !    This routine calls the LAPACK routine DGESVD to compute the
    !    factorization.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
    !    decomposition we are investigating.
    !
    !    Output, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
    !    that form the singular value decomposition of A.
    !

        implicit none

        integer ( kind = 4 ) m
        integer ( kind = 4 ) n

        real ( kind = 8 ) a(m,n)
        real ( kind = 8 ) a_copy(m,n)
        integer ( kind = 4 ) i
        integer ( kind = 4 ) info
        integer ( kind = 4 ) lda
        integer ( kind = 4 ) ldu
        integer ( kind = 4 ) ldv
        character jobu
        character jobv
        integer ( kind = 4 ) lwork
        real ( kind = 8 ) sdiag(min(m,n))
        real ( kind = 8 ) s(m,n)
        real ( kind = 8 ) u(m,m)
        real ( kind = 8 ) v(n,n)
        real ( kind = 8 ), allocatable, dimension ( : ) :: work

        lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

        allocate ( work(1:lwork) )
    !
    !  Compute the eigenvalues and eigenvectors.
    !
        jobu = 'A'
        jobv = 'A'
        lda = m
        ldu = m
        ldv = n
    !
    !  The input matrix is destroyed by the routine.  Since we need to keep
    !  it around, we only pass a copy to the routine.
    !
        a_copy(1:m,1:n) = a(1:m,1:n)

        call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, lwork, info )

        if ( info /= 0 ) then
!            write ( *, '(a)' ) ' '
!            write ( *, '(a)' ) 'R8MAT_SVD_LAPACK - Failure!'
!            write ( *, '(a)' ) '  The SVD could not be calculated.'
!            write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
!            write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
            return
        end if
    !
    !  Make the MxN matrix S from the diagonal values in SDIAG.
    !
        s(1:m,1:n) = 0.0D+00
        do i = 1, min ( m, n )
            s(i,i) = sdiag(i)
        end do
    !
    !  Transpose V.
    !
        v = transpose ( v )

        deallocate ( work )

        return


    end subroutine get_SVD

    subroutine pseudo_inverse(m, n, u, s, v, a_pseudo)

    !*****************************************************************************
    !
    !! PSEUDO_INVERSE computes the pseudoinverse.
    !
    !  Discussion:
    !
    !    Given the singular value decomposition of a real MxN matrix A:
    !
    !      A = U * S * V'
    !
    !    where
    !
    !      U is MxM orthogonal,
    !      S is MxN, and entirely zero except for the diagonal;
    !      V is NxN orthogonal.
    !
    !    the pseudo inverse is the NxM matrix A+ with the form
    !
    !      A+ = V * S+ * U'
    !
    !    where
    !
    !      S+ is the NxM matrix whose nonzero diagonal elements are
    !      the inverses of the corresponding diagonal elements of S.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
    !    decomposition we are investigating.
    !
    !    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
    !    that form the singular value decomposition of A.
    !
    !    Output, real ( kind = 8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
    !
        implicit none

        integer ( kind = 4 ) m
        integer ( kind = 4 ) n

        real ( kind = 8 ) a_pseudo(n,m)
        integer ( kind = 4 ) i
        real ( kind = 8 ) s(m,n)
        real ( kind = 8 ) sp(n,m)
        real ( kind = 8 ) u(m,m)
        real ( kind = 8 ) v(n,n)

        sp(1:n,1:m) = 0.0D+00
        do i = 1, min ( m, n )
            if ( s(i,i) /= 0.0D+00 ) then
                sp(i,i) = 1.0D+00 / s(i,i)
            end if
        end do

        a_pseudo(1:n,1:m) = matmul ( v(1:n,1:n), matmul ( sp(1:n,1:m), transpose ( u(1:m,1:m) ) ) )

        return

    end subroutine pseudo_inverse


    subroutine print_matrix ( m, n, a, title )

    !*****************************************************************************
    !
    !! PRINT_MATRIX prints an R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is a matrix of real ( kind = 8 ) values.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, the number of rows in A.
    !
    !    Input, integer ( kind = 4 ) N, the number of columns in A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix.
    !
    !    Input, character ( len = * ) TITLE, a title.
    !
        implicit none

        integer ( kind = 4 ) m
        integer ( kind = 4 ) n

        real ( kind = 8 ) a(m,n)
        character ( len = * ) title

        call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

        return
    end subroutine print_matrix

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

    subroutine r8mat_uniform_01 ( m, n, seed, r )

    !*****************************************************************************
    !
    !! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
    !
    !  Discussion:
    !
    !    An R8MAT is a matrix of real ( kind = 8 ) values.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the array.
    !
    !    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
    !    should NOT be 0.  On output, SEED has been updated.
    !
    !    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
    !
        implicit none

        integer ( kind = 4 ) m
        integer ( kind = 4 ) n

        integer ( kind = 4 ) i
        integer ( kind = 4 ) j
        integer ( kind = 4 ) k
        integer ( kind = 4 ) seed
        real ( kind = 8 ) r(m,n)

        do j = 1, n

            do i = 1, m

                k = seed / 127773

                seed = 16807 * ( seed - k * 127773 ) - k * 2836
  
                    if ( seed < 0 ) then
                        seed = seed + 2147483647
                    end if

                r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10

            end do
        end do

        return
    end subroutine r8mat_uniform_01

    subroutine svd_product_test ( m, n, a, u, s, v )
    
    !*****************************************************************************
    !
    !! SVD_PRODUCT_TEST tests that A = U * S * V'.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
    !    decomposition we are investigating.
    !
    !    Input, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
    !    that form the singular value decomposition of A.
    !
        implicit none
      
        integer ( kind = 4 ) m
        integer ( kind = 4 ) n
      
        real ( kind = 8 ) a(m,n)
        real ( kind = 8 ) a_norm
        real ( kind = 8 ) dif_norm
        integer ( kind = 4 ) i
        real ( kind = 8 ) s(m,n)
        real ( kind = 8 ) u(m,m)
        real ( kind = 8 ) usv(m,n)
        real ( kind = 8 ) v(n,n)
      
        a_norm = sqrt ( sum ( a(1:m,1:n)**2 ) )
        usv(1:m,1:n) = matmul ( u(1:m,1:m), &
                       matmul ( s(1:m,1:n), transpose ( v(1:n,1:n) ) ) )
      
        call print_matrix ( m, n, usv, '  The product U * S * V'':' )
      
        dif_norm = sqrt ( sum ( ( a(1:m,1:n) - usv(1:m,1:n) )**2 ) )
!        write ( *, '(a)' ) ' '
!        write ( *, '(a,g14.6)' ) '  Frobenius Norm of A, A_NORM = ', a_norm
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '  ABSOLUTE ERROR for A = U*S*V'':'
!        write ( *, '(a,g14.6)' ) '  Frobenius norm of difference A-U*S*V'' = ', dif_norm
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '  RELATIVE ERROR for A = U*S*V'':'
!        write ( *, '(a,g14.6)' ) '  Ratio of DIF_NORM / A_NORM = ', dif_norm / a_norm
      
        return
    end subroutine svd_product_test

    subroutine pseudo_product_test ( m, n, a, a_pseudo )
    
    !*****************************************************************************
    !
    !! PSEUDO_PRODUCT_TEST examines pseudoinverse products.
    !
    !  Discussion:
    !
    !    Given an MxN matrix A, and its pseudoinverse A+, we must have
    !
    !      A+ * A * A+ = A+
    !      A * A+ * A = A
    !      ( A * A+ )' = A * A+ (MxM symmetry)
    !      ( A+ * A )' = A+ * A (NxN symmetry)
    !
    !    If M <= N, A * A+ may be "interesting" (equal to or "like" the identity),
    !    if N <= M, A+ * A may be "interesting" (equal to or "like" the identity).
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns
    !    in the matrix A.
    !
    !    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
    !    decomposition we are investigating.
    !
    !    Input, real ( kind = 8 ) A_PSEUDO(N,M), the pseudo_inverse of A.
    !
        implicit none
      
        integer ( kind = 4 ) m
        integer ( kind = 4 ) n
      
        real ( kind = 8 ) a(m,n)
        real ( kind = 8 ) a_pseudo(n,m)
        real ( kind = 8 ) bmm(m,m)
        real ( kind = 8 ) tbmm(m,m)
        real ( kind = 8 ) bmn(m,n)
        real ( kind = 8 ) bnm(n,m)
        real ( kind = 8 ) bnn(n,n)
        real ( kind = 8 ) tbnn(n,n)
        real ( kind = 8 ) dif1
        real ( kind = 8 ) dif2
        real ( kind = 8 ) dif3
        real ( kind = 8 ) dif4
        integer ( kind = 4 ) i
      
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) 'PSEUDO_PRODUCT_TEST'
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '  The following relations MUST hold:'
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '   A  * A+ * A  = A'
!        write ( *, '(a)' ) '   A+ * A  * A+ = A+'
!        write ( *, '(a)' ) ' ( A  * A+ ) is MxM symmetric;'
!        write ( *, '(a)' ) ' ( A+ * A  ) is NxN symmetric'
      
        bmn(1:m,1:n) = matmul ( a(1:m,1:n), &
                       matmul ( a_pseudo(1:n,1:m), a(1:m,1:n) ) )
      
        dif1 = sqrt ( sum ( ( a(1:m,1:n) - bmn(1:m,1:n) )**2 ) )
      
        bnm(1:n,1:m) = matmul ( a_pseudo(1:n,1:m), &
                       matmul ( a(1:m,1:n), a_pseudo(1:n,1:m) ) )
      
        dif2  =sqrt ( sum ( ( a_pseudo(1:n,1:m) - bnm(1:n,1:m) )**2 ) )

        bmm(1:m,1:m) = matmul ( a(1:m,1:n), a_pseudo(1:n,1:m) )
        tbmm = transpose(bmm)

        dif3 = sqrt ( sum ( ( bmm(1:m,1:m) - tbmm(1:m,1:m) )**2 ) )

        bnn(1:n,1:n) = matmul ( a_pseudo(1:n,1:m), a(1:m,1:n) )
        tbnn = transpose(bnn)

        dif4 = sqrt ( sum ( ( bnn(1:n,1:n) - tbnn(1:n,1:n) )**2 ) )

!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '  Here are the Frobenius norms of the errors'
!        write ( *, '(a)' ) '  in these relationships:'
!        write ( *, '(a)' ) ' '
!        write ( *, '(a,g14.6)' ) '   A  * A+ * A  = A            ', dif1
!        write ( *, '(a,g14.6)' ) '   A+ * A  * A+ = A+           ', dif2
!        write ( *, '(a,g14.6)' ) ' ( A  * A+ ) is MxM symmetric; ', dif3
!        write ( *, '(a,g14.6)' ) ' ( A+ * A  ) is NxN symmetric; ', dif4
      
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '  In some cases, the matrix A * A+'
!        write ( *, '(a)' ) '  may be interesting (if M <= N, then'
!        write ( *, '(a)' ) '  it MIGHT look like the identity.)'
!        write ( *, '(a)' ) ' '
      
        call print_matrix ( m, m, bmm, '  A * A+:' )
      
!        write ( *, '(a)' ) ' '
!        write ( *, '(a)' ) '  In some cases, the matrix A+ * A'
!        write ( *, '(a)' ) '  may be interesting (if N <= M, then'
!        write ( *, '(a)' ) '  it MIGHT look like the identity.)'
!        write ( *, '(a)' ) ' '
      
        call print_matrix ( n, n, bnn, '  A+ * A:' )
      
        return
    end subroutine pseudo_product_test

end module get_pseudo_inverse_mod
