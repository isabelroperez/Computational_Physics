program Project_1
    implicit none
    real*8 :: alpha, N, T, V, A, l, beta, W, WORK, O, suma, a_0
    integer :: m, i, j, k, q, ITYPE, LDA, LDB, LWORK, INFO, primeros, p
    character :: JOBZ, UPLO
    parameter (m = 12, LDA = m, LDB = m, LWORK = 3 * m**2 - 1, &
        JOBZ = 'V', UPLO = 'U', ITYPE = 1, a_0 = 1.d0, primeros = 3)
    dimension :: alpha(m), T(m, m), V(m, m), A(LDA, m), N(LDB, m), W(m), WORK(LWORK), O(m, m)
    
    !--------------------------------------------------MAIN PROGRAM START
    l = 0.d0
    call CoefficientsAlpha(alpha, m)
    
    do q = 0, 2
        print*, 'For l equal to', q
        
        beta = -1.d0
        call MatrixN(N, alpha, m, l)
        call MatrixT(N, alpha, m, T)
        call MatrixV(N, alpha, m, V, beta, l)
        call SumMatrixTMatrixV(m, T, V, A)
        call DSYGV(ITYPE, JOBZ, UPLO, m, A, LDA, N, LDB, W, WORK, LWORK, INFO)
        
        print*, 'The 3 eigenvalues of minimum energy are'
        write(*, '(20f10.3)') (W(i), i = 1, primeros)
        
        if (l == 1.d0) then
            call RadialFunction(alpha, m, l, A)
        end if
        
        beta = 1.d0
        
        do while (beta /= 3.d0)
            call MatrixN(N, alpha, m, l)
            call MatrixO(N, alpha, m, O, beta, l)
            do p = 1, primeros
                suma = 0.d0
                do j = 1, m
                    do i = 1, m
                        suma = suma + A(i, p) * O(j, i) * A(j, p)
                    end do
                end do
                if (beta == 1.d0) then
                    print*, '<r> ::', suma
                else
                    print*, '<r**2> ::', suma
                end if
            end do
            beta = beta + 1.d0
        end do
        
        l = l + 1.d0
    end do
    
    stop
end program

!-----------------------------------------------------------------SUBROUTINES    
!---------------------------------------------------------ALPHA COEFFICIENTS      
subroutine CoefficientsAlpha(alpha, m)
    implicit none
    real*8 :: alpha
    dimension alpha(m)
    integer :: i, m
    
    !print*,'Assigning coefficients alpha(i)...'
    
    open(10, file = 'alpha.txt') 
    do i = 1, m
        read(10, *) alpha(i)
    end do
    close(10)
    
    return
end subroutine

!--------------------------------------------------------------MATRIX N(j,i)
subroutine MatrixN(N, alpha, m, l)
    implicit none
    real*8 :: N, alpha, l
    dimension N(m, m), alpha(m)
    integer :: i, j, m
    
    !print*,'Generating coefficients N(j,i)...'
    
    do j = 1, m
        do i = 1, m
            N(j, i) = ((4.d0 * alpha(i) * alpha(j))**(l + 3.d0 / 2.d0)) / &
                      ((alpha(i) + alpha(j))**(2.d0 * l + 3.d0))
        end do
    end do
    
    return
end subroutine

!--------------------------------------------------------------MATRIX T(j,i)
subroutine MatrixT(N, alpha, m, T)
    implicit none
    real*8 :: N, alpha, T
    dimension N(m, m), T(m, m), alpha(m)
    integer :: i, j, m
    
    !print*,'Generating coefficients T(j,i)...'
    
    do j = 1, m
        do i = 1, m
            T(j, i) = N(j, i) * alpha(i) * alpha(j)
        end do
    end do     
    
    return
end subroutine

!--------------------------------------------------------------MATRIX V(j,i)
subroutine MatrixV(N, alpha, m, V, beta, l)
    implicit none
    real*8 :: N, alpha, V, beta, l
    dimension N(m, m), alpha(m), V(m, m)
    integer :: i, j, m
    
    !print*,'Generating coefficients V(j,i)...'
    
    do j = 1, m
        do i = 1, m
            V(j, i) = (-2.d0 * N(j, i) * gamma(3.d0 + 2.d0 * l + beta)) / &
                      (((alpha(i) + alpha(j))**beta) * gamma(3.d0 + 2.d0 * l))
        end do
    end do

    return
end subroutine

!--------------------------------------------------------------MATRIX O(j,i)
subroutine MatrixO(N, alpha, m, O, beta, l)
    implicit none
    real*8 :: N, alpha, O, beta, l
    dimension N(m, m), alpha(m), O(m, m)
    integer :: i, j, m

    !print*,'Generating coefficients O(j,i)...'
    
    do j = 1, m
        do i = 1, m
            O(j, i) = (N(j, i) * gamma(3.d0 + 2.d0 * l + beta)) / &
                      (((alpha(i) + alpha(j))**beta) * gamma(3.d0 + 2.d0 * l))
        end do
    end do
    
    return
end subroutine

!------------------------------------------------------MATRIX A(j,i)(LAPACK)              
subroutine SumMatrixTMatrixV(m, T, V, A)
    implicit none
    real*8 :: T, V, A
    dimension V(m, m), T(m, m), A(m, m)
    integer :: i, j, m

    !print*,'Generating coefficients A(j,i)=V(j,i)+T(j,i)...'
   
    j = 1
    do while (j /= m + 1)
        do i = 1, m 
            A(j, i) = V(j, i) + T(j, i)
        end do
        j = j + 1
    end do

    return
end subroutine

!-------------------------------------------------------------RADIAL FUNCTION
subroutine RadialFunction(alpha, m, l, A)
    real*8 :: l, alpha, x, funphis, barphi, A, h
    integer :: i, m, j, n_points, x_max
    dimension alpha(m), barphi(m), A(m, m)

    n_points = 300
    x_max = 30
    h = x_max / real(n_points)
    x = 0.d0
    
    open(45, file = 'radial_data.txt')
    do i = 0, n_points
        funphis = 0.d0
        do j = 1, m
            barphi(j) = ((((2.d0 * alpha(j))**3.d0) / (gamma(2.d0 * l + 2.d0)))**
 &                (1.d0 / 2.d0)) * (2.d0 * alpha(j) * x)**l * exp(-alpha(j) * x)
            funphis = funphis + A(j, 1) * barphi(j) 
        end do
        write(45, '(20f10.4)') x, funphis
        x = i * h 
    end do
    close(45)
    
    return
end subroutine
!-------------------------------------------------------------END SUBROUTINES