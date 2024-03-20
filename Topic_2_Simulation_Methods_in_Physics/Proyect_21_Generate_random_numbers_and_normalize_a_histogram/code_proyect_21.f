program RANDOM_NUMBER_GENERATOR
    implicit none
    integer :: i, m, j, n
    integer*8 :: px
    real*8 :: z, x, f, pi, a, b, h, x_i, x_medio, px_normalized
    parameter (m = 10**4, pi = 3.14159, n = 100)
    dimension :: x(m), z(m), px(m), x_i(m), x_medio(m), px_normalized(m)

    b = 5.d-1 * pi
    a = -5.d-1 * pi  
    h = (b - a) / n

    do i = 1, n + 1
        x_i(i) = -pi / 2.d0 + (i - 1) * h 
    end do

    do i = 1, n
        x_medio(i) = (x_i(i + 1) + x_i(i)) / 2.d0
    end do

    call init_random_seed()
    call RANDOM_NUMBER(z)

    do i = 1, n
        px(i) = 0
        px_normalized(i) = 0
    end do  
    
    do i = 1, m
        x(i) = f(z(i))
        j = int((x(i) - a) / h)
        px(j + 1) = px(j + 1) + 1
        px_normalized(j + 1) = px(j + 1) / (m * h)
    end do

    open(45, file = 'bin_data_10.txt')
    do i = 1, n
        write(45, *) x_medio(i), px_normalized(i)
    end do     
    close(45)

    suma = 0.d0
    do i = 1, n
        suma = suma + px_normalized(i) * h
    end do
    print*, 'if it equals 1 it is correct', suma
    
    stop
end program

! Subroutine to initialize the random number generator
subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = n)
    allocate(seed(n))
    
    ! First try if the OS provides a random number generator
    open(newunit = un, file = "/dev/urandom", access = "stream", &
        form = "unformatted", action = "read", status = "old", iostat = istat)
    if (istat == 0) then
        read(un) seed
        close(un)
    else
        ! Fallback to XOR:ing the current time and pid.
        call system_clock(count)
        if (count /= 0) then
            t = transfer(count, t)
        else
            call date_and_time(values = dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                + dt(5) * 60 * 60 * 1000 &
                + dt(6) * 60 * 1000 + dt(7) * 1000 &
                + dt(8)
            t = transfer(tms, t)
        end if
        s = ieor(t(1), t(2))
        pid = getpid() + 1099279 ! Add a prime
        s = ieor(s, pid)
        if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
                seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
        else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
        end if
    end if
    call random_seed(put = seed)
end subroutine

! Function definition for f(z)
function f(z)
    real*8 :: f, z
    f = asin(2.d0 * z - 1.d0)
    return
end function

