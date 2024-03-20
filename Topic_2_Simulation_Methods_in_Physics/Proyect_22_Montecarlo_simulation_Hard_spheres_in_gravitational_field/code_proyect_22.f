program Slit
    ! Objective: Simulate a bulk of ideal gas in a three-dimensional space with a gravitational field.
    implicit none
    real*8 :: xold, yold, zold, tonto, bmg
    integer*8 :: idum, val, nmed, ncm
    real*8 :: paso, x(5000), y(5000), z(5000), rm, nal, box, boxz
    integer*8 :: i, j, k, intentos, nbolas, validos, valtot, nescr, npe
    integer*8 :: nbinz, binx(200), binz(200), nbinx, cero
    integer*8 :: binx2(200), binz2(200)
    character(3) :: cini

    ! Reading initial data
    call inicia(nbolas, box, boxz, cini, intentos, paso, nbinz, nbinx, nmed, nescr, bmg)
    npe = nmed / nescr
    print*, 'Writing every', npe, 'measurements'
    cero = 0

    ! Initializing random seeds
    CALL init_random_seed()

    ! Initializing positions
    if (cini == 'red') then
        call posini(x, y, z, box, boxz, nbolas)
    else if (cini == 'rnd') then
        call posinir(x, y, z, box, boxz, nbolas)
    else
        print*, 'Unknown initial method'
        stop
    end if
    call escribepos(x, y, z, nbolas, cero)

    ! Thermalization without measurement
    print*, 'Preheating...'
    call mueve(intentos / 10, validos, nbolas, x, y, z, box, boxz, paso, bmg)
    print*, 'Acceptance rate = ', float(validos) / (nbolas * intentos / 10)

    ! Running and measuring
    print*, 'Let"s go'

    ncm = intentos / nmed
    valtot = 0
    do j = 1, nmed
        call mueve(ncm, validos, nbolas, x, y, z, box, boxz, paso, bmg)
        valtot = valtot + validos
        call midez(z, binz, binz2, nbolas, nbinz, box, boxz, j)
        if (mod(j, npe) == 0) then
            print*, 'Writing files', j / npe, 'after measurement', j
            call escribez(binz, binz2, nbinz, box, boxz, j / npe, j, nbolas)
        end if
    end do
    cero = 9999
    call escribepos(x, y, z, nbolas, cero)
    print*, 'Acceptance rate = ', float(valtot) / (nbolas * intentos)

    stop
end program

! Subroutine to read initial data
subroutine inicia(nbolas, box, boxz, cini, intentos, paso, nbinz, nbinx, nmed, nescr, bmg)
    implicit none
    integer*8 :: nbolas, nmed, intentos, nescr, nbinz, nbinx
    real*8 :: box, boxz, paso, bmg
    character(3) :: cini

    open(34, file = 'entrada.txt')
    read(34, *) nbolas
    write(*, '(a24,i5)') 'Nparts = ', nbolas
    read(34, *) box
    write(*, '(a24,f8.3)') 'Box = ', box
    read(34, *) boxz
    write(*, '(a24,f8.3)') 'Height = ', boxz
    write(*, '(a24,f8.5)') 'Average Density = ', nbolas / box / box / boxz
    read(34, *) cini
    write(*, '(a24,a3)') 'Initial Placement = ', cini
    read(34, *) intentos
    write(*, '(a24,1p,e8.2)') 'No. of Attempts = ', float(intentos)
    read(34, *) paso
    write(*, '(a24,f8.5)') 'Step = ', paso
    read(34, *) nbinz
    write(*, '(a24,i5)') 'nbinz = ', nbinz
    read(34, *) nmed
    write(*, '(a24,1p,e8.2)') 'Nmed = ', float(nmed)
    read(34, *) nescr
    write(*, '(a24,1p,e8.2)') 'Nescr = ', float(nescr)
    read(34, *) bmg
    write(*, '(a24,1p,e8.2)') 'mg/kT = ', bmg
    close(34)

    return
end

! Subroutine for initial positions
subroutine posini(x, y, z, box, boxz, nbolas)
    implicit none
    real*8 :: x(5000), y(5000), z(5000), box, boxz, dist, distz
    integer*8 :: nbolas, n3, i, j, k, part

    n3 = int(nbolas ** 0.33333333d0) + 1
    dist = box / n3
    distz = (boxz - 1.d0) / n3

    if ((dist < 1.d0) .or. (distz < 1.d0)) then
        print*, 'Cannot fit in the box!'
        stop
    end if

    do i = 0, n3 - 1
        do j = 0, n3 - 1
            do k = 0, n3 - 1
                part = i * n3 * n3 + j * n3 + k + 1
                x(part) = k * dist
                y(part) = j * dist
                z(part) = 0.5d0 + i * distz
            end do
        end do
    end do

    return
end

! Subroutine for random initial positions
subroutine posinir(x, y, z, box, boxz, nbolas)
    implicit none
    real*8 :: x(5000), y(5000), z(5000), r2, r, box, boxz, tonto(3)
    real*8 :: xr, yr, zr
    integer*8 :: i, j, nbolas
    logical :: solap
    i = 1
    do while (i < nbolas + 1)
        call random_number(tonto)
        x(i) = tonto(1) * box
        y(i) = tonto(2) * box
        z(i) = 0.5d0 + tonto(3) * (boxz - 1.d0)
        solap = .false.
        do j = 1, i - 1
            xr = x(i) - x(j)
            xr = xr - box * nint(xr / box)
            yr = y(i) - y(j)
            yr = yr - box * nint(yr / box)
            zr = z(i) - z(j)
            r2 = xr * xr + yr * yr + zr * zr
            if (r2 < 1.d0) solap = .true.
        end do
        if (.not. solap) i = i + 1
    end do
    return
end

! Subroutine to write positions
subroutine escribepos(x, y, z, nbolas, j)
    implicit none
    real*8 :: x(5000), y(5000), z(5000)
    integer*8 :: i, j, nbolas
    character(4) :: anumer
    write(anumer, '(i4.4)') j
    open(34, file = 'pos' // anumer)
    do i = 1, nbolas
        write(34, *) x(i), y(i), z(i)
    end do
    close(34)
    return
end

! Subroutine for movement
subroutine mueve(intentos, validos, nbolas, x, y, z, box, boxz, paso, bmg)
    implicit none
    integer*8 :: intentos, validos, nbolas, i, j, k, sol_wall
    real*8 :: x(5000), y(5000), z(5000), box, boxz, paso, bmg
    real*8 :: xn, yn, zn, xr, yr, zr, r2, hi, hf, aleat
    logical :: acp_mov

    ! Random particle selection
    validos = 0
    do k = 1, nbolas * intentos
        call random_number(aleat)
        i = int(aleat * nbolas) + 1

        call random_number(aleat)
        zn = z(i) + (aleat - 5.d-1) * paso

        call random_number(aleat)
        xn = x(i) + (aleat - 5.d-1) * paso

        call random_number(aleat)
        yn = y(i) + (aleat - 5.d-1) * paso

        ! Applying Periodic Boundary Conditions (PBC) for movement
        xn = xn - box * floor(xn / box)
        yn = yn - box * floor(yn / box)

        ! Initial and final energy calculation
        hi = bmg * z(i)
        hf = bmg * zn

        if ((zn < (boxz - 5.d-1)) .and. (zn > (5.d-1))) then
            call random_number(aleat)
            if (aleat < exp(-(hf - hi))) then
                sol_wall = 1
                j = 1
                do while (j < nbolas)
                    if (i /= j) then
                        xr = xn - x(j)
                        xr = xr - box * nint(xr / box)
                        yr = yn - y(j)
                        yr = yr - box * nint(yr / box)
                        zr = zn - z(j)
                        r2 = xr**2 + yr**2 + zr**2

                        ! Applying PBC for distance
                        xr = xr - box * nint(xr / box)
                        yr = yr - box * nint(yr / box)
                        zr = zr - boxz * nint(zr / boxz)

                        if (r2 < 1.d0) then
                            sol_wall = 0
                        end if
                    end if
                    j = j + 1
                end do
                if (sol_wall == 1) then
                    x(i) = xn
                    y(i) = yn
                    z(i) = zn
                    validos = validos + 1
                end if
            end if
        end if
    end do

    return
end

! Subroutine for measuring in Z
subroutine midez(z, binz, binz2, nbolas, nbinz, box, boxz, k)
    implicit none
    real*8 :: x(5000), y(5000), z(5000), box, boxz, tbin, suma
    integer*8 :: binz(200), nbolas, nbinz, k
    integer*8 :: i, part(200), ib
    integer*8 :: binz2(200)

    ! bin size (tbin) = boxz / nbinz
    ! ib iterates over the bins and assigns balls when appropriate
    ! binz2 is used to calculate dispersion, but not utilized here

    tbin = (boxz / nbinz)
    if ((k == 0) .or. (k == 1)) then
        do i = 1, nbinz
            binz(i) = 0
        end do
    end if
    do i = 1, nbolas
        ib = int(z(i) / tbin) + 1
        binz(ib) = binz(ib) + 1
    end do

    return
end

! Subroutine for writing histograms on the Z axis
subroutine escribez(binz, binz2, nbinz, box, boxz, numer, k, nbolas)
    implicit none
    real*8 :: box, boxz, bin, distz, tbin, suma
    integer*8 :: nbinz, binz(200), binz2(200), k, nbolas, numer, i
    character(4) :: anumer

    tbin = (boxz / nbinz)

    ! File name followed by a four-digit integer (i4)
    write(anumer, '(i4)') numer
    open(34, file = 'binz' // anumer)
    suma = 0.d0

    do i = 1, nbinz
        suma = (suma + binz(i))

        ! Writing histogram: bin reference point and normalized number of particles inside
        write(34, *) (i - (5.d-1)) * tbin, binz(i) / (tbin * box * box * k)
    end do

    close(34)

    return
end

! Subroutine to initialize the random number generator
subroutine init_random_seed()
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid, t(2), s
    integer(8) :: count, tms

    call random_seed(size = n)
    allocate(seed(n))

    ! First try if the OS provides a random number generator
    open(newunit = un, file = "/dev/urandom", access = "stream", form = "unformatted", action = "read", status = "old", iostat = istat)
    if (istat == 0) then
        read(un) seed
        close(un)
    else
        ! Fallback to XORing the current time and pid.
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
end subroutine init_random_seed

      

