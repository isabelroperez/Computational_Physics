! Solving the diffusion equation using Crank-Nicolson
PROGRAM main
  implicit none
  double precision, parameter :: pi = 4.0 * datan(1.0d0)
  integer, parameter :: N = 100  ! Number of cells
  integer :: i, j, k, nitmax, nskip, info
  double precision :: L, dx, dt, diff, u1, u2, sigma, tmax, alpha, xmin
  double precision :: u(0:N), uold(0:N), A(0:N,0:N), b(0:N,1), ipiv(N+1)
  double precision :: media(0:N), maxima(0:N), minima(0:N)
  ! uold(i) = u_i^n, u(i) = u_i^(n+1)

  ! Problem parameters
  L = 4.0d0    ! Depth
  dt = 2.5e-4  ! Year (time step)
  diff = 24.0d-1 ! m^2/year (diffusion constant)
  u1 = 10.0d0 ! Dirichlet condition (homogeneous surface temperature at x=0)
  tmax = 4.0d0 ! Year
  nskip = 100 ! Write solution every nskip iterations

  dx = L / N ! m (spatial step)
  nitmax = int(tmax / dt)
  alpha = diff * dt / dx**2 ! Convenience
  xmin = 0.0d0   ! Minimum depth

  ! Initial condition u(x, t=0)
  do i = 0, N
    uold(i) = 10.0d0 ! Initial temperature at all depths
    media(i) = 0.0d0 ! Initialization for average temperature calculation (Tmed)
    maxima(i) = 0.0d0 ! Initialization for maximum temperature (Tmax)
    minima(i) = 100.0d0 ! Initialization for minimum temperature (Tmin, above 10 minimum)
  end do

  ! Write initial state
  open(100, FILE='T1.dat')
  open(101, FILE='T2.dat')
  open(102, FILE='T3.dat')
  do i = 0, N
    write(100, *) 0.0, i * dx, uold(i)
  end do

  ! Main loop
  do k = 1, nitmax
    ! Initialize coefficient matrix A_ij = 0
    do j = 0, N
      do i = 0, N
        A(i, j) = 0.0d0
      end do
    end do

    ! Boundary conditions
    b(0, 1) = 10.0d0 + 25.0d0 * sin(2.0d0 * pi * dt * k)
    b(N, 1) = 0.0d0
    A(0, 0) = 1.0d0
    A(N, N) = 1.0d0
    A(N, N-1) = -1.0d0

    ! Non-zero elements of A and b
    do i = 1, N-1
      A(i, i) = 1.0d0 + alpha
      A(i, i-1) = -alpha / 2.0d0
      A(i, i+1) = -alpha / 2.0d0
      b(i, 1) = (alpha / 2.0d0) * (uold(i-1) + uold(i+1)) + (1.0d0 - alpha) * uold(i)
    end do

    ! Use LAPACK
    call DGESV(N+1, 1, A, N+1, ipiv, b, N+1, info)

    ! Update uold(i) and u(i) for the next iteration
    do i = 0, N
      uold(i) = b(i, 1)
    end do

    ! Write to file if (k/nskip) is an integer
    if (mod(k, nskip) .eq. 0) then
      do i = 0, N
        write(100, *) k * dt, i * dx, uold(i) ! 3 columns: t, x, u
        ! Calculate average temperature (Tmed)
        media(i) = media(i) + uold(i)
        ! Conditions to write Tmax, Tmin in T2.dat
        ! Compare maxima with uold for each i, store the last one
        if (maxima(i) < uold(i)) then
          maxima(i) = uold(i)
        end if
        if (minima(i) > uold(i)) then
          minima(i) = uold(i)
          ! Check if Tmin is less than 0 (freezing)
          ! The last minima(i) stored is the lowest negative temperature
          ! The corresponding x will be (i+1)
          if (minima(i) < 0.0d0) then
            xmin = (i + 1) * dx
          end if
        end if
        ! For x=0.6=i*dx, dx fixed then i=(0.6/dx)+1, analogously for the rest
        write(101, *) k * dt, uold(1), uold(16), uold(26), uold(51)
      end do
    end if
  end do

  do i = 0, N
    write(102, *) i * dx, maxima(i), minima(i), (media(i) * dt * N) / tmax
  end do

  close(100)
  close(101)
  close(102)
  print*, 'Minimum depth before freezing:'
  print*, xmin

END













