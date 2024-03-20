c     Program to obtain the solution of an ODE with mixed boundary conditions
c     linf, lsup: lower and upper limits of the interval
c     Differential equation: y''(x) + 9*y(x) = sin(x)
c     Nodes at the ends of the cells
c     (i=0,N) -> (x_0=linf, x_N=lsup): boundary points
c     A(i,j) -> Coefficient matrix

program MAIN
  integer, parameter :: N = 200
  double precision, parameter :: pi = 4 * atan(1.0d0)
  integer :: i, info
  double precision :: linf, lsup, x, h, yexact, error
  double precision :: A(0:N, 0:N), b(0:N, 1), y(0:N), ipiv(N+1)

  linf = 0.0d0
  lsup = 2.0d0
  h = (lsup - linf) / N ! cell size

c     Initialize coefficient matrix A_ij = 0
  do j = 0, N
    do i = 0, N
      A(i, j) = 0.0d0
    end do
  end do

c     Specify the non-zero coefficients of A and the vector b
  do i = 1, N-1
    x = linf + i * h ! relationship between i (discrete variable) and x (continuous variable)
    A(i, i) = ... ! complete
    A(i, i-1) = ... ! complete
    A(i, i+1) = ... ! complete
    b(i, 1) = ... ! complete
  end do

c     Coefficients to impose boundary conditions
  A(0, 0) = ... ! complete
  b(0, 1) = ... ! complete
  A(N, N) = ... ! complete
  b(N, 1) = ... ! complete

  call DGESV(N+1, 1, A, N+1, ipiv, b, N+1, info)

c     Write the solution
  open(1, file='sol1.dat')

  do i = 0, N
    x = linf + i * h
    write(1, *) x, b(i, 1)
  end do

  close(1)

end program


