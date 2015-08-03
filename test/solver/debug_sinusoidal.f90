program debug_sinusoidal
use elliptic_tools
use field_tools
implicit none

integer, parameter :: nx = 200, ny = 200, max_iter = 100000, check_step = 100, &
    &                 converge_time = 10, lost_rate = 5, debug = 2

real(4), parameter :: res_abs=1e-5, res_rel=3e-2, alpha = 1.0, &
    &                 PI = acos(-1.0), PII = 2.0*PI, &
    &                 Lx=600000.0, Ly=600000.0, dx = Lx/(nx-1), dy=Ly/(ny-1)

real(4) :: psi(nx,ny), f(nx,ny), coe(9, nx, ny), &
    &      a(nx-1, ny-2), b(nx-1, ny-1), c(nx-2, ny-1), workspace(nx,ny), &
    &      ideal_psi(nx,ny), ideal_f(nx,ny)

integer :: max_iter_now
real(4) :: res_abs_now, res_rel_now

integer :: i,j, err
real(4) :: x, y, kx, ky, d1, d2, d3

a(:,:) = 1.0
b(:,:) = 0.0
c(:,:) = 1.0

kx =  2.0 * PI / Lx
ky =  3.0 * PI / Ly
do i=1,nx
    do j=1,ny
        x = (i-1) * dx
        y = (j-1) * dy
    
        ideal_psi(i,j) = sin(x * kx) * sin(y * ky)
        ideal_f(i,j)   = - (kx**2.0 + ky**2.0) * sin(x * kx) * sin(y * ky)
        
    end do
end do

call write_2Dfield(11, "ideal_psi.bin", ideal_psi, nx, ny)
call write_2Dfield(11, "ideal_f.bin", ideal_f, nx, ny)

call cal_coe(a, b, c, coe, dx, dy, nx, ny, err)
call do_elliptic(ideal_psi, coe, f, nx, ny, err)
call write_2Dfield(11, "cal_f_from_ideal_psi.bin", f, nx, ny)

d1 = cal_max_dist(f, ideal_f, nx*ny)

psi=0
print *, "Solve cal_f"
res_abs_now = res_abs
res_rel_now = res_rel
max_iter_now = max_iter
call solve_elliptic(max_iter_now, check_step, converge_time, lost_rate,     &
    &               res_abs_now, res_rel_now, alpha, psi, coe, f,       &
    &               workspace, nx, ny, err, debug)
print *, "max_iter:", max_iter_now, ", res_abs:", res_abs, ", res_rel:", res_rel
call write_2Dfield(11, "cal_psi_from_cal_f.bin", psi, nx, ny)
d2 = cal_max_dist(psi, ideal_psi, nx*ny)

print *, "|ideal_f - cal_f_from_ideal_psi| = ", d1 
print *, "|ideal_psi - cal_psi_from_cal_f| = ", d2

contains

real(4) function cal_max_dist(dat1, dat2, n)
implicit none
real(4) :: dat1(n), dat2(n)
integer :: n

cal_max_dist = maxval( abs(dat1 - dat2))

end function


real(4) function cal_avg_dist(dat1, dat2, n)
implicit none
real(4) :: dat1(n), dat2(n)
integer :: n

cal_avg_dist = sum( (dat1 - dat2)**2) / n

end function

end program
