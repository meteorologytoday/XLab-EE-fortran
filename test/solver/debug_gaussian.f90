program debug
use elliptic_tools
use field_tools
implicit none


integer, parameter :: nx = 201, ny = 201, max_iter = 5000
real(4), parameter :: res_abs=1e-5, res_rel=1e-5, &
    &                 PI = acos(-1.0), PII = 2*PI,
    &                 Lx=600000.0, Ly=600000.0, dx = Lx/(nx-1), dy=Ly/(ny-1)

real(4) :: psi(nx,ny), f(nx,ny), coe(9, nx, ny), &
    &      a(nx-1, ny-2), b(nx-1, ny-1), c(nx-2, ny-1), workspace(nx,ny)
real(4) :: ideal_psi(nx,ny), ideal_f(nx,ny)

real(4) :: res_abs_now, res_rel_now

integer :: i,j, err
real(4) :: x, y, rr1, rr2, cntx1, cntx2, cnty1, cnty2, d1, d2, d3

a(:,:) = 1.0
b(:,:) = 0.0
c(:,:) = 1.0

! psi = 2 gaussians
! A = 1.0
! B = sin(x) * sin(y)
! C = 1.0
cntx1 = 0.3 * Lx
cnty1 = 0.5 * Ly
cntx2 = 0.7 * Lx
cnty2 = 0.5 * Ly

do i=1,nx
    do j=1,ny
        x = (i-1) * dx
        y = (j-1) * dy
        rr1 = (x - cntx1)**2 + (y - cnty1)**2
        rr2 = (x - cntx2)**2 + (y - cnty2)**2
    
        ideal_psi(i,j) = 100.0 * exp(-rr1 / (0.05*Ly)**2) &
&                       + 80.0 * exp(-rr2 / (0.05*Ly)**2)

    end do
end do

call write_2Dfield(11, "debug_ideal_psi.bin", ideal_psi, nx, ny)

call cal_coe(a, b, c, coe, dx, dy, nx, ny, err)
call do_elliptic(ideal_psi, coe, f, nx, ny, err)
call write_2Dfield(11, "debug_cal_f_from_ideal_psi.bin", f, nx, ny)

d1 = cal_max_dist(f, ideal_f, nx*ny)

psi=0
print *, "Solve cal_f"
strategy_now = strategy
strategy_r_now = strategy_r
call solve_elliptic(max_iter, strategy_now,  strategy_r_now, psi, coe, f, workspace, nx, ny, err, 1)

call write_2Dfield(11, "debug_cal_psi_from_cal_f.bin", psi, nx, ny)
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
