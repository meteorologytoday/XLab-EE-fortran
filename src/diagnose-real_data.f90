program diagnose
use elliptic_tools
use field_tools
use constants
implicit none

integer        :: nr, nz
character(256) :: A_file, B_file, C_file, Q_file, input_folder, output_folder, output_file

real(4), pointer   :: psi(:,:), chi(:,:), f(:,:), Q(:,:), coe(:, :, :),    &
&                     a(:,:), b(:,:), c(:,:), workspace(:,:),              &
&                     aa(:,:), bb(:,:), cc(:,:), JJ(:,:),                  &
&                     w_mom(:,:), theta(:,:), eta(:,:), eta_tmp(:,:),          &
&                     ra(:), za(:), exner(:), sigma(:)

real(4)            :: testing_dt, Lr, Lz, dr, dz, eta_avg_b, eta_avg_nob


integer :: p_strategy, strategy, max_iter
real(4) :: p_strategy_r, strategy_r

integer :: i,j, m, n, err
real(4) :: r, z, eta_avg, sigmaa, tmp1, tmp2, tmp3

read(*,*) Lr, Lz;    read(*,*) nr, nz;
read(*,'(a)') input_folder
read(*,'(a)') output_folder
read(*,'(a)') A_file;
read(*,'(a)') B_file;
read(*,'(a)') C_file;
read(*,'(a)') Q_file;
read(*,*) p_strategy, p_strategy_r, max_iter;

print *, "Lr:", Lr, ", Lz:", Lz
print *, "nr:", nr, ", nz:", nz
print *, "Input folder: ", trim(input_folder)
print *, "Output folder: ", trim(output_folder)
print *, "A file: ", trim(A_file)
print *, "B file: ", trim(B_file)
print *, "C file: ", trim(C_file)
print *, "Q file: ", trim(Q_file)
print *, "Strategy wanted: ", p_strategy
print *, "Strategy specific value: ", p_strategy_r
print *, "Strategy max iteration: ", max_iter

print *, "Input done. Doing my job..."
allocate(psi(nr,nz));     allocate(chi(nr,nz));     allocate(eta(nr-1,nz));
allocate(f(nr,nz));       allocate(Q(nr-1,nz));     allocate(JJ(nr-1,nz));
allocate(workspace(nr,nz));
allocate(eta_tmp(nr-1,nz));                         allocate(coe(9,nr,nz));
allocate(a(nr-1, nz-2));  allocate(b(nr-1, nz-1));  allocate(c(nr-2, nz-1));
allocate(aa(nr-1, nz-2)); allocate(bb(nr-1, nz-1)); allocate(cc(nr-2, nz-1));
allocate(w_mom(nr-1,nz)); allocate(theta(nr-1,nz)); allocate(eta(nr-1,nz));
allocate(ra(nr));         allocate(za(nz));         allocate(sigma(nz));
allocate(exner(nz))

dr = Lr / (nr-1); dz = Lz / (nz-1);
! assign x and y
do i=1,nr
    ra(i) = (i-1.0) * dr
end do

! calculate sigma (pseudo-density)
do j=1,nz
    za(j) = (j-1.0) * dz
    exner(j) = 1.0 - za(j) / h0
    sigma(j) = p0 / (theta0 * Rd) * exner(j)**(1.0 / kappa - 1.0)
end do


call read_2Dfield(15, trim(input_folder)//"/"//A_file, a, nr-1, nz-2)
call read_2Dfield(15, trim(input_folder)//"/"//B_file, b, nr-1, nz-1)
call read_2Dfield(15, trim(input_folder)//"/"//C_file, c, nr-2, nz-1)
call read_2Dfield(15, trim(input_folder)//"/"//Q_file, Q, nr-1, nz  )

! normalize coefficient aa = a / (sigma * r), bb = b / (sigma * r)
!                     , cc = c / (sigma * r)
do i=1,nr-1
    do j=1,nz-2
        r = (i-0.5) * dr ! r
        sigmaa = sigma(j+1)
        aa(i,j) = a(i,j) / (r*sigmaa)
    end do
end do

do i=1,nr-1
    do j=1,nz-1
        r = (i-0.5) * dr
        sigmaa = (sigma(j) + sigma(j+1)) / 2.0
        bb(i,j) = b(i,j) / (r*sigmaa)
    end do
end do

do i=1,nr-2
    do j=1,nz-1
        r = i * dr ! r
        sigmaa = (sigma(j) + sigma(j+1)) / 2.0
        cc(i,j) = c(i,j) / (r*sigmaa)
    end do
end do

! assign JJ = Q / (Cp * Pi)
do i=1,nr-1
    do j=1,nz
        JJ(i,j) = Q(i,j) / (Cp * exner(j))
    end do
end do

output_file = trim(output_folder)//"/j.bin"
call write_2Dfield(11, output_file, JJ, nr-1, nz)



f=0.0
! assign f= g0/theta0 * dtheta/dr --- using b
do i=2,nr-1
    do j=2,nz-1
        f(i,j) = - (b(i,j) + b(i-1,j) + b(i-1,j-1) + b(i,j-1)) / 4.0
    end do
end do

! ===== [STAGE  I] ===== !

!!!!! Solve raw eta with bb /= 0
print *, "Baroclinic Experiment"
call cal_coe(aa, bb, cc, coe, dr, dz, nr, nz, err)
strategy = p_strategy; strategy_r = p_strategy_r;
chi=0.0 ! Note we arbitrary decide the boundary condition of chi
call solve_elliptic(max_iter, strategy, strategy_r, chi, coe, f, workspace, nr, nz, err, 0)
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(nr, nz, chi, eta, ra, sigma, exner)
call write_2Dfield(11, trim(output_folder)//"/eta-b.bin",eta,nr-1,nz)
eta_avg_b = cal_eta_avg(sigma, Q, eta, ra, za, nr, nz)

!!!!! Solve raw eta with bb = 0
print *, "Barotropic Experiment"
bb=0.0
call cal_coe(aa, bb, cc, coe, dr, dz, nr, nz, err)
strategy = p_strategy; strategy_r = p_strategy_r;
chi=0.0 ! Note we arbitrary decide the boundary condition of chi
call solve_elliptic(max_iter, strategy, strategy_r, chi, coe, f, workspace, nr, nz, err, 0)
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(nr, nz, chi, eta, ra, sigma, exner)
call write_2Dfield(11, trim(output_folder)//"/eta-nob.bin",eta,nr-1,nz)
eta_avg_nob = cal_eta_avg(sigma, Q, eta, ra, za, nr, nz)

open (unit=15,file=trim(output_folder)//"/efficiency.txt",action="write",status="replace")
write (15,*) "Average Efficiency without b (%): ", eta_avg_nob
write (15,*) "Average Efficiency with    b (%): ", eta_avg_b
close (15)


print *, "Average Efficiency without b (%): ", eta_avg_nob
print *, "Average Efficiency with    b (%): ", eta_avg_b

contains

real(4) function cal_eta_avg(sigma, Q, eta, ra, za, nr, nz)
implicit none
real(4) :: Q(nr-1, nz), eta(nr-1, nz), ra(nr), za(nz), sigma(nz)
integer :: nr, nz

real(4) :: sum_q, sum_eta_q, r, weight
integer :: i,j

sum_q = 0; sum_eta_q = 0;
do i = 1,nr-1
    do j = 1, nz
        r = (ra(i) + ra(i+1))/2
        weight = MERGE(0.5,1.0,j==1 .or. j==nz)
        sum_eta_q = sum_eta_q + eta(i,j) * Q(i,j) * sigma(j) * r * weight
        sum_q = sum_q + Q(i,j) * sigma(j) * r * weight
    end do
end do

print *, "sum_eta_q:", sum_eta_q, ", sum_q:", sum_q

cal_eta_avg = sum_eta_q / sum_q
end function

subroutine cal_eta(nr, nz, chi, eta, ra, sigma, exner)
implicit none
integer, intent(in) :: nr, nz
real(4), intent(in) :: chi(nr,nz), ra(nr), sigma(j), exner(j)
real(4), intent(out):: eta(nr-1,nz)

real(4) :: r,dr
integer :: i,j

do i=1,nr-1
    do j=1,nz
        r = (ra(i) + ra(i+1))/2
        dr = ra(i+1)-ra(i)
        eta(i,j) = (chi(i+1,j) - chi(i,j)) / (r*dr) * g0 / (sigma(j) * Cp * exner(j) * theta0)
    end do
end do

eta = eta * 100.0 ! in percent

end subroutine


end program
