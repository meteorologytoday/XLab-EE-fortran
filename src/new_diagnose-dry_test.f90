program diagnose
use elliptic_tools
use field_tools
use constants
implicit none

integer        :: nr, nz
character(256) :: A_file, B_file, C_file, Q_file, input_folder, output_folder, &
&                 output_file, yes_or_no, previous_path, chi_bc_file, psi_bc_file

real(4), pointer   :: psi(:,:), chi(:,:), f(:,:), Q(:,:), coe(:, :, :),    &
&                     a(:,:), b(:,:), c(:,:), workspace(:,:),              &
&                     solver_a(:,:), solver_b(:,:), solver_c(:,:), JJ(:,:),&
&                     w_mom(:,:), u_mom(:,:), theta(:,:), eta(:,:), eta_tmp(:,:),&
&                     ra(:), za(:), exner(:), sigma(:), dthetadr(:,:),     &
&                     psi_bc(:,:), chi_bc(:,:)

real(4)            :: testing_dt, Lr, Lz, dr, dz, eta_avg_b, eta_avg_nob

integer :: p_strategy, strategy, max_iter
real(4) :: p_strategy_r, strategy_r

integer :: i,j, m, n, err
real(4) :: r, z, eta_avg, tmp1, tmp2, tmp3
logical :: use_previous, file_exists, use_chi_bc, use_psi_bc

read(*,'(a)') yes_or_no
if(yes_or_no == 'yes') then
    read(*,'(a)') previous_path
    use_previous = .true.
else
    use_previous = .false.
end if
read(*,*) testing_dt
read(*,*) Lr, Lz;    read(*,*) nr, nz;
read(*,'(a)') input_folder
read(*,'(a)') output_folder
read(*,'(a)') A_file;
read(*,'(a)') B_file;
read(*,'(a)') C_file;
read(*,'(a)') Q_file;
read(*,*) p_strategy, p_strategy_r, max_iter;
read(*, '(a)') yes_or_no
if(yes_or_no == 'yes') then
    read(*,'(a)') psi_bc_file;
    use_psi_bc = .true.
else
    use_psi_bc = .false.
end if
read(*, '(a)') yes_or_no
if(yes_or_no == 'yes') then
    read(*,'(a)') chi_bc_file
    use_chi_bc = .true.
else
    use_chi_bc = .false.
end if



print *, "Use previous relaxation: ", merge('Y'//'('// trim(previous_path) //')', 'N', use_previous .eqv. .true.)
print *, "Testing time: ", testing_dt
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
print *, "Use PSI boundary condition: ", merge('Y'//'('//trim(psi_bc_file)//')', 'N', use_psi_bc .eqv. .true.)
print *, "Use CHI boundary condition: ", merge('Y'//'('//trim(chi_bc_file)//')', 'N', use_chi_bc .eqv. .true.)


allocate(psi(nr,nz));   allocate(chi(nr,nz));      allocate(eta(nr,nz));
allocate(f(nr,nz));     allocate(Q(nr,nz));        allocate(JJ(nr,nz));
allocate(workspace(nr,nz));

allocate(eta_tmp(nr-1,nz));                        allocate(coe(9,nr,nz));
allocate(a(nr, nz));allocate(b(nr, nz));   allocate(c(nr, nz));
allocate(solver_a(nr-1, nz-2));allocate(solver_b(nr-1, nz-1));allocate(solver_c(nr-2, nz-1));

allocate(w_mom(nr,nz)); allocate(u_mom(nr,nz));
allocate(theta(nr,nz));  allocate(eta(nr,nz));
allocate(ra(nr));       allocate(za(nz));          allocate(sigma(nz));
allocate(exner(nz));    allocate(dthetadr(nr,nz));

call read_2Dfield(15, trim(input_folder)//"/"//A_file, a, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//B_file, b, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//C_file, c, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//Q_file, Q, nr, nz)

if(use_psi_bc .eqv. .true.) then
    allocate(psi_bc(nr,nz));
    call read_2Dfield(15, trim(input_folder)//"/"//psi_bc_file, psi_bc, nr, nz)
end if


if(use_chi_bc .eqv. .true.) then
    allocate(chi_bc(nr,nz));
    call read_2Dfield(15, trim(input_folder)//"/"//chi_bc_file, chi_bc, nr, nz)
end if

dr = Lr / (nr-1); dz = Lz / (nz-1);
! assign x and y
do i=1,nr
    ra(i) = (i-1) * dr
end do

! calculate sigma (pseudo-density)
do j=1,nz
    za(j) = (j-1) * dz
    exner(j) = 1.0 - za(j) / h0
    sigma(j) = p0 / (theta0 * Rd) * exner(j)**(1.0 / kappa - 1.0)
end do



! normalize coefficient solver_a = a / (sigma * r), solver_b = b / (sigma * r)
!                     , solver_c = c / (sigma * r)

call A2solverA
call B2solverB
call C2solverC

! assign JJ = Q / (Cp * Pi)
do i=1,nr
    do j=1,nz
        JJ(i,j) = Q(i,j) / (Cp * exner(j))
    end do
end do

output_file = trim(output_folder)//"/j.bin"
call write_2Dfield(11, output_file, JJ, nr-1, nz)


f=0.0
! assign f= g0/theta0 * dJJ/dr
call d_dr(JJ, f)
f = f * (g0 * theta0)

call write_2Dfield(11, trim(output_folder)//"/djdr.bin", f, nr, nz)

call cal_coe(solver_a, solver_b, solver_c, coe, dr, dz, nr, nz, err)

! ===== [STAGE   I] ===== !
inquire(file=trim(previous_path)//'/psi.bin', exist=file_exists)
if(file_exists .eqv. .true.) then
    print *, "Using previous psi.bin..."
    call read_2Dfield(15, trim(previous_path)//'/psi.bin', psi, nr, nz)
else
    psi = 0.0
end if

if(use_psi_bc .eqv. .true.) then
    psi = psi_bc
end if

print *, "Solving psi..."
strategy = p_strategy; strategy_r = p_strategy_r;
call solve_elliptic(max_iter, strategy, strategy_r, psi, coe, f, workspace, nr, nz, err, 0)
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

output_file = trim(output_folder)//"/psi.bin"
call write_2Dfield(11, output_file, psi, nr, nz)


! ===== [STAGE  II] ===== !

! calculate w and dtheta_dt
call d_rdr(psi, w_mom)
call d_dz(psi, u_mom); u_mom = -u_mom
do i=1,nr
    do j=1,nz
        ! dtheta/dt = J - w dtheta/dz - u dtheta/dr
        theta(i,j) = JJ(i,j) - theta0/g0 * a(i,j) * w_mom(i,j)/sigma(j) &
            &                + theta0/g0 * b(i,j) * u_mom(i,j)/sigma(j)
        
        w_mom(i,j) = w_mom(i,j)/sigma(j)
        u_mom(i,j) = u_mom(i,j)/sigma(j)
    end do
end do

call write_2Dfield(11, trim(output_folder)//"/w.bin",w_mom,nr,nz)
call write_2Dfield(11, trim(output_folder)//"/u.bin",u_mom,nr,nz)
call write_2Dfield(11, trim(output_folder)//"/dtheta_dt.bin",theta,nr,nz)

! Calculate forced B field, A and C assumed to be unperturbed for now
! tmp1 --- [center] --- tmp2
! i+m         i         i+n
call d_dr(theta, workspace)
theta = workspace * testing_dt
b = b - g0 / theta0 * theta

call B2solverB
call write_2Dfield(11, trim(output_folder)//"/coe-new-b.bin", b, nr, nz)

! ===== [STAGE III : Solve] ===== !

f = - theta0/g0 * b

print *, "Solving CHI with B=0"
workspace(1:nr-1,1:nz-1) = solver_b;solver_b = 0
call cal_coe(solver_a, solver_b, solver_c, coe, dr, dz, nr, nz, err)
solver_b = workspace(1:nr-1,1:nz-1) ! restore

inquire(file=trim(previous_path)//'/CHI-nob.bin', exist=file_exists)
if(file_exists .eqv. .true.) then
    print *, "Using previous CHI-nob.bin..."
    call read_2Dfield(15, trim(previous_path)//'/CHI-nob.bin', chi, nr, nz)
else
    chi=0.0
end if

if(use_chi_bc .eqv. .true.) then
    chi = chi_bc
end if


strategy = p_strategy; strategy_r = p_strategy_r;
call solve_elliptic(max_iter, strategy, strategy_r, chi, coe, f, workspace, nr, nz, err, 0)
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(nr, nz, chi, eta, ra, sigma, exner);
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta-nob.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/CHI-nob.bin",chi,nr,nz)

eta_avg_nob = cal_eta_avg(sigma, Q, eta, ra, za, nr, nz)

print *, "Solving CHI with B!=0"
call cal_coe(solver_a, solver_b, solver_c, coe, dr, dz, nr, nz, err)
strategy = p_strategy; strategy_r = p_strategy_r;
call solve_elliptic(max_iter, strategy, strategy_r, chi, coe, f, workspace, nr, nz, err, 0)
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(nr, nz, chi, eta, ra, sigma, exner);
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta-b.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/CHI-b.bin",chi,nr,nz)


eta_avg_b = cal_eta_avg(sigma, Q, eta, ra, za, nr, nz)

print *, "Average Efficiency without b (%): ", eta_avg_nob
print *, "Average Efficiency with    b (%): ", eta_avg_b

open (unit=15,file=trim(output_folder)//"/efficiency.txt",action="write",status="replace")
write (15,*) "Average Efficiency witout b (%): ", eta_avg_nob
write (15,*) "Average Efficiency with   b (%): ", eta_avg_b
close (15)


contains

subroutine d_dz(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr,nz)
integer :: i,j,m,n
do i = 1, nr
    do j = 1, nz
        if(j == 1) then
            m = 0; n = 1
        else if(j == nz) then
            m = -1; n = 0
        else
            m = -1; n = 1
        end if
        to_dat(i,j) = (from_dat(i,j+m) - from_dat(i,j+n)) / (za(j+m) - za(j+n)) 
    end do
end do
end subroutine

subroutine d_dr(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr,nz)
integer :: i,j,m,n
do i = 1, nr
    do j = 1, nz
        if(i == 1) then
            m = 0; n = 1
        else if(i == nr) then
            m = -1; n = 0
        else
            m = -1; n = 1
        end if
        to_dat(i,j) = (from_dat(i+m,j) - from_dat(i+n,j)) / (ra(i+m) - ra(i+n)) 
    end do
end do
end subroutine

subroutine d_rdr(from_dat, to_dat)
! notice that all the points with r=0 become inifinite
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr,nz)
integer :: i,j,m,n

call d_dr(from_dat, to_dat)
do i = 1, nr
    do j = 1, nz
        to_dat(i,j) = to_dat(i,j) / ra(i) 
    end do
end do
end subroutine


real(4) function cal_eta_avg(sigma, Q, eta, ra, za, nr, nz)
implicit none
real(4) :: Q(nr, nz), eta(nr, nz), ra(nr), za(nz), sigma(nz)
integer :: nr, nz

real(4) :: sum_q, sum_eta_q, r, weight
integer :: i,j

sum_q = 0; sum_eta_q = 0;
do i = 2,nr
    do j = 1, nz
        r = ra(i)
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

call d_rdr(chi, eta)

do i=1,nr
    do j=1,nz
        eta(i,j) = eta(i,j) * g0 / (sigma(j) * Cp * exner(j) * theta0)
    end do
end do

eta = eta * 100.0 ! in percent

end subroutine

subroutine A2solverA
implicit none
integer :: i,j

do i=1,nr-1
    do j=1,nz-2
        solver_a(i,j) = (a(i  ,j+1)/ra(i  )/sigma(j+1)  &
           &           + a(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
    end do
end do

do i=1,nr-1
    do j=1,nz-1
        solver_b(i,j) = ( b(i  ,j  )/ra(i  )/sigma(j  )  &
           &            + b(i+1,j  )/ra(i+1)/sigma(j  )  &
           &            + b(i  ,j+1)/ra(i  )/sigma(j+1)  &
           &            + b(i+1,j+1)/ra(i+1)/sigma(j+1)) / 4.0
    end do
end do

do i=1,nr-2
    do j=1,nz-1
        r = i * dr ! r
        solver_c(i,j) = (c(i+1,j  )/ra(i+1)/sigma(j  )  &
           &           + c(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
    end do
end do
end subroutine


subroutine B2solverB
implicit none
integer :: i,j

do i=1,nr-1
    do j=1,nz-2
        solver_a(i,j) = (a(i  ,j+1)/ra(i  )/sigma(j+1)  &
           &           + a(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
    end do
end do

do i=1,nr-1
    do j=1,nz-1
        solver_b(i,j) = ( b(i  ,j  )/ra(i  )/sigma(j  )  &
           &            + b(i+1,j  )/ra(i+1)/sigma(j  )  &
           &            + b(i  ,j+1)/ra(i  )/sigma(j+1)  &
           &            + b(i+1,j+1)/ra(i+1)/sigma(j+1)) / 4.0
    end do
end do

do i=1,nr-2
    do j=1,nz-1
        r = i * dr ! r
        solver_c(i,j) = (c(i+1,j  )/ra(i+1)/sigma(j  )  &
           &           + c(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
    end do
end do
end subroutine


subroutine C2solverC
implicit none
integer :: i,j

do i=1,nr-1
    do j=1,nz-2
        solver_a(i,j) = (a(i  ,j+1)/ra(i  )/sigma(j+1)  &
           &           + a(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
    end do
end do

do i=1,nr-1
    do j=1,nz-1
        solver_b(i,j) = ( b(i  ,j  )/ra(i  )/sigma(j  )  &
           &            + b(i+1,j  )/ra(i+1)/sigma(j  )  &
           &            + b(i  ,j+1)/ra(i  )/sigma(j+1)  &
           &            + b(i+1,j+1)/ra(i+1)/sigma(j+1)) / 4.0
    end do
end do

do i=1,nr-2
    do j=1,nz-1
        r = i * dr ! r
        solver_c(i,j) = (c(i+1,j  )/ra(i+1)/sigma(j  )  &
           &           + c(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
    end do
end do
end subroutine

end program
