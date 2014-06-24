program diagnose
use elliptic_tools
use field_tools
use constants
implicit none

integer        :: nr, nz
character(256) :: A_file, B_file, C_file, Q_file, input_folder, output_folder, &
&                 output_file, yes_or_no, chi_bc_file, psi_bc_file

real(4), pointer   :: psi(:,:), chi(:,:), f(:,:), Q(:,:), coe(:, :, :),    &
&                     a(:,:), b(:,:), c(:,:), workspace(:,:),              &
&                     solver_a(:,:), solver_b(:,:), solver_c(:,:), JJ(:,:),&
&                     w_mom(:,:), u_mom(:,:), theta(:,:), eta(:,:),        &
&                     ra(:), za(:), exner(:), sigma(:), dthetadr(:,:),     &
&                     psi_bc(:,:), chi_bc(:,:)

real(4)            :: testing_dt, Lr, Lz, dr, dz, eta_avg_b, eta_avg_nob

integer :: saved_strategy, strategy, max_iter
real(4) :: saved_strategy_r, strategy_r

integer :: i,j, m, n, err
real(4) :: r, z, eta_avg, tmp1, tmp2, tmp3
logical :: file_exists, use_chi_bc, use_psi_bc

read(*,*) testing_dt
read(*,*) Lr, Lz;    read(*,*) nr, nz;
read(*,'(a)') input_folder
read(*,'(a)') output_folder
read(*,'(a)') A_file;
read(*,'(a)') B_file;
read(*,'(a)') C_file;
read(*,'(a)') Q_file;
read(*,*) saved_strategy, saved_strategy_r, max_iter;
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



print *, "Testing time: ", testing_dt
print *, "Lr:", Lr, ", Lz:", Lz
print *, "nr:", nr, ", nz:", nz
print *, "Input folder: ", trim(input_folder)
print *, "Output folder: ", trim(output_folder)
print *, "A file: ", trim(A_file)
print *, "B file: ", trim(B_file)
print *, "C file: ", trim(C_file)
print *, "Q file: ", trim(Q_file)
print *, "Strategy wanted: ", saved_strategy
print *, "Strategy specific value: ", saved_strategy_r
print *, "Strategy max iteration: ", max_iter
!print *, "Use PSI boundary condition: ", merge('Yes'//'('//trim(psi_bc_file)//')', 'No', use_psi_bc .eqv. .true.)
!print *, "Use CHI boundary condition: ", merge('Yes'//'('//trim(chi_bc_file)//')', 'No', use_chi_bc .eqv. .true.)


allocate(psi(nr,nz));   allocate(chi(nr,nz));      allocate(eta(nr,nz));
allocate(f(nr,nz));     allocate(Q(nr,nz));        allocate(JJ(nr,nz));
allocate(workspace(nr,nz));

allocate(coe(9,nr,nz));
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

print *, "sum(A): ", sum(a)
print *, "sum(B): ", sum(b)
print *, "sum(C): ", sum(c)
print *, "sum(Q): ", sum(Q)

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

call write_2Dfield(11, trim(output_folder)//"/j.bin", JJ, nr, nz)
call write_2Dfield(11, trim(output_folder)//"/a.bin", solver_a, nr-1, nz-2)
call write_2Dfield(11, trim(output_folder)//"/b.bin", solver_b, nr-1, nz-1)
call write_2Dfield(11, trim(output_folder)//"/c.bin", solver_c, nr-2, nz-1)


f=0.0
! assign f= g0/theta0 * dJJ/dr
call d_dr(JJ, f)
f = f * (g0 * theta0)

call write_2Dfield(11, trim(output_folder)//"/djdr.bin", f, nr, nz)

call cal_coe(solver_a, solver_b, solver_c, coe, dr, dz, nr, nz, err)

! ===== [STAGE   I] ===== !
psi = 0.0
if(use_psi_bc .eqv. .true.) then
    psi = psi_bc
end if

print *, "Solving psi..."
strategy = saved_strategy; strategy_r = saved_strategy_r;
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
!workspace(1:nr-1,1:nz-1) = solver_b(1:nr-1,1:nz-1); solver_b = 0.0
call cal_coe(solver_a, solver_b, solver_c, coe, dr, dz, nr, nz, err)
if(hasNan3(coe)) then
    print *, "coe has Nan"
end if
!solver_b(1:nr-1,1:nz-1) = workspace(1:nr-1,1:nz-1) ! restore

chi = 0.0
if(use_chi_bc .eqv. .true.) then
    chi = chi_bc
end if


strategy = saved_strategy; strategy_r = saved_strategy_r;
call solve_elliptic(max_iter, strategy, strategy_r, chi, coe, f, workspace, nr, nz, err, 1)
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(chi, eta);
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta-nob.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/CHI-nob.bin",chi,nr,nz)

eta_avg_nob = cal_eta_avg(Q, eta)

print *, "Solving CHI with B!=0"
if(hasNan2(chi)) then
    print *, "chi has Nan"
    stop
end if
call cal_coe(solver_a, solver_b, solver_c, coe, dr, dz, nr, nz, err)
strategy = saved_strategy; strategy_r = saved_strategy_r;
call solve_elliptic(max_iter, strategy, strategy_r, chi, coe, f, workspace, nr, nz, err, 1)
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(chi, eta);
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta-b.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/CHI-b.bin",chi,nr,nz)


eta_avg_b = cal_eta_avg(Q, eta)

print *, "Average Efficiency without b (%): ", eta_avg_nob
print *, "Average Efficiency with    b (%): ", eta_avg_b

open (unit=15,file=trim(output_folder)//"/efficiency.txt",action="write",status="replace")
write (15,*) "Average Efficiency witout b (%): ", eta_avg_nob
write (15,*) "Average Efficiency with   b (%): ", eta_avg_b
close (15)


contains

logical function hasNan2(mtx)
implicit none
real(4) :: mtx(:,:)

integer :: n1,n2
integer :: i,j

n1 = size(mtx,1)
n2 = size(mtx,2)
hasNan2 = .false. 
do i=1,n1
    do j=1,n2
        if(isnan(mtx(i,j))) then
            hasNan2 = .true.
            exit
        end if
    end do
end do

end function


logical function hasNan3(mtx)
implicit none
real(4) :: mtx(:,:,:)

integer :: n1,n2,n3
integer :: i,j,k

n1 = size(mtx,1)
n2 = size(mtx,2)
n3 = size(mtx,3)
hasNan3 = .false. 
do i=1,n1
    do j=1,n2
        do k=1,n3
            if(isnan(mtx(i,j,k))) then
                hasNan3 = .true.
                exit
            end if
        end do
    end do
end do

end function

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

do i = 1, nr
    if(ra(i) == 0) then
        do j = 1, nz
            ! using Taylor expansion in 1st order to interpolate value at 
            ! ra(i) from ra(i+1) and ra(i+2) if ra(i)=0
            to_dat(i,j) = 2.0 * to_dat(i+1,j) - to_dat(i+2,j)
        end do
    end if
end do


end subroutine


real(4) function cal_eta_avg(Q, eta)
implicit none
real(4) :: Q(nr, nz), eta(nr, nz)

real(4) :: sum_q, sum_eta_q, r, weight
integer :: i,j

sum_q = 0; sum_eta_q = 0;
do i = 1,nr
    do j = 1, nz
        weight = MERGE(0.5,1.0,j==1 .or. j==nz)
        sum_eta_q = sum_eta_q + eta(i,j) * Q(i,j) * sigma(j) * ra(i) * weight
        sum_q = sum_q + Q(i,j) * sigma(j) * ra(i) * weight
    end do
end do

print *, "sum_eta_q:", sum_eta_q, ", sum_q:", sum_q
print *, "efficiency(%): ", sum_eta_q / sum_q

cal_eta_avg = sum_eta_q / sum_q
end function

subroutine cal_eta(chi, eta)
implicit none
real(4) :: chi(nr,nz), eta(nr,nz)
integer :: i,j

call d_rdr(chi, eta)

do i=1,nr
    do j=1,nz
        eta(i,j) = eta(i,j) * g0 / (sigma(j) * Cp * exner(j) * theta0)
    end do
end do

end subroutine

subroutine A2solverA
implicit none
integer :: i,j

do i=1,nr-1
    do j=1,nz-2
        if(ra(i) == 0) then
            ! using Taylor expansion in 1st order to interpolate value at 
            ! ra(i+0.5) from ra(i+1) and ra(i+2) if ra(i)=0
            solver_a(i,j) =   3.0/2.0 * a(i+1,j+1)/ra(i+1)/sigma(j+1)    &
                &           - 1.0/2.0 * a(i+2,j+1)/ra(i+2)/sigma(j+1)
        else
            solver_a(i,j) = (a(i  ,j+1)/ra(i  )/sigma(j+1)   &
                &          + a(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
        end if
    end do
end do

end subroutine


subroutine B2solverB
implicit none
integer :: i,j

do i=1,nr-1
    do j=1,nz-1
        if(ra(i) == 0) then
            ! using Taylor expansion in 1st order to interpolate value at 
            ! ra(i+0.5) from ra(i+1) and ra(i+2) if ra(i)=0
            solver_b(i,j) = (  3.0/2.0 * b(i+1,j  )/ra(i+1)/sigma(j  )    &
                &            - 1.0/2.0 * b(i+2,j  )/ra(i+2)/sigma(j  )    &
                &            + 3.0/2.0 * b(i+1,j+1)/ra(i+1)/sigma(j+1)    &
                &            - 1.0/2.0 * b(i+2,j+1)/ra(i+2)/sigma(j+1) )/2.0
 
        else
            solver_b(i,j) = ( b(i  ,j  )/ra(i  )/sigma(j  )  &
               &            + b(i+1,j  )/ra(i+1)/sigma(j  )  &
               &            + b(i  ,j+1)/ra(i  )/sigma(j+1)  &
               &            + b(i+1,j+1)/ra(i+1)/sigma(j+1)) / 4.0
        end if
    end do
end do

end subroutine


subroutine C2solverC
implicit none
integer :: i,j

do i=1,nr-2
    do j=1,nz-1
        r = i * dr ! r
        solver_c(i,j) = (c(i+1,j  )/ra(i+1)/sigma(j  )  &
           &           + c(i+1,j+1)/ra(i+1)/sigma(j+1)) / 2.0
    end do
end do
end subroutine

end program
