program diagnose
use elliptic_tools
use field_tools
use constants
use read_input_tools
implicit none

integer, parameter :: stdin=5, fd=15
integer        :: nr, nz
character(256) :: A_file, B_file, C_file, Q_file, F_file, input_folder, output_folder, &
&                 output_file, yes_or_no, rchi_bc_file, rpsi_bc_file, mode_str,          &
&                 word(4), buffer



! Grid point are designed as follows: 
! O A O
! C B C
! O A O
! 
! O  : (nr   , nz  ) rpsi rchi a_in b_in c_in Q_in f dthetadr
! A  : (nr-1 , nz  ) w eta
! B  : (nr-1 , nz-1) m theta F Q solver_B
! C  : (nr   , nz-1) u
! sA : (nr-1 , nz-2) solver_A
! sC : (nr-2 , nz-1) solver_C
!


real(4), pointer   :: rpsi(:,:), rchi(:,:), f(:,:), Q_in(:,:), coe(:, :, :), &
&                     F_in(:,:),                                             &
&                     a_in(:,:), b_in(:,:), c_in(:,:),                       &
&                     wksp_O(:,:), wksp_A(:,:), wksp_B(:,:), wksp_C(:,:),    &
&                     solver_a_A(:,:), solver_b_B(:,:), solver_c_C(:,:),     &
&                     solver_b_basic_B(:,:), solver_B_anomaly_B(:,:),        &
&                     JJ_B(:,:), RHS_rpsi_thm(:,:), RHS_rpsi_mom(:,:),       &
&                     a_A(:,:), b_B(:,:), b_C(:,:), c_C(:,:),                &
&                     w_A(:,:), u_C(:,:), theta(:,:), eta(:,:), m2(:,:),     &
&                     ra(:), za(:), exner(:), sigma(:),                      & 
&                     rpsi_bc(:,:), rchi_bc(:,:), wtheta(:,:), f_basic(:,:), &
&                     f_anomaly(:,:), b_basic_B(:,:), b_anomaly_B(:,:)

real(4)            :: testing_dt, Lr, Lz, dr, dz, eta_avg_b, eta_avg_nob,  &
                      time_beg, time_end

integer :: saved_strategy_rpsi, max_iter_rpsi, &
&          saved_strategy_rchi, max_iter_rchi, strategy
real(4) :: saved_strategy_rpsi_r, saved_strategy_rchi_r, strategy_r, alpha,  &
&          alpha_rpsi, alpha_rchi

integer :: i,j, m, n, err, mode(4)
real(4) :: r, z, eta_avg, tmp1, tmp2, tmp3, eta_avg_b_wtheta
logical :: file_exists, use_rchi_bc, use_rpsi_bc, debug_mode



inquire(file="./debug_mode", exist=debug_mode);

call cpu_time(time_beg)
call read_input(stdin, mode_str)
i = 1; n = 1
do n=1, size(mode)
    j = INDEX(mode_str(i:), "-")
    if (j == 0) then
        word(n) = mode_str(i:)
        exit
    end if
    word(n) = mode_str(i:i+j-2)
    i = i + j
end do

mode = 0
if(word(1) == 'TENDENCY') then
    mode(1) = 0
else if(word(1) == 'INSTANT') then
    mode(1) = 1
else
    print *, "Unknown Mode : [", trim(word(1)), "]", word(1)=='TENDENCY'
    stop
end if

if(word(2) == 'DENSITY_NORMAL') then
    mode(2) = 0
else if(word(2) == 'DENSITY_BOUSSINESQ') then
    mode(2) = 1
else
    print *, "Unknown Mode :", trim(word(2))
    stop
end if

! ABC : [ABC] flag
if(word(3) == 'ABC_UPDATE') then
    mode(3) = 1
else if(word(3) == 'ABC_NOUPDATE') then
    mode(3) = 0
else
    print *, "Unknown Mode :", trim(word(3))
    stop
end if

if(word(4) == 'INTEGRAL_CHECK') then
    mode(4) = 1
else if(word(4) == 'INTEGRAL_NOCHECK') then
    mode(4) = 0
else
    print *, "Unknown Mode :", trim(word(4))
    stop
end if
if(mode(1) == 0) then
    call read_input(stdin, buffer); read(buffer, *) testing_dt
end if
call read_input(stdin, buffer); read(buffer, *) Lr, Lz;
call read_input(stdin, buffer); read(buffer, *) nr, nz;
call read_input(stdin, input_folder);
call read_input(stdin, output_folder);
call read_input(stdin, A_file);
call read_input(stdin, B_file);
call read_input(stdin, C_file);
call read_input(stdin, Q_file);
call read_input(stdin, F_file);

call read_input(stdin, buffer);
read(buffer, *) saved_strategy_rpsi, saved_strategy_rpsi_r, max_iter_rpsi, alpha_rpsi;

call read_input(stdin, buffer);
read(buffer, *) saved_strategy_rchi, saved_strategy_rchi_r, max_iter_rchi, alpha_rchi;

call read_input(stdin, yes_or_no)
if(yes_or_no == 'yes') then
    call read_input(stdin, rpsi_bc_file);
    use_rpsi_bc = .true.
else
    use_rpsi_bc = .false.
end if
call read_input(stdin, yes_or_no)
if(yes_or_no == 'yes') then
    call read_input(stdin, rchi_bc_file);
    use_rchi_bc = .true.
else
    use_rchi_bc = .false.
end if


print *, "mode: ", mode(1),",",mode(2),",",mode(3)
if(mode(1) == 0) then
    print *, "Testing time: ", testing_dt
end if
print *, "Lr:", Lr, ", Lz:", Lz
print *, "nr:", nr, ", nz:", nz
print *, "Input folder: ", trim(input_folder)
print *, "Output folder: ", trim(output_folder)
print *, "A file: ", trim(A_file)
print *, "B file: ", trim(B_file)
print *, "C file: ", trim(C_file)
print *, "Q file: ", trim(Q_file)
print *, "F file: ", trim(F_file)
print *, "rpsi's strategy, residue, iter: ", saved_strategy_rpsi, &
&           saved_strategy_rpsi_r, max_iter_rpsi, alpha_rpsi
print *, "rchi's strategy, residue, iter: ", saved_strategy_rchi, &
&           saved_strategy_rchi_r, max_iter_rchi, alpha_rchi

if(use_rpsi_bc .eqv. .true.) then
    print *, "Use rpsi boundary condition: Yes (", trim(rpsi_bc_file), ")"
else
    print *, "Use rpsi boundary condition: No "
end if

if(use_rchi_bc .eqv. .true.) then
    print *, "Use CHI boundary condition: Yes (", trim(rchi_bc_file), ")"
else
    print *, "Use rpsi boundary condition: No "
end if


allocate(rpsi(nr,nz));   allocate(rchi(nr,nz));   
allocate(f(nr,nz));     allocate(Q_in(nr-1,nz-1)); allocate(JJ_B(nr-1,nz-1));
allocate(F_in(nr-1,nz-1));
allocate(RHS_rpsi_thm(nr,nz)); allocate(RHS_rpsi_mom(nr,nz));
allocate(wksp_O(nr,nz));
allocate(wksp_A(nr-1,nz));
allocate(wksp_B(nr-1,nz-1));
allocate(wksp_C(nr,nz-1));
allocate(f_basic(nr,nz));
allocate(f_anomaly(nr,nz));

allocate(coe(9,nr,nz));
allocate(a_in(nr, nz));  allocate(b_in(nr, nz));   allocate(c_in(nr, nz));
allocate(a_A(nr-1, nz)); allocate(b_C(nr, nz-1));  allocate(b_B(nr-1,nz-1));
allocate(c_C(nr, nz-1)); 
allocate(b_basic_B(nr-1,nz-1));                    allocate(b_anomaly_B(nr-1,nz-1));

allocate(solver_a_A(nr-1, nz-2)); allocate(solver_b_B(nr-1, nz-1));allocate(solver_c_C(nr-2, nz-1));
allocate(solver_b_basic_B(nr-1, nz-1));
allocate(solver_b_anomaly_B(nr-1, nz-1));
allocate(wtheta(nr-1,nz-1));
allocate(w_A(nr-1,nz)); allocate(u_C(nr,nz-1));
allocate(theta(nr-1,nz-1)); allocate(eta(nr-1,nz)); allocate(m2(nr-1,nz-1));
allocate(ra(nr));       allocate(za(nz));          allocate(sigma(nz));
allocate(exner(nz));    

call read_2Dfield(15, trim(input_folder)//"/"//A_file, a_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//B_file, b_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//C_file, c_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//Q_file, Q_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//F_file, F_in, nr, nz)

if(use_rpsi_bc .eqv. .true.) then
    allocate(rpsi_bc(nr,nz));
    call read_2Dfield(15, trim(input_folder)//"/"//rpsi_bc_file, rpsi_bc, nr, nz)
end if


if(use_rchi_bc .eqv. .true.) then
    allocate(rchi_bc(nr,nz));
    call read_2Dfield(15, trim(input_folder)//"/"//rchi_bc_file, rchi_bc, nr, nz)
end if


! ### Calculate dr, dz, ra, za, exner, sigma (pseudo-density)
dr = Lr / (nr-1); dz = Lz / (nz-1);
do i=1,nr
    ra(i) = (i-1) * dr
end do
do j=1,nz
    za(j) = (j-1) * dz
    exner(j) = merge(1.0 - za(j) / h0, 1.0, mode(2) == 0)
    sigma(j) = merge(p0 / (theta0 * Rd) * exner(j)**(1.0 / kappa - 1.0), 1.0, &
&                   mode(2) == 0)
end do



! ### Normalize coefficient solver_a_A = a / (sigma * r), solver_b_B = b / (sigma * r)
!                         , solver_c_C = c / (sigma * r)

do i=1,nr-1
    do j=1,nz-2
        solver_a_A(i,j) = (a_in(i,j+1) + a_in(i+1,j+1))         &
            &           / (ra(i) + ra(i+1)) / sigma(j+1)
    end do
end do


do i=1,nr-1
    do j=1,nz-1
        solver_b_B(i,j) = ( b_in(i  ,j  ) + b_in(i+1,j  )       &
            &             + b_in(i  ,j+1) + b_in(i+1,j+1) )     &
            &             /(ra(i)+ra(i+1))/(sigma(j)+sigma(j+1))
    end do
end do

solver_b_basic_B = solver_b_B

do i=1,nr-2
    do j=1,nz-1
        solver_c_C(i,j) = (c_in(i+1,j) + c_in(i+1,j+1))         &
            &           / ra(i+1) / (sigma(j)+sigma(j+1))  
    end do
end do

if(hasNan2(solver_a_A)) then
    print *, "SOLVER a has NAN"
end if
if(hasNan2(solver_b_B)) then
    print *, "SOLVER b has NAN"
end if
if(hasNan2(solver_c_C)) then
    print *, "SOLVER c has NAN"
end if



! ### Calculate a_A, b_C, b_B, c_C

do i=1,nr-1
    do j=1,nz
        a_A(i,j) = (a_in(i,j) + a_in(i+1,j)) / 2.0
    end do
end do

do i=1,nr
    do j=1,nz-1
        b_C(i,j) = (b_in(i,j) + b_in(i,j+1)) / 2.0
    end do
end do

do i=1,nr-1
    do j=1,nz-1
        b_B(i,j) = ( b_in(i  ,j  ) + b_in(i+1,j  )       &
            &      + b_in(i  ,j+1) + b_in(i+1,j+1) ) /4.0
    end do
end do

b_basic_B = b_B ! Copy basic state

do i=1,nr
    do j=1,nz-1
        c_C(i,j) = (c_in(i,j) + c_in(i,j+1)) / 2.0
    end do
end do

! ### Derive angular momentum square from C term (c_C)
! The integration order matters! Please be aware of it.
do i=1, nr-1
    do j=1, nz-1
        if(i==1) then
            m2(i,j) = 0
        else
            m2(i,j) = m2(i-1,j) + (ra(i)**3.0) * c_C(i,j) * (ra(i+1) - ra(i-1)) / 2.0
        end if
    end do
end do


! ### Assign JJ_in JJ_B
do i=1,nr-1
    do j=1,nz-1
        JJ_B(i,j) = Q_in(i,j) / (Cp * exner(j))
    end do
end do

call write_2Dfield(11, trim(output_folder)//"/J-B.bin",     JJ_B, nr-1, nz-1)
call write_2Dfield(11, trim(output_folder)//"/solver_a-sA.bin", solver_a_A, nr-1, nz-2)
call write_2Dfield(11, trim(output_folder)//"/solver_b-B.bin", solver_b_B, nr-1, nz-1)
call write_2Dfield(11, trim(output_folder)//"/solver_c-sC.bin", solver_c_C, nr-2, nz-1)


! ### Assign 1st part of RHS = g0/theta0 * dJJ/dr
RHS_rpsi_thm = 0.0
call d_dr_B2C(JJ_B, wksp_C)
!print *, "djdr max: " , maxval(abs(wksp_C))
! wksp_C to f
do i=2,nr-1
    do j=2,nz-1
        RHS_rpsi_thm(i,j) = (wksp_C(i,j)+wksp_C(i,j-1))/2.0
    end do
end do
!print *, "djdr max: " , maxval(abs(djdr))
RHS_rpsi_thm = RHS_rpsi_thm * g0 / theta0
!print *, "djdr max: " , maxval(abs(djdr))
call write_2Dfield(11, trim(output_folder)//"/RHS_rpsi_thm-O.bin", RHS_rpsi_thm, nr, nz)


! ### Assign 2nt part of RHS = - 2 * r^(-2) * d(mf)/dz
RHS_rpsi_mom = 0.0
do i=1,nr-1
    do j=1,nz-1
         wksp_B(i,j) = sqrt(m2(i,j)) * F_in(i,j)
    end do
end do
call d_dz_B2A(wksp_B, wksp_A)
! wksp_A to f
do i=2,nr-1
    do j=2,nz-1
         ! Notice that 2.0 factor of averge wksp_A cancelled by another in 2*m*F
         RHS_rpsi_mom(i,j) = - (wksp_A(i,j) + wksp_A(i-1,j))/(ra(i)**2.0)
    end do
end do

    
!print *, "RHS_rpsi_mom max: " , maxval(abs(djdr))
call write_2Dfield(11, trim(output_folder)//"/RHS_rpsi_mom-O.bin", RHS_rpsi_mom, nr, nz)



print *, "Initialization complete."

! TENDENCY MODE
if(mode(1) == 0) then
    ! ### STAGE I : invert to get rpsi  !
    call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)
    rpsi = 0.0
    if(use_rpsi_bc .eqv. .true.) then
        rpsi = rpsi_bc
    end if

    print *, "Solving rpsi..."
    f = RHS_rpsi_thm + RHS_rpsi_mom
    strategy = saved_strategy_rpsi; strategy_r = saved_strategy_rpsi_r
    alpha = alpha_rpsi
    call solve_elliptic(max_iter_rpsi, strategy, strategy_r, alpha, rpsi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

    call write_2Dfield(11, trim(output_folder)//"/rpsi_before-O.bin", rpsi, nr, nz)

    ! ### STAGE II : forcast thermal field to generate f  !

    ! calculate u, w and dtheta_dt
    call rpsiToUW(rpsi, u_C, w_A)

    theta = 0.0
    do i=1,nr-1
        do j=1,nz-1
            ! dtheta/dt = J - w dtheta/dz - u dtheta/dr
            theta(i,j) = JJ_B(i,j) - theta0/g0 * ( a_A(i,j  ) * w_A(i,j  )         &
            &                                + a_A(i,j+1) * w_A(i,j+1) ) / 2.0 &
            &                  + theta0/g0 * ( b_C(i  ,j) * u_C(i  ,j)         &
            &                                + b_C(i+1,j) * u_C(i+1,j) ) / 2.0 
        end do
    end do

    print *, "Max dtheta_dt: ", maxval(theta)

    call write_2Dfield(11, trim(output_folder)//"/w_before-A.bin",w_A,nr-1,nz)
    call write_2Dfield(11, trim(output_folder)//"/u_before-C.bin",u_C,nr,nz-1)
    call write_2Dfield(11, trim(output_folder)//"/dtheta_dt-B.bin", theta,nr-1,nz-1)
    
    theta = theta * testing_dt
    ! Calculate forced B field, A and C assumed to be unperturbed for now
    ! dtheta/dr for f
    call d_dr_B2B(theta, wksp_B)
    b_anomaly_B = -g0/theta0 * wksp_B ! Copy anomaly
    b_B = b_B + b_anomaly_B
    
    ! Notice that solver only needs 2:nz-1,
    ! so the value of a_A on top and boundary is not important.
    ! recalculate a_A is for the checking step
    call d_dz_B2A(theta, wksp_A)
    a_A(:,2:nz-1) = a_A(:,2:nz-1) + g0 / theta0 * wksp_A(:,2:nz-1)

    do i=1,nr-1
        do j=1,nz-1
            solver_b_anomaly_B(i,j) = b_anomaly_B(i,j)          &
            &  / ((ra(i)+ra(i+1))/2.0)/((sigma(j)+sigma(j+1))/2.0)
        end do
    end do

end if
! ### STAGE III : Solve !

! doing right hand side
f_basic = 0; f_anomaly = 0;
do i=2,nr-1
    do j=2,nz-1
        f_basic(i,j) = - ( b_basic_B(i-1,j-1) &
            &      + b_basic_B(i-1,j  ) &
            &      + b_basic_B(i  ,j  ) &
            &      + b_basic_B(i  ,j-1) )/4.0

        f_anomaly(i,j) = - ( b_anomaly_B(i-1,j-1) &
            &      + b_anomaly_B(i-1,j  ) &
            &      + b_anomaly_B(i  ,j  ) &
            &      + b_anomaly_B(i  ,j-1) )/4.0
 
    end do
end do

call write_2Dfield(11, trim(output_folder)//"/RHS_rchi-O.bin", f, nr, nz)


print *, "Solving CHI with L(A,B=B0   ,C) = -dB"
print *, "Solving CHI with L(A,B=B0   ,C) = -B0"
print *, "Solving CHI with L(A,B=B0+dB,C) = -dB"
print *, "Solving CHI with L(A,B=B0+dB,C) = -B0"



rchi = 0.0
if(use_rchi_bc .eqv. .true.) then
    rchi = rchi_bc
end if

!==================================================================================
print *, "Solving CHI with L(A,B=0,C) = -B0"
solver_b_B = 0.0; f = f_basic;
call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta);
eta_avg_nob = cal_eta_avg(Q_in, eta)
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta_[0_B0]-O.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi_[0_B0]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=0,C) = -dB"
solver_b_B = 0.0; f = f_anomaly;
call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta);
eta_avg_nob = cal_eta_avg(Q_in, eta)
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta_[0_dB]-O.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi_[0_dB]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=B0,C) = -B0"
solver_b_B = solver_b_basic_B; f = f_basic;
call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta);
eta_avg_nob = cal_eta_avg(Q_in, eta)
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta_[B0_B0]-O.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi_[B0_B0]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=B0,C) = -dB"
solver_b_B = solver_b_basic_B; f = f_anomaly;
call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta);
eta_avg_nob = cal_eta_avg(Q_in, eta)
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta_[B0_dB]-O.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi_[B0_dB]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=B0+dB,C) = -B0"
solver_b_B = solver_b_basic_B + solver_b_anomaly_B; f = f_basic;
call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta);
eta_avg_nob = cal_eta_avg(Q_in, eta)
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta_[B0dB_dB]-O.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi_[B0dB_dB]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=B0+dB,C) = -dB"
solver_b_B = solver_b_basic_B + solver_b_anomaly_B; f = f_anomaly;
call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta);
eta_avg_nob = cal_eta_avg(Q_in, eta)
eta = eta * 100.0 ! in percent

call write_2Dfield(11,trim(output_folder)//"/eta_[B0dB_dB]-O.bin",eta,nr,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi_[B0dB_dB]-O.bin",rchi,nr,nz)

!==================================================================================

if(mode(4) == 1) then
    print *, "Integral check..."
    if(use_rpsi_bc .eqv. .true.) then
        rpsi = rpsi_bc
    end if

    print *, "Solving rpsi..."
    f = RHS_rpsi_thm + RHS_rpsi_mom;
    call cal_coe(solver_a_A, solver_b_B, solver_c_C, coe, dr, dz, nr, nz, err)
    strategy = saved_strategy_rpsi; strategy_r = saved_strategy_rpsi_r;
    alpha = alpha_rpsi;
    call solve_elliptic(max_iter_rpsi, strategy, strategy_r, alpha, rpsi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."
    call write_2Dfield(11, trim(output_folder)//"/rpsi-O.bin", rpsi, nr, nz)

    call rpsiToUW(rpsi, u_C, w_A)

    ! Update b_C for checking mode
    do i=2, nr-1
        do j=1, nz-1
            b_C(i,j) = ( b_B(i-1, j  ) &
                &      + b_B(i  , j  ) )/2.0
        end do
    end do
    call relativeTheta(theta, a_A * (theta0 / g0), b_C * (- theta0 / g0))
    eta_avg_b_wtheta = cal_eta_avg_wtheta(Q_in, w_A, theta)
    print *, eta_avg_b_wtheta
    eta_avg_b_wtheta = cal_eta_avg_wtheta2(Q_in, rpsi, theta)
    print *, eta_avg_b_wtheta
    call write_2Dfield(11,trim(output_folder)//"/wtheta-B.bin",wtheta,nr-1,nz-1)

    call write_2Dfield(11,trim(output_folder)//"/w-A.bin",w_A,nr-1,nz)
    call write_2Dfield(11,trim(output_folder)//"/u-C.bin",u_C,nr,nz-1)
    call write_2Dfield(11,trim(output_folder)//"/theta-B.bin",theta,nr-1,nz-1)

end if

print *, "Average Efficiency without b using Q eta   (%): ", eta_avg_nob
print *, "Average Efficiency with    b using Q eta   (%): ", eta_avg_b
if(mode(3)==0) then
    print *, "Average Efficiency with    b using w theta (%): ", eta_avg_b_wtheta
end if

open (unit=15,file=trim(output_folder)//"/efficiency.txt",action="write",status="replace")
write (15,*) "Average Efficiency without b using Q eta   (%): ", eta_avg_nob
write (15,*) "Average Efficiency with    b using Q eta   (%): ", eta_avg_b

! ABC_UPDATE
if(mode(3)==0) then
    write (15,*) "Average Efficiency with    b using w theta (%): ", eta_avg_b_wtheta
end if
close (15)

call cpu_time(time_end)

print *, "Time elapsed (sec): ", (time_end - time_beg)

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

subroutine relativeTheta(theta_B, dtheta_dz_A, dtheta_dr_C)
implicit none
real(4) :: theta_B(nr-1,nz-1), dtheta_dz_A(nr-1,nz), dtheta_dr_C(nr,nz-1)

real(4) :: dist
integer :: i,j
theta_B = theta0;
do i=2,nr-1
    dist = (ra(i+1) - ra(i-1)) / 2.0
    theta_B(i,1) = theta_B(i-1,1) + dist * dtheta_dr_C(i,1)
end do

do i=1,nr-1
    do j=2,nz-1
        dist = (za(j+1) - za(j-1)) / 2.0
        theta_B(i,j) = theta_B(i,j-1) + dist * dtheta_dz_A(i,j)
    end do
end do

end subroutine


subroutine rpsiToUW(from_rpsi, to_u, to_w)
implicit none
real(4) :: from_rpsi(nr,nz), to_w(nr-1,nz), to_u(nr,nz-1)
real(4) :: r
integer :: i,j

call d_rdr_O2A(from_rpsi, to_w)
call d_dz_O2C(from_rpsi, to_u); to_u = -to_u
! notice gradient of rpsi gets momentum flux, so we need to divide it by sigma.
do i=1,nr-1
    do j=1,nz
        to_w(i,j) = to_w(i,j)/sigma(j)
    end do
end do

do i=1,nr
    do j=1,nz-1
        r = ra(i)
        if(r /= 0) then 
            to_u(i,j) = to_u(i,j)/(r*(sigma(j)+sigma(j+1))/2.0)
        else
            to_u(i,j) = 0.0
        end if
    end do
end do

end subroutine

subroutine d_dz_B2A(from_dat, to_dat)
! Notice that the top and bottom will not
! be altered since grid points' arrangement is not matched
implicit none
real(4) :: from_dat(nr-1,nz-1), to_dat(nr-1,nz)
integer :: i,j,m,n
do i = 1, nr-1
    do j = 2, nz-2
        to_dat(i,j) = (from_dat(i,j) - from_dat(i,j-1)) &
            &       / ((za(j+1) - za(j-1))/2.0)
    end do
end do
end subroutine


subroutine d_dz_O2C(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr,nz-1)
integer :: i,j,m,n
do i = 1, nr
    do j = 1, nz-1
        to_dat(i,j) = (from_dat(i,j+1) - from_dat(i,j)) / (za(j+1) - za(j))
    end do
end do
end subroutine

subroutine d_dr_B2B(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr-1,nz-1), to_dat(nr-1,nz-1)
integer :: i,j,m,n
do i = 1, nr-1
    do j = 1, nz-1
        if(i == 1) then
            m = 0; n = 1
        else if(i == nr-1) then
            m = -1; n = 0
        else
            m = -1; n = 1
        end if
        to_dat(i,j) = (from_dat(i+m,j) - from_dat(i+n,j)) / (ra(i+m) - ra(i+n)) 
    end do
end do
end subroutine

subroutine d_dr_B2C(from_dat, to_dat)
! Notice that the first and last columns will not
! be altered since grid points' arrangement is not matched
implicit none
real(4) :: from_dat(nr-1,nz-1), to_dat(nr,nz-1)
integer :: i,j,m,n
do i = 2, nr-1
    do j = 1, nz-1
        to_dat(i,j) = (from_dat(i,j) - from_dat(i-1,j)) &
            &       / ((ra(i+1) - ra(i-1))/2.0)
    end do
end do
end subroutine

subroutine d_dr_O2A(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr-1,nz)
integer :: i,j
do i = 1, nr-1
    do j = 1, nz
        to_dat(i,j) = (from_dat(i+1,j) - from_dat(i,j)) / (ra(i+1) - ra(i)) 
    end do
end do
end subroutine


subroutine d_rdr_O2A(from_dat, to_dat)
! notice that all the points with r=0 become inifinite
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr-1,nz)
real(4) :: r
integer :: i,j

call d_dr_O2A(from_dat, to_dat)

do i = 1, nr-1
    do j = 1, nz
        r = (ra(i)+ra(i+1))/2.0
        to_dat(i,j) = to_dat(i,j) / r
    end do
end do

end subroutine

real(4) function cal_eta_avg(Q, eta)
implicit none
real(4) :: Q(nr-1, nz-1), eta(nr-1, nz)
real(4) :: sum_q, sum_eta_q, r, dr, dz
integer :: i,j

sum_q = 0; sum_eta_q = 0;
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)

        sum_eta_q = sum_eta_q + ((eta(i,j) + eta(i,j+1)) /2.0) &
           &                    * Q(i,j) * sigma(j) * r * dr * dz
        sum_q = sum_q + Q(i,j) * sigma(j) * r * dr * dz
    end do
end do

print *, "sum_eta_q:", sum_eta_q, ", sum_q:", sum_q
print *, "efficiency(%): ", sum_eta_q / sum_q * 100

cal_eta_avg = sum_eta_q / sum_q * 100
end function


real(4) function cal_eta_avg_wtheta(Q_B, w_A, theta_B)
implicit none
real(4) :: Q_B(nr-1, nz-1), w_A(nr-1, nz), theta_B(nr-1, nz-1)
real(4) :: sum_q, sum_wtheta, r, tmp, sum_pos, sum_neg, dr, dz
integer :: i,j

sum_q = 0; sum_wtheta = 0;
sum_pos = 0; sum_neg = 0;
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)

        wtheta(i,j) =  ((w_A(i,j) * sigma(j) + w_A(i,j+1) * sigma(j+1)) /2.0) &
           &          * theta_B(i,j)

        tmp = g0 / theta0 &
           &          * ((w_A(i,j) * sigma(j) + w_A(i,j+1) * sigma(j+1)) /2.0) &
           &          * theta_B(i,j) * r * dr * dz

        sum_wtheta = sum_wtheta + g0 / theta0 &
           &          * ((w_A(i,j) * sigma(j) + w_A(i,j+1) * sigma(j+1)) /2.0) &
           &          * theta_B(i,j) * r * dr * dz
        if(tmp > 0) then
            sum_pos = sum_pos + tmp
        else
            sum_neg = sum_neg + tmp
        end if
        sum_q = sum_q + Q_B(i,j) * sigma(j) * r * dr * dz 
    end do
end do

print *, "sum_pos:", sum_pos, ", sum_neg:", sum_neg
print *, "sum_wtheta:", sum_wtheta, ", sum_q:", sum_q
print *, "efficiency(%): ", sum_wtheta / sum_q * 100.0

cal_eta_avg_wtheta = sum_wtheta / sum_q * 100.0
end function

real(4) function cal_eta_avg_wtheta2(Q_B, rpsi_O, theta_B)
implicit none
real(4) :: Q_B(nr-1, nz-1), rpsi_O(nr, nz), theta_B(nr-1, nz-1)
real(4) :: sum_q, sum_wtheta, r, tmp, sum_pos, sum_neg, dr, dz
integer :: i,j

sum_q = 0; sum_wtheta = 0;
sum_pos = 0; sum_neg = 0;
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)

        tmp = g0 / theta0 &
           &    * ((rpsi_O(i+1,j) + rpsi_O(i+1,j+1) - rpsi_O(i,j) - rpsi_O(i,j+1)) /2.0) &
           &    * theta_B(i,j) * dz

        sum_wtheta = sum_wtheta + tmp

        if(tmp > 0) then
            sum_pos = sum_pos + tmp
        else
            sum_neg = sum_neg + tmp
        end if
        sum_q = sum_q + Q_B(i,j) * sigma(j) * r * dr * dz 
    end do
end do

print *, "sum_pos:", sum_pos, ", sum_neg:", sum_neg
print *, "sum_wtheta:", sum_wtheta, ", sum_q:", sum_q
print *, "efficiency(%): ", sum_wtheta / sum_q * 100.0

cal_eta_avg_wtheta2 = sum_wtheta / sum_q * 100.0
end function


subroutine cal_eta(rchi, eta)
implicit none
real(4) :: rchi(nr,nz), eta(nr-1,nz)
integer :: i,j

call d_rdr_O2A(rchi, eta)

do i=1,nr-1
    do j=1,nz
        eta(i,j) = eta(i,j) * g0 / (sigma(j) * Cp * exner(j) * theta0)
    end do
end do

end subroutine


end program
