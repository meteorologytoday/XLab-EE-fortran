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
&                 word(3), buffer



! Grid point are designed as follows: 
! O A O
! C B C
! O A O
! 
! O  : (nr   , nz  ) rpsi rchi rhoA_in rhoB_in rhoC_in Q_in f dthetadr
! A  : (nr-1 , nz  ) w eta
! B  : (nr-1 , nz-1) m theta F Q solver_B
! C  : (nr   , nz-1) u
! sA : (nr-1 , nz-2) solver_A
! sC : (nr-2 , nz-1) solver_C
!


real(4), pointer   :: rpsi(:,:), rchi(:,:), f(:,:), Q_in(:,:), coe(:, :, :), &
&                     F_in(:,:),                                             &
&                     rhoA_in(:,:), rhoB_in(:,:), rhoC_in(:,:),                       &
&                     wksp_O(:,:), wksp_A(:,:), wksp_B(:,:), wksp_C(:,:),    &
&                     solverA_A(:,:), solverB_B(:,:), solverC_C(:,:),     &
&                     solver_b_basic_B(:,:), solver_B_anomaly_B(:,:),        &
&                     JJ_B(:,:), RHS_rpsi_thm(:,:), RHS_rpsi_mom(:,:),       &
&                     rhoA_A(:,:), rhoB_B(:,:), rhoB_C(:,:), rhoC_C(:,:),                &
&                     w_A(:,:), u_C(:,:), theta(:,:), eta(:,:), m2(:,:),     &
&                     ra(:), za(:), exner(:), rho(:),                      & 
&                     rpsi_bc(:,:), rchi_bc(:,:), wtheta_B(:,:), f_basic(:,:), &
&                     f_anomaly(:,:), b_basic_B(:,:), b_anomaly_B(:,:),      &
&                     bndconv(:,:)

real(4)            :: testing_dt, Lr, Lz, dr, dz, eta_avg_b, eta_avg_nob,  &
                      time_beg, time_end

integer :: saved_strategy_rpsi, max_iter_rpsi, &
&          saved_strategy_rchi, max_iter_rchi, strategy
real(4) :: saved_strategy_rpsi_r, saved_strategy_rchi_r, strategy_r, alpha,  &
&          alpha_rpsi, alpha_rchi

integer :: i,j, m, n, err, mode(4)
real(4) :: r, z, tmp1, tmp2, tmp3, &
&          sum_Q, &
&          sum_Qeta_0_0  , sum_Qeta_dB_0  , sum_Qeta_B0_0  , sum_Qeta_B0dB_0 , &
&          sum_Qeta_0_dB , sum_Qeta_dB_dB , sum_Qeta_B0_dB , sum_Qeta_B0dB_dB, &
&          sum_Qeta_0_B0 , sum_Qeta_dB_B0 , sum_Qeta_B0_B0 , sum_Qeta_B0dB_B0, &
&          sum_wtheta_0_JF, sum_wtheta_dB_JF, sum_wtheta_B0_JF, sum_wtheta_B0dB_JF, &
&          sum_bndconv_0, sum_bndconv_dB, sum_bndconv_B0, sum_bndconv_B0dB, &
&          sum_bndconv2_0,sum_bndconv2_dB,sum_bndconv2_B0,sum_bndconv2_B0dB

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


print *, "mode: ", mode(1),",",mode(2),",",mode(3),",",mode(4)
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
allocate(rhoA_in(nr, nz));  allocate(rhoB_in(nr, nz));   allocate(rhoC_in(nr, nz));
allocate(rhoA_A(nr-1, nz)); allocate(rhoB_C(nr, nz-1));  allocate(rhoB_B(nr-1,nz-1));
allocate(rhoC_C(nr, nz-1)); 
allocate(b_basic_B(nr-1,nz-1));                    allocate(b_anomaly_B(nr-1,nz-1));

allocate(solverA_A(nr-1, nz-2)); allocate(solverB_B(nr-1, nz-1));allocate(solverC_C(nr-2, nz-1));
allocate(solver_b_basic_B(nr-1, nz-1));
allocate(solver_b_anomaly_B(nr-1, nz-1));
allocate(wtheta_B(nr-1,nz-1));
allocate(w_A(nr-1,nz)); allocate(u_C(nr,nz-1));
allocate(theta(nr-1,nz-1)); allocate(eta(nr-1,nz)); allocate(m2(nr-1,nz-1));
allocate(ra(nr));       allocate(za(nz));          allocate(rho(nz));
allocate(exner(nz));    allocate(bndconv(nr-1,2));

call read_2Dfield(15, trim(input_folder)//"/"//A_file, rhoA_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//B_file, rhoB_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//C_file, rhoC_in, nr, nz)
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


! ### Calculate dr, dz, ra, za, exner, rho (pseudo-density)
dr = Lr / (nr-1); dz = Lz / (nz-1);
do i=1,nr
    ra(i) = (i-1) * dr
end do
do j=1,nz
    za(j) = (j-1) * dz
    exner(j) = merge(1.0 - za(j) / h0, 1.0, mode(2) == 0)
    rho(j) = merge(p0 / (theta0 * Rd) * exner(j)**(1.0 / kappa - 1.0), 1.0, &
&                   mode(2) == 0)
end do


! ### sum Q ### !
sum_Q = cal_sum_Q(Q_in);



! ### Normalize coefficient solverA_A = a / (rho * r), solverB_B = b / (rho * r)
!                         , solverC_C = c / (rho * r)

do i=1,nr-1
    do j=1,nz-2
        solverA_A(i,j) = (rhoA_in(i,j+1) + rhoA_in(i+1,j+1))         &
            &           / (ra(i) + ra(i+1)) / rho(j+1)
    end do
end do


do i=1,nr-1
    do j=1,nz-1
        solverB_B(i,j) = ( rhoB_in(i  ,j  ) + rhoB_in(i+1,j  )       &
            &             + rhoB_in(i  ,j+1) + rhoB_in(i+1,j+1) )     &
            &             /(ra(i)+ra(i+1))/(rho(j)+rho(j+1))
    end do
end do

solver_b_basic_B = solverB_B

do i=1,nr-2
    do j=1,nz-1
        solverC_C(i,j) = (rhoC_in(i+1,j) + rhoC_in(i+1,j+1))         &
            &           / ra(i+1) / (rho(j)+rho(j+1))  
    end do
end do

if(hasNan2(solverA_A)) then
    print *, "SOLVER a has NAN"
end if
if(hasNan2(solverB_B)) then
    print *, "SOLVER b has NAN"
end if
if(hasNan2(solverC_C)) then
    print *, "SOLVER c has NAN"
end if



! ### Calculate rhoA_A, rhoB_C, rhoB_B, rhoC_C

do i=1,nr-1
    do j=1,nz
        rhoA_A(i,j) = (rhoA_in(i,j) + rhoA_in(i+1,j)) / 2.0
    end do
end do

do i=1,nr
    do j=1,nz-1
        rhoB_C(i,j) = (rhoB_in(i,j) + rhoB_in(i,j+1)) / 2.0
    end do
end do

do i=1,nr-1
    do j=1,nz-1
        rhoB_B(i,j) = ( rhoB_in(i  ,j  ) + rhoB_in(i+1,j  )       &
            &      + rhoB_in(i  ,j+1) + rhoB_in(i+1,j+1) ) /4.0
    end do
end do

b_basic_B = rhoB_B ! Copy basic state

do i=1,nr
    do j=1,nz-1
        rhoC_C(i,j) = (rhoC_in(i,j) + rhoC_in(i,j+1)) / 2.0
    end do
end do

! ### Derive angular momentum square from C term (rhoC_C)
! The integration order matters! Please be aware of it.
do i=1, nr-1
    do j=1, nz-1
        if(i==1) then
            m2(i,j) = 0
        else
            m2(i,j) = m2(i-1,j) + (ra(i)**3.0) * rhoC_C(i,j) * (ra(i+1) - ra(i-1)) / 2.0
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
call write_2Dfield(11, trim(output_folder)//"/solver_a-sA.bin", solverA_A, nr-1, nz-2)
call write_2Dfield(11, trim(output_folder)//"/solver_b-B.bin", solverB_B, nr-1, nz-1)
call write_2Dfield(11, trim(output_folder)//"/solver_c-sC.bin", solverC_C, nr-2, nz-1)


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
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)
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
            theta(i,j) = JJ_B(i,j) - theta0/g0 * ( rhoA_A(i,j  ) * w_A(i,j  )         &
            &                                + rhoA_A(i,j+1) * w_A(i,j+1) ) / 2.0 &
            &                  + theta0/g0 * ( rhoB_C(i  ,j) * u_C(i  ,j)         &
            &                                + rhoB_C(i+1,j) * u_C(i+1,j) ) / 2.0 
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
    rhoB_B = rhoB_B + b_anomaly_B
    
    ! Notice that solver only needs 2:nz-1,
    ! so the value of rhoA_A on top and bottom boundary is not important.
    ! recalculate rhoA_A is for the checking step
    call d_dz_B2A(theta, wksp_A)
    rhoA_A(:,2:nz-1) = rhoA_A(:,2:nz-1) + g0 / theta0 * wksp_A(:,2:nz-1)

    ! Update rhoB_C 
    do i=2, nr-1
        do j=1, nz-1
            rhoB_C(i,j) = ( rhoB_B(i-1, j  ) &
                   &      + rhoB_B(i  , j  ) )/2.0
        end do
    end do
    call relativeTheta(theta, rhoA_A * (theta0 / g0), rhoB_C * (- theta0 / g0))
    call write_2Dfield(11,trim(output_folder)//"/theta_after-B.bin",theta,nr-1,nz-1)


    do i=1,nr-1
        do j=1,nz-1
            solver_b_anomaly_B(i,j) = b_anomaly_B(i,j)          &
            &  / ((ra(i)+ra(i+1))/2.0)/((rho(j)+rho(j+1))/2.0)
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

f = f_basic + f_anomaly;
call write_2Dfield(11, trim(output_folder)//"/RHS_rchi-O.bin", f, nr, nz)



! The order of the following is such that the initial guessing field is better.
!==================================================================================

if(use_rchi_bc .eqv. .true.) then
    rchi = rchi_bc
    print *, "Solving CHI with L(A,B=0,C) = 0 with boundary condition"
    solverB_B = 0.0; f = 0;
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
    alpha = alpha_rchi;
    call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

    call cal_eta(rchi, eta)
    sum_Qeta_0_0 = cal_sum_Qeta(Q_in, eta)

    call write_2Dfield(11,trim(output_folder)//"/eta-[0_0]-A.bin",eta,nr-1,nz)
    call write_2Dfield(11,trim(output_folder)//"/rchi-[0_0]-O.bin",rchi,nr,nz)

    print *, "Solving CHI with L(A,B=dB,C) = 0 with boundary condition"
    solverB_B = solver_b_anomaly_B; f = 0;
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
    alpha = alpha_rchi;
    call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

    call cal_eta(rchi, eta)
    sum_Qeta_dB_0 = cal_sum_Qeta(Q_in, eta)

    call write_2Dfield(11,trim(output_folder)//"/eta-[dB_0]-A.bin",eta,nr-1,nz)
    call write_2Dfield(11,trim(output_folder)//"/rchi-[dB_0]-O.bin",rchi,nr,nz)

    print *, "Solving CHI with L(A,B=B0,C) = 0 with boundary condition"
    solverB_B = solver_b_basic_B; f = 0;
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
    alpha = alpha_rchi;
    call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

    call cal_eta(rchi, eta)
    sum_Qeta_B0_0 = cal_sum_Qeta(Q_in, eta)

    call write_2Dfield(11,trim(output_folder)//"/eta-[B0_0]-A.bin",eta,nr-1,nz)
    call write_2Dfield(11,trim(output_folder)//"/rchi-[B0_0]-O.bin",rchi,nr,nz)

    print *, "Solving CHI with L(A,B=B0+dB,C) = 0 with boundary condition"
    solverB_B = solver_b_basic_B + solver_b_anomaly_B; f = 0;
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
    alpha = alpha_rchi;
    call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

    call cal_eta(rchi, eta)
    sum_Qeta_B0dB_0 = cal_sum_Qeta(Q_in, eta)

    call write_2Dfield(11,trim(output_folder)//"/eta-[B0dB_0]-A.bin",eta,nr-1,nz)
    call write_2Dfield(11,trim(output_folder)//"/rchi-[B0dB_0]-O.bin",rchi,nr,nz)


end if

rchi = 0.0
!==================================================================================
print *, "Solving CHI with L(A,B=0,C) = -dB"
solverB_B = 0.0; f = f_anomaly;
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_0_dB = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[0_dB]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[0_dB]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=B0,C) = -dB"
solverB_B = solver_b_basic_B; f = f_anomaly;
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_B0_dB = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[B0_dB]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[B0_dB]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=dB,C) = -dB"
solverB_B = solver_b_anomaly_B; f = f_anomaly;
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_dB_dB = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[dB_dB]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[dB_dB]-O.bin",rchi,nr,nz)

rchi = 0.0;
!
!==================================================================================
print *, "Solving CHI with L(A,B=B0+dB,C) = -dB"
solverB_B = solver_b_basic_B + solver_b_anomaly_B; f = f_anomaly;
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_B0dB_dB = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[B0dB_dB]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[B0dB_dB]-O.bin",rchi,nr,nz)

rchi = 0.0;
!==================================================================================
print *, "Solving CHI with L(A,B=0,C) = -B0"
solverB_B = 0.0; f = f_basic;
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_0_B0 = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[0_B0]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[0_B0]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=B0,C) = -B0"
solverB_B = solver_b_basic_B; f = f_basic;
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r;
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_B0_B0 = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[B0_B0]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[B0_B0]-O.bin",rchi,nr,nz)

!==================================================================================
print *, "Solving CHI with L(A,B=dB,C) = -B0"
solverB_B = solver_b_anomaly_B; f = f_basic
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_dB_B0 = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[dB_B0]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[dB_B0]-O.bin",rchi,nr,nz)


!==================================================================================
print *, "Solving CHI with L(A,B=B0+dB,C) = -B0"
solverB_B = solver_b_basic_B + solver_b_anomaly_B; f = f_basic
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

strategy = saved_strategy_rchi; strategy_r = saved_strategy_rchi_r
alpha = alpha_rchi;
call solve_elliptic(max_iter_rchi, strategy, strategy_r, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."

call cal_eta(rchi, eta)
sum_Qeta_B0dB_B0 = cal_sum_Qeta(Q_in, eta)

call write_2Dfield(11,trim(output_folder)//"/eta-[B0dB_B0]-A.bin",eta,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/rchi-[B0dB_B0]-O.bin",rchi,nr,nz)

!==================================================================================


! Integral Check

print *, "Integral check..."

! Apply boundary condition
rpsi = 0.0;
if(use_rpsi_bc .eqv. .true.) then
    rpsi = rpsi_bc
end if

print *, "Solving rpsi... L(A, B=0, C) = dJ/dr + dF/dz"
f = RHS_rpsi_thm + RHS_rpsi_mom
solverB_B = 0.0
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)
strategy = saved_strategy_rpsi; strategy_r = saved_strategy_rpsi_r;
alpha = alpha_rpsi
call solve_elliptic(max_iter_rpsi, strategy, strategy_r, alpha, rpsi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."
call rpsiToUW(rpsi, u_C, w_A)
call write_2Dfield(11, trim(output_folder)//"/rpsi_after-[0]-O.bin", rpsi, nr, nz)
call write_2Dfield(11,trim(output_folder)//"/w_after-[0]-A.bin",w_A,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/u_after-[0]-C.bin",u_C,nr,nz-1)

call cal_wtheta(w_A, theta, wtheta_B)
sum_wtheta_0_JF = cal_sum_wtheta(wtheta_B) * (g0/theta0)
call write_2Dfield(11,trim(output_folder)//"/wtheta_JF_after-[0]-B.bin",wtheta_B,nr-1,nz-1)


print *, "Solving rpsi... L(A, B=dB, C) = dJ/dr + dF/dz"
f = RHS_rpsi_thm + RHS_rpsi_mom
solverB_B = solver_b_anomaly_B
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)
strategy = saved_strategy_rpsi; strategy_r = saved_strategy_rpsi_r;
alpha = alpha_rpsi
call solve_elliptic(max_iter_rpsi, strategy, strategy_r, alpha, rpsi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."
call rpsiToUW(rpsi, u_C, w_A)
call write_2Dfield(11, trim(output_folder)//"/rpsi_after-[dB]-O.bin", rpsi, nr, nz)
call write_2Dfield(11,trim(output_folder)//"/w_after-[dB]-A.bin",w_A,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/u_after-[dB]-C.bin",u_C,nr,nz-1)

call cal_wtheta(w_A, theta, wtheta_B)
sum_wtheta_dB_JF = cal_sum_wtheta(wtheta_B) * (g0/theta0)
call write_2Dfield(11,trim(output_folder)//"/wtheta_JF_after-[dB]-B.bin",wtheta_B,nr-1,nz-1)


print *, "Solving rpsi... L(A, B=B0, C) = dJ/dr + dF/dz"
f = RHS_rpsi_thm + RHS_rpsi_mom
solverB_B = solver_b_basic_B
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)
strategy = saved_strategy_rpsi; strategy_r = saved_strategy_rpsi_r;
alpha = alpha_rpsi
call solve_elliptic(max_iter_rpsi, strategy, strategy_r, alpha, rpsi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."
call rpsiToUW(rpsi, u_C, w_A)
call write_2Dfield(11, trim(output_folder)//"/rpsi_after-[B0]-O.bin", rpsi, nr, nz)
call write_2Dfield(11,trim(output_folder)//"/w_after-[B0]-A.bin",w_A,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/u_after-[B0]-C.bin",u_C,nr,nz-1)

call cal_wtheta(w_A, theta, wtheta_B)
sum_wtheta_B0_JF = cal_sum_wtheta(wtheta_B) * (g0/theta0)
call write_2Dfield(11,trim(output_folder)//"/wtheta_JF_after-[B0]-B.bin",wtheta_B,nr-1,nz-1)

print *, "Solving rpsi... L(A, B=B0dB, C) = dJ/dr + dF/dz"
f = RHS_rpsi_thm + RHS_rpsi_mom
solverB_B = solver_b_basic_B + solver_b_anomaly_B
call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)
strategy = saved_strategy_rpsi; strategy_r = saved_strategy_rpsi_r;
alpha = alpha_rpsi
call solve_elliptic(max_iter_rpsi, strategy, strategy_r, alpha, rpsi, coe, f, wksp_O, nr, nz, err, &
    & merge(1,0,debug_mode .eqv. .true.))
print *, "Relaxation uses ", strategy, " steps. Final residue is ", strategy_r, "."
call rpsiToUW(rpsi, u_C, w_A)
call write_2Dfield(11, trim(output_folder)//"/rpsi_after-[B0dB]-O.bin", rpsi, nr, nz)
call write_2Dfield(11,trim(output_folder)//"/w_after-[B0dB]-A.bin",w_A,nr-1,nz)
call write_2Dfield(11,trim(output_folder)//"/u_after-[B0dB]-C.bin",u_C,nr,nz-1)

call cal_wtheta(w_A, theta, wtheta_B)
sum_wtheta_B0dB_JF = cal_sum_wtheta(wtheta_B) * (g0/theta0)
call write_2Dfield(11,trim(output_folder)//"/wtheta_JF_after-[B0dB]-B.bin",wtheta_B,nr-1,nz-1)


! Exchange conversion term check !


if(use_rchi_bc .eqv. .true.) then

    print *, "Exchange conversion term check..."
    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[0]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[0_0]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[0_dB]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call read_2Dfield(15, trim(output_folder)//"/rchi-[0_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv_0)
    call write_2Dfield(11,trim(output_folder)//"/bndconv-[0].bin", bndconv,nr-1,2)

    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[B0]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0_0]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0_dB]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv_B0)
    call write_2Dfield(11,trim(output_folder)//"/bndconv-[B0].bin", bndconv,nr-1,2)

    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[dB]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[dB_0]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[dB_dB]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call read_2Dfield(15, trim(output_folder)//"/rchi-[dB_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv_dB)
    call write_2Dfield(11,trim(output_folder)//"/bndconv-[dB].bin", bndconv,nr-1,2)

    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[B0dB]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_0]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_dB]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv_B0dB)
    call write_2Dfield(11,trim(output_folder)//"/bndconv-[B0dB].bin", bndconv,nr-1,2)

end if


! Boundary conversion method 2 !
if(use_rchi_bc .eqv. .true.) then

    print *, "Exchange conversion term check..."
    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[0]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[0_dB]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[0_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv2_0)
    call write_2Dfield(11,trim(output_folder)//"/bndconv2-[0].bin", bndconv,nr-1,2)

    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[B0]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0_dB]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv2_B0)
    call write_2Dfield(11,trim(output_folder)//"/bndconv2-[B0].bin", bndconv,nr-1,2)

    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[dB]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[dB_dB]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[dB_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv2_dB)
    call write_2Dfield(11,trim(output_folder)//"/bndconv2-[dB].bin", bndconv,nr-1,2)

    call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[B0dB]-O.bin", rpsi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_dB]-O.bin", rchi, nr, nz)
    call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
    call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv2_B0dB)
    call write_2Dfield(11,trim(output_folder)//"/bndconv2-[B0dB].bin", bndconv,nr-1,2)

end if




call cpu_time(time_end)

open (unit=15,file=trim(output_folder)//"/efficiency.txt",action="write",status="replace")
write (15,*) "Time elapsed (sec): ", (time_end - time_beg)
write (15,*) "sum Q                 : ", sum_Q

write (15,*) "# Boundary efficiency"
if(use_rchi_bc .eqv. .true.) then
    write (15,*) "eta [L(B=0)    = 0]      w/  boundary : ", sum_Qeta_0_0, ", ", (sum_Qeta_0_0 / sum_Q)
    write (15,*) "eta [L(B=dB)   = 0]      w/  boundary : ", sum_Qeta_dB_0, ", ", (sum_Qeta_dB_0 / sum_Q)
    write (15,*) "eta [L(B=B0)   = 0]      w/  boundary : ", sum_Qeta_B0_0, ", ", (sum_Qeta_B0_0 / sum_Q)
    write (15,*) "eta [L(B=B0dB) = 0]      w/  boundary : ", sum_Qeta_B0dB_0, ", ", (sum_Qeta_B0dB_0 / sum_Q)
end if

write (15,*) "# Internal efficiency"
write (15,*) "eta [L(B=0)    = dB]     wo/ boundary : ", sum_Qeta_0_dB, ", ", (sum_Qeta_0_dB / sum_Q)
write (15,*) "eta [L(B=dB)   = dB]     wo/ boundary : ", sum_Qeta_dB_dB, ", ", (sum_Qeta_dB_dB / sum_Q)
write (15,*) "eta [L(B=B0)   = dB]     wo/ boundary : ", sum_Qeta_B0_dB, ", ", (sum_Qeta_B0_dB / sum_Q)
write (15,*) "eta [L(B=B0dB) = dB]     wo/ boundary : ", sum_Qeta_B0dB_dB, ", ", (sum_Qeta_B0dB_dB / sum_Q)

write (15,*) "eta [L(B=0)    = B0]     wo/ boundary : ", sum_Qeta_0_B0, ", ", (sum_Qeta_0_B0 / sum_Q)
write (15,*) "eta [L(B=dB)   = B0]     wo/ boundary : ", sum_Qeta_dB_B0, ", ", (sum_Qeta_dB_B0 / sum_Q)
write (15,*) "eta [L(B=B0)   = B0]     wo/ boundary : ", sum_Qeta_B0_B0, ", ", (sum_Qeta_B0_B0 / sum_Q)
write (15,*) "eta [L(B=B0dB) = B0]     wo/ boundary : ", sum_Qeta_B0dB_B0, ", ", (sum_Qeta_B0dB_B0 / sum_Q)

if(use_rchi_bc .eqv. .true.) then
    write (15,*) "# Boundary conversion (Method 1)"
    write (15,*) "bndconv [L(B=0) = B0dB]   w/ boundary : ", sum_bndconv_0, ", ", (sum_bndconv_0 / sum_Q)
    write (15,*) "bndconv [L(B=dB) = B0dB]  w/ boundary : ", sum_bndconv_dB, ", ", (sum_bndconv_dB / sum_Q)
    write (15,*) "bndconv [L(B=B0) = B0dB]  w/ boundary : ", sum_bndconv_B0, ", ", (sum_bndconv_B0 / sum_Q)
    write (15,*) "bndconv [L(B=B0dB) = B0dB]w/ boundary : ", sum_bndconv_B0dB, ", ", (sum_bndconv_B0dB / sum_Q)

    write (15,*) "# Boundary conversion (Method 2)"
    write (15,*) "bndconv2 [L(B=0) = B0dB]   w/ boundary : ", sum_bndconv2_0, ", ", (sum_bndconv2_0 / sum_Q)
    write (15,*) "bndconv2 [L(B=dB) = B0dB]  w/ boundary : ", sum_bndconv2_dB, ", ", (sum_bndconv2_dB / sum_Q)
    write (15,*) "bndconv2 [L(B=B0) = B0dB]  w/ boundary : ", sum_bndconv2_B0, ", ", (sum_bndconv2_B0 / sum_Q)
    write (15,*) "bndconv2 [L(B=B0dB) = B0dB]w/ boundary : ", sum_bndconv2_B0dB, ", ", (sum_bndconv2_B0dB / sum_Q)

end if

write (15,*) "# Decomposition sum"

tmp1 = sum_Qeta_0_0 + sum_Qeta_0_dB + sum_Qeta_0_B0
if(use_rchi_bc .eqv. .true.) then; tmp1 = tmp1 + sum_bndconv_0; end if
write (15,*) "etaQ [L(B=0)    = J F] w/  boundary : ", tmp1, ", ", (tmp1 / sum_Q)
tmp1 = sum_Qeta_dB_0 + sum_Qeta_dB_dB + sum_Qeta_dB_B0
if(use_rchi_bc .eqv. .true.) then; tmp1 = tmp1 + sum_bndconv_dB; end if
write (15,*) "etaQ [L(B=B0)   = J F] w/  boundary : ", tmp1, ", ", (tmp1 / sum_Q)
tmp1 = sum_Qeta_B0_0 + sum_Qeta_B0_dB + sum_Qeta_B0_B0
if(use_rchi_bc .eqv. .true.) then; tmp1 = tmp1 + sum_bndconv_B0; end if
write (15,*) "etaQ [L(B=dB)   = J F] w/  boundary : ", tmp1, ", ", (tmp1 / sum_Q)
tmp1 = sum_Qeta_B0dB_0 + sum_Qeta_B0dB_dB + sum_Qeta_B0dB_B0
if(use_rchi_bc .eqv. .true.) then; tmp1 = tmp1 + sum_bndconv_B0dB; end if
write (15,*) "etaQ [L(B=B0dB) = J F] w/  boundary : ", tmp1, ", ", (tmp1 / sum_Q)


write (15,*) "# wtheta integral"

write (15,*) "wtheta [L(B=0)    = J F] w/  boundary : ", sum_wtheta_0_JF, ", ", (sum_wtheta_0_JF / sum_Q)
write (15,*) "wtheta [L(B=B0)   = J F] w/  boundary : ", sum_wtheta_B0_JF, ", ", (sum_wtheta_B0_JF / sum_Q)
write (15,*) "wtheta [L(B=dB)   = J F] w/  boundary : ", sum_wtheta_dB_JF, ", ", (sum_wtheta_dB_JF / sum_Q)
write (15,*) "wtheta [L(B=B0dB) = J F] w/  boundary : ", sum_wtheta_B0dB_JF, ", ", (sum_wtheta_B0dB_JF / sum_Q)

close (15)

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
! notice gradient of rpsi gets momentum flux, so we need to divide it by rho.
do i=1,nr-1
    do j=1,nz
        to_w(i,j) = to_w(i,j)/rho(j)
    end do
end do

do i=1,nr
    do j=1,nz-1
        r = ra(i)
        if(r /= 0) then 
            to_u(i,j) = to_u(i,j)/(r*(rho(j)+rho(j+1))/2.0)
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

real(4) function cal_sum_Q(Q)
implicit none
real(4) :: Q(nr-1, nz-1)
real(4) :: r, dr, dz, rho_
integer :: i,j

cal_sum_Q = 0.0;
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_Q = cal_sum_Q + Q(i,j) * rho_ * r * dr * dz
    end do
end do

end function

real(4) function cal_sum_Qeta(Q, eta)
implicit none
real(4) :: Q(nr-1, nz-1), eta(nr-1, nz)
real(4) :: r, dr, dz, rho_
integer :: i,j

cal_sum_Qeta = 0.0
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_Qeta = cal_sum_Qeta + ((eta(i,j) + eta(i,j+1)) /2.0) * Q(i,j) * rho_ * r * dr * dz
    end do
end do

end function

real(4) function cal_sum_wtheta(wtheta_B)
implicit none
real(4) :: wtheta_B(nr-1, nz-1)
real(4) :: r, dr, dz, rho_
integer :: i,j

cal_sum_wtheta = 0.0
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_wtheta = cal_sum_wtheta + wtheta_B(i,j) * rho_ * r * dr * dz
    end do
end do

end function



subroutine cal_wtheta(w_A, theta_B, wtheta_B)
implicit none
real(4) :: w_A(nr-1, nz), theta_B(nr-1, nz-1), wtheta_B(nr-1, nz-1)
integer :: i,j

do i = 1,nr-1
    do j = 1, nz-1
        wtheta_B(i,j) = ((w_A(i,j) + w_A(i,j+1)) /2.0) * theta_B(i,j)
    end do
end do
end subroutine

subroutine cal_eta(rchi, eta)
implicit none
real(4) :: rchi(nr,nz), eta(nr-1,nz)
integer :: i,j

call d_rdr_O2A(rchi, eta)

do i=1,nr-1
    do j=1,nz
        eta(i,j) = eta(i,j) * g0 / (rho(j) * Cp * exner(j) * theta0)
    end do
end do
end subroutine

subroutine cal_exchange_conversion(rpsi, rchi, rhoC, bndconv, sum_bndconv)
implicit none
real(4) :: rpsi(nr,nz), rchi(nr, nz), rhoC(nr,nz), bndconv(nr-1,2), sum_bndconv
integer :: i, r, dr, dz

sum_bndconv = 0.0
dz = za(2) - za(1)
dr = ra(2) - ra(1)
do i=1, nr-1
    r = (ra(i)+ra(i+1))/2.0
    
    ! Bottom
    bndconv(i, 1) = ((rhoC(i, 1) + rhoC(i+1, 1))/(2.0*rho(1))) * (&
    &    ( (rpsi(i,1) + rpsi(i+1,1))/2.0 ) * &
    &     ( (rchi(i,2) + rchi(i+1,2) - rchi(i,1) - rchi(i+1,1))/(2.0*dz)) - &
    &    ( (rchi(i,1) + rchi(i+1,1))/2.0 ) * &
    &     ( (rpsi(i,2) + rpsi(i+1,2) - rpsi(i,1) - rpsi(i+1,1))/(2.0*dz)) ) &
    &    / r**2.0

    ! Top
    bndconv(i, 2) = ((rhoC(i, nz) + rhoC(i+1, nz))/(2.0*rho(nz))) * (&
    &    ( (rpsi(i,nz) + rpsi(i+1,nz))/2.0 ) * &
    &     ( (rchi(i,nz) + rchi(i+1,nz) - rchi(i,nz-1) - rchi(i+1,nz-1))/(2.0*dz)) - &
    &    ( (rchi(i,nz) + rchi(i+1,nz))/2.0 ) * &
    &     ( (rpsi(i,nz) + rpsi(i+1,nz) - rpsi(i,nz-1) - rpsi(i+1,nz-1))/(2.0*dz)) ) &
    &    / r**2.0

    sum_bndconv = sum_bndconv - (bndconv(i,2) - bndconv(i,1)) * r * dr
end do


end subroutine


end program
