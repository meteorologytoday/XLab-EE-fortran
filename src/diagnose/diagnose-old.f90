program diagnose
use elliptic_tools
use field_tools
use constants
use read_input_tools
use message_tools
implicit none

include "variables.f90"


inquire(file="./debug_mode", exist=debug_mode);

call cpu_time(time_beg)

include "read-control-file.f90"
include "initialize-variables.f90"
print *, "Initialization complete."

! TENDENCY MODE
if(mode(2) == 0) then
    include "tendency.f90"
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

    if(mode(4) == 0 .or. mode(4) == 2) then
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

    end if

    if(mode(4) == 1 .or. mode(4) == 2) then
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

end if

rchi = 0.0


!==================================================================================
if(mode(4) == 0 .or. mode(4) == 2) then

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
end if


!==================================================================================
if(mode(4) == 1 .or. mode(4) == 2) then
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
end if


!==================================================================================
if(mode(4) == 0 .or. mode(4) == 2) then

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

end if


!==================================================================================
if(mode(4) == 1 .or. mode(4) == 2) then
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
end if
!==================================================================================


! Integral Check

print *, "Integral check..."

! Apply boundary condition
rpsi = 0.0;
if(use_rpsi_bc .eqv. .true.) then
    rpsi = rpsi_bc
end if

if(mode(4) == 0 .or. mode(4) == 2) then
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
end if

if(mode(4) == 1 .or. mode(4) == 2) then
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
end if

! Exchange conversion term check !


if(use_rchi_bc .eqv. .true.) then

    print *, "Exchange conversion term check..."
    if(mode(4) == 0 .or. mode(4) == 2) then
        call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[0]-O.bin", rpsi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[0_0]-O.bin", rchi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[0_dB]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
        call read_2Dfield(15, trim(output_folder)//"/rchi-[0_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
        call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv_0)
        call write_2Dfield(11,trim(output_folder)//"/bndconv-[0].bin", bndconv,nr-1,2)
    end if

    if(mode(4) == 1 .or. mode(4) == 2) then
        call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[B0dB]-O.bin", rpsi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_0]-O.bin", rchi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_dB]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
        call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
        call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv_B0dB)
        call write_2Dfield(11,trim(output_folder)//"/bndconv-[B0dB].bin", bndconv,nr-1,2)
    end if
end if


! Boundary conversion method 2 !
if(use_rchi_bc .eqv. .true.) then

    print *, "Exchange conversion term check..."
    if(mode(4) == 0 .or. mode(4) == 2) then
        call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[0]-O.bin", rpsi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[0_dB]-O.bin", rchi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[0_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
        call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv2_0)
        call write_2Dfield(11,trim(output_folder)//"/bndconv2-[0].bin", bndconv,nr-1,2)
    end if

    if(mode(4) == 1 .or. mode(4) == 2) then
        call read_2Dfield(15, trim(output_folder)//"/rpsi_after-[B0dB]-O.bin", rpsi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_dB]-O.bin", rchi, nr, nz)
        call read_2Dfield(15, trim(output_folder)//"/rchi-[B0dB_B0]-O.bin", wksp_O, nr, nz); rchi = rchi + wksp_O
        call cal_exchange_conversion(rpsi, rchi, rhoC_in, bndconv, sum_bndconv2_B0dB)
        call write_2Dfield(11,trim(output_folder)//"/bndconv2-[B0dB].bin", bndconv,nr-1,2)
    end if
end if




call cpu_time(time_end)

if(mode(4) == 0 .or. mode(4) == 2) then
    open (unit=15,file=trim(output_folder)//"/efficiency.txt",action="write",status="replace")
    write (15,*) "Time elapsed (sec)                          : ", (time_end - time_beg)
    write (15,*) "sum Q                                       : ", sum_Q
    write (15,*) "sum dtheta_dt                               : ", sum_dtheta_dt
    write (15,*) "Local heat response (sum Q / sum dtheta_dt) : ", sum_dtheta_dt / sum_Q

    write (15,*) "# Boundary efficiency"
    if(use_rchi_bc .eqv. .true.) then
        write (15,*) "eta [L(B=0)    = 0]      w/  boundary : ", sum_Qeta_0_0, ", ", (sum_Qeta_0_0 / sum_Q)
    end if

    write (15,*) "# Internal efficiency"
    write (15,*) "eta [L(B=0)    = dB]     wo/ boundary : ", sum_Qeta_0_dB, ", ", (sum_Qeta_0_dB / sum_Q)
    write (15,*) "eta [L(B=0)    = B0]     wo/ boundary : ", sum_Qeta_0_B0, ", ", (sum_Qeta_0_B0 / sum_Q)

    if(use_rchi_bc .eqv. .true.) then
        write (15,*) "# Boundary conversion (Method 1)"
        write (15,*) "bndconv [L(B=0) = B0dB]   w/ boundary : ", sum_bndconv_0, ", ", (sum_bndconv_0 / sum_Q)

        write (15,*) "# Boundary conversion (Method 2)"
        write (15,*) "bndconv2 [L(B=0) = B0dB]   w/ boundary : ", sum_bndconv2_0, ", ", (sum_bndconv2_0 / sum_Q)
    end if

    write (15,*) "# Decomposition sum"

    tmp1 = sum_Qeta_0_0 + sum_Qeta_0_dB + sum_Qeta_0_B0
    if(use_rchi_bc .eqv. .true.) then; tmp1 = tmp1 + sum_bndconv_0; end if
    write (15,*) "etaQ [L(B=0)    = J F] w/  boundary : ", tmp1, ", ", (tmp1 / sum_Q)

    write (15,*) "# wtheta integral"
    write (15,*) "wtheta [L(B=0)    = J F] w/  boundary : ", sum_wtheta_0_JF, ", ", (sum_wtheta_0_JF / sum_Q)
end if

if(mode(4) == 1 .or. mode(4) == 2) then
    write (15,*) "# Boundary efficiency"
    if(use_rchi_bc .eqv. .true.) then
        write (15,*) "eta [L(B=B0dB) = 0]      w/  boundary : ", sum_Qeta_B0dB_0, ", ", (sum_Qeta_B0dB_0 / sum_Q)
    end if

    write (15,*) "# Internal efficiency"
    write (15,*) "eta [L(B=B0dB) = dB]     wo/ boundary : ", sum_Qeta_B0dB_dB, ", ", (sum_Qeta_B0dB_dB / sum_Q)
    write (15,*) "eta [L(B=B0dB) = B0]     wo/ boundary : ", sum_Qeta_B0dB_B0, ", ", (sum_Qeta_B0dB_B0 / sum_Q)

    if(use_rchi_bc .eqv. .true.) then
        write (15,*) "# Boundary conversion (Method 1)"
        write (15,*) "bndconv [L(B=B0dB) = B0dB]w/ boundary : ", sum_bndconv_B0dB, ", ", (sum_bndconv_B0dB / sum_Q)

        write (15,*) "# Boundary conversion (Method 2)"
        write (15,*) "bndconv2 [L(B=B0dB) = B0dB]w/ boundary : ", sum_bndconv2_B0dB, ", ", (sum_bndconv2_B0dB / sum_Q)
    end if

    write (15,*) "# Decomposition sum"

    tmp1 = sum_Qeta_B0dB_0 + sum_Qeta_B0dB_dB + sum_Qeta_B0dB_B0
    if(use_rchi_bc .eqv. .true.) then; tmp1 = tmp1 + sum_bndconv_B0dB; end if
    write (15,*) "etaQ [L(B=B0dB) = J F] w/  boundary : ", tmp1, ", ", (tmp1 / sum_Q)

    write (15,*) "# wtheta integral"
    write (15,*) "wtheta [L(B=B0dB) = J F] w/  boundary : ", sum_wtheta_B0dB_JF, ", ", (sum_wtheta_B0dB_JF / sum_Q)
end if

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

call d_rcuvdr_O2A(from_rpsi, to_w)
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
            to_u(i,j) = to_u(i,j)/(rcuva(i)*(rho(j)+rho(j+1))/2.0)
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


subroutine d_rcuvdr_O2A(from_dat, to_dat)
! notice that all the points with r=0 become inifinite
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr-1,nz)
integer :: i,j

call d_dr_O2A(from_dat, to_dat)

do i = 1, nr-1
    do j = 1, nz
        to_dat(i,j) = to_dat(i,j) / ((rcuva(i) + rcuva(i+1)) / 2.0)
    end do
end do

end subroutine

real(4) function integrate_weight_B(weight)
implicit none
real(4) :: weight(nr-1, nz-1)
real(4) :: r, rcuv, dr, dz, rho_
integer :: i,j

integrate_weight_B = 0.0;
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        rcuv = (rcuva(i)+rcuva(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        integrate_weight_B = integrate_weight_B + weight(i,j) * rho_ * rcuv * dr * dz
    end do
end do

end function

real(4) function cal_sum_Q(Q)
implicit none
real(4) :: Q(nr-1, nz-1)
real(4) :: r, rcuv, dr, dz, rho_
integer :: i,j

cal_sum_Q = 0.0;
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        rcuv = (rcuva(i)+rcuva(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_Q = cal_sum_Q + Q(i,j) * rho_ * rcuv * dr * dz
    end do
end do

cal_sum_Q = integrate_weight_B(Q)

end function

real(4) function cal_sum_Qeta(Q, eta)
implicit none
real(4) :: Q(nr-1, nz-1), eta(nr-1, nz)
real(4) :: r, rcuv, dr, dz, rho_
integer :: i,j

cal_sum_Qeta = 0.0
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        rcuv = (rcuva(i)+rcuva(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_Qeta = cal_sum_Qeta + ((eta(i,j) + eta(i,j+1)) /2.0) * Q(i,j) * rho_ * rcuv * dr * dz
    end do
end do

end function

real(4) function cal_sum_wtheta(wtheta_B)
implicit none
real(4) :: wtheta_B(nr-1, nz-1)
real(4) :: r, rcuv, dr, dz, rho_
integer :: i,j

cal_sum_wtheta = 0.0
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        rcuv = (rcuva(i)+rcuva(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_wtheta = cal_sum_wtheta + wtheta_B(i,j) * rho_ * rcuv * dr * dz
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

call d_rcuvdr_O2A(rchi, eta)

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
