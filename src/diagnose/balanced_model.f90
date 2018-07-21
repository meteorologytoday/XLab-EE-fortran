program simulation
use elliptic_tools
use field_tools
use constants
use read_input_tools
use message_tools

print *, "Simulation of balanced vortex model"
real(4), parameter :: MATH_PI = acos(-1.0), RAD2DEG = 180.0 / MATH_PI, DEG2RAD = MATH_PI / 180.0
integer, parameter :: stdin=5, fd=15, CYLINDRICAL_MODE=0, SPHERICAL_MODE=1
integer        :: nr, nz


character(256) :: A_file, B_file, C_file, Q_file, F_file, input_folder, output_folder, &
&                 output_file, yes_or_no, rchi_bc_file, rpsi_bc_file, mode_str,          &
&                 word(4), buffer


real(4), pointer   :: rpsi(:,:), f(:,:), Q_in(:,:), coe(:, :, :),   &

real(4)            :: testing_dt, Lr(2), Lz(2), Lat(2), dr, dz, dlat,          &
                      eta_avg_b, eta_avg_nob,  &
                      time_beg, time_end

real(4)            :: m_ref
integer            :: m_ref_r, m_ref_z

integer :: saved_strategy_rpsi, max_iter_rpsi, &
&          saved_strategy_rchi, max_iter_rchi, strategy
real(4) :: saved_strategy_rpsi_r, saved_strategy_rchi_r, strategy_r, alpha,  &
&          alpha_rpsi, alpha_rchi

integer :: i,j, m, n, err, mode(size(word))
real(4) :: r, z, tmp1, tmp2, tmp3, planet_radius,&
&          sum_Q, sum_dtheta_dt, &
&          sum_Qeta_0_0  , sum_Qeta_dB_0  , sum_Qeta_B0_0  , sum_Qeta_B0dB_0 , &
&          sum_Qeta_0_dB , sum_Qeta_dB_dB , sum_Qeta_B0_dB , sum_Qeta_B0dB_dB, &
&          sum_Qeta_0_B0 , sum_Qeta_dB_B0 , sum_Qeta_B0_B0 , sum_Qeta_B0dB_B0, &
&          sum_wtheta_0_JF, sum_wtheta_dB_JF, sum_wtheta_B0_JF, sum_wtheta_B0dB_JF, &
&          sum_bndconv_0, sum_bndconv_dB, sum_bndconv_B0, sum_bndconv_B0dB, &
&          sum_bndconv2_0,sum_bndconv2_dB,sum_bndconv2_B0,sum_bndconv2_B0dB

logical :: file_exists, use_rchi_bc, use_rpsi_bc, debug_mode



inquire(file="./debug_mode", exist=debug_mode)

call cpu_time(time_beg)



! initialize A, B, C, Q

! Read setting
call read_input(stdin, buffer); read(buffer, *) Lr(1), Lr(2), Lz(1), Lz(2);
call read_input(stdin, buffer); read(buffer, *) dt
call read_input(stdin, buffer); read(buffer, *) timesteps


print *, "Testing time: ", testing_dt
print *, "Lr:", Lr(1), Lr(2)
print *, "Lz:", Lz(1), Lz(2)

call read_input(stdin, theta_file);
call read_input(stdin,     v_file);
call read_input(stdin,     Q_file);

call read_input(stdin, buffer);
read(buffer, *) rpsi_err_sum_threshold_abs, rpsi_err_sum_threshold_rel, max_iter_rpsi;



! turn v into m


! Loop

! for timestep
!     for rk4
!         backtrace point
!         determine which four points are responsible
!         interpolate
!
!     get next time step











! output






end subroutine

end program
