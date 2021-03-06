integer, parameter :: ERROR_INPUT = 1


real(4), parameter :: MATH_PI = acos(-1.0), RAD2DEG = 180.0 / MATH_PI, DEG2RAD = MATH_PI / 180.0

integer, parameter :: stdin=5, fd=15, &
&                     CYLINDRICAL_MODE=0, SPHERICAL_MODE=1, &
&                     DIAGPARAM_DYNAMIC_EFFICIENCY = 0, &
&                     DIAGPARAM_SECONDARY_CIRCULATION = 1, &
&                     DIAGPARAM_NONE = 2, &
&                     DENSITY_NORMAL = 0, &
&                     DENSITY_BOUSSINESQ = 1, &
&                     BARO_ALL = 2, &
&                     BAROCLINIC = 1, &
&                     BAROTROPIC = 0
                    

integer        :: nr, nz, diag_param
character(256) :: A_file, B_file, C_file, forcing_file, input_folder, output_folder, &
&                 output_file, bc_init_file, mode_str,          &
&                 word(4), buffer, format_str



! Grid point are designed as follows: 
! O A O
! C B C
! O A O
! 
! O  : (nr   , nz  ) rpsi strf rhoA_in rhoB_in rhoC_in Q_in f dthetadr
! A  : (nr-1 , nz  ) w eta
! B  : (nr-1 , nz-1) m theta F Q solver_B
! C  : (nr   , nz-1) u
! sA : (nr-1 , nz-2) solver_A
! sC : (nr-2 , nz-1) solver_C
!


real(4), pointer   :: strf(:,:), f(:,:),coe(:, :, :),   &
&                     sin_table(:),                                 &
&                     rhoA_in(:,:), rhoB_in(:,:), rhoC_in(:,:),                &
&                     forcing_in(:,:), &
&                     wksp_O(:,:), wksp_A(:,:), wksp_B(:,:), wksp_C(:,:),      &
&                     solverA_A(:,:), solverB_B(:,:), solverC_C(:,:),          &
&                     saved_solverB_B(:,:),          &
&                     rhoA_A(:,:), rhoB_B(:,:), rhoB_C(:,:), rhoC_C(:,:),      &
&                     eta(:,:), u_C(:,:), w_A(:,:),      &
&                     ra(:), rcuva(:), za(:), exner(:), rho(:),                &
&                     bc_init_in(:,:)

real(4)            :: Lr(2), Lz(2), Lat(2), dr, dz, dlat,          &
                      eta_avg_b, eta_avg_nob,  &
                      time_beg, time_end

integer :: saved_strategy_strf, max_iter_strf, saved_max_iter_strf, strategy
real(4) :: saved_strategy_strf_r1, saved_strategy_strf_r2, strategy_r1, strategy_r2, &
&          alpha, alpha_strf

integer :: i,j, m, n, err, geometry, density_mode, operator_complexity
real(4) :: r, z, tmp1, tmp2, tmp3, planet_radius

logical :: file_exists
integer :: debug_mode

