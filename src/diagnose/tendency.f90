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


! ### Assign 2nt part of RHS = - 2 * rcuv^(-2) * d(mf)/dz
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
         RHS_rpsi_mom(i,j) = - (wksp_A(i,j) + wksp_A(i-1,j))/(rcuva(i)**2.0)
        
    end do
!    print *, "rcuva(",i,") = ", rcuva(i)
end do

    
!print *, "RHS_rpsi_mom max: " , maxval(abs(djdr))
call write_2Dfield(11, trim(output_folder)//"/RHS_rpsi_mom-O.bin", RHS_rpsi_mom, nr, nz)

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
    
    sum_dtheta_dt = integrate_weight_B(theta);
     
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
            &  / ((rcuva(i)+rcuva(i+1))/2.0)/((rho(j)+rho(j+1))/2.0)
        end do
    end do


