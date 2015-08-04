
f = forcing_in

if(operator_complexity == 0 .or. operator_complexity == 2) then
    print *, "Solving CHI with L(A,B=0,C) = -B"
    solverB_B = 0.0
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    max_iter_strf = saved_max_iter_strf
    strategy_r1 = saved_strategy_strf_r1
    strategy_r2 = saved_strategy_strf_r2
    alpha = alpha_strf
    strf = bc_init_in

    call solve_elliptic(max_iter_strf, 100, 10, 5, strategy_r1, strategy_r2, alpha, strf, coe, f, wksp_O, nr, nz, err, debug_mode)
    print *, "Relaxation uses ", max_iter_strf, " steps. Final residue is ", strategy_r1, ",", strategy_r2
    
    if(diag_param == DIAGPARAM_DYNAMIC_EFFICIENCY) then
        call cal_eta(strf, eta)
        call write_2Dfield(11,trim(output_folder)//"/eta-[BAROTROPIC]-A.bin",eta,nr-1,nz)
    else if(diag_param == DIAGPARAM_SECONDARY_CIRCULATION) then
        call cal_uw(strf, u_C, w_A)
        call write_2Dfield(11,trim(output_folder)//"/w-[BAROTROPIC]-A.bin",w_A,nr-1,nz)
        call write_2Dfield(11,trim(output_folder)//"/u-[BAROTROPIC]-C.bin",u_C,nr,nz-1)
    end if
    call write_2Dfield(11,trim(output_folder)//"/strf-[BAROTROPIC]-O.bin",strf,nr,nz)
end if

if(operator_complexity == 1 .or. operator_complexity == 2) then
    print *, "Solving CHI with L(A,B,C) = -B"
    solverB_B = saved_solverB_B
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    max_iter_strf = saved_max_iter_strf
    strategy_r1 = saved_strategy_strf_r1
    strategy_r2 = saved_strategy_strf_r2
    alpha = alpha_strf
    strf = bc_init_in

    call solve_elliptic(max_iter_strf, 100, 10, 5, strategy_r1, strategy_r2, alpha, strf, coe, f, wksp_O, nr, nz, err,debug_mode)
    print *, "Relaxation uses ", max_iter_strf, " steps. Final residue is ", strategy_r1, ",", strategy_r2

    if(diag_param == DIAGPARAM_DYNAMIC_EFFICIENCY) then
        call cal_eta(strf, eta)
        call write_2Dfield(11,trim(output_folder)//"/eta-[BAROCLINIC]-A.bin",eta,nr-1,nz)
    else if(diag_param == DIAGPARAM_SECONDARY_CIRCULATION) then
        call cal_uw(strf, u_C, w_A)
        call write_2Dfield(11,trim(output_folder)//"/w-[BAROCLINIC]-A.bin",w_A,nr-1,nz)
        call write_2Dfield(11,trim(output_folder)//"/u-[BAROCLINIC]-C.bin",u_C,nr,nz-1)
    end if

    call write_2Dfield(11,trim(output_folder)//"/strf-[BAROCLINIC]-O.bin",strf,nr,nz)

end if
