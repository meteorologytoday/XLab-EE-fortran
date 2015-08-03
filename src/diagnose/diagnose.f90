if(operator_complexity == 0 .or. operator_complexity == 2) then
    print *, "Solving CHI with L(A,B=0,C) = -B"
    solverB_B = 0.0
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    max_iter_rchi = saved_max_iter_rchi
    strategy_r1 = saved_strategy_rchi_r1
    strategy_r2 = saved_strategy_rchi_r2
    alpha = alpha_rchi
    rchi = 0.0

    call solve_elliptic(max_iter_rchi, 100, 10, 5, strategy_r1, strategy_r2, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", max_iter_rchi, " steps. Final residue is ", strategy_r1, ",", strategy_r2

    call cal_eta(rchi, eta)
    call write_2Dfield(11,trim(output_folder)//"/eta-[0]-A.bin",eta,nr-1,nz)
    call write_2Dfield(11,trim(output_folder)//"/rchi-[0]-O.bin",rchi,nr,nz)
end if

if(operator_complexity == 1 .or. operator_complexity == 2) then
    print *, "Solving CHI with L(A,B,C) = -B"
    solverB_B = saved_solverB_B
    call cal_coe(solverA_A, solverB_B, solverC_C, coe, dr, dz, nr, nz, err)

    max_iter_rchi = saved_max_iter_rchi
    strategy_r1 = saved_strategy_rchi_r1
    strategy_r2 = saved_strategy_rchi_r2
    alpha = alpha_rchi
    rchi = 0.0

    call solve_elliptic(max_iter_rchi, 100, 10, 5, strategy_r1, strategy_r2, alpha, rchi, coe, f, wksp_O, nr, nz, err, &
        & merge(1,0,debug_mode .eqv. .true.))
    print *, "Relaxation uses ", max_iter_rchi, " steps. Final residue is ", strategy_r1, ",", strategy_r2

    call cal_eta(rchi, eta)
    call write_2Dfield(11,trim(output_folder)//"/eta-[B]-A.bin",eta,nr-1,nz)
    call write_2Dfield(11,trim(output_folder)//"/rchi-[B]-O.bin",rchi,nr,nz)
end if
