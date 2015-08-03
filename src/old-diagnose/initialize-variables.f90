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
allocate(sin_table(nr));

allocate(solverA_A(nr-1, nz-2)); allocate(solverB_B(nr-1, nz-1));allocate(solverC_C(nr-2, nz-1));
allocate(solver_b_basic_B(nr-1, nz-1));
allocate(solver_b_anomaly_B(nr-1, nz-1));
allocate(wtheta_B(nr-1,nz-1));
allocate(w_A(nr-1,nz)); allocate(u_C(nr,nz-1));
allocate(theta(nr-1,nz-1)); allocate(eta(nr-1,nz)); allocate(m2(nr-1,nz-1));
allocate(ra(nr));       allocate(za(nz));          allocate(rho(nz));
allocate(rcuva(nr));      
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


! ### Calculate dr, dz, dlat, ra, rcuva, za, exner, rho (pseudo-density)
dr = (Lr(2) - Lr(1)) / (nr-1); dz = (Lz(2) - Lz(1)) / (nz-1);

do i=1,nr
    ra(i) = Lr(1) + (i-1) * dr
end do

! Notice the mode(3) specifies rho to be constant or not
do j=1,nz
    za(j) = Lz(1) + (j-1) * dz
    exner(j) = merge(1.0 - za(j) / h0, 1.0, mode(3) == 0)
    rho(j) = merge(p0 / (theta0 * Rd) * exner(j)**(1.0 / kappa - 1.0), 1.0, &
&                   mode(3) == 0)
end do


if(mode(1) == 0) then
    rcuva = ra 
else if(mode(1) == 1) then
    dlat = (Lat(2) - Lat(1)) / (nr-1)
    do i=1, nr
        rcuva(i) = planet_radius * cos(Lat(1) + (i-1) * dlat)
        sin_table(i) = sin(Lat(1) + (i-1) * dlat)
    end do
end if


! ### sum Q ### !
sum_Q = cal_sum_Q(Q_in);

! ### Normalize coefficient solverA_A = a / (rho * r), solverB_B = b / (rho * r)
!                         , solverC_C = c / (rho * r)

do i=1,nr-1
    do j=1,nz-2
        solverA_A(i,j) = (rhoA_in(i,j+1) + rhoA_in(i+1,j+1))         &
            &           / (rcuva(i) + rcuva(i+1)) / rho(j+1)
    end do
end do


do i=1,nr-1
    do j=1,nz-1
        solverB_B(i,j) = ( rhoB_in(i  ,j  ) + rhoB_in(i+1,j  )       &
            &             + rhoB_in(i  ,j+1) + rhoB_in(i+1,j+1) )     &
            &             /(rcuva(i)+rcuva(i+1))/(rho(j)+rho(j+1))
    end do
end do

solver_b_basic_B = solverB_B

do i=1,nr-2
    do j=1,nz-1
        solverC_C(i,j) = (rhoC_in(i+1,j) + rhoC_in(i+1,j+1))         &
            &           / rcuva(i+1) / (rho(j)+rho(j+1))  
    end do
end do

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

if(mode(1) == CYLINDRICAL_MODE) then

    m2(1,:) = (((rcuva(2) - rcuva(1)) / 4.0)**3.0) * rhoC_C(i,j) * (ra(2) - ra(1)) / 2.0
    
    do i=1, nr-1
        do j=1, nz-1
            m2(i,j) = m2(i-1,j) + (rcuva(i)**3.0) * rhoC_C(i,j) * (ra(i+1) - ra(i-1)) / 2.0
        end do
    end do

else if(mode(1) == SPHERICAL_MODE) then

    m2(1,:) = (((rcuva(2) - rcuva(1)) / 4.0)**3.0) * rhoC_C(i,j) * (ra(2) - ra(1)) / 2.0 &
        &     / ((sin_table(2) + 3.0 * sin_table(1)) / 4.0)

    do i = 2, nr-1
        do j = 1, nz-1
            m2(i,j) = m2(i-1,j) + (rcuva(i)**3.0) * rhoC_C(i,j) * (ra(i+1) - ra(i-1)) / 2.0 &
                &    / sin_table(i)
        end do
    end do
end if

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



