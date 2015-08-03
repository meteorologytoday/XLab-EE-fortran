allocate(rchi(nr,nz));   
allocate(f(nr,nz));
allocate(wksp_O(nr,nz));
allocate(wksp_A(nr-1,nz));
allocate(wksp_B(nr-1,nz-1));
allocate(wksp_C(nr,nz-1));
allocate(f_basic(nr,nz));

allocate(coe(9,nr,nz));
allocate(rhoA_in(nr, nz));  allocate(rhoB_in(nr, nz));   allocate(rhoC_in(nr, nz));
allocate(rhoA_A(nr-1, nz)); allocate(rhoB_C(nr, nz-1));  allocate(rhoB_B(nr-1,nz-1));
allocate(rhoC_C(nr, nz-1)); 
allocate(b_basic_B(nr-1,nz-1));                    allocate(b_anomaly_B(nr-1,nz-1));
allocate(sin_table(nr));

allocate(solverA_A(nr-1, nz-2)); allocate(solverB_B(nr-1, nz-1));allocate(solverC_C(nr-2, nz-1));
allocate(saved_solverB_B(nr-1, nz-1));
allocate(theta(nr-1,nz-1)); allocate(eta(nr-1,nz));
allocate(ra(nr));       allocate(za(nz));          allocate(rho(nz));
allocate(rcuva(nr));
allocate(exner(nz));    allocate(bndconv(nr-1,2));

call read_2Dfield(15, trim(input_folder)//"/"//A_file, rhoA_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//B_file, rhoB_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//C_file, rhoC_in, nr, nz)

if(use_rchi_bc .eqv. .true.) then
    allocate(rchi_bc(nr,nz));
    call read_2Dfield(15, trim(input_folder)//"/"//rchi_bc_file, rchi_bc, nr, nz)
end if

! ### Calculate dr, dz, dlat, ra, rcuva, za, exner, rho (pseudo-density)
dr = (Lr(2) - Lr(1)) / (nr-1); dz = (Lz(2) - Lz(1)) / (nz-1);

do i=1,nr
    ra(i) = Lr(1) + (i-1) * dr
end do

! Notice the density_mode specifies rho to be constant or not
do j=1,nz
    za(j) = Lz(1) + (j-1) * dz
    exner(j) = merge(1.0 - za(j) / h0, 1.0, density_mode == 0)
    rho(j) = merge(p0 / (theta0 * Rd) * exner(j)**(1.0 / kappa - 1.0), 1.0, &
&                   density_mode == 0)
end do

if(geometry == 0) then
    rcuva = ra 
else if(geometry == 1) then
    dlat = (Lat(2) - Lat(1)) / (nr-1)
    do i=1, nr
        rcuva(i) = planet_radius * cos(Lat(1) + (i-1) * dlat)
        sin_table(i) = sin(Lat(1) + (i-1) * dlat)
    end do
end if

! ### Calculate solverA_A, solverB_B, solverC_C

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

call write_2Dfield(11, trim(output_folder)//"/solver_a-sA.bin", solverA_A, nr-1, nz-2)
call write_2Dfield(11, trim(output_folder)//"/solver_b-B.bin", solverB_B, nr-1, nz-1)
call write_2Dfield(11, trim(output_folder)//"/solver_c-sC.bin", solverC_C, nr-2, nz-1)

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
