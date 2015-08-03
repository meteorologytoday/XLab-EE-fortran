allocate(strf(nr,nz));   
allocate(f(nr,nz));
allocate(wksp_O(nr,nz));
allocate(wksp_A(nr-1,nz));
allocate(wksp_B(nr-1,nz-1));
allocate(wksp_C(nr,nz-1));

allocate(coe(9,nr,nz));
allocate(rhoA_in(nr, nz));  allocate(rhoB_in(nr, nz));   allocate(rhoC_in(nr, nz));
allocate(forcing_in(nr, nz)); allocate(bc_init_in(nr, nz));
allocate(rhoA_A(nr-1, nz)); allocate(rhoB_C(nr, nz-1));  allocate(rhoB_B(nr-1,nz-1));
allocate(rhoC_C(nr, nz-1)); 
allocate(sin_table(nr));

allocate(solverA_A(nr-1, nz-2)); allocate(solverB_B(nr-1, nz-1));allocate(solverC_C(nr-2, nz-1));
allocate(saved_solverB_B(nr-1, nz-1));
allocate(eta(nr-1,nz));
allocate(ra(nr));       allocate(za(nz));          allocate(rho(nz));
allocate(rcuva(nr));
allocate(exner(nz));
print *, "Allocation complete."

call read_2Dfield(15, trim(input_folder)//"/"//A_file, rhoA_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//B_file, rhoB_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//C_file, rhoC_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//forcing_file, forcing_in, nr, nz)
call read_2Dfield(15, trim(input_folder)//"/"//bc_init_file, bc_init_in, nr, nz)

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

print *, "Geometry complete."
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

saved_solverB_B = solverB_B

do i=1,nr-2
    do j=1,nz-1
        solverC_C(i,j) = (rhoC_in(i+1,j) + rhoC_in(i+1,j+1))         &
            &           / rcuva(i+1) / (rho(j)+rho(j+1))  
    end do
end do

print *, "Solver coe part I complete."
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

print *, "Solver coe part II complete."

do i=1,nr
    do j=1,nz-1
        rhoC_C(i,j) = (rhoC_in(i,j) + rhoC_in(i,j+1)) / 2.0
    end do
end do

call write_2Dfield(11, trim(output_folder)//"/solver_a-sA.bin", solverA_A, nr-1, nz-2)
call write_2Dfield(11, trim(output_folder)//"/solver_b-B.bin", solverB_B, nr-1, nz-1)
call write_2Dfield(11, trim(output_folder)//"/solver_c-sC.bin", solverC_C, nr-2, nz-1)

print *, "Solver complete."
! doing right hand side
!f = 0.0
!do i=2,nr-1
!    do j=2,nz-1
!        f(i,j) = - ( rhoB_B(i-1,j-1) &
!            &      + rhoB_B(i-1,j  ) &
!            &      + rhoB_B(i  ,j  ) &
!            &      + rhoB_B(i  ,j-1) )/4.0
!    end do
!end do
