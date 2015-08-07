allocate(strf(nr,nz));   
allocate(f(nr,nz));
allocate(wksp_O(nr,nz));
allocate(wksp_A(nr-1,nz));
allocate(wksp_B(nr-1,nz-1));
allocate(wksp_C(nr,nz-1));

allocate(coe(9,nr,nz));
allocate(rhoA_in(nr-1, nz-2));  allocate(rhoB_in(nr-1, nz-1));   allocate(rhoC_in(nr-2, nz-1));
allocate(Q_in(nr-1, nz-1)); allocate(F_in(nr-1, nz-1)); 
allocate(rhoA_A(nr-1, nz)); allocate(rhoB_C(nr, nz-1));  allocate(rhoB_B(nr-1,nz-1));
allocate(rhoC_C(nr, nz-1)); 
allocate(sin_table(nr));

allocate(solverA_A(nr-1, nz-2)); allocate(solverB_B(nr-1, nz-1));allocate(solverC_C(nr-2, nz-1));
allocate(saved_solverB_B(nr-1, nz-1));

allocate(ra(nr));       allocate(za(nz));          allocate(rho(nz));
allocate(rcuva(nr));
allocate(exner(nz));

allocate(w_A(nr-1,nz));
allocate(u_C(nr,nz-1));

print *, "Allocation complete."

call read_2Dfield(15, trim(input_folder)//"/"//A_file, rhoA_in, nr-1, nz-2)
call read_2Dfield(15, trim(input_folder)//"/"//B_file, rhoB_in, nr-1, nz-1)
call read_2Dfield(15, trim(input_folder)//"/"//C_file, rhoC_in, nr-2, nz-1)
call read_2Dfield(15, trim(input_folder)//"/"//Q_file, rhoB_in, nr-1, nz-1)
call read_2Dfield(15, trim(input_folder)//"/"//F_file, rhoB_in, nr-1, nz-1)

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

if(geometry == CYLINDRICAL_MODE) then
    rcuva = ra 
else if(geometry == SPHERICAL_MODE) then
    dlat = (Lat(2) - Lat(1)) / (nr-1)
    do i=1, nr
        rcuva(i) = planet_radius * cos(Lat(1) + (i-1) * dlat)
        sin_table(i) = sin(Lat(1) + (i-1) * dlat)
    end do
end if

print *, "Geometry complete."

