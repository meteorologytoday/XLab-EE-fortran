! Read setting
call read_input(stdin, mode_str)
do i=1, size(word)
    if(split_line(mode_str, word(i), "-") /= 0 .and. i /= size(word)) then
        write (format_str, "(A, I1, A)") "There should be ", size(word), " inputs&
                                        & separated by dashes."
        call error_msg("INIT", ERROR_INPUT, format_str)
        stop
    end if
end do

mode(:) = 0
if(word(1) == 'CYLINDRICAL') then
    mode(1) = CYLINDRICAL_MODE
else if(word(1) == 'SPHERICAL') then
    mode(1) = SPHERICAL_MODE
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(1)) // "]")
    stop
end if

if(word(2) == 'TENDENCY') then
    mode(2) = 0
else if(word(2) == 'INSTANT') then
    mode(2) = 1
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(2)) // "]")
    stop
end if

if(word(3) == 'DENSITY_NORMAL') then
    mode(3) = 0
else if(word(3) == 'DENSITY_BOUSSINESQ') then
    mode(3) = 1
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(3)) // "]")
    stop
end if
if(word(4) == 'BARO_ALL') then
    mode(4) = 2
else if(word(4) == 'BAROCLINIC') then
    mode(4) = 1
else if(word(4) == 'BAROTROPIC') then
    mode(4) = 0
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(4)) // "]" )
    stop
end if


! TENDENCY MODE
if(mode(2) == 0) then
    call read_input(stdin, buffer); read(buffer, *) testing_dt
end if

! Domain size. In spherical mode, domain is by default from south pole
! to north pole
if(mode(1) == CYLINDRICAL_MODE) then
    call read_input(stdin, buffer); read(buffer, *) Lr(1), Lr(2), Lz(1), Lz(2);

    if(Lr(2) .le. Lr(1)) then
        call error_msg("INIT", ERROR_INPUT, "Domain size in radial direction must be positive.")
    end if

    if(Lz(2) .le. Lz(1)) then
        call error_msg("INIT", ERROR_INPUT, "Domain size in z direction must be positive.")
    end if

else if(mode(1) == SPHERICAL_MODE) then
    call read_input(stdin, buffer); read(buffer, *) planet_radius, Lz(1), Lz(2);
    Lat(1) = -90.0; Lat(2) = 90.0;
    Lr(1) = Lat(1) * DEG2RAD * planet_radius
    Lr(2) = Lat(2) * DEG2RAD * planet_radius

    if(Lz(2) .le. Lz(1)) then
        call error_msg("INIT", ERROR_INPUT, "Domain size in z direction must be positive.")
    end if

end if

call read_input(stdin, buffer); read(buffer, *) nr, nz;
call read_input(stdin, input_folder);
call read_input(stdin, output_folder);
call read_input(stdin, A_file);
call read_input(stdin, B_file);
call read_input(stdin, C_file);
call read_input(stdin, Q_file);
call read_input(stdin, F_file);

call read_input(stdin, buffer);
read(buffer, *) saved_strategy_rpsi, saved_strategy_rpsi_r, max_iter_rpsi, alpha_rpsi;

call read_input(stdin, buffer);
read(buffer, *) saved_strategy_rchi, saved_strategy_rchi_r, max_iter_rchi, alpha_rchi;

call read_input(stdin, yes_or_no)
if(yes_or_no == 'yes') then
    call read_input(stdin, rpsi_bc_file);
    use_rpsi_bc = .true.
else
    use_rpsi_bc = .false.
end if
call read_input(stdin, yes_or_no)
if(yes_or_no == 'yes') then
    call read_input(stdin, rchi_bc_file);
    use_rchi_bc = .true.
else
    use_rchi_bc = .false.
end if


print *, "mode: ", mode(1), ",", mode(2),",",mode(3),",",mode(4)
if(mode(2) == 0) then
    print *, "Testing time: ", testing_dt
end if
if(mode(1) == CYLINDRICAL_MODE) then 
    print *, "Lr:", Lr(1), Lr(2)
    print *, "Lz:", Lz(1), Lz(2)
else if(mode(1) == SPHERICAL_MODE) then
    print *, "Using spherical mode, domain is forced to be global."
    print *, "Planet Radius: ", planet_radius
    print *, "Lat:", Lat(1), Lat(2)
    print *, "Lz:", Lz(1), Lz(2)
end if

print *, "nr:", nr, ", nz:", nz
print *, "Input folder: ", trim(input_folder)
print *, "Output folder: ", trim(output_folder)
print *, "A file: ", trim(A_file)
print *, "B file: ", trim(B_file)
print *, "C file: ", trim(C_file)
print *, "Q file: ", trim(Q_file)
print *, "F file: ", trim(F_file)
print *, "rpsi's strategy, residue, iter: ", saved_strategy_rpsi, &
&           saved_strategy_rpsi_r, max_iter_rpsi, alpha_rpsi
print *, "rchi's strategy, residue, iter: ", saved_strategy_rchi, &
&           saved_strategy_rchi_r, max_iter_rchi, alpha_rchi

if(use_rpsi_bc .eqv. .true.) then
    print *, "Use rpsi boundary condition: Yes (", trim(rpsi_bc_file), ")"
else
    print *, "Use rpsi boundary condition: No "
end if

if(use_rchi_bc .eqv. .true.) then
    print *, "Use rchi boundary condition: Yes (", trim(rchi_bc_file), ")"
else
    print *, "Use rpsi boundary condition: No "
end if

