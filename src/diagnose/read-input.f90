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

if(word(1) == 'DYNAMIC_EFFICIENCY') then
    diag_param = DIAGPARAM_DYNAMIC_EFFICIENCY
else if(word(1) == 'SECONDARY_CIRCULATION') then
    diag_param = DIAGPARAM_SECONDARY_CIRCULATION
else if(word(1) == 'NONE') then
    diag_param = DIAGPARAM_NONE
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(1)) // "]")
    stop
end if


if(word(2) == 'CYLINDRICAL') then
    geometry = CYLINDRICAL_MODE
else if(word(2) == 'SPHERICAL') then
    geometry = SPHERICAL_MODE
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(2)) // "]")
    stop
end if

if(word(3) == 'DENSITY_NORMAL') then
    density_mode = 0
else if(word(3) == 'DENSITY_BOUSSINESQ') then
    density_mode = 1
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(3)) // "]")
    stop
end if

if(word(4) == 'BARO_ALL') then
    operator_complexity = 2
else if(word(4) == 'BAROCLINIC') then
    operator_complexity = 1
else if(word(4) == 'BAROTROPIC') then
    operator_complexity = 0
else
    call error_msg("INIT", ERROR_INPUT, "Unknown Mode [" // trim(word(4)) // "]" )
    stop
end if

! Domain size. In spherical mode, domain is by default from south pole
! to north pole
if(geometry == CYLINDRICAL_MODE) then
    call read_input(stdin, buffer); read(buffer, *) Lr(1), Lr(2), Lz(1), Lz(2);

    if(Lr(2) .le. Lr(1)) then
        call error_msg("INIT", ERROR_INPUT, "Domain size in radial direction must be positive.")
    end if

    if(Lz(2) .le. Lz(1)) then
        call error_msg("INIT", ERROR_INPUT, "Domain size in z direction must be positive.")
    end if

else if(geometry == SPHERICAL_MODE) then
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
if (diag_param == DIAGPARAM_SECONDARY_CIRCULATION) then
    call read_input(stdin, forcing_file);
end if
call read_input(stdin, bc_init_file);
call read_input(stdin, buffer);
read(buffer, *) saved_strategy_strf_r1, saved_strategy_strf_r2, saved_max_iter_strf, alpha_strf;

print *, "----- Diagnose Input -----"
print *, "Diagnose parameter: ", diag_param
print *, "Geometry: ", geometry
print *, "Density distribution: ", density_mode
print *, "Operator complexity: ", operator_complexity
if(geometry == CYLINDRICAL_MODE) then 
    print *, "Lr:", Lr(1), Lr(2)
    print *, "Lz:", Lz(1), Lz(2)
else if(geometry == SPHERICAL_MODE) then
    print *, "Using spherical mode, domain is forced to be global."
    print *, "Planet Radius: ", planet_radius
    print *, "Lat:", Lat(1), Lat(2)
    print *, "Lz:", Lz(1), Lz(2)
end if

print *, "nr:", nr, ", nz:", nz
print *, "Input folder:  ", trim(input_folder)
print *, "Output folder: ", trim(output_folder)
print *, "A file:        ", trim(A_file)
print *, "B file:        ", trim(B_file)
print *, "C file:        ", trim(C_file)
if (diag_param == DIAGPARAM_SECONDARY_CIRCULATION) then
    print *, "forcing file:  ", trim(forcing_file)
end if
print *, "bc_init file:  ", trim(bc_init_file)
print *, "absolute, relative residue, iter: ", saved_strategy_strf_r1, &
&       saved_strategy_strf_r1, saved_max_iter_strf, alpha_strf
print *, "--------------------------"
