program main
use elliptic_tools
use field_tools
use constants
use read_input_tools
use message_tools
implicit none

include "variables.f90"

print *, "Balanced Vortex Prediction Program"

debug_mode = 0
inquire(file="./debug_mode_1", exist=file_exists);
if(file_exists .eqv. .true.) then
    debug_mode = 1
end if
inquire(file="./debug_mode_2", exist=file_exists);
if(file_exists .eqv. .true.) then
    debug_mode = 2
end if

include "read-input.f90"
print *, "Read input complete."
include "initialize-variables.f90"
print *, "Initialization complete."
call cpu_time(time_beg)
include "evolution.f90"
call cpu_time(time_end)
print *, "Diagnose complete."

include "write-output.f90"

contains

include "quick-tools1.f90"
include "quick-tools2.f90"

end program main
