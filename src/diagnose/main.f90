program main
use elliptic_tools
use field_tools
use constants
use read_input_tools
use message_tools
implicit none

include "variables.f90"

inquire(file="./debug_mode", exist=debug_mode);

include "read-input.f90"
include "initialize-variables.f90"

call cpu_time(time_beg)
include "diagnose.f90"
call cpu_time(time_end)

include "write-output.f90"

contains

include "quick-tools.f90"

end program main
