module message_tools
implicit none

contains

subroutine error_msg(err_type, err_code, err_msg)
implicit none
character(*) :: err_type, err_msg
integer :: err_code

write (*, "(A,A,A,I3,A,A)") "ERROR: [", trim(err_type), ",", err_code, "] : " , trim(err_msg)
end subroutine


subroutine system_msg(sys_type, sys_msg)
implicit none
character(*) :: sys_type, sys_msg
write (*, "(A,A,A,A)") "[", trim(sys_type), "] : " , trim(sys_msg)

end subroutine

end module message_tools
