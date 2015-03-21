module read_input_tools
contains

! read_input skips comments following the javascript tradition:
! // line comment
! Notice that it reads 256 characters most in a line
subroutine read_input(fd, line)
implicit none
integer :: fd
character(256)   :: line

character(256) :: buffer
logical :: loop_end
integer :: idx



loop_end = .false.
do while(loop_end .eqv. .false.)
    
    read (fd, '(a)') buffer

    ! Search line comment
    idx = index(buffer, '//')
    if(idx /= 0) then
        ! Truncate the line
        buffer = buffer(1:idx-1)
    end if
    buffer = trim(buffer)
    if(buffer /= '') then
        loop_end = .true.
    end if
    
end do
!print *, 'result: ', trim(buffer)
line = trim(buffer)

end subroutine


integer function split_line(line, output, delimiter)
implicit none
character(*), intent(inout)   :: line, output
character(*)                  :: delimiter


integer :: i


i = INDEX(line, delimiter)

if (i == 0) then ! No match and length is not zero
    output = line
    line = ''
    split_line = 1
else
    output = line(1:i-1)
    line = line(i+1:)
    split_line = 0
end if

end function

end module
    
