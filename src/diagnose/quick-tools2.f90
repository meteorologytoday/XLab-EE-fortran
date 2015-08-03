subroutine d_rcuvdr_O2A(from_dat, to_dat)
! notice that all the points with r=0 become inifinite
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr-1,nz)
integer :: i,j

call d_dr_O2A(from_dat, to_dat)

do i = 1, nr-1
    do j = 1, nz
        to_dat(i,j) = to_dat(i,j) / ((rcuva(i) + rcuva(i+1)) / 2.0)
    end do
end do

end subroutine

