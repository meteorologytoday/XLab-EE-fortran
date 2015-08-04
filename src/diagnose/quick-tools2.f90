subroutine d_dz_B2A(from_dat, to_dat)
! Notice that the top and bottom will not
! be altered since grid points' arrangement is not matched
implicit none
real(4) :: from_dat(nr-1,nz-1), to_dat(nr-1,nz)
integer :: i,j,m,n
do i = 1, nr-1
    do j = 2, nz-2
        to_dat(i,j) = (from_dat(i,j) - from_dat(i,j-1)) &
            &       / ((za(j+1) - za(j-1))/2.0)
    end do
end do
end subroutine


subroutine d_dz_O2C(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr,nz-1)
integer :: i,j,m,n
do i = 1, nr
    do j = 1, nz-1
        to_dat(i,j) = (from_dat(i,j+1) - from_dat(i,j)) / (za(j+1) - za(j))
    end do
end do
end subroutine

subroutine d_dr_B2B(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr-1,nz-1), to_dat(nr-1,nz-1)
integer :: i,j,m,n
do i = 1, nr-1
    do j = 1, nz-1
        if(i == 1) then
            m = 0; n = 1
        else if(i == nr-1) then
            m = -1; n = 0
        else
            m = -1; n = 1
        end if
        to_dat(i,j) = (from_dat(i+m,j) - from_dat(i+n,j)) / (ra(i+m) - ra(i+n)) 
    end do
end do
end subroutine

subroutine d_dr_B2C(from_dat, to_dat)
! Notice that the first and last columns will not
! be altered since grid points' arrangement is not matched
implicit none
real(4) :: from_dat(nr-1,nz-1), to_dat(nr,nz-1)
integer :: i,j,m,n
do i = 2, nr-1
    do j = 1, nz-1
        to_dat(i,j) = (from_dat(i,j) - from_dat(i-1,j)) &
            &       / ((ra(i+1) - ra(i-1))/2.0)
    end do
end do
end subroutine

subroutine d_dr_O2A(from_dat, to_dat)
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr-1,nz)
integer :: i,j
do i = 1, nr-1
    do j = 1, nz
        to_dat(i,j) = (from_dat(i+1,j) - from_dat(i,j)) / (ra(i+1) - ra(i)) 
    end do
end do
end subroutine


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