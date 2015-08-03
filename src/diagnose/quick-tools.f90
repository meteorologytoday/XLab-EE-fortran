subroutine cal_eta(rchi, eta)
implicit none
real(4) :: rchi(nr,nz), eta(nr-1,nz)
integer :: i,j

call d_rcuvdr_O2A(rchi, eta)

do i=1,nr-1
    do j=1,nz
        eta(i,j) = eta(i,j) * g0 / (rho(j) * Cp * exner(j) * theta0)
    end do
end do
end subroutine

