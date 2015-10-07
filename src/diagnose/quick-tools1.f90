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

subroutine cal_uw(from_rpsi, to_u, to_w)
implicit none
real(4) :: from_rpsi(nr,nz), to_w(nr-1,nz), to_u(nr,nz-1)
real(4) :: r
integer :: i,j

call d_rcuvdr_O2A(from_rpsi, to_w)
call d_dz_O2C(from_rpsi, to_u); to_u = -to_u
! notice gradient of rpsi gets momentum flux, so we need to divide it by rho.
do i=1,nr-1
    do j=1,nz
        to_w(i,j) = to_w(i,j)/rho(j)
    end do
end do

do i=1,nr
    do j=1,nz-1
        r = ra(i)
        if(r /= 0) then 
            to_u(i,j) = to_u(i,j)/(rcuva(i)*(rho(j)+rho(j+1))/2.0)
        else
            to_u(i,j) = 0.0
        end if
    end do
end do

end subroutine