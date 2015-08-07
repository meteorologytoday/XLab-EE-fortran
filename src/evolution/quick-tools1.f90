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

subroutine cyl_ABC_to_m2theta(rhoA_sA, rhoB_B, rhoC_sC, m2, theta)
implicit none
real(4) :: rhoA_sA(nr-1, nz-2), rhoB_B(nr-1, nz-1), rhoC_sC(nr-2, nz-1), &
    &      m2_B(nr-1, nz-1), theta_B(nr-1, nz-1)

integer :: i, j

! invert m2 
do i=1, nr-1
    do j=1, nz-1
        if(i == 1 .and. j == 1) then
            m2(i,j) = 0.0
        end if
        else
    end do
end do


end subroutine


