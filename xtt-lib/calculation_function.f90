module calculation_function

contains

subroutine relativeMM(theta_B, dtheta_dz_A, dtheta_dr_C)
implicit none
real(4) :: theta_B(nr-1,nz-1), dtheta_dz_A(nr-1,nz), dtheta_dr_C(nr,nz-1)

real(4) :: dist
integer :: i,j
theta_B = theta0;

m2(1,1) = 0

do j = 0

do i=1, nr-1
    do j=1, nz-1
            m2(i,j) = m2(i-1,j) + (rcuva(i)**3.0) * rhoC_C(i,j) * (ra(i+1) - ra(i-1)) / 2.0
    end do
end do

do i=2,nr-1
    dist = (ra(i+1) - ra(i-1)) / 2.0
    theta_B(i,1) = theta_B(i-1,1) + dist * dtheta_dr_C(i,1)
end do

do i=1,nr-1
    do j=2,nz-1
        dist = (za(j+1) - za(j-1)) / 2.0
        theta_B(i,j) = theta_B(i,j-1) + dist * dtheta_dz_A(i,j)
    end do
end do

end subroutine




subroutine relativeTheta(theta_B, dtheta_dz_A, dtheta_dr_C)
implicit none
real(4) :: theta_B(nr-1,nz-1), dtheta_dz_A(nr-1,nz), dtheta_dr_C(nr,nz-1)

real(4) :: dist
integer :: i,j
theta_B = theta0;
do i=2,nr-1
    dist = (ra(i+1) - ra(i-1)) / 2.0
    theta_B(i,1) = theta_B(i-1,1) + dist * dtheta_dr_C(i,1)
end do

do i=1,nr-1
    do j=2,nz-1
        dist = (za(j+1) - za(j-1)) / 2.0
        theta_B(i,j) = theta_B(i,j-1) + dist * dtheta_dz_A(i,j)
    end do
end do

end subroutine


subroutine rpsiToUW(from_rpsi, to_u, to_w)
implicit none
real(4) :: from_rpsi(nr,nz), to_w(nr-1,nz), to_u(nr,nz-1)
real(4) :: r
integer :: i,j

call d_rdr_O2A(from_rpsi, to_w)
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
            to_u(i,j) = to_u(i,j)/(r*(rho(j)+rho(j+1))/2.0)
        else
            to_u(i,j) = 0.0
        end if
    end do
end do

end subroutine

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


subroutine d_rdr_O2A(from_dat, to_dat)
! notice that all the points with r=0 become inifinite
implicit none
real(4) :: from_dat(nr,nz), to_dat(nr-1,nz)
real(4) :: r
integer :: i,j

call d_dr_O2A(from_dat, to_dat)

do i = 1, nr-1
    do j = 1, nz
        r = (ra(i)+ra(i+1))/2.0
        to_dat(i,j) = to_dat(i,j) / r
    end do
end do

end subroutine

real(4) function cal_sum_Q(Q)
implicit none
real(4) :: Q(nr-1, nz-1)
real(4) :: r, dr, dz, rho_
integer :: i,j

cal_sum_Q = 0.0;
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_Q = cal_sum_Q + Q(i,j) * rho_ * r * dr * dz
    end do
end do

end function

real(4) function cal_sum_Qeta(Q, eta)
implicit none
real(4) :: Q(nr-1, nz-1), eta(nr-1, nz)
real(4) :: r, dr, dz, rho_
integer :: i,j

cal_sum_Qeta = 0.0
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_Qeta = cal_sum_Qeta + ((eta(i,j) + eta(i,j+1)) /2.0) * Q(i,j) * rho_ * r * dr * dz
    end do
end do

end function

real(4) function cal_sum_wtheta(wtheta_B)
implicit none
real(4) :: wtheta_B(nr-1, nz-1)
real(4) :: r, dr, dz, rho_
integer :: i,j

cal_sum_wtheta = 0.0
do i = 1,nr-1
    do j = 1, nz-1
        r = (ra(i)+ra(i+1))/2.0
        dr = ra(i+1) - ra(i)
        dz = za(j+1) - za(j)
        rho_ = (rho(j+1) + rho(j))/2.0

        cal_sum_wtheta = cal_sum_wtheta + wtheta_B(i,j) * rho_ * r * dr * dz
    end do
end do

end function



subroutine cal_wtheta(w_A, theta_B, wtheta_B)
implicit none
real(4) :: w_A(nr-1, nz), theta_B(nr-1, nz-1), wtheta_B(nr-1, nz-1)
integer :: i,j

do i = 1,nr-1
    do j = 1, nz-1
        wtheta_B(i,j) = ((w_A(i,j) + w_A(i,j+1)) /2.0) * theta_B(i,j)
    end do
end do
end subroutine

subroutine cal_eta(rchi, eta)
implicit none
real(4) :: rchi(nr,nz), eta(nr-1,nz)
integer :: i,j

call d_rdr_O2A(rchi, eta)

do i=1,nr-1
    do j=1,nz
        eta(i,j) = eta(i,j) * g0 / (rho(j) * Cp * exner(j) * theta0)
    end do
end do
end subroutine

subroutine cal_exchange_conversion(rpsi, rchi, rhoC, bndconv, sum_bndconv)
implicit none
real(4) :: rpsi(nr,nz), rchi(nr, nz), rhoC(nr,nz), bndconv(nr-1,2), sum_bndconv
integer :: i, r, dr, dz

sum_bndconv = 0.0
dz = za(2) - za(1)
dr = ra(2) - ra(1)
do i=1, nr-1
    r = (ra(i)+ra(i+1))/2.0
    
    ! Bottom
    bndconv(i, 1) = ((rhoC(i, 1) + rhoC(i+1, 1))/(2.0*rho(1))) * (&
    &    ( (rpsi(i,1) + rpsi(i+1,1))/2.0 ) * &
    &     ( (rchi(i,2) + rchi(i+1,2) - rchi(i,1) - rchi(i+1,1))/(2.0*dz)) - &
    &    ( (rchi(i,1) + rchi(i+1,1))/2.0 ) * &
    &     ( (rpsi(i,2) + rpsi(i+1,2) - rpsi(i,1) - rpsi(i+1,1))/(2.0*dz)) ) &
    &    / r**2.0

    ! Top
    bndconv(i, 2) = ((rhoC(i, nz) + rhoC(i+1, nz))/(2.0*rho(nz))) * (&
    &    ( (rpsi(i,nz) + rpsi(i+1,nz))/2.0 ) * &
    &     ( (rchi(i,nz) + rchi(i+1,nz) - rchi(i,nz-1) - rchi(i+1,nz-1))/(2.0*dz)) - &
    &    ( (rchi(i,nz) + rchi(i+1,nz))/2.0 ) * &
    &     ( (rpsi(i,nz) + rpsi(i+1,nz) - rpsi(i,nz-1) - rpsi(i+1,nz-1))/(2.0*dz)) ) &
    &    / r**2.0

    sum_bndconv = sum_bndconv - (bndconv(i,2) - bndconv(i,1)) * r * dr
end do



end module calculation_function
