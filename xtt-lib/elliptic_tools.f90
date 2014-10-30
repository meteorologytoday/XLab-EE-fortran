module elliptic_tools

contains

subroutine cal_coe(a, b, c, workspace, dx, dy, nx, ny, err)
implicit none
real(4), intent(in) :: a(nx-1, ny-2), b(nx-1, ny-1), c(nx-2, ny-1), &
&                      dx, dy
real(4), intent(inout) :: workspace(9, nx, ny)
integer, intent(in)    :: nx, ny
integer, intent(inout) :: err

integer :: i,j
real(4) :: PP, QQ, PQ4, Ap, Am, Cp, Cm, BXp, BXm, BYp, BYm

! 1 2 3
! 4 5 6
! 7 8 9
!
! +-A-+-A-+-A-+
! | B C B C B |
! +-A-+-A-+-A-+
! | B C B C B |
! +-A-+-A-+-A-+

PP = dx**2
QQ = dy**2
PQ4 = 4*dx*dy

err = 1

do i=2, nx-1
    do j=2, ny-1 
        Ap  = a(i  , j-1) / PP
        Am  = a(i-1, j-1) / PP
        Cp  = c(i-1, j  ) / QQ
        Cm  = c(i-1, j-1) / QQ
        BXp = (b(i  ,j  ) + b(i  ,j-1))/(2.0 * PQ4)
        BXm = (b(i-1,j  ) + b(i-1,j-1))/(2.0 * PQ4)
        BYp = (b(i-1,j  ) + b(i  ,j  ))/(2.0 * PQ4)
        BYm = (b(i-1,j-1) + b(i  ,j-1))/(2.0 * PQ4)

        workspace(1, i, j) = - (BXm + BYp)
        workspace(2, i, j) = Cp + (BXp - BXm)
        workspace(3, i, j) = BXp + BYp
        workspace(4, i, j) = Am - (BYp - BYm)
        workspace(5, i, j) = - (Am + Ap + Cm + Cp)
        workspace(6, i, j) = Ap + (BYp - BYm)
        workspace(7, i, j) = BXm + BYm
        workspace(8, i, j) = Cm - (BXp - BXm)
        workspace(9, i, j) = - (BXp + BYm)
    end do
end do

err = 0

end subroutine



subroutine do_elliptic(psi, coe, outdat, nx, ny, err)
implicit none
real(4), intent(in)    :: psi(nx, ny), coe(9, nx,ny)
real(4), intent(inout) :: outdat(nx,ny)
integer, intent(in)    :: nx, ny
integer, intent(inout) :: err

integer :: i,j,cnt

err = 1

    do i=2, nx-1
        do j=2, ny-1 
            outdat(i,j)  =    coe(1, i, j) * psi(i-1, j+1)  &
&                           + coe(2, i, j) * psi(i  , j+1)  &
&                           + coe(3, i, j) * psi(i+1, j+1)  &
&                           + coe(4, i, j) * psi(i-1, j  )  &
&                           + coe(5, i, j) * psi(i  , j  )  &
&                           + coe(6, i, j) * psi(i+1, j  )  &
&                           + coe(7, i, j) * psi(i-1, j-1)  &
&                           + coe(8, i, j) * psi(i  , j-1)  &
&                           + coe(9, i, j) * psi(i+1, j-1)

        end do
end do

end subroutine


subroutine solve_elliptic(max_iter, strategy, strategy_r, alpha, dat, coe, f, &
&                         workspace, nx, ny, err, debug)
!
! max_iter{Integer}
!     Maxima iterations. The relaxation stops If pmax_iter]
!     is reached but the convergence condition is not met.
!
! strategy{Integer}
! strategy_r{Real(4)}
!     [strategy] specifies which convergence condition is used. 
!     
!     strategy "1":
!     convergence by comparing the error shown below
!     => Err := (N^{-1} \sum_{i,j} ( f_{cal}(i,j) - f(i,j))^2 )^{1/2}
!        ( N = (nx-2) * (ny-2), 2 <= i <= nx-1, 2 <= j <= ny-1 )
!     the condition is met when [Err] < [strategy_r].
!
!     When solve_elliptic finishes, [strategy] is set to the iteration
!     times it takes, and [strategy_r] is set to the [Err] the solution has.
!
!     strategy "2":
!     This strategy uses the same error estimation as strategy "1", while
!     its breaking point is when the "Variation" of error relative to the
!     previous iteration is smaller than [strategy_r].
!
!     strategy "3":
!     Error is the maximaum
!
!
! dat{real(:,:)}
!     [dat] is the field wanted to be reverted to.
!     Boundary condition can be specified on this field.
!     Initial guess can also be specified here.
!
! f{real(:,:)}
!     [f] is the right hand side.
!
! err{Integer}
!     [err] is the error code.
!
! alpha{Real}
!     [alpha] is for overrelaxation.
!     phi^{i+1}_{j} = phi^{i} + alpha * (residue)
!
!     notice that [alpha] = 1 means non-oveerrelaxation
! 
! debug{Integer}
!     Output debug message when [debug] == 1, blank otherwise.

implicit none
real(4), intent(inout), target :: dat(nx,ny), workspace(nx, ny)
real(4), intent(inout) :: strategy_r
real(4), intent(in)    :: coe(9, nx, ny), f(nx, ny), alpha
integer, intent(in)    :: max_iter, nx, ny
integer, intent(inout) :: strategy, err
integer                :: debug, converge_cnt

integer :: i,j,cnt, tmp_err
integer :: check_step
real(4) :: tmp_r, err_now, err_before, stg3_max_err, ratio 

real(4), pointer :: fr_dat(:,:), to_dat(:,:), tmp_ptr(:,:)
logical :: flag ! after every [check_step] iterations this flag will be turned on. 

print *, "alpha: ", alpha
print *, "debug: ", debug
converge_cnt = 0

err_before=Huge(err_before)
check_step = 100
err = 1

workspace(1, 1:ny) = dat(1,1:ny)
workspace(nx,1:ny) = dat(nx,1:ny)
workspace(1:nx, 1) = dat(1:nx,1)
workspace(1:nx,ny) = dat(1:nx,ny)

workspace = dat

fr_dat => workspace
to_dat => dat

do cnt=1, max_iter

    if(mod(cnt,check_step) == 0) then
        flag = .true.
    else
        flag = .false.
    end if

    tmp_ptr => fr_dat
    fr_dat => to_dat
    to_dat => tmp_ptr

    call do_elliptic(fr_dat, coe, to_dat, nx, ny, tmp_err)
    to_dat(2:nx-1, 2:ny-1) = to_dat(2:nx-1,2:ny-1) - f(2:nx-1,2:ny-1)

    if(flag .eqv. .true.) then
        if(strategy == 1 .or. strategy == 2) then
            err_now = 0
            do i = 2, nx-1
                do j = 2, ny-1
                    err_now = err_now + to_dat(i,j)**2.0
                end do
            end do
            err_now = sqrt(err_now/((nx-2)*(ny-2)))
        else if(strategy == 3 .or. strategy == 4) then
            err_now = maxval(abs(to_dat))
        end if
    end if
    if(flag .eqv. .true.) then
        ratio = (err_before - err_now) / err_before;
        if(debug == 1) then
            print *, "Iter: ",cnt, "; err_now: ", err_now, "; ratio: ", ratio
        end if
    end if
    
    do i = 2,nx-1
        do j = 2,ny-1
            if((isnan(coe(5,i,j)) .eqv. .true.)) then
                print *, "coe is nan at (",i,",",j,")"
                stop
            end if
            if(isnan(fr_dat(i,j)) .eqv. .true.) then
                print *, "fr_dat is nan at (",i,",",j,")"
                stop
            end if
            if(isnan(to_dat(i,j)) .eqv. .true.) then
                print *, "to_dat is nan at (",i,",",j,")"
                stop
            end if

            to_dat(i,j) = fr_dat(i,j) + alpha * to_dat(i,j) / (- coe(5,i,j))
        end do
    end do

    if(flag .eqv. .true.) then
        if(strategy == 1 .or. strategy == 3) then
            if(err_now < strategy_r) then
                err = 0
            end if
        else if(strategy == 2 .or. strategy == 4) then
            if(err_before == 0) then
                err = 0
                print *, "Error = 0, hardly to see this!"
            else if(abs(ratio) < strategy_r) then
                converge_cnt = converge_cnt + 1
                print *, "converge_cnt: ", converge_cnt
                if(converge_cnt >= 10)  then
                    err = 0
                end if
            end if
            err_before = err_now
        end if

        if(err == 0 .or. cnt == max_iter) then
            if(debug == 1) then
                print *, "iteration : ", cnt, ", err_avg = ", err_now, "fr_dat: "
            end if
            strategy = cnt; strategy_r = err_now
            exit
        end if
    end if
end do
if(associated(to_dat, dat) .eqv. .false.) then
    if(debug == 1) then
        print *,"to_dat is not associated with dat"
    end if
    dat(:,:) = workspace(:,:)
end if


end subroutine

subroutine solve_elliptic_AC(max_iter, max_r, psi, a, c, f, workspace, dx, dy, nx, ny, err)
implicit none
real(4), intent(inout), target :: psi(nx,ny), workspace(nx, ny)
real(4), intent(in)    :: max_r, dx, dy
real(4), intent(in)    :: a(nx-1, ny-2), c(nx-2, ny-1), f(nx, ny)
integer, intent(in)    :: max_iter, nx, ny
integer, intent(inout) :: err

integer :: i,j,cnt
real(4) :: P, Q, tmp_r, err_now, ax1, ax2, cy1, cy2

real(4), pointer :: fr_dat(:,:), to_dat(:,:), tmp_ptr(:,:)


P = dx**2
Q = dy**2

err = 1

workspace(1, 1:ny) = psi(1,1:ny)
workspace(nx,1:ny) = psi(nx,1:ny)
workspace(1:nx, 1) = psi(1:nx,1)
workspace(1:nx,ny) = psi(1:nx,ny)


fr_dat => workspace
to_dat => psi 

do cnt=1, max_iter

    tmp_ptr => fr_dat
    fr_dat => to_dat
    to_dat => tmp_ptr

    err_now = 0
    do i=2, nx-1
        do j=2, ny-1 
            tmp_r = ( fr_dat(i+1, j)          &
&                   + fr_dat(i-1, j)          &
&                   + fr_dat(i, j+1)          &
&                   + fr_dat(i, j-1)          &
&                   - 4*fr_dat(i,j)) / dx**2  &
&                   - f(i,j)
            err_now = err_now + tmp_r ** 2
            
            to_dat(i,j) = fr_dat(i,j) + tmp_r / (ax1 + ax2 + cy1 + cy2)
        end do
    end do
    err_now = sqrt(err_now / ((nx-2)*(ny-2)))
    print *, "iteration : ", cnt, ", err_avg = ", err_now, "fr_dat: ",   &
&       associated(fr_dat, psi), "to_dat: ", associated(to_dat, workspace)


    if(err_now < err) then
       err = 0
        exit
    end if
    
end do
print *, "out"
if(associated(to_dat, psi) .eqv. .false.) then
    print *,"to_dat is not associated with psi"
    psi(:,:) = workspace(:,:)
end if
end subroutine

end module elliptic_tools
