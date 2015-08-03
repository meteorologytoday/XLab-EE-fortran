module elliptic_tools
implicit none
integer, parameter :: err_over_max_iteration = ISHFT(1,0), &
    &                 err_explode = ISHFT(1,1)

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


subroutine solve_elliptic(max_iter, check_step, converge_time, lost_rate, &
&                         strategy_r1, strategy_r2, alpha, dat, coe, f, &
&                         workspace, nx, ny, err, debug)
implicit none
real(4), intent(inout), target :: dat(nx,ny), workspace(nx, ny)
real(4), intent(inout) :: strategy_r1, strategy_r2
real(4), intent(in)    :: coe(9, nx, ny), f(nx, ny), alpha
integer, intent(in)    :: nx, ny, check_step, converge_time, lost_rate, debug
integer, intent(inout) :: max_iter, err
integer                :: lose_chance_cnt, check_step_use, converge_cnt, converge_time_use, lost_rate_use

integer :: i,j,cnt, tmp_err
real(4) :: err_now, err_before, ratio 

real(4), pointer :: fr_dat(:,:), to_dat(:,:), tmp_ptr(:,:)
logical :: stop_iteration, flag ! after every [check_step] iterations this flag will be turned on.
logical :: check_abs_err, check_rel_err


if (strategy_r1 > 0) then
    check_abs_err = .true.
else
    check_abs_err = .false.
    strategy_r1 = HUGE(strategy_r1)
end if

if (strategy_r2 > 0) then
    check_rel_err = .true.
else
    check_rel_err = .false.
    strategy_r2 = HUGE(strategy_r2)
end if

if((check_abs_err .eqv. .false.) .and. (check_rel_err .eqv. .false.)) then
    print *, "ERROR: [check_abs_err] and [check_rel_err] cannot both be non-positive."
    stop
end if

check_step_use = 100
if(check_step > 0) then
    check_step_use = check_step
end if

converge_time_use = 10
if(converge_time > 0) then
    converge_time_use = converge_time
end if

lost_rate_use = 5
if(lost_rate > 0) then
    lost_rate_use = lost_rate
end if

if(debug == 1 .or. debug == 2) then
    print *, "----- Solve Elliptic Inputs -----"
    print *, "  max_iter       : ", max_iter
    print *, "  strategy_r1    : ", strategy_r1
    print *, "  strategy_r2    : ", strategy_r2
    print *, "  alpha          : ", alpha
    print *, "  (nx, ny)       : (", nx, ", ", ny, ")"
    print *, "  alpha          : ", alpha
    print *, "  debug          : ", debug
    print *, "  check step     : ", check_step_use
    print *, "  converge time : ", converge_time_use
    print *, "---------------------------------"
end if

converge_cnt = 0
lose_chance_cnt = 0

err_before=Huge(err_before)
err = 0

workspace(1, 1:ny) = dat(1,1:ny)
workspace(nx,1:ny) = dat(nx,1:ny)
workspace(1:nx, 1) = dat(1:nx,1)
workspace(1:nx,ny) = dat(1:nx,ny)

workspace = dat

fr_dat => workspace
to_dat => dat

stop_iteration = .false.
do cnt=1, max_iter

    if(mod(cnt, check_step_use) == 0) then
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
        err_now = 0
        do i = 2, nx-1
            do j = 2, ny-1
                err_now = err_now + to_dat(i,j)**2.0
            end do
        end do
        err_now = sqrt(err_now/((nx-2)*(ny-2)))

        ratio = (err_before - err_now) / err_before;
        if(debug == 2) then
            write(*,'(A,I8,A,ES12.3E2,A,ES12.3E2)') "Iter: ",cnt, ", err_now: ", err_now, ", ratio: ", ratio
        end if
        ratio = abs(ratio)
        if(err_before == 0) then
            stop_iteration = .true.
            if(debug == 2) then
                print *, "Error = 0, hardly to see this!"
            end if
        else if((err_now < strategy_r1) .and. (ratio < strategy_r2)) then
            converge_cnt = converge_cnt + 1
            lose_chance_cnt = 0
            if(debug == 2) then
                print *, "converge_cnt: ", converge_cnt
            end if
            if(converge_cnt >= converge_time_use)  then
                stop_iteration = .true.
            end if
        else
            if(converge_cnt > 0) then
                lose_chance_cnt = lose_chance_cnt + 1
            
                if(lose_chance_cnt >= lost_rate_use) then
                    converge_cnt = converge_cnt - 1
                    lose_chance_cnt = 0
                    if(debug == 2) then
                        print *, "Lose one count! converge_cnt now is: ", converge_cnt
                    end if
                end if
            end if
        end if
        err_before = err_now
    end if

    do i = 2,nx-1
        do j = 2,ny-1
            to_dat(i,j) = fr_dat(i,j) + alpha * to_dat(i,j) / (- coe(5,i,j))
        end do
    end do        
    
    if(cnt == max_iter) then
        stop_iteration = .true.
        err = IOR(err, err_over_max_iteration)
        if(debug == 2) then
            print *, "Max iteration reached. Exit iteration."
        end if
    end if
    if(stop_iteration .eqv. .true.) then
        if(debug == 2) then
            print *, "iter : ", cnt, ", err_avg = ", err_now
        end if
        max_iter = cnt; strategy_r1 = err_now; strategy_r2 = ratio
        call judge_error(err)
        exit
    end if
end do

if(associated(to_dat, dat) .eqv. .false.) then
    if(debug == 2) then
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

subroutine judge_error(err)
implicit none
integer, intent(in) :: err
logical :: known_err

known_err = .false.

if(err == 0) then
    print *, "Elliptic Tools: Iteration success."
    known_err = .true.
end if

if(IAND(err, err_over_max_iteration) /= 0) then
    print *, "Elliptic Tools: [Error] Max iteration reached."
    known_err = .true.
end if

if(IAND(err, err_explode) /= 0) then
    print *, "Elliptic Tools: [Error] Iteration explodes."
    known_err = .true.
end if

if(known_err .eqv. .false.) then
    print *, "Elliptic Tools: Unknown error code ", err
end if
end subroutine

end module elliptic_tools
