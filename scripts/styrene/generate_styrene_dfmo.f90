! ###########################################
! # STYRENE multiobjective problem
! ###########################################

subroutine setdim(n,m,q)
    implicit none
    integer :: n,m,q

    n = 8
    m = 9
    q = 3
    return
end subroutine setdim

! NB: The DFMO algorithm always calls functs then
! fconstriq after
subroutine functs(n,x,q,f)
    implicit none
    integer :: n, q
    integer :: i
    real*8 :: x(n), f(q)
    real*8 :: outputs(12)
    real*8 :: l(n), u(n)
    integer :: info

    l = 0
    l = l + 1e-7
    u = 100
    u = u - 1e-7

    ! Open temporary file to put inputs of the blackbox
    ! To avoid being outside bounds (most of the time, with small precision)
    ! force to be inside bounds
    open(1, file='styrene_tmp_x.txt', status='replace')
    do i=1,8
        write(1,"(es16.9)") min(max(x(i), l(i)), u(i))
    enddo

    ! Call blackbox and save inputs into a temporary file
    call execute_command_line('./truth.exe styrene_tmp_x.txt > styrene_tmp_outputs.txt', wait=.true.)

    ! Close and delete inputs tmp file
    close(1)
    ! call execute_command_line('rm solar8_tmp_x.txt', wait=.true.)

    ! Detect if has failed
    call execute_command_line('grep -q ERROR styrene_tmp_outputs.txt', wait=.true., exitstat=info)
    if (info == 0) then
        f = 1e27 ! take a big value
    else
        ! Read file outputs
        open(2, file='styrene_tmp_outputs.txt', status='old')
        read(2,*) outputs
        close(2)

        write(*,*) outputs

        f(1) = outputs(12) ! net present value of the project (f1)
        f(2) = outputs(5) ! minimal purity of produced styrene (f2)
        f(3) = outputs(7) ! overall ethylbenzene conversion into styrene (f3)

    endif

    return
end subroutine functs

subroutine fconstriq(n,m,x,ciq)
    implicit none
    integer :: n,m
    real*8 :: x(n), ciq(m)
    real*8 :: outputs(12)
    integer :: info

    call execute_command_line('grep -q ERROR styrene_tmp_outputs.txt', wait=.true., exitstat=info)
    if (info == 0) then
        ciq = 1e27 ! big values
    else
        ! Read file outputs
        open(3, file='styrene_tmp_outputs.txt', status='old')
        read(3,*) outputs
        close(3)

        ciq(1:4) = outputs(1:4) ! boolean constraints: set to 1 if simulation fails, 0 otherwise.
        ciq(5) = outputs(6) ! minimal purity of produced benzene
        ciq(6:9) = outputs(8:11) ! Four constraints relating to payout time, cashflow, investment and annual costs
    endif

    ! Delete it
    call execute_command_line('rm styrene_tmp_outputs.txt', wait=.true.)

    return
end subroutine fconstriq

subroutine setbounds(n,l,u)
    implicit none
    integer :: n
    real*8 :: l(n), u(n)

    l = 0
    u = 100
    return
end subroutine setbounds

subroutine startp(n,x)
    implicit none
    integer :: n
    real*8 :: x(n)

    x = (/54.0, 66.0, 86.0, 8.0, 29.0, 51.0, 32.0, 15.0/)
    return
end subroutine startp
