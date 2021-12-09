! ###########################################
! # Solar8 multiobjective problem
! ###########################################

subroutine setdim(n,m,q)
    implicit none
    integer :: n,m,q

    n = 11
    m = 9
    q = 2
    return
end subroutine setdim

! NB: The DFMO algorithm always calls functs then
! fconstriq after
subroutine functs(n,x,q,f)
    implicit none
    integer :: n, q
    integer :: i
    real*8 :: x(n), f(q)
    real*8 :: outputs(11)
    real*8 :: l(n), u(n)

    l = (/1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1e-1, 5e-3, 6e-3/)
    l = l + 1e-7
    u = (/40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 5.00, 1e-1, 1e-1/)
    u = u - 1e-7

    ! Open temporary file to put inputs of the blackbox
    ! To avoid being outside bounds (most of the time, with small precision)
    ! force to be inside bounds
    open(1, file='solar8_tmp_x.txt', status='replace')
    do i=1,5
        write(1,"(es16.9)") min(max(x(i), l(i)), u(i))
    enddo
    ! Maximum number of heliostats
    write(1,*) 2650
    do i=6,8
        write(1,"(es16.9)") min(max(x(i), l(i)), u(i))
    enddo
    ! Receiver number of tubes
    write(1,*) 36
    do i=9,11
        write(1,"(es16.9)") min(max(x(i), l(i)), u(i))
    enddo

    ! Call blackbox and save inputs into a temporary file
    call execute_command_line('./solar_bb.exe 8 solar8_tmp_x.txt > solar8_tmp_outputs.txt', wait=.true.)

    ! Close and delete inputs tmp file
    close(1)

    ! Read file outputs
    open(2, file='solar8_tmp_outputs.txt', status='old')
    read(2,*) outputs
    close(2)

    f = outputs(1:2)

    return
end subroutine functs

subroutine fconstriq(n,m,x,ciq)
    implicit none
    integer :: n,m
    real*8 :: x(n), ciq(m)
    real*8 :: outputs(11)

    ! Read file outputs
    open(3, file='solar8_tmp_outputs.txt', status='old')
    read(3,*) outputs
    close(3)

    ! Delete it
    call execute_command_line('rm solar8_tmp_outputs.txt', wait=.true.)

    ciq = outputs(3:11)

    return
end subroutine fconstriq

subroutine setbounds(n,l,u)
    implicit none
    integer :: n
    real*8 :: l(n), u(n)

    l = (/1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.01, 0.005, 0.0060/)
    u = (/40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 5.00, 0.100, 0.1000/)
    return
end subroutine setbounds

subroutine startp(n,x)
    implicit none
    integer :: n
    real*8 :: x(n)

    x = (/11.0, 11.0, 200.0, 10.0, 10.0, 89.0, 0.5, 8.0, 0.30, 0.020, 0.0216/)

    return
end subroutine startp
