! ###########################################
! # Solar9 multiobjective problem
! ###########################################
subroutine setdim(n,m,q)
    implicit none
    integer :: n,m,q

    n = 22
    m = 17
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
    real*8 :: outputs(19)
    real*8 :: l(n), u(n)

    l =  (/1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 793.0, 1.0, 1.0, 0.01, &
        0.01, 495.0, 0.01, 0.0050, 0.006, 0.007, 0.5, 0.0050, 0.006, 0.15/)
    l = l + 1e-7
    u = (/40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 995.0, 50.0, 30.0, &
        5.00, 5.00, 650.0, 5.00, 0.1000, 0.100, 0.200, 10.0, 0.1000, 0.100, 0.40/)
    u = u - 1e-7

    ! Open temporary file to put inputs of the blackbox
    ! To avoid being outside bounds (most of the time, with small precision)
    ! force to be inside bounds
    open(1, file='solar9_tmp_x.txt', status='replace')
    do i=1,5
        write(1,"(es16.9)") min(max(x(i), l(i)), u(i))
    enddo
    ! Maximum number of heliostats
    write(1,*) 1000
    do i=6,14
        write(1,"(es16.9)") min(max(x(i), l(i)), u(i))
    enddo
    ! Receiver number of tubes
    write(1,*) 500
    do i=15,22
        write(1,"(es16.9)") min(max(x(i), l(i)), u(i))
    enddo
    ! Exchanger number of baffles
    write(1,*) 3
    ! Exchanger number of tubes
    write(1,*) 12000
    ! Exchanger number of shells
    write(1,*) 1
    ! Exchanger number of passes per shell
    write(1,*) 2
    ! Type of turbine
    write(1,*) 2

    ! Call blackbox and save inputs into a temporary file
    call execute_command_line('./solar_bb.exe 9 solar9_tmp_x.txt > solar9_tmp_outputs.txt', wait=.true.)

    ! Close and delete inputs tmp file
    close(1)

    ! Read file outputs
    open(2, file='solar9_tmp_outputs.txt', status='old')
    read(2,*) outputs
    close(2)

    write(*,*) outputs

    f = outputs(1:2)

    return
end subroutine functs

subroutine fconstriq(n,m,x,ciq)
    implicit none
    integer :: n,m
    real*8 :: x(n), ciq(m)
    real*8 :: outputs(19)

    ! Read file outputs
    open(3, file='solar9_tmp_outputs.txt', status='old')
    read(3,*) outputs
    close(3)

    ! Delete it
    call execute_command_line('rm solar9_tmp_outputs.txt', wait=.true.)

    ciq = outputs(3:19)

    return
end subroutine fconstriq

subroutine setbounds(n,l,u)
    implicit none
    integer :: n
    real*8 :: l(n), u(n)

    l =  (/1.0, 1.0, 20.0, 1.0, 1.0, 1.0, 0.0, 1.0, 793.0, 1.0, 1.0, 0.01, &
        0.01, 495.0, 0.01, 0.0050, 0.006, 0.007, 0.5, 0.0050, 0.006, 0.15/)
    u = (/40.0, 40.0, 250.0, 30.0, 30.0, 89.0, 20.0, 20.0, 995.0, 50.0, 30.0, &
        5.00, 5.00, 650.0, 5.00, 0.1000, 0.100, 0.200, 10.0, 0.1000, 0.100, 0.40/)
    return
end subroutine setbounds

subroutine startp(n,x)
    implicit none
    integer :: n
    real*8 :: x(n)

    x = (/9.0, 9.0, 150.0, 6.0, 8.0, 45.0, 0.5, 5.0, 900.0, 9.0, 9.0, 0.30, 0.20, 560.0, &
        0.30, 0.0165, 0.018, 0.017, 10.0, 0.0155, 0.016, 0.20/)

    return
end subroutine startp
