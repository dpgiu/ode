program test
    use precision
    use functions
    use solvers
    implicit none

    real(dp), parameter        :: pi = acos(-1.0_dp)
    real(dp)                   :: err, t, dt, x0, v0, Qmin, Qmax, dQ, AA, BB, xex
    real(dp), allocatable      :: u(:), u0(:), param(:)
    integer                    :: i, j, k, NQ, Nstep, Ncicli, num_args, funit            ! funit = file unit
    procedure(solver), pointer :: psolver
    character(50), allocatable :: args(:)   ! stringa
    character(15) :: solname
    ! logical :: boleano = .true. || .false.

    if (command_argument_count() < 10) then
        write(*,*) "test Nstep Ncicli method Qmin Qmax dQ A w x0 v0"
        stop
    endif

    num_args = command_argument_count()
    allocate(args(num_args))
    allocate(param(num_args-1))
    allocate(u(2))
    allocate(u0(2))

    do i = 1, num_args
        call get_command_argument(i, args(i))
        if (i == 3) then
            read(args(i),*) solname
        else
            read(args(i),*) param(i)
        end if
    end do

    Nstep  = int(param(1))
    Ncicli = int(param(2))
    Qmin   = param(4)
    Qmax   = param(5)
    dQ     = param(6)
    A      = param(7)
    w      = param(8)                            ! per mostrare i diagrammi di biforcazione in genere si sceglie il valore 2/3
    u0     = (/param(9), param(10)/)
    NQ     = int((Qmax - Qmin) / dQ)      
    t      = 0.0_dp
    dt     = real((2.0_dp * pi)/(w * Nstep), dp) 
    err    = 0.0_dp

    !arg = 'rk2  ' <- trailing white spaces -->'
    ! if (arg == 'rk2') then
    !  ...
    ! else if (trim(arg) = 'rk4') then
    ! ...
    ! end if
    ! select case è un costrutto equivalente

    open(newunit=funit, file="solution.dat")       ! pbm arcode c'è quando si usa con numeri grandi e si apre più volte file, 
    write(funit,*) "time  ", "theta  ", "speed"    ! così prende valore che non viene usato sicuramente da nessuno

    select case(trim(solname))
    case('rk2')
        do i = 1, Nstep * Ncicli
            write(funit,*) t, u0(1), u0(2)
            call rk2(pendolo, t, dt, u0, u)
            t = t+dt
            u0 = u
        end do
    case('rk4_notrans')
        do k = 0, NQ
            Q = Qmin + k * dQ
            do j = 1, Ncicli
                if (j > 3 * Ncicli / 4) then        ! cambiato da 2/3 a 3/4, aumentando iterazioni relative al transiente 
                    write(funit,*) Q, u0(2)         ! cambiato da t, u0(1), u0(2), per graficare la sezione Poincaré v vs Q
                end if
                do i = 1, Nstep
                    call rk4(pendolo, t, dt, u0, u)
                    t = t+dt
                    u0 = u
                end do
            end do
        end do
    case('rk4')
        do i = 1, Nstep * Ncicli
            write(funit,*) t, u0(1), u0(2)
            call rk4(pendolo, t, dt, u0, u)
            t = t+dt
            u0 = u
        end do
    case('dp54')
        do i = 1, Nstep * Ncicli
            call dopri54(pendolo, t, dt, u0, u, err)
            t = t+dt
            u0 = u
        end do
    case('dp87')
        do i = 1, Nstep * Ncicli
            call dopri87(pendolo, t, dt, u0, u, err)
            t = t+dt
            u0 = u
        end do
    case default
        write(funit,*) 'error: solver ',trim(solname),' not found'
        stop
    end select

    close(funit)
    ! write(*,*) '   time         computed            exact            error           step-error'
    ! write(*,'(f10.4,4(es18.8))') t, u0(1), xex, abs(u0(1)-xex), err

end program test

