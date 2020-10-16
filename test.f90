program test
  use precision
  use functions
  use solvers 
  implicit none

  real(dp) :: err, t, dt, x0, v0, AA, BB, xex
  real(dp), allocatable :: u(:), u0(:) 
  integer :: i, Nstep, funit
  procedure(solver), pointer :: psolver 
  character(50) :: arg
  character(10) :: solname
  logical :: booleano ! = .true. || .false.

  if (command_argument_count() < 8) then
     write(*,*) "test dt Nstep method Q A w x0 v0"
     stop
  endif

  call get_command_argument(1,arg)
  read(arg,*) dt 

  call get_command_argument(2,arg)
  read(arg,*) Nstep

  call get_command_argument(3,solname)

  call get_command_argument(4,arg)
  read(arg,*) Q

  call get_command_argument(5,arg)
  read(arg,*) A

  call get_command_argument(6,arg)
  read(arg,*) w

  call get_command_argument(7,arg)
  read(arg,*) x0

  call get_command_argument(8,arg)
  read(arg,*) v0
  

  allocate(u(2))
  allocate(u0(2))

  ! u0(1) = u0
  ! u0(2) = v0
  
  u0 = (/x0, v0/)
  t = 0.0_dp
  err = 0.0_dp

  ! arg ='rk2  <-- trailing white spaces -->'
  !if (trim(arg) == 'rk2') then
  !   ...
  !else if (trim(arg) = 'rk4') then
  !  ...
  !end if

  open(newunit=funit, file="solution.dat")
     
  select case(trim(solname))
  case('rk2')
     do i = 1, Nstep
      write(funit,*) t, u0(1), u0(2)   
      call rk2(harmonic, t, dt, u0, u)
      t = t+dt
      u0 = u
    end do  
  case('rk4')
    do i = 1, Nstep
      write(funit,*) t, u0(1), u0(2)  
      call rk4(harmonic, t, dt, u0, u)
      t = t+dt
      u0 = u
    end do  
  case('dp54')    
    do i = 1, Nstep
      call dopri54(ff, t, dt, u0, u, err)
      t = t+dt
      xex = AA*cos(-k1*t)
      u0 = u
    end do  
  case('dp87')
    do i = 1, Nstep
      call dopri87(ff, t, dt, u0, u, err)
      t = t+dt
      xex = AA*cos(-k1*t)
      u0 = u
    end do
  case default
     write(*,*) 'error: solver ',trim(arg),' not found'
     stop
  end select

  close(funit)
  
  !write(*,*) '   time         computed            exact            error           step-error'
  !write(*,'(f10.4,4(es18.8))') t, u0(1), xex, abs(u0(1)-xex), err
    
end program test

