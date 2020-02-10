program test
  use precision
  use functions
  use solvers 
  implicit none

  interface 
    subroutine  solver(ff, t, dt, u0, u, err)
       use precision
       interface 
         function ff(t,u) result(up)   
           use precision    
           real(dp), intent(in) :: t    
           real(dp), intent(in) :: u(:)    
           real(dp), allocatable :: up(:)    
         end function
       end interface
       real(dp), intent(in) :: t
       real(dp), intent(in) :: dt
       real(dp), intent(in) :: u0(:)
       real(dp), intent(inout) :: u(:)
       real(dp), intent(inout) :: err
    end subroutine solver
  end interface 

  real(dp) :: err, t, dt, x0, v0, AA, BB, xex
  real(dp), allocatable :: u(:), u0(:) 
  integer :: i, Nstep
  procedure(solver), pointer :: psolver 
  character(50) :: arg

  if (command_argument_count() < 3) then
     write(*,*) "test dt Nstep method"
     stop
  endif

  call get_command_argument(1,arg)
  read(arg,*) dt 

  call get_command_argument(2,arg)
  read(arg,*) Nstep

  call get_command_argument(3,arg)

  allocate(u(2))
  allocate(u0(2))

  k1 = -1.0_dp
  k2 = 0.0_dp
  x0 = 10.0_dp
  v0 = 0.0_dp
  u0 = (/x0, v0/)
  t = 0.0_dp
  err = 0.0_dp

  AA = x0
  !AA = (v0 - k2*x0) / (k1-k2)  
  !BB = (x0*k1 - v0) / (k1-k2)   
  select case(trim(arg))
  case('rk2')
    do i = 1, Nstep
      call rk2(ff, t, dt, u0, u)
      t = t+dt
      xex = AA*cos(-k1*t)
      u0 = u
    end do  
  case('rk4')
    do i = 1, Nstep
      call rk4(ff, t, dt, u0, u)
      t = t+dt
      xex = AA*cos(-k1*t)
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
  end select
  write(*,*) '   time         computed            exact            error           step-error'
  write(*,'(f10.4,4(es18.8))') t, u0(1), xex, abs(u0(1)-xex), err
    
end program test

