module caotic_pendulum
  use precision
  use functions
  use solvers
  implicit none

  real(dp), parameter :: pi = 3.14159265358979323846264338327_dp 
  real(dp), parameter :: eps = 1.d-10

contains

  subroutine solve(Nstep, Nperiods, x0, v0)
    integer, intent(in) :: Nstep, Nperiods    
    real(dp), intent(in) :: x0, v0

    real(dp) :: t, tfin, ttrans, dt, err
    real(dp), allocatable :: u(:), u0(:), uex(:)
    integer :: Ntot, Ntrans, id, id2
    integer :: i
    
    dt = 2.0_dp*pi/(w*real(Nstep, dp))

    Ntot = Nstep*Nperiods  
    tfin = Ntot*dt
    Ntrans = int(Nperiods*0.0)
    ttrans = Nstep*Ntrans*dt

    write(*,*) 'dt=',dt

    allocate(u0(2))
    allocate(uex(2))
    u0 = (/x0, v0/)
    u = u0
    t = 0.0_dp

    do while (t<ttrans)
       call dopri54(pendolo, t, dt, u0, u, err)
       !call rk4(pendolo, t, dt, u0, u)
       t = t + dt
       u0 = u
    end do

    
    open(newunit=id, file='sol.dat')
    !open(newunit=id2, file='poincare.dat', position='APPEND')

    i = 0
    do while (t<tfin)
       write(id,'(F20.8,4(ES20.8))') t, u, u(2)*u(2)/2.0_dp + 1.0_dp - cos(u(1)), err
       write(*,'(F20.8,3(ES20.8))') t, u, err

       !if (mod(i,Nstep) == 0 ) then
       !   write(id2,*) Q, u(2)
       !end if

       call dopri54(pendolo, t, dt, u0, u, err)
       !call rk4(pendolo, t, dt, u0, u)
       t = t + dt
       i = i + 1
       u0 = u
    end do

    close(id)
    !close(id2)

  end subroutine solve

end module caotic_pendulum

