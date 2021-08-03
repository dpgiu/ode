program main

    use solvers
  
    implicit none
    
    integer :: N, err, ii
    real(8) :: dt, t
    real(8), dimension(2) :: u, u0, v, v0


    dt = 0.1
    t = 0.0
    u0 = 1.0   !mette a 0 tutte le componenti di u0
    v0 = u0
    N = 628

    do ii=1,N
       call eulero(g, t, dt, u0, u, err)
       call rk2(g, t, dt, v0, v, err)
        print*, t, " ", u, " ", v
        u0 = u
        v0 = v
        t = t + dt
    end do


  contains

    function g(t,u) result(up)
      real(8) :: k = 1.0
      real(8), intent(in) :: t    
      real(8), intent(in) :: u(:)    
      real(8), allocatable :: up(:)    
      allocate(up(size(u)))

      up(1) = u(2)
      up(2) = -k**2. * u(1)

    end function g
    
  
end program main
