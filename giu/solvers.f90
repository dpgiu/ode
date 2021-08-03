module solvers

  implicit none

  interface 
         function ff(t,u) result(up)   
          !use precision    
           real(8), intent(in) :: t    
           real(8), intent(in) :: u(:)    
           real(8), allocatable :: up(:)    
         end function
  end interface      


  contains
  
        subroutine eulero(f, t, dt, u0, u, err)
            procedure(ff) :: f
            real(8), intent(in) :: t, dt, u0(:)
            real(8), intent(inout) :: u(:)
            integer :: err
        
            u = u0 + f(t, u0)*dt
        end subroutine eulero


        subroutine rk2(f, t, dt, u0, u, err)
          procedure(ff) :: f    
          real(8), intent(in) :: t
          real(8), intent(in) :: dt
          real(8), intent(in) :: u0(:)
          real(8), intent(inout) :: u(:)
          integer :: err
 
          real(8), allocatable :: k1(:), k2(:)
          allocate(k1(size(u0)), k2(size(u0)))
 
          k1 = f(t, u0)          !f(t,u0)
          k2 = f(t+dt, u0 + dt*k1)
          u = u0 + (k1+k2)*dt*0.5 

        end subroutine rk2






































 end module 
