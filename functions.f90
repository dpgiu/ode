module functions
  use precision
  implicit none
  private

  public :: ff
  public :: pendolo
  public :: harmonic
  public :: func
  
  real(dp), public :: k1 
  real(dp), public :: k2
  real(dp), public :: Q
  real(dp), public :: A
  real(dp), public :: w

  interface 
    function func(t,u) result(up)   
      use precision    
      real(dp), intent(in) :: t    
      real(dp), intent(in) :: u(:)    
      real(dp), allocatable :: up(:)    
    end function
  end interface
  

  contains

  subroutine fsub(t,u,up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)     
    real(dp), intent(out) :: up(:)       

    up(1) = u(2)
    up(2) = k1*u(1) + k2*u(2)
      
  end subroutine fsub
      
  function f1(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))

    up(1) = k1*u(1) 

  end function f1

  
  function ff(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)
    up(2) = k1*u(1) + k2*u(2)

  end function ff

  function pendolo(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)
    up(2) = -sin(u(1)) - u(2)/Q + A*cos(w*t)

  end function pendolo

    function harmonic(t,u) result(up)
    real(dp), intent(in) :: t    
    real(dp), intent(in) :: u(:)
    real(dp), allocatable :: up(:)
    
    allocate(up(size(u)))
    
    up(1) = u(2)
    up(2) = -u(1) - u(2)/Q + A*cos(w*t)

  end function harmonic

 

  
end module functions


