module rec_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_1D_mod
  implicit none

  integer(i4b),                        private :: n                 ! Number of grid points
  real(dp), allocatable, dimension(:), private :: x_rec             ! Grid
  real(dp), allocatable, dimension(:), private :: tau, tau2, tau22  ! Splined tau and second derivatives
  real(dp), allocatable, dimension(:), private :: n_e, n_e2         ! Splined (log of) electron density, n_e
  real(dp), allocatable, dimension(:), private :: logn_e, logn_e2   ! Splined (log of) electron density, n_e                      
  real(dp), allocatable, dimension(:), private :: g, g2, g22        ! Splined visibility function

contains
  subroutine initialize_rec_mod
    implicit none
    
    integer(i4b)                        :: i, j, k
    real(dp)                            :: saha_limit, xstart, xstop !y, dydx, xmin, xmax, dx, f, n_e0, X_e0, xstart, xstop
    real(dp)                            :: incrementx, first, second, third, sahastuff
    logical(lgt)                        :: use_saha
    real(dp), allocatable, dimension(:) :: X_e ! Fractional electron density, n_e / n_H
    real(dp), parameter                 :: out_unit=20
    real(dp), parameter                 :: out_unit2=20
    real(dp), parameter                 :: out_unit3=20

    saha_limit = 0.99d0                 ! Switch from Saha to Peebles when X_e < 0.99
    xstart     = log(1.d-10)            ! Start grids at a = 10^-10
    xstop      = 0.d0                   ! Stop  grids at a = 1
    n          = 1000                   ! Number of grid points between xstart and xstopo
    
    epsilon    = 1.d-10                 ! Tolerance used in odeint                                               
    stepmin    = 0.d0                   ! Minimal stepsize, odeint reduces stepsize towards this value                 
    yp1        = 1.d30                  ! “natural spline”, zero second derivatives at endpoints                              
    ypn        = 1.d30

    allocate(x_rec(n))
    allocate(X_e(n))
    allocate(n_e(n))
    allocate(n_e2(n))
    allocate(logn_e(n))
    allocate(logn_e2(n))
    allocate(tau(n))
    allocate(tau2(n))
    allocate(tau22(n))
    allocate(g(n))
    allocate(g2(n))
    allocate(g22(n))
    
    ! Task: Fill in x (rec) grid
    ! Making a uniform xgrid, 1000 points between a=10^-10 and 1                           
    incrementx = (-xstart)/(n-1)   !Distance between gridpoints in uniform xgrid, 1000 points  
    do i=1,n
       x_rec(i)=xstart+((i-1)*incrementx);
    enddo

    ! Task: Compute X_e and n_e at all grid times
    step = abs(1.d-2*(x_rec(1)-x_rec(2)))    ! steplength for ode integration                

    use_saha = .true.
    do i = 1, n
       if (use_saha) then
          ! Use the Saha equation
          first     = (m_H*(exp(x_rec(i)))**3.d0)/(Omega_b*rho_c)
          second    = (m_e*T_0*k_b/(exp(x_rec(i))*2.d0*pi*(hbar**2)))**(3.d0/2.d0)
          third     = exp(-epsilon_0*exp(x_rec(i))/(T_0*k_b))
          sahastuff = first*second*third
          X_e(i)    = (-sahastuff+sqrt(sahastuff*(sahastuff+4.d0)))/2.d0 
          ! Used only plus because we want a positive density
          if (X_e(i) < saha_limit) use_saha = .false.
       else
          ! Use the Peebles equation
          X_e(i) = X_e(i-1)
          call odeint(X_e(i:i), x_rec(i-1), x_rec(i), epsilon, step, stepmin, Xe_derivs, bsstep, output)
       end if
    end do

    !open (unit=out_unit, file="n_eN.txt", action="write", status="replace")        
    do i = 1,n
       n_e(i) = ((Omega_b*rho_c)/(m_H*(exp(x_rec(i))**3.d0)))*X_e(i)
       !write (out_unit,*), x_rec(i), X_e(i), n_e(i)
    enddo
    !close (out_unit)

    ! Task: Compute splined (log of) electron density function    
    !n_e = log(n_e)
    call spline(x_rec, n_e, yp1, ypn, n_e2)

    !n_e2 = exp(logn_e2)
    !open (unit=out_unit, file="n_eN.txt", action="write", status="replace") 
    do i = 1,n
       !n_e(i) = get_n_e(x_rec(i))
       if (n_e(i) /= get_n_e(x_rec(i))) then
          write(*,*) 'get_n_e is WRONG'
       endif
       !write (out_unit,*), x_rec(i), X_e(i), n_e(i)                                      
    enddo
    !close (out_unit)
    
    ! Task: Compute optical depth at all grid points
    tau(n) = 0.d0
    ! Looping backwards from tau=0
    do i = n-1, 1, -1
       tau(i) = tau(i+1) !Setting next point in eta to be the previous, this gets updated in odeint      
       call odeint(tau(i:i), x_rec(i+1), x_rec(i), epsilon, step, stepmin, tau_derivs, bsstep, output)
    end do

    ! Task: Compute splined (log of) optical depth
    ! Did not use log as it did not work... log(0).
    call spline(x_rec, (tau), yp1, ypn, tau2)

    ! Task: Compute splined second derivative of (log of) optical depth                   
    call spline(x_rec, (tau2), yp1, ypn, tau22) 

    !open (unit=out_unit2, file="tauN.txt", action="write", status="replace")
    do i = 1,n
       if (tau(i) /= get_tau(x_rec(i))) then
          write(*,*) 'get_tau is WRONG'
       endif
       !write (out_unit2,*), get_tau(x_rec(i)), get_dtau(x_rec(i)), get_ddtau(x_rec(i))
    enddo
    !close (out_unit2)

    ! Task: Compute splined visibility function
    do i = 1,n
       g(i) = -get_dtau(x_rec(i))*exp(-get_tau(x_rec(i)))
    enddo
    call spline(x_rec, (g), yp1, ypn, g2)

    ! Task: Compute splined second derivative of visibility function
    call spline(x_rec, (g2), yp1, ypn, g22)

    !open (unit=out_unit3, file="gN.txt", action="write", status="replace") 
    do i = 1,n
       if (g(i) /= get_g(x_rec(i))) then
          write(*,*) 'get_g is WRONG'
       endif
       !write (out_unit3,*), get_g(x_rec(i)), get_dg(x_rec(i)), get_ddg(x_rec(i))                  
    enddo
    !close (out_unit3)                                                                                                             
  end subroutine initialize_rec_mod

  ! Finds the derivative of Xe wrt x at a given x, Peebles equation                          
  subroutine Xe_derivs(x, X_e, derivative)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: X_e
      real(dp), dimension(:), intent(out) :: derivative
      real(dp)                            :: a, T_b
      real(dp) :: Lambda2s, n_H, n1s, Lambda_a, phi2, alpha2, beta, beta2, Cr, Xe
      a            = exp(x)
      T_b          = T_0/a
      Xe           = X_e(1)
      Lambda2s     = 8.227d0                         ! s-1     
      n_H          = (Omega_b*rho_c)/(m_H*(a**3.d0)) ! m-3
      n1s          = (1.d0-Xe)*n_H                   ! m-3
      Lambda_a     = get_H(x)*((3.d0*epsilon_0)**3.d0)/(((8.d0*pi)**2.d0)*n1s*((hbar*c)**3))  ! s-1           
      phi2         = 0.448d0*log(epsilon_0/(T_b*k_b))! no unit
      alpha2       = ((64.d0*pi*alpha**2.d0)/(sqrt(27.d0)*pi*m_e**2.d0))*sqrt(epsilon_0/(T_b*k_b))*phi2*((hbar**2)/c) ! m3/s
      beta         = alpha2*((m_e*T_b*k_b/(2.d0*pi*(hbar**2)))**(3.d0/2.d0))*exp(-epsilon_0/(T_b*k_b)) ! s-1
      ! Setting a temperature limit to change 10e-200 to zero to avoid NaNs
      if(T_b <= 169.d0) then
            beta2  = 0.d0
        else
            beta2  = beta*exp((3.d0*epsilon_0)/(4.d0*k_b*T_b))          ! s-1
        end if
      Cr           = (Lambda2s+Lambda_a)/(Lambda2s+Lambda_a+beta2)      ! s-1      
      derivative   = (Cr/get_H(x))*(beta*(1.d0-Xe)-n_H*alpha2*Xe**2.d0) ! s-1          
    end subroutine Xe_derivs

  ! Finds the derivative of the optical depth wrt x at a given x                                                      
  subroutine tau_derivs(x, tau, dtau)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: tau
      real(dp), dimension(:), intent(out) :: dtau
      dtau = -get_n_e(x)*sigma_T*c/get_H(x)
    end subroutine tau_derivs

  ! Task: Complete routine for computing n_e at arbitrary x, using precomputed information
  ! Hint: Remember to exponentiate...
  function get_n_e(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_n_e
    get_n_e = (splint(x_rec, n_e, n_e2, x))
  end function get_n_e

  ! Task: Complete routine for computing tau at arbitrary x, using precomputed information
  function get_tau(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_tau
    get_tau = (splint(x_rec, tau, tau2, x))
  end function get_tau

  ! Task: Complete routine for computing the derivative of tau at arbitrary x, using precomputed information
  function get_dtau(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_dtau
    get_dtau = splint_deriv(x_rec, tau, tau2, x)
  end function get_dtau

  ! Task: Complete routine for computing the second derivative of tau at arbitrary x, 
  ! using precomputed information
  function get_ddtau(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_ddtau
    get_ddtau = splint(x_rec, tau2, tau22, x) ! right??? should it be squared or /T
  end function get_ddtau

  ! Task: Complete routine for computing the visibility function, g, at arbitray x
  function get_g(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_g
    get_g = splint(x_rec, g, g2, x) 
  end function get_g

  ! Task: Complete routine for computing the derivative of the visibility function, g, at arbitray x
  function get_dg(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_dg
    get_dg = splint_deriv(x_rec, g, g2, x)
  end function get_dg

  ! Task: Complete routine for computing the second derivative of the visibility function, g, at arbitray x
  function get_ddg(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_ddg
    get_ddg = splint(x_rec, g2, g22, x)
  end function get_ddg

end module rec_mod
