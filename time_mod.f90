module time_mod
  use healpix_types
  use params
  use spline_1D_mod
  use ode_solver
  implicit none
  
  real(dp), allocatable, dimension(:) :: zgrid, agrid, xgrid   ! My version of x_t and a_t (one dim, if two then :,: in brackets)
  real(dp), allocatable, dimension(:) :: xgriduniform, agriduniform, zgriduniform ! My version of x_eta and a_eta
  real(dp), allocatable, dimension(:) :: xgridold, xgridnew    ! Two xgrids that gets put together to the non-uniform one
  real(dp), allocatable, dimension(:) :: H                     ! The Hubble parameter, array of values at different times
  real(dp), allocatable, dimension(:) :: Omega_m_array, Omega_r_array, Omega_b_array, Omega_v_array ! Arrays of densities

  integer(i4b)                        :: n_t                ! Number of x-values, 500                                                
  integer(i4b)                        :: n_eta              ! Number of eta grid poins, 1000                                         
  real(dp), allocatable, dimension(:) :: eta, eta2          ! Eta and eta''(second derivative) at each grid point   

  !real(dp), allocatable, dimension(:) :: x_t               ! Grid of relevant x-values, HKEs xgrid                                  
  !real(dp), allocatable, dimension(:) :: a_t               ! Grid of relevant a-values, HKEs agrid
  !real(dp), allocatable, dimension(:) :: x_eta             ! Grid points for eta, HKEs xgriduniform  

  real(dp) :: yp1, ypn, epsilon, stepmin, eta_init, step    ! Other variables used, defined below
    
contains    
  subroutine initialize_time_mod
    implicit none

    integer(i4b) :: i, n, n1, n2                               ! Integers used
    real(dp)     :: z_start_rec, z_end_rec, z_0, x_start_rec, x_end_rec, x_0, dx, x_eta1, x_eta2, a_init  ! Defined below
    real(dp)     :: xincrementold, xincrementnew, incrementxuniform, etad1(1), etadn(1)                   ! Defined below

    ! Define two epochs, 1) during and 2) after recombination.
    n1          = 200                       ! Number of grid points during recombination ("old" array)
    n2          = 300                       ! Number of grid points after recombination  ("new" array)
    n_t         = n1 + n2                   ! Total number of grid points in non uniform grids
    z_start_rec = 1630.4d0                  ! Redshift of start of recombination
    z_end_rec   = 614.2d0                   ! Redshift of end of recombination
    z_0         = 0.d0                      ! Redshift today
    x_start_rec = -log(1.d0 + z_start_rec)  ! x of start of recombination
    x_end_rec   = -log(1.d0 + z_end_rec)    ! x of end of recombination
    x_0         = 0.d0                      ! x today
    
    n_eta       = 1000                      ! Number of eta grid points (for spline)
    a_init      = 1.d-10                    ! Start value of a for eta evaluation
    x_eta1      = log(a_init)               ! Start value of x for eta evaluation
    x_eta2      = 0.d0                      ! End value of x for eta evaluation

    yp1         = 1.d30                     ! “natural spline”, zero second derivatives at endpoints
    ypn         = 1.d30

    eta_init    = a_init/(H_0*sqrt(Omega_r))! Initial conformal time
    epsilon     = 1.d-10                    ! Tolerance used in odeint
    stepmin     = 0.d0                      ! Minimal stepsize, odeint reduces stepsize towards this value

    ! Task: Fill in x and a grids
    ! 200 points between z = 1630 and z = 614  OLD
    ! 300 points between z = 614 and z = 0     NEW

    allocate(zgrid(n_t))
    allocate(agrid(n_t))
    allocate(xgrid(n_t))
    allocate(xgridnew(n2))
    allocate(xgridold(n1))
    allocate(xgriduniform(n_eta))
    allocate(agriduniform(n_eta))
    allocate(zgriduniform(n_eta))
    !allocate(x_t(n_t))     !HKEs version of xgrid                                                                                      
    !allocate(a_t(n_t))     !HKEs version of agrid  
    !allocate(x_eta(n_eta)) !HKEs xgriduniform                                                                       
    allocate(eta(n_eta))
    allocate(eta2(n_eta))

    ! Setting distance between gridpoints in old and new part of non-uniform xgrid
    xincrementold = (x_end_rec - x_start_rec)/(n1-1)
    xincrementnew = (x_0-x_end_rec)/n2

    !New, grid after recombination                                                                              
    do i=1,n2
       xgridnew(i)=x_end_rec+(i*xincrementnew);
    enddo

    !Old, grid during recombination
    do i=1,n1
       xgridold(i)=x_start_rec+((i-1)*xincrementold);
    enddo

    xgrid = (/xgridold,xgridnew/) !Joining old and new grids to one non-uniform xgrid with 500 points
    agrid = exp(xgrid)            !Making a corresponding grid of the scale factor a
    zgrid = (1.d0/agrid)-1.d0     !Making a corresponding grid of the redshift z

    ! Making a uniform xgrid, 1000 points between a=10^-10 and 1                             
    incrementxuniform = (-log(a_init))/n_eta   !Distance between gridpoints in uniform xgrid, 1000 points      
    do i=1,n_eta
       xgriduniform(i)=log(a_init)+((i-1)*incrementxuniform);
    enddo

    agriduniform = exp(xgriduniform)        !Making a corresponding grid of the scale factor a, 1000 points
    zgriduniform = (1.d0/agriduniform)-1.d0 !Making a corresponding grid of the redshift z, 1000 points

    ! Task: 1) Compute the conformal time at each eta time step
    !       2) Spline the resulting function, using the provided "spline" routine in spline_1D_mod.f90
    eta(1) = eta_init                                      ! Initial conformal time
    step = abs(1.d-2*(agriduniform(1)-agriduniform(2)))    ! steplength for ode integration

    ! Making a grid of conformal time eta grid, corresponding to the uniform x grid w 1000 values
    do i = 2, n_eta
       eta(i) = eta(i-1) !Setting next point in eta to be the previous, this gets updated in odeint                              
       call odeint(eta(i:i), agriduniform(i-1), agriduniform(i), epsilon, step, stepmin, eta_derivs, bsstep, output)
    end do

    ! Using the spline function to make conformal time continous between gridpoints in x
    call eta_derivs(agriduniform(1), eta(1:1), etad1)             ! Finds the derivative of eta at first value of a
    call eta_derivs(agriduniform(n_eta), eta(n_eta:n_eta), etadn) ! Finds the derivative of eta at last value of a 
    ! (i:i gives an array of length 1 rather than a number) COULD HAVE JUST USED YP1 AND YPN!!! MAYBE SHOULD HAVE
    call spline(xgriduniform, eta, etad1(1)*exp(xgriduniform(1)), etadn(1)**exp(xgriduniform(1)), eta2)
    ! spline needs y'', y', y, x, have x and y, spline finds others (etad is first derivative, eta2 is second)

    ! Make array of H values and densities, for plots in Milestone 1
    allocate(Omega_m_array(n_eta))
    allocate(Omega_r_array(n_eta))
    allocate(Omega_b_array(n_eta))
    allocate(Omega_v_array(n_eta))
    allocate(H(n_eta))

    do i = 1,n_eta
       H(i)=get_H(xgriduniform(i))
       Omega_m_array(i)=Omega_m*(agriduniform(i)**-3)
       Omega_b_array(i)=Omega_b*(agriduniform(i)**-3)
       Omega_r_array(i)=Omega_r*(agriduniform(i)**-4)
       Omega_v_array(i)=Omega_lambda
    enddo

    end subroutine initialize_time_mod

  ! Finds the derivative of the conformal time wrt scale factor at a given a
  subroutine eta_derivs(a, eta, derivative) 
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: a
      real(dp), dimension(:), intent(in)  :: eta
      real(dp), dimension(:), intent(out) :: derivative
      real(dp)                            :: x
      x = log(a)
      derivative = c/(a*a*get_H(x)) !d(eta)/da, eq (11) in Milestone 1, used in odeint
    end subroutine eta_derivs

  ! This subroutine needed to be copy-pasted from odesolver in order to make odeint work
  subroutine output(x, y)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
    end subroutine output

  ! Task: Write a function that computes H at given x
  ! Computes the Hubble constant for a given x=log(scale factor)
 function get_H(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_H
    real(dp)             :: a
    a = exp(x)
    get_H = H_0*sqrt((Omega_b+Omega_m)*a**-3.d0 + (Omega_r+Omega_nu)*a**-4.d0 + Omega_lambda)
  end function get_H

  ! Task: Write a function that computes H' = a*H  at given x
  function get_H_p(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_H_p
    real(dp)             :: a
    a = exp(x)
    get_H_p = a*get_H(x) ! Hubble constant prime is just H*a
  end function get_H_p

  ! Task: Write a function that computes dH'/dx at given x
  ! The derivative of H primed wrt x.
  function get_dH_p(x)
    implicit none
    real(dp), intent(in) :: x
    real(dp)             :: get_dH_p
    real(dp)             :: top
    real(dp)             :: bottom
    top = -(Omega_m + Omega_b)*exp(-x) - 2*(Omega_r + Omega_nu)*exp(-2*x) + 2*Omega_lambda*exp(2*x);
    bottom = sqrt((Omega_m + Omega_b)*exp(-x) + (Omega_r + Omega_nu)*exp(-2*x) + Omega_lambda*exp(2*x));
    get_dH_p = (H_0/2)*(top/bottom)
  end function get_dH_p

  ! Task: Write a function that computes eta(x), using the previously precomputed splined function
  function get_eta(x_in)
    implicit none
    real(dp), intent(in) :: x_in
    real(dp)             :: get_eta
    real(dp)             :: a_in
    a_in = exp(x_in)
    get_eta = splint(agriduniform, eta, eta2, a_in)
  end function get_eta
end module time_mod

