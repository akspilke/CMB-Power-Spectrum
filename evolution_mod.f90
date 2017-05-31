module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use spline_2D_mod
  use rec_mod
  implicit none

  ! Accuracy parameters
  real(dp),     parameter, private        :: a_init   = 1.d-8
  real(dp),     parameter, private        :: x_init   = log(a_init)
  real(dp),     parameter, private        :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private        :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter                 :: n_k      = 100
  integer(i4b), parameter                 :: n_k_hi   = 5000
  integer(i4b), parameter                 :: n_x_hi   = 5000
  integer(i4b), parameter, private        :: lmax_int = 6

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  ! Fourier mode list
  real(dp), allocatable, dimension(:)     :: ks
  ! Book-keeping variables
  real(dp),     private                   :: k_current
  integer(i4b), private                   :: npar = 6+lmax_int

contains
  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(k, x, S)
    implicit none
    real(dp), pointer,     dimension(:),   intent(out) :: k, x
    real(dp), pointer,     dimension(:,:), intent(out) :: S
    real(dp), allocatable, dimension(:,:)              :: S_lores
    real(dp), allocatable, dimension(:,:,:,:)          :: S_coeff
    real(dp), allocatable, dimension(:)                :: x_hi, k_hi
    integer(i4b) :: ki, xi, i
    real(dp)     :: g, dg, ddg, tau, dtau, ddtau, H_p, dH_p, dofHdH, Pi, dPi, ddPi, xval
    real(dp)     :: incrementx, incrementk
    real(dp)     :: x_start_rec, ck_cur
    x_start_rec = -log(1.d0 + 1630.4d0)

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 
    ! Substeps:                                                                                                       
      ! 1) First compute the source function over the existing k and x grids
    allocate(S_lores(n_k,n_t))
    allocate(S_coeff(4,4,n_k,n_t))
    allocate(x_hi(n_x_hi))
    allocate(k_hi(n_k_hi))
    allocate(S(n_k_hi,n_x_hi))

    do ki = 1, n_k
       k_current = ks(ki)
       ck_cur    = k_current*c
       do xi     = 1, n_t
          Pi     = Theta(xi,2,ki) ! without pol
          dPi    = dTheta(xi,2,ki)
          xval   = xgrid(xi)
          g      = get_g(xval)
          dg     = get_dg(xval)
          ddg    = get_ddg(xval)
          dtau   = get_dtau(xval)
          ddtau  = get_ddtau(xval)
          H_p    = get_H_p(xval)
          dH_p   = get_dH_p(xval)
          dofHdH = H_0**2*((Omega_b+Omega_m)*exp(-xval)+4.d0*(Omega_r+Omega_nu)*exp(-2.d0*xval)+4.d0*Omega_lambda*exp(2.d0*xval))/2.d0 
          ddPi   = (ck_cur/(5.d0*H_p))*(2.d0*(dH_p*Theta(xi,1,ki)/H_p+dTheta(xi,1,ki))&
                   -3.d0*(dH_p*Theta(xi,3,ki)/H_p+dTheta(xi,3,ki)))&
                   +0.9d0*(ddtau*Pi+dtau*dPi)
          !S_lores(ki,xi) = g*(Theta(xi,0,ki)+Psi(xi,ki)+(Pi/4.d0))+exp(get_tau(xgrid(xi)))*(dPsi(xi,ki)-dPhi(xi,ki))&
          !                 -(H_p*g*dv_b(xi,ki)+H_p*dg*v_b(xi,ki)+dH_p*g*v_b(xi,ki))+(3.d0/(4.d0*k_current))&
          !                 *(dH_p)**2*g*Pi+H_p*ddH_p*g*Pi+H_p*dH_p*dg*Pi+H_p*dH_p*g*dPi&
          !                 +2.d0*H_p*dH_p*dg*Pi+H_p**2*ddg*Pi+H_p**2*dg*dPi&
          !                 +2.d0*H_p*dH_p*g*dPi+H_p**2*dg*dPi+H_p**2*g*ddPi
          ! More like Callin
          !Pi = 0.d0
          S_lores(ki,xi) = g*(Theta(xi,0,ki)+Psi(xi,ki)+(Pi/4.d0))+exp(-get_tau(xval))*(dPsi(xi,ki)-dPhi(xi,ki))&
                           -(H_p*g*dv_b(xi,ki)+H_p*dg*v_b(xi,ki)+dH_p*g*v_b(xi,ki))/ck_cur+(3.d0/(4.d0*ck_cur**2))&
                           *((g*Pi*dofHdH)+3.d0*H_p*dH_p*(dg*Pi+g*dPi)+H_p**2*(ddg*Pi+2.d0*dg*dPi+g*dPi))
          write(*,*) xi, S_lores(ki,xi)
       enddo
       !write(*,*), S_lores
       stop
    enddo
    !write(*,*)
    write(*,*) 'S_lores done'
      ! 2) Then spline this function with a 2D spline
    call splie2_full_precomp(ks, xgrid, S_lores, S_coeff)
    write(*,*) 'S_lores splined'
      ! 3) Finally, resample the source function on a high-resolution uniform
      !    5000 x 5000 grid and return this, together with corresponding
      !    high-resolution k and x arrays
    ! Making x and k arrays of length 5000
    incrementx = (-x_start_rec)/(n_x_hi-1)  
    incrementk = (k_max-k_min)/(n_k_hi-1)
    do i = 1,n_x_hi
       x_hi(i) = x_start_rec+((i-1)*incrementx);
       k_hi(i) = k_min+((i-1)*incrementk);
    enddo
    x_hi(n_x_hi)=0.d0 ! Making sure last value in x_hi is zero
    ! Computing hires source function
    do ki = 1, n_k_hi
       do xi = 1, n_x_hi
          !write(*,*) minval(ks), k_hi(ki), maxval(ks)
          !write(*,*) minval(xgrid), x_hi(xi), maxval(xgrid)
          S(ki,xi) = splin2_full_precomp(ks, xgrid, S_coeff, k_hi(ki), x_hi(xi))
       enddo
    enddo
    write(*,*) 'Done with Source func!'
  end subroutine get_hires_source_function

  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none
    integer(i4b) :: l, i
    real(dp)     :: Hp0, dtau0

    ! Task: Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do i = 1, n_k ! n_k is 100
       ks(i) = k_min + (k_max-k_min)*((i-1)/99.d0)**2
    enddo

    ! Allocate arrays for perturbation quantities 
    allocate(Theta(  0:n_t, 0:lmax_int, n_k))   ! 0:n_t makes array index start at 0 w length n_t=500
    allocate(delta(  0:n_t, n_k))
    allocate(delta_b(0:n_t, n_k))
    allocate(v(      0:n_t, n_k))
    allocate(v_b(    0:n_t, n_k))
    allocate(Phi(    0:n_t, n_k))
    allocate(Psi(    0:n_t, n_k))
    allocate(dPhi(   0:n_t, n_k))
    allocate(dPsi(   0:n_t, n_k))
    allocate(dv_b(   0:n_t, n_k))
    allocate(dTheta( 0:n_t, 0:lmax_int, n_k))

    Hp0   = get_H_p(x_init)
    dtau0 = get_dtau(x_init)

    ! Task: Set up initial conditions for the Boltzmann and Einstein equations
    Phi(0,:)        = 1.d0                 ! for all ks these quantities start here at x=x_init
    delta(0,:)      = (3.d0/2.d0)*Phi(0,:)
    delta_b(0,:)    = delta(0,:)
    Theta(0,0,:)    = Phi(0,:)/2.d0 
    do i = 1, n_k   ! These initial conditions are different for different ks
       v(0,i)       = (c*ks(i)*Phi(0,i))/(2.d0*Hp0)
       v_b(0,i)     = v(0,i)
       Theta(0,0,i) = Phi(0,i)/2.d0
       Theta(0,1,i) = -(c*ks(i)*Phi(0,i))/(6.d0*Hp0)
       Theta(0,2,i) = -(20.d0*c*ks(i)*Theta(0,1,i))/(45.d0*Hp0*dtau0)! without polarization  
       do l = 3, lmax_int
          Theta(0,l,i) = -l*c*ks(i)*Theta(0,l-1,i)/((2.d0*l+1)*Hp0*dtau0)
       end do
    end do
  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none
    integer(i4b)                        :: i, j, k, l
    real(dp)                            :: x1, x2, R
    real(dp)                            :: eps, hmin, h1, x_tc, H_p, dt, t1, t2, xincrement1, xincrement2
    real(dp), allocatable, dimension(:) :: y, y_tight_coupling, dydx
    real(dp), allocatable, dimension(:) :: x_ev, x1_ev, x2_ev

    eps    = 1.d-8
    hmin   = 0.d0

    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))    
    allocate(x1_ev(n_k))
    allocate(x_ev(n_k+n_t))

    write(*,*) 'running over 10 not 100 ks! (ev_mod)'
    ! Propagate each k-mode independently
    do k = 1, n_k
       k_current = ks(k)  ! Store k_current as a global module variable
       h1        = 1.d-5
       ! Initialize equation set for tight coupling
       y_tight_coupling(1) = delta(0,k)
       y_tight_coupling(2) = delta_b(0,k)
       y_tight_coupling(3) = v(0,k)
       y_tight_coupling(4) = v_b(0,k)
       y_tight_coupling(5) = Phi(0,k)
       y_tight_coupling(6) = Theta(0,0,k) !Theta0 (l=0)
       y_tight_coupling(7) = Theta(0,1,k) !Theta1 (l=1)

       ! Find the time to which tight coupling is assumed, 
       ! and integrate equations to that time
       x_tc = get_tight_coupling_time(k_current)

       ! Making x grid, use if plot from x_init
       !xincrement1 = (x_tc-x_init)/(n_k-1)
       !do i = 1, n_k
       !   x1_ev(i) = x_init+(i-1)*xincrement1
       !enddo
       !x_ev = (/x1_ev,xgrid/) ! grid from a_init to today
 
       ! Task: Integrate from x_init until the end of tight coupling, using
       !       the tight coupling equations
       call odeint(y_tight_coupling, x_init, x_tc, eps, h1, hmin, y_tc_derivs, bsstep, output)
       ! LOOP IF YA WANNA PLOT FROM START
       !do i = 2, 100                                             
       !   call odeint(y_tight_coupling, x1_ev(i-1), x1_ev(i), eps, h1, hmin, y_tc_derivs, bsstep, output)
       !enddo

       ! Task: Set up variables for integration from the end of tight coupling 
       ! until today
       y(1:7)    = y_tight_coupling(1:7)
       y(8)      = -20.d0*c*ks(k)*y(7)/(45.d0*get_H_p(x_tc)*get_dtau(x_tc)) !Theta2
       do l = 3, lmax_int! thetal
          y(6+l) = -(l*c*ks(k)*y(5+l))/((2.d0*l+1)*get_H_p(x_tc)*get_dtau(x_tc))
       end do

       call odeint(y, x_tc, xgrid(1), eps, h1, hmin, y_derivs, bsstep, output)

       write(*,*) 'k = ',k, 'of 100'
      
       do i = 2, n_t
          ! Task: Integrate equations from tight coupling to today
          call odeint(y, xgrid(i-1), xgrid(i), eps, h1, hmin, y_derivs, bsstep, output)
          ! Store derivatives
          call y_derivs(xgrid(i), y, dydx) 
          ! Task: Store variables at time step i in global variables
          delta(i,k)      = y(1)
          delta_b(i,k)    = y(2)
          v(i,k)          = y(3)
          v_b(i,k)        = y(4)
          Phi(i,k)        = y(5)
          do l = 0, lmax_int
             Theta(i,l,k) = y(6+l)
          end do
          Psi(i,k)        = -y(5)-12.d0*H_0**2*Omega_r*y(8)/((c*ks(k)*exp(xgrid(i)))**2)
          ! Task: Store derivatives that are required for C_l estimation
          dv_b(i,k)       = dydx(4)
          dPhi(i,k)       = dydx(5)
          dTheta(i,:,k)   = dydx(6)
          dPsi(i,k)       = -dPhi(i,k)-12.d0*(H_0/(c*ks(k)*exp(xgrid(i))))**2*Omega_r*(dTheta(i,2,k)-2.d0*Theta(i,2,k))
       end do
    end do

    deallocate(y_tight_coupling)
    deallocate(y)
    deallocate(dydx)
  end subroutine integrate_perturbation_eqns

  ! Task: Complete the following routine, such that it returns the time at which
  !       tight coupling ends. In this project, we define this as either when
  !       dtau < 10 or c*k/(H_p*dtau) > 0.1 or x > x(start of recombination) 
  function get_tight_coupling_time(k)
    implicit none
    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    integer(i4b)          :: i,n
    real(dp)              :: x, x_start_rec
    n = 1d4
    x_start_rec = -log(1.d0 + 1630.4d0)
    do i = 0,n
       x = x_init + i*(0.d0-x_init)/n
       if (abs(get_dtau(x) < 10.d0) .and. abs(c*k/(get_H_p(x)*get_dtau(x)) <= 0.1d0) .and. x < x_start_rec) then
           get_tight_coupling_time = x
       endif
    enddo
  end function get_tight_coupling_time
 
  ! Finds the derivative of the Einstein and Boltzmann eqs after tight coupling                          
  subroutine y_derivs(x, y, dydx)
      use healpix_types
      implicit none
      real(dp),                 intent(in)  :: x
      real(dp), dimension(:),   intent(in)  :: y
      real(dp), dimension(:),   intent(out) :: dydx
      real(dp), dimension(:),   allocatable :: Theta, d_Theta
      real(dp)     :: a, q, Hp, R, dtau, Psi, k
      real(dp)     :: d_Phi, d_delta, d_deltab, d_v, d_vb
      real(dp)     :: delta, delta_b, v, v_b, Phi, ckoH
      integer(i4b) :: l !, k
      allocate(Theta(0:lmax_int), d_Theta(0:lmax_int))
      
      delta      = y(1)
      delta_b    = y(2)
      v          = y(3)
      v_b        = y(4)
      Phi        = y(5)
      !Theta(0)   = y(6)
      !Theta(1)   = y(7)
      !Theta(2)   = y(8)
      do l = 0, lmax_int
       Theta(l) = y(6+l)
      end do

      k          = k_current
      a          = exp(x)
      dtau       = get_dtau(x)
      Hp         = get_H_p(x)
      R          = (4.d0*Omega_r)/(3.d0*Omega_b*a)
      ckoH       = c*k/Hp

      Psi        = -Phi-12.d0*H_0**2*Omega_r*Theta(2)/((c*k*a)**2)
      d_Phi      = Psi-((ckoH**2)*Phi/(3.d0))+(H_0**2/(2.d0*Hp**2))*((Omega_m*delta/a)&
                   +(Omega_b*delta_b/a)+(4*Omega_r*Theta(0)/a**2))
      d_delta    = (ckoH*v)-3.d0*d_Phi
      d_deltab   = (ckoH*v_b)-3.d0*d_Phi
      d_Theta(0) = -ckoH*Theta(1)-d_Phi
      d_v        = -v-(ckoH*Psi)
      d_vb       = -v_b-(ckoH*Psi)+dtau*R*(3.d0*Theta(1)+v_b)
      d_Theta(1) = (ckoH/3.d0)*(Theta(0)-(2.d0*Theta(2))+Psi)+dtau*(Theta(1)+v_b/3.d0)
     !do l = 2,lmax_int-1
     !   Theta(l)   = y(l+6)
     !   d_Theta(l) = (ckoH*l*Theta(l-1)/(2.d0*l+1))-(l+1)*ckoH*Theta(l+1)/(2.d0*l+1)&
     !                +dtau*(Theta(l)-Theta(l)*abs(l==2)/10.d0) ! delta l,2 KROENEKER!!!!!!    
     !enddo
      do l = 2, lmax_int-1
      d_Theta(l) = real(l,dp)*c*k/(real(2*l+1,dp)*Hp)*Theta(l-1) - &
                   &real(l+1,dp)*c*k/(real(2*l+1,dp)*Hp)*Theta(l+1) + dtau*Theta(l)
      if (l == 2)  d_Theta(l) = d_Theta(l) - 0.1d0*dtau*Theta(l)
      end do
      d_Theta(lmax_int) = (ckoH*Theta(lmax_int-1))-c*(lmax_int+1)*Theta(lmax_int)/(Hp*get_eta(x))&
                          +dtau*Theta(lmax_int)  
      dydx(1)    = d_delta
      dydx(2)    = d_deltab
      dydx(3)    = d_v
      dydx(4)    = d_vb
      dydx(5)    = d_Phi
      do l = 0, lmax_int
         dydx(6+l) = d_Theta(l)
      enddo

      deallocate(Theta, d_Theta)
    end subroutine y_derivs

  ! Finds the derivative of the Einstein and Boltzmann eqs before/during tight coupling        
  subroutine y_tc_derivs(x, y_tc, dydx_tc)
      use healpix_types
      implicit none
      real(dp),                 intent(in)  :: x
      real(dp), dimension(:),   intent(in)  :: y_tc
      real(dp), dimension(:),   intent(out) :: dydx_tc
      real(dp), dimension(:),   allocatable :: Theta, d_Theta
      real(dp)     :: a, q, Hp, dHp, R, dtau, Psi, ckoH, k1
      real(dp)     :: d_Phi, d_delta, d_deltab, d_v, d_vb
      real(dp)     :: delta, delta_b, v, v_b, Phi
      allocate(Theta(0:lmax_int), d_Theta(0:lmax_int))

      delta      = y_tc(1)
      delta_b    = y_tc(2)
      v          = y_tc(3)
      v_b        = y_tc(4)
      Phi        = y_tc(5)
      Theta(0)   = y_tc(6)
      Theta(1)   = y_tc(7)

      k1         = k_current
      a          = exp(x)
      dtau       = get_dtau(x)
      Hp         = get_H_p(x)
      dHp        = get_dH_p(x)
      R          = (4.d0*Omega_r)/(3.d0*Omega_b*a)   
      ckoH       = c*k1/Hp
      Theta(2)   = (-20.d0*ckoH*Theta(1))/(45.d0*dtau)
      Psi        = -Phi-12.d0*H_0**2.d0*Omega_r*Theta(2)/((c*k1*a)**2.d0)
      
      d_Phi      = Psi-((ckoH**2.d0)*Phi/(3.d0))+(H_0**2.d0/(2.d0*Hp**2.d0))*((Omega_m*delta/a)&
                   +(Omega_b*delta_b/a)+(4.d0*Omega_r*Theta(0)/a**2.d0))
      d_Theta(0) = -ckoH*Theta(1)-d_Phi
      !q         = -(((1.d0-2.d0*R)*dtau+((1.d0+R)*get_ddtau(x))*(3.d0*Theta(1)+v_b))-(ckoH*Psi)&
      !            +(1.d0-(get_dH_p(x)/Hp))*ckoH*(-Theta(0)+2.d0*Theta(2))-(ckoH*d_Theta(0)))&
      !            /((1.d0+R)*dtau+(get_dH_p(x)/Hp)-1.d0)
      q          = 1.d0/((1.d0+R)*dtau + dHp/Hp - 1.d0) * (-((1.d0-2.d0*R)*dtau + (1.d0+R)&
                   *get_ddtau(x))*(3.d0*Theta(1)+v_b) - c*k1/Hp*Psi + (1.d0 - dHp/Hp)*c*k1&
                   /Hp*(-Theta(0)+2.d0*Theta(2)) + c*k1/Hp*(-d_Theta(0))) ! Setting dTheta(2) = 0 here

      d_delta    = (ckoH*v)-3.d0*d_Phi
      d_deltab   = (ckoH*v_b)-3.d0*d_Phi
      d_v        = -v-(ckoH*Psi)
      d_vb       = (1.d0/(1.d0+R))*(-v_b-(ckoH*Psi)+(R*(q+ckoH*(-Theta(0)+2.d0*Theta(2))-(ckoH*Psi))))
      d_Theta(1) = (q-d_vb)/3.d0

      dydx_tc(1) = d_delta
      dydx_tc(2) = d_deltab
      dydx_tc(3) = d_v
      dydx_tc(4) = d_vb
      dydx_tc(5) = d_Phi
      dydx_tc(6) = d_Theta(0)
      dydx_tc(7) = d_Theta(1)

      deallocate(Theta, d_Theta)
    end subroutine y_tc_derivs

end module evolution_mod
