module cl_mod
  use healpix_types
  use evolution_mod
  use sphbess_mod
  implicit none

  real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
  real(dp),     allocatable, dimension(:)       :: z_spline
  real(dp),     pointer,     dimension(:)       :: x_hi, k_hi ! Changed to pointers, right???    

contains
  ! Driver routine for (finally!) computing the CMB power spectrum
  subroutine compute_cls
    implicit none

    integer(i4b) :: i, j, l, l_num, x_num, n_spline, ki, xi, ni
    real(dp)     :: dx, S_func, j_func, z, eta, eta0, x0, x_min, x_max, d, e
    integer(i4b), allocatable, dimension(:)       :: ls
    real(dp),     allocatable, dimension(:)       :: Cl, Cl_hi, Cl2 !Power spectrum!!!
    real(dp),     allocatable, dimension(:)       :: integrand_tr
    real(dp),     allocatable, dimension(:,:)     :: integrand, integrand2, integrand_hi
  ! real(dp),     pointer,     dimension(:,:)     :: j_l, j_l2
    real(dp),     pointer,     dimension(:)       :: x_arg, int_arg, lsd, l_hi
    real(dp),     pointer,     dimension(:)       :: k, x
    real(dp),     pointer,     dimension(:,:,:,:) :: S_coeff
    real(dp),     pointer,     dimension(:,:)     :: S, S2
    real(dp),     allocatable, dimension(:,:)     :: Theta, Theta2, Theta_hi
    real(dp),     allocatable, dimension(:)       :: j_l_spline, j_l_spline2, etax
 !  real(dp),     pointer,     dimension(:)       :: x_hi, k_hi ! Changed to pointers, right???
    real(dp), parameter                           :: out_unit=20
    real(dp), parameter                           :: out_unitlk=44
    real(dp), parameter                           :: out_unitl1=43
    real(dp), parameter                           :: out_unitl5=42
    real(dp), parameter                           :: out_unitl10=41
    real(dp), parameter                           :: out_unitl20=40
    real(dp), parameter                           :: out_unitl30=39
    real(dp), parameter                           :: out_unitl44=38
    real(dp), parameter                           :: out_unitk1=11
    real(dp), parameter                           :: out_unitk2=12
    real(dp), parameter                           :: out_unitk3=13
    real(dp), parameter                           :: out_unitk4=14
    real(dp), parameter                           :: out_unitk5=15
    real(dp), parameter                           :: out_unitk6=16

    real(dp)                                      :: t1, t2, integral
    logical(lgt)                                  :: exist
    character(len=128)                            :: filename
    real(dp), allocatable, dimension(:)           :: y, y2
    real(dp)                                      :: n_PS = 0.96d0
    integer(i4b)                                  :: nlhi = 1200
    real(dp)                                      :: incrementx, incrementk
    real(dp)                                      :: k_min    = 0.1d0 * H_0 / c
    real(dp)                                      :: k_max    = 1.d3  * H_0 / c
    real(dp)                                      :: x_start_rec, h_trapez, h_trapez_tr, etanull, jl_splined, Cl_temp
    x_start_rec = -log(1.d0 + 1630.4d0)

    ! Set up which l's to compute, these are the ones used in Callin
    l_num = 44
    allocate(ls(l_num))
    ls = (/ 2, 3, 4, 6, 8, 10, 12, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, &
         & 120, 140, 160, 180, 200, 225, 250, 275, 300, 350, 400, 450, 500, 550, &
         & 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200 /)

    !write(*,*) ls(1), ls(5), ls(10), ls(20), ls(30), ls(44)

    ! Task: Get source function from evolution_mod
    allocate(S(n_x_hi,n_k_hi))
    allocate(x_hi(n_x_hi))
    allocate(k_hi(n_k_hi))

    !Making x and k arrays of length 5000                                                    
    incrementx = (-x_start_rec)/(n_x_hi-1)
    incrementk = (k_max-k_min)/(n_k_hi-1)
    do i = 1,n_x_hi
       x_hi(i) = x_start_rec+((i-1)*incrementx);
       k_hi(i) = k_min+((i-1)*incrementk);
    enddo
    x_hi(n_x_hi)=0.d0 ! Making sure last value in x_hi is zero 

    call initialize_perturbation_eqns
    call integrate_perturbation_eqns    
    call get_hires_source_function(k_hi,x_hi,S)
    !write(*,*) 'SOURCE FUNC S = ', S
    
    ! Computing Source function only if it has not already been done                                 
    inquire(file="Shires.unf", exist=exist)
    if (exist) then
       open(68,file="Shires.unf",form="unformatted")
       read(68) S
       read(68) x_hi
       read(68) k_hi
       close(68)
       write(*,*) 'S_hires exists'
    else
       write(*,*) 'Computing Shires'
       call initialize_perturbation_eqns                                                           
       call integrate_perturbation_eqns 
       call get_hires_source_function(k_hi,x_hi,S)
       open(68,file="Shires.unf",form="unformatted")
       write(68) S
       write(68) x_hi
       write(68) k_hi
       close(68)
    end if

    ! Task: Initialize spherical Bessel functions for each l; use 5400 sampled points between 
    !       z = 0 and 3500. Each function must be properly splined
    ! Hint: It may be useful for speed to store the splined objects on disk in an unformatted
    !       Fortran (= binary) file, so that these only has to be computed once. Then, if your
    !       cache file exists, read from that; if not, generate the j_l's on the fly.
    n_spline = 5400
    allocate(z_spline(n_spline))         ! Note: z is *not* redshift, but simply the dummy argument of j_l(z)
    allocate(j_l(n_spline, l_num))       ! Bessel functions
    allocate(j_l2(n_spline, l_num))      ! 2nd deriv of bessel funcs, for splint
    allocate(integrand_tr(n_x_hi))       ! integrand in transfer func
    allocate(Theta(l_num,n_k_hi))        ! transfer function (plot)
    allocate(Theta2(l_num,n_k_hi))
    allocate(Theta_hi(nlhi,n_k_hi))
    allocate(integrand(l_num,n_k_hi))    ! integrand in power spectrum (plot)
    allocate(integrand2(l_num,n_k_hi))
    allocate(integrand_hi(nlhi,n_k_hi))
    allocate(Cl(l_num))                  ! POWER SPECTRUM <3. NOT L_NUM, MAYBE Z_SPLINE???
    allocate(Cl2(l_num))
    allocate(Cl_hi(nlhi))
    allocate(etax(n_x_hi))

    ! Computing bessel functions only if it has not already been done
    inquire(file="besselfuncs.unf", exist=exist)
    if (exist) then
       open(58,file="besselfuncs.unf",form="unformatted")
       read(58) z_spline
       read(58) j_l
       read(58) j_l2
       close(58)
       write(*,*) 'Besselfuncs exist'
    else
       write(*,*) 'Computing Besselfuncs j_l'
       do i = 1, n_spline
          z_spline(i) = (i-1)*3500.d0/(n_spline-1.d0)
          do l = 1, l_num
             if (z_spline(i) > 2.d0) then
                call sphbes(ls(l),z_spline(i),j_l(i,l))
             endif
          enddo
       enddo
       write(*,*) 'Splining Besselfuncs j_l'
       do l = 1, l_num
          call spline(z_spline, j_l(:,l), yp1, ypn, j_l2(:,l))
       enddo
       open(58,file="besselfuncs.unf",form="unformatted")
       write(58) z_spline
       write(58) j_l
       write(58) j_l2
       close(58)
    end if
    !write(*,*) 'BESSEL j_l =', j_l

    ! Overall task: Compute the C_l's for each given l
    write(*,*) 'Computing integrand and transfer func'
    h_trapez_tr = abs((x_hi(n_x_hi)-x_hi(1))/n_x_hi) ! stepsize for trapezoidal method
    !h_trapez_tr = abs((xgrid(n_t)-xgrid(1))/n_t)
    h_trapez    = abs((k_hi(n_k_hi)-k_hi(1))/n_k_hi)
    !write(*,*) 'Running with n_k instead of n_k_hi in cl_mod'
    write(*,*) 'Running over n_x_hi! (cl_mod)'
    open (unit=out_unitlk, file="integrand_tr1001.txt", action="write", status="replace")
    open (unit=out_unitl1,   file="integrandl41.txt",   action="write", status="replace")
    open (unit=out_unitl5,   file="integrandl501.txt",  action="write", status="replace")
    open (unit=out_unitl10,  file="integrandl2001.txt", action="write", status="replace")
    open (unit=out_unitl20,  file="integrandl5001.txt", action="write", status="replace")
    open (unit=out_unitl30,  file="integrandl8001.txt", action="write", status="replace")
    open (unit=out_unitl44, file="integrandl12001.txt", action="write", status="replace")
    etanull = get_eta(0.d0)
    do i = 1, n_x_hi
       etax(i) = get_eta(x_hi(i)) !hi
    enddo
    do l = 1, l_num
       write(*,*) 'l = ', l
       do ki = 1, n_k_hi !SHOULD BE HI
          ! Computing integrand for transfer function
          do xi = 1, n_x_hi! t ! n_x_hi
             jl_splined = splint(z_spline, j_l(:,l), j_l2(:,l), k_hi(ki)*(etanull-etax(xi)))
             !integrand_tr(xi) = S(ki,xi)*Bessellkx(l,k_hi(ki),x_hi(xi)) !JUST THIS TAKES 20 MINUTES!!!! ca 30 sec per k
             integrand_tr(xi) = S(ki,xi)*jl_splined
             if ((l == 17) .and. (ki == 2000)) then
                write (out_unitlk,*), x_hi(xi) , integrand_tr(xi)
             endif
          enddo 
          !write(*,*) 'INTEGRAND TR = ', integrand_tr
          !stop
       ! Task: Compute the transfer function, Theta_l(k)
          ! Using trapezoidal method CHECK IN NUMERICAL RECIPES, NOT SURE IF RIGHT
          !Theta(l,ki) = (integrand_tr(1)+integrand_tr(n_x_hi))/2.d0
          Theta(l,ki) = integrand_tr(1)+integrand_tr(n_x_hi) !n_x_hi
          do xi = 2, n_x_hi !t !n_x_hi-1
             !Theta(l,ki) = Theta(l,ki)+integrand_tr(xi)
             !Theta(l,ki) = (x_hi(xi+1)-x_hi(xi))*(integrand_tr(xi+1)+integrand_tr(xi))/2.d0 ! NON UNIFORM
             Theta(l,ki) = Theta(l,ki)+2.d0*integrand_tr(xi)
          enddo
          Theta(l,ki) = h_trapez_tr*Theta(l,ki)!/2.d0
          !write(*,*) Theta(17,ki)
          ! Computing integrand for Power spectrum
          integrand(l,ki) = (c*k_hi(ki)/H_0)**(n_PS-1.d0)*Theta(l,ki)**2/k_hi(ki)
          !integrand(l,ki) = Theta(l,ki)**2/k_hi(ki) !WITHOUT P
          ! writing integrand to file for 6 different ls 
          if (ls(l) == 4) then
             write (out_unitl1,*), k_hi(ki), integrand(l,ki)
          endif
          if (ls(l) == 50) then
             write (out_unitl5,*), k_hi(ki), integrand(l,ki)
          endif
          if (ls(l) == 200) then
             write (out_unitl10,*), k_hi(ki), integrand(l,ki)
          endif
          if (ls(l) == 500) then                                                                
             write (out_unitl20,*), k_hi(ki), integrand(l,ki)
          endif
          if (ls(l) == 800) then
             write (out_unitl30,*), k_hi(ki), integrand(l,ki)
          endif
          if (l == 44) then
             write (out_unitl44,*), k_hi(ki), integrand(l,ki)
          endif
       enddo
       !write(*,*) integrand
       ! Task: Integrate P(k) * (Theta_l^2 / k) over k to find un-normalized C_l's
       ! Using trapezoidal method
       Cl_temp = (integrand(l,1)+integrand(l,n_k_hi))!/2.d0
       !Cl(l) = integrand(l,1)+integrand(l,n_k_hi)
       do ki = 2, n_k_hi-1
          Cl_temp = Cl_temp+2.d0*integrand(l,ki)
          !Cl(l) = (k_hi(ki+1)-k_hi(ki))*(integrand(l+1,ki)+integrand(l,ki))/2.d0 NON-UNIFORM
          !Cl_temp = Cl_temp+2.d0*integrand(l,ki)
       enddo
       ! Task: Store C_l in an array. Optionally output to file
       !Cl(l) = h_trapez*Cl(l)*ls(l)*(ls(l)+1.d0)/(4.d0*pi)
       Cl(l) = h_trapez*Cl_temp*ls(l)*(ls(l)+1.d0)/(2.d0*pi) ! units for plot
    end do
    close (out_unitlk)
    close (out_unitl1)
    close (out_unitl5)
    close (out_unitl10)
    close (out_unitl20)
    close (out_unitl30)
    close (out_unitl44)
    write(*,*) 'Done with the massive loop!!!'

    !write(*,*) integrand_tr!, integrand, Theta
    !write(*,*)' CL NOT HIRES = ', Cl

    ! Task: Spline C_l's found above, and output smooth C_l curve for each integer l
    allocate(lsd(l_num))
    do l = 1, l_num
       lsd(l) = ls(l) ! Need double precision in spline function!!!
    enddo
    call spline(lsd, Cl, yp1, ypn, Cl2)
    call spline(lsd, Theta(:,l), yp1, ypn, Theta2(:,l))

    ! New lgrid with unit stepsize, also dp 
    allocate(l_hi(nlhi))
    do l = 1, nlhi
       l_hi(l) = l
    enddo
    write(*,*) 'Making hires transferfunc'
    ! Making and writing hires transferfunc to file
    open (unit=out_unitk1,  file="thetak1.txt", action="write", status="replace")
    open (unit=out_unitk2,  file="thetak2.txt", action="write", status="replace")
    open (unit=out_unitk3,  file="thetak3.txt", action="write", status="replace")
    open (unit=out_unitk4,  file="thetak4.txt", action="write", status="replace")
    open (unit=out_unitk5,  file="thetak5.txt", action="write", status="replace")
    open (unit=out_unitk6,  file="thetak6.txt", action="write", status="replace")
    do l = 1, nlhi
       Theta_hi(l,k) = splint(lsd, Theta(:,l), Theta2(:,l), l_hi(l))
       write (out_unitk1,*), l, Theta_hi(l,2)
       write (out_unitk2,*), l, Theta_hi(l,100)
       write (out_unitk3,*), l, Theta_hi(l,500)
       write (out_unitk4,*), l, Theta_hi(l,1000)
       write (out_unitk5,*), l, Theta_hi(l,3000)
       write (out_unitk6,*), l, Theta_hi(l,4999)
    enddo
    close (out_unitk1)
    close (out_unitk2)
    close (out_unitk3)
    close (out_unitk4)
    close (out_unitk5)
    close (out_unitk6)

    write(*,*) 'Making the final POWERSPECTRUM!!!!'
    ! Making hires array for final power spectrum!
    open (unit=out_unit, file="POWERSPECTRUM1.txt", action="write", status="replace")
    do l = 1, nlhi
       Cl_hi(l) = splint(lsd, Cl, Cl2, l_hi(l))
       write (out_unit,*), l, Cl_hi(l)
    enddo
    !write(*,*) Cl_hi
    !do l = 1, l_num
    !   write(out_unit,*) l, Cl(l)
    !enddo
    close (out_unit)
    write(*,*) 'Done, donedonedone done'
  end subroutine compute_cls
  
  !function Bessellkx(l,k,x)
  !    implicit none
  !    integer(i4b), intent(in) :: l
  !    real(dp),     intent(in) :: x, k
  !    real(dp)                 :: Bessellkx
  !    real(dp)                 :: etanull, etax
  !    etanull = get_eta(0.d0)
  !    etax    = get_eta(x)
  !    Bessellkx = splint(z_spline, j_l(:,l), j_l2(:,l), k*(etanull-etax))
  !end function Bessellkx

end module cl_mod
