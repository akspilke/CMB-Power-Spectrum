program cmbspec
  use healpix_types
  use params
  use time_mod
  use rec_mod
  use evolution_mod
  use cl_mod
  implicit none

  !real(dp), parameter :: out_unit=20
  !real(dp), parameter :: out_unit2=20
  !real(dp), parameter :: out_unit3=20
  !real(dp), parameter :: out_unit4=20
  !real(dp), parameter :: out_unit5=20
  !real(dp), parameter :: out_unit6=20
  !real(dp), parameter :: out_unit7=20
  integer  :: i
  real(dp) :: k1, k2, k3, k4, k5, k6

  ! Initialize time grids
  call initialize_time_mod
  
  write(*,*) 'CMB program works, wohoo!'
  
   call initialize_rec_mod
   ! initialize evolution_mod                                               
   !call initialize_perturbation_eqns
   !call integrate_perturbation_eqns
   !call get_hires_souce_function
   call compute_cls

  !open (unit=out_unit, file="densities.txt", action="write", status="replace")
  !do i=1,1000
  !write (out_unit,*), Omega_b_array(i), Omega_m_array(i), Omega_r_array(i)
  !enddo
  !close (out_unit)
  
  !open (unit=out_unit2, file="etax1000take2.txt", action="write", status="replace")
  !do i=1,1000
  !write (out_unit2,*) eta(i), xgriduniform(i)
  !enddo
  !close (out_unit2)

  !open (unit=out_unit3, file="Hz.txt", action="write", status="replace")
  !do i=1,1000
  !write (out_unit3,*) H(i), zgriduniform(i)
  !enddo
  !close (out_unit3)

  !call initialize_rec_mod
  ! initialize evolution_mod
  !call initialize_perturbation_eqns
  !call integrate_perturbation_eqns
  
  k1 = 1
  k2 = 10
  k3 = 30
  k4 = 50
  k5 = 75
  k6 = 100

  !open (unit=out_unit, file="Phi.txt", action="write", status="replace")
  !open (unit=out_unit2, file="Psi.txt", action="write", status="replace")
  !open (unit=out_unit3, file="delta.txt", action="write", status="replace") 
  !open (unit=out_unit4, file="deltab.txt", action="write", status="replace") 
  !open (unit=out_unit5, file="v.txt", action="write", status="replace")
  !open (unit=out_unit6, file="vb.txt", action="write", status="replace")
  !open (unit=out_unit7, file="Theta0.txt", action="write", status="replace")

  !do i = 1, n_t
  !write(out_unit,'(7F20.7)')  xgrid(i), Phi(i,k1), Phi(i,k2), Phi(i,k3), Phi(i,k4), Phi(i,k5), Phi(i,k6)
  !write(out_unit2,'(6F20.7)') Psi(i,k1), Psi(i,k2), Psi(i,k3), Psi(i,k4), Psi(i,k5), Psi(i,k6)
  !write(out_unit3,'(6F20.7)') delta(i,k1), delta(i,k2), delta(i,k3), delta(i,k4), delta(i,k5), delta(i,k6)
  !write(out_unit4,'(6F20.7)') delta_b(i,k1), delta_b(i,k2), delta_b(i,k3), delta_b(i,k4), delta_b(i,k5), delta_b(i,k6)
  !write(out_unit5,'(6F20.7)') v(i,k1), v(i,k2), v(i,k3), v(i,k4), v(i,k5), v(i,k6)
  !write(out_unit6,'(6F20.7)') v_b(i,k1), v_b(i,k2), v_b(i,k3), v_b(i,k4), v_b(i,k5), v_b(i,k6)
  !write(out_unit7,'(6F20.7)') Theta(i,0,k1), Theta(i,0,k2), Theta(i,0,k3), Theta(i,0,k4), Theta(i,0,k5), Theta(i,0,k6)
  !enddo

  !close (out_unit)
  !close (out_unit2)
  !close (out_unit3)
  !close (out_unit4)
  !close (out_unit5)
  !close (out_unit6)
  !close (out_unit7)

  !open (unit=out_unit, file="Phi.txt", action="write", status="replace")              
  !do i=1,n_t                                                                        
  !write (out_unit,'(7F20.7)')  xgrid(i), Phi(i,k1), Phi(i,k2), Phi(i,k3), Phi(i,k4), Phi(i,k5), Phi(i,k6)
  !enddo                                                                           
  !close (out_unit)

  !open (unit=out_unit2, file="Psi.txt", action="write", status="replace")
  !do i=1,n_t
  !write (out_unit2,'(6F20.7)') Psi(i,k1), Psi(i,k2), Psi(i,k3), Psi(i,k4), Psi(i,k5), Psi(i,k6) 
  !enddo
  !close (out_unit2)

  !open (unit=out_unit3, file="delta.txt", action="write", status="replace")
  !do i=1,n_t
  !write (out_unit3,'(6F20.7)') delta(i,k1), delta(i,k2), delta(i,k3), delta(i,k4), delta(i,k5), delta(i,k6)
  !enddo
  !close (out_unit3)

  !open (unit=out_unit4, file="deltab.txt", action="write", status="replace")
  !do i=1,n_t
  !write (out_unit4,'(6F20.7)') delta_b(i,k1), delta_b(i,k2), delta_b(i,k3), delta_b(i,k4), delta_b(i,k5), delta_b(i,k6)
  !enddo
  !close (out_unit4)

  !open (unit=out_unit5, file="v.txt", action="write", status="replace")
  !do i=1,n_t
  !write (out_unit5,'(6F20.7)') v(i,k1), v(i,k2), v(i,k3), v(i,k4), v(i,k5), v(i,k6)
  !enddo
  !close (out_unit5)

  !open (unit=out_unit6, file="vb.txt", action="write", status="replace")
  !do i=1,n_t
  !write (out_unit6,'(6F20.7)') v_b(i,k1), v_b(i,k2), v_b(i,k3), v_b(i,k4), v_b(i,k5), v_b(i,k6)
  !enddo
  !close (out_unit6)

  !open (unit=out_unit7, file="theta0.txt", action="write", status="replace")
  !do i=1,n_t
  !write (out_unit7,'(6F20.7)') Theta(i,0,k1), Theta(i,0,k2), Theta(i,0,k3), Theta(i,0,k4), Theta(i,0,k5), Theta(i,0,k6)
  !enddo
  !close (out_unit7)

  !call get_hires_souce_function
  !call compute_cls

end program cmbspec
