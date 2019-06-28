!******************************************************************************
! MODULE: massgrwoth
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Modules that planet mass increases by pebble accretion 
!
!******************************************************************************

module massgrowth

  use orbital_elements
  use types_numeriques
  use physical_constant, only : PI, TWOPI, MSUN, K2, EARTH_MASS, AU
!  implicit none
  
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> Beibei  Liu
!
!> @date 12 Jan 2018
!
! DESCRIPTION: 
!> @brief Calculates the mass increase of the planet by pebble accretion
!!\n\n based on the formulas of Liu & Ormel 2018, Ormel & Liu2018
!! The routine also gives the compostion X, planet radius Rp of each object.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine pebbleaccretion (t,hstep,nbod,m,x,v,rphys,rho,m_x)
!  t          = current epoch [days]
!  hstep          = timestep [days]
!  x      = coordinates (x,y,z) with respect to the central body [AU]
!  v      = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  m          = planet mass (in solar masses * K2)
!  nbod          = number of bodies 
!  rphys          = planet radius [AU]
!  rho          = planet density
!  m_x          =  mass of element X in planet (in solar masses * K2)

!  use types_numeriques
!  use physical_constant, only : PI, TWOPI, MSUN, K2, EARTH_MASS, AU
  use mercury_globals, only : opt,alpha_t, kap, mdot_gas0, L_s, Rdisk0, Mdot_peb, alpha, tau_s
  use user_module, only: gasdisk 
 
  implicit none

  ! Input
  integer, intent(in) :: nbod
  real(double_precision),intent(in) :: t, hstep
  real(double_precision),intent(inout) :: x(3,nbod), v(3,nbod)
  real(double_precision),intent(inout) :: m(nbod)
  real(double_precision),intent(inout) :: m_x(nbod)
  real(double_precision),intent(inout) :: rphys(nbod)
  real(double_precision),intent(inout) :: rho(nbod)
  real(double_precision) :: ts(nbod) ! pebble's stokes number
  real(double_precision) :: f_ice ! ice fraction in pebbles
  real(double_precision) :: asemi(nbod), ecc(nbod), inc(nbod), rplanet(nbod), rdisk(nbod)
  real(double_precision) :: prob(nbod), prob2d(nbod), prob3d(nbod)
  real(double_precision) ::  siggas, hgas, etagas, rhogas
  real(double_precision) :: q, rp_to_a, mdot
  real(double_precision) :: rtemp(nbod), filter(nbod),temp, newprob(nbod), detm(nbod)
  real(double_precision) :: rho_ice, rho_sil, rho_old, vol_old, vol_new_ice, vol_new_sil
  real(double_precision) ::  m_iso, vol_tot 
  real(double_precision) :: rtrans, siggas_vis, siggas_irr
  integer :: ilist(nbod)
  integer i, k, j

   ! ### material density of water ice and silicate ###
   rho_ice = 1.d0* ( AU * AU * AU * K2 / MSUN)
   rho_sil = 3.d0* ( AU * AU * AU * K2 / MSUN)
  
  !  ### pebble mass flux (Mdot_peb): in earth mass/year ###
  !  ### mdot: in solar mass/day  ###
  mdot  = Mdot_peb*EARTH_MASS/365.25d0 



  !#######################################################################
  !### calculate from xyz to semimajor axis, ecc & inclination(radian) ###
  !#######################################################################
  call  xyztoaei(nbod,m,x,v,rdisk,asemi,ecc,inc)
  
  do i = 2, nbod
    !##########################################
    !### get disk infomation: etagas, hgas ####
    !##########################################
    call gasdisk (t,rdisk(i),x(3,i),m(1),alpha,kap,hgas,etagas,siggas,siggas_vis,siggas_irr,rhogas,rtrans)

    !###########################################
    !### get the pebble accretion efficiency ###
    !########################################### 
    ts(i) = tau_s
    rp_to_a = rphys(i)/asemi(i)  ! ratio of planet radius to semimajor axis
    if (rp_to_a > 0.d0) then ! still bound to the system
      q = m(i)/m(1) ! mass ratio between the planet and the central star 
      prob(i) = eps_PA (ts(i), q, ecc(i), inc(i), alpha_t, etagas, hgas, rp_to_a)
    else ! planet escape from the system (a<0), then no accretion 
      prob(i) = 0.d0
    end if 
  end do 


 
  ! ### initial set-up sepeartion rtemp ###
  do i = 1, nbod
    rdisk(i) = (x(1, i)**2 + x(2, i)**2)**0.5d0
    rtemp(i) = rdisk(i)
  end do 
   
  !### define a radii array that modst distant body's position ranks the first ###
  do i =1, nbod-1
    do j = i+1, nbod
        if (rtemp(i) < rtemp(j)) then
          temp = rtemp(i)
          rtemp(i) = rtemp(j)
          rtemp(j) = temp 
        end if 
    end do 
  end do 




  !### to calculate the filing factor of i planet ###
  prob(1) = 0.0d0 ! no mass increase of central star 
  !#### give the ilist where ilist(i)=1 means rdisk(i) = maximum #### 
  do i =1, nbod
    do j = 1, nbod
        if (rdisk(i) ==rtemp(j)) then
          ilist(i) = j
          newprob(j) = prob(i)
        end if
    end do
  end do


! ### filter is the flux ratio = M_real(i)/M_disk ###
! ### i=1,2,3 countes from the distant body ### 
! ### newprob(1,2,3) is the efficiency from the distant body while prob(1,2,3) os the efficiency from the list i ###
  do i =1, nbod
    if ( i ==1 ) then
      filter(i) = 1.d0
    else
      filter(i) = filter(i-1)*(1-newprob(i-1))
    end if
  end do


  ! ####################################################################################
  ! #### add the mass/radius increase dm = prob* Mdot_real = prob* Mdot_disk*filter ###
  ! ####################################################################################
  do i =2, nbod
    ! #### mass increase by pebble accretion #### 
    detm(i) = newprob(ilist(i))*hstep*mdot*K2*filter(ilist(i))
    f_ice = 0.5d0 ! right now assume outside of the water ice line 


    ! update the mass increase, consider the pebble isolation 
    m_iso = pebble_iso(hgas,alpha,etagas,ts(i),m(1))  ! in earth mass 
    m_iso = m_iso*EARTH_MASS*K2 !!!! unsolved, think about isolation mass around low-mass star

    if (m(i) > m_iso .and. opt(13)==1) then  ! mass growth quenched by pebble isolation
      m(i) = m(i) 
      m_x(i) = m_x(i) 
    else
      m(i) = m(i) + (0.5d0 + f_ice)*detm(i)
      m_x(i) = m_x(i) + f_ice*detm(i)
    end if 

    !### calculate the new density based on pebble accretion ###
    vol_tot = m_x(i)/rho_ice +(m(i)-m_x(i))/rho_sil
    rho(i) = m(i) /vol_tot
    ! ### update the radius due to mass increase ###
    rphys(i)=(3d0*m(i)/(4d0*rho(i)*PI))**(1.d0/3.d0)

    ! ###############################################################
    ! ###stop the integration if planet orbit is less than 0.1 AU ###
    ! ###############################################################
    !if (rdisk(i) <0.1d0 .and. asemi(i) < 0.1d0) stop
  end do


end subroutine pebbleaccretion 



!###############################################
! #### calculate the pebble isolation mass  ####
!###############################################
function pebble_iso(hgas,alpha,etagas,ts,mstar)  ! in Earth mass  
!  hgas       =  gas disk aspect ratio 
!  alpha      =  turbulent viscos alpha
!  etagas     =  headwind prefactor 
!  ts         =  stokes number 
  use physical_constant, only : K2
  implicit none
  real(double_precision),intent(in) :: hgas, alpha, etagas,ts, mstar
  real(double_precision) :: f_miso, pebble_iso
  f_miso = (hgas/5d-2)**3*(0.34d0*(log10(alpha)/log10(1d-3))**4+ 0.66d0)!* &
!     (1d0 - (2.5d0 +(-2d0*etagas/hgas**2) )/6d0)
  pebble_iso = 25d0*f_miso !+ alpha/(2d0*ts)/(4.76d-3/f_miso)
  pebble_iso = pebble_iso*(mstar/K2) ! consider the stellar mass dependence  

end function pebble_iso



!######################################################
! #### calculate from x, v to a, e, i  ####
!######################################################
subroutine xyztoaei(nbod,m,x,v,rdisk,asemi,ecc,inc)
!  nbod          = number of bodies
!  m          = planet mass (in solar masses * K2)
!  x      = coordinates (x,y,z) with respect to the central body [AU]
!  v      = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  asemi          =  semimajor axis
!  rdisk         =  distance
!  ecc          =  eccentricity
!  inc          =  inclination (radian)

  implicit none

  ! Input
  integer :: k, i
  integer, intent(in) :: nbod
  real(double_precision),intent(in) :: x(3,nbod), v(3,nbod), m(nbod)
  real(double_precision), intent(inout) :: rdisk(nbod), asemi(nbod), ecc(nbod), inc(nbod)
  real(double_precision) :: r2, v2, rv, r, zmb, semi, ecc2, h2, hz
  asemi(1) = 0.d0
  ecc(1) = 0.d0
  inc(1) = 0.d0
  do i = 2, nbod
    r2 = 0.d0
    v2 = 0.d0
    rv = 0.d0
    do  k = 1,3
      r2 = r2 + x(k, i)**2
      v2 = v2 + v(k, i)**2
      rv = rv + x(k, i)*v(k, i)
    end do
    r = r2**0.5d0 ! 3 dimensional distance 
    rdisk(i) = (x(1, i)**2 + x(2, i)**2)**0.5d0 ! xy dimensional distance 
    zmb = m(1) + m(nbod)
    semi = 2.0d0/r - v2/zmb
    semi = 1.0d0/semi
    asemi(i) = semi
    ecc2 = (1.0d0 - r/semi)**2 + rv**2/(semi*zmb)
    ecc(i) = ecc2**0.5d0
    hz=x(1, i)*v(2, i)-x(2, i)*v(1, i)
    h2=((x(2, i)*v(3, i)-x(3, i)*v(2, i))**2 &
  + (x(3, i)*v(1, i)-x(1, i)*v(3, i))**2+hz**2)**0.5d0
    if (abs(hz/h2) < 1.d0) then
      inc(i) = acos(hz/h2)
    else
      if (hz/h2 > 0d0) inc(i) = 0.d0
      if (hz/h2 < 0d0) inc(i) = PI ! in radian unit
    end if
  end do


  end subroutine xyztoaei





function eps_PA (tau, qp, ecc, inc, alpha, etagas, hgas, rp_to_a)
  integer, parameter      :: dp = kind(1.d0)
  real(dp), intent(in)    :: tau, qp, ecc, inc, alpha, etagas, hgas, rp_to_a
  real(dp)                :: vast, vturb, vcir, vecc, vz, detv
  real(dp)                :: arg_xy, arg_z, ftrans, q_c, vxy, fturb
  real(dp), parameter     :: PI = 4.d0*DATAN(1.d0)
  real(dp)                :: A2 = 0.32d0, A3 =0.39d0, ash=0.52d0, acir=5.7d0,&
                                    aturb = 0.33d0, aset = 0.5d0, ai =0.68d0, ae=0.76d0
  real(dp)                :: prob2d_set, prob2d_bal, prob2d, hpeb, heff, prob3d_set, prob3d_bal, prob3d
  real(dp)                :: eps_PA
  
  q_c = etagas**3/tau ! critical mass ratio separating the headwind and the shear regimes
  ! #### fitting factor to combine the settling and the balistic regime/settling factor ####
  vast = (qp/tau)**(1.d0/3.d0)
  vturb =  hgas*alpha**0.5/(1.0d0 + tau)**0.5d0
  ! turbulent velocity of particle = sqrt(D_gas/t_stop + t_corr)
  fturb = vast/(vast**2 + aturb*vturb**2)**0.5d0
  vcir = 1.0d0/(1.0d0 + acir*(qp/q_c))*etagas + ash*(qp*tau)**(1.d0/3d0)
  vecc = ae*ecc
  vz = abs(ai*inc)
  ! #### xy plane relative velocity ####
  vxy = max(vcir, vecc)
  ! #### relative velocity between the planet and the pebble ####
  detv = max(vcir, vecc + vz)
  ! #### ratio between systematic and random motion ####
  ! #### in our case, random velocity is only applied in z-direction ####
  !!!arg_xy = vxy**2/vast**2
  arg_xy = vxy**2/(vast**2 + aturb*vturb**2) !include vturb in xy plane
  arg_z = vz**2/(vast**2 + aturb*vturb**2)
  if(aset*(arg_xy  + arg_z) >40d0) then !too large velocity ratio, then f=0
    ftrans = 0d0
  else
    ftrans = (exp(-aset*arg_xy )*fturb) * (exp(-aset*arg_z )*fturb)
  end if

  ! #### 2D accretion probability in different regimes ####
  prob2d_set = A2*qp**0.5d0/etagas*tau**(-0.5d0)*detv**0.5d0
  prob2d_bal = rp_to_a/(2*PI*etagas*tau)*(detv**2 +2.d0*qp/rp_to_a)**0.5d0
  ! ####   2D accretion probability combined two regimes ####
  prob2d = ftrans*prob2d_set + (1.0d0 -ftrans)*prob2d_bal
  ! prob2d = 1.d0 - exp(-prob2d)
  ! #### 3D accretion probability in different regimes ####
  hpeb = (alpha/(alpha+tau))**0.5d0 *hgas
  heff = (hpeb**2 + PI*inc**2/2.0d0*(1.0d0-exp(-inc/2d0/hpeb)))**0.5d0
  prob3d_set = A3*qp/(etagas*heff)
  !prob3d_set = A3*q/(etagas*hpeb)!CWO debug
  prob3d_bal = 1.0d0/(4.d0*(2*PI)**0.5d0*etagas*heff*tau)*(2d0*qp/detv*rp_to_a + rp_to_a**2*detv)
  prob3d = prob3d_set*ftrans**2 + prob3d_bal*(1.0d0 - ftrans**2)
  eps_PA =  (prob2d**(-2) + prob3d**(-2))**(-0.5d0)

end function eps_PA


end module massgrowth
