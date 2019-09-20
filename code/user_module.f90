!******************************************************************************
! MODULE: user_module
!******************************************************************************
!
! DESCRIPTION: 
!> @brief Module that contain user defined function. This module can call other module and subroutine. 
!! The only public routine is mfo_user, that return an acceleration that 
!! mimic a random force that depend on what the user want to model.
!
!******************************************************************************

module user_module

  use types_numeriques
  use physical_constant

  implicit none
  
  private
  
  public :: mfo_user,  gasdisk
  
  contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!> @author 
!> John E. Chambers
!
!> @date 2 March 2001
!
! DESCRIPTION: 
!> @brief Applies an arbitrary force, defined by the user.
!!\n\n
!! If using with the symplectic algorithm MAL_MVS, the force should be
!! small compared with the force from the central object.
!! If using with the conservative Bulirsch-Stoer algorithm MAL_BS2, the
!! force should not be a function of the velocities.
!! \n\n
!! Code Units are in AU, days and solar mass * K2 (mass are in solar mass, but multiplied by K2 earlier in the code).
!
!> @note All coordinates and velocities must be with respect to central body
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
subroutine mfo_user (time,jcen,n_bodies,n_big_bodies,mass, phyradius, position, velocity,acceleration)
!  mass          = mass (in solar masses * K2)
!  position      = coordinates (x,y,z) with respect to the central body [AU]
!  velocity      = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  n_bodies      = current number of bodies (INCLUDING the central object)
!  n_big_bodies  =    "       "    " big bodies (ones that perturb everything else)
!  time          = current epoch [days]

  use physical_constant, only : PI, TWOPI, AU, MSUN, K2
  
  implicit none

  
  ! Input
  integer, intent(in) :: n_bodies !< [in] current number of bodies (1: star; 2-nbig: big bodies; nbig+1-nbod: small bodies)
  integer, intent(in) :: n_big_bodies !< [in] current number of big bodies (ones that perturb everything else)
  real(double_precision), intent(in) :: time !< [in] current epoch (days)
  real(double_precision), intent(in) :: jcen(3) !< [in] J2,J4,J6 for central body (units of RCEN^i for Ji)
  real(double_precision), intent(in) :: mass(n_bodies) !< [in] mass (in solar masses * K2)
  real(double_precision), intent(in) :: phyradius(n_bodies) !< planet radius 
  real(double_precision), intent(in) :: position(3,n_bodies)
  real(double_precision), intent(in) :: velocity(3,n_bodies)
  
  ! Output
  real(double_precision),intent(out) :: acceleration(3,n_bodies)
  
  !------------------------------------------------------------------------------ 
  !------Local-------

  ! loop integers
  integer :: planet
  real(double_precision) :: radii ! distace of the planet ( in AU)
  real(double_precision) :: acc_mig(3), acc_drag(3)
  real(double_precision) :: msmall
  
  do planet = 2, n_bodies
! planet's location in disk
    radii = (position(1,planet)**2 + position(2,planet)**2)**0.5 
! truncate the disk at 0.2 AU
    if (radii >= 0.2d0) then 
      if (planet <= n_big_bodies) then  ! for big bodies 
! calculate the acc from the migration 
        call  type12 (time,mass(1),planet,mass(planet),position(:,planet),velocity(:,planet),acc_mig(:))
        acceleration(1,planet) =  + acc_mig(1)
        acceleration(2,planet) =  + acc_mig(2)
        acceleration(3,planet) =  + acc_mig(3)
      else 
        acceleration(1,planet) =  0.0d0
        acceleration(2,planet) =  0.0d0
        acceleration(3,planet) =  0.0d0 
      end if
     else ! inside the cavity 
        acceleration(1,planet) =  0.0d0
        acceleration(2,planet) =  0.0d0
        acceleration(3,planet) =  0.0d0 
    end if
  end do
  
  !------------------------------------------------------------------------------
  return
end subroutine mfo_user



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!%%%   type I miration subroutine %%%        
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! Author: Beibei Liu
! Calculates the type I migration acceleration on a planet.
! Use the migration torque describled by Paarderkooper2011
! Use the eccentricity and inclination damping from Cresswell & Nelson 2008


subroutine type12 (t,mstar,num,mpl,x,v,acc)
!  mpl    = planet mass (in solar masses * K2)
!  mstar  = star mass (in solar masses * K2)
!  x      = coordinates (x,y,z) with respect to the central body [AU]
!  v      = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  num    = current number of bodies
!  t      = current epoch [days]

  use physical_constant, only : PI, TWOPI, AU, MSUN, K2, EARTH_MASS
  use mercury_globals, only : opt,alpha_t, kap, mdot_gas0, L_s, Rdisk0, alpha, tau_s,t0_dep,tau_dep
  use orbital_elements

  implicit none

  
  ! Input
  integer, intent(in) :: num
  real(double_precision), intent(in) :: t, mpl, mstar
  real(double_precision), intent(in) :: x(3), v(3)
  real(double_precision), intent(out) :: acc(3)
  real(double_precision) :: asemi(3), aecc(3), ainc
  real(double_precision) :: taumig, tauecc, tauinc
  real(double_precision) :: pecc, hgas, etagas, rhogas
  real(double_precision) :: rtrans, siggas_vis, siggas_irr,f1,f2,f_tot,fs
  real(double_precision) :: tauwave, rdisk, siggas, Omega, m_gap, m_planet,delta_m
  real(double_precision) :: r2, v2, rv, r
  real(double_precision) :: ecc, inc
  real(double_precision) ::peri,long,node,lmean
  integer k
! planet's location in the disk 
  rdisk = (x(1)**2 + x(2)**2 )**0.5
!  ###########################################
!  ### get disk properties:surface density ###
!  ###########################################
  call gasdisk (t,rdisk,x(3),mstar,alpha,kap,hgas,etagas,siggas,siggas_vis,siggas_irr,rhogas,rtrans)

! angular velocity at planet's location
  Omega = (mstar/rdisk**3)**0.5d0 ! in day^{-1}

  r2 = 0.d0
  v2 = 0.d0
  rv = 0.d0
  do  k = 1,3
    r2 = r2 + x(k)**2
    v2 = v2 + v(k)**2
    rv = rv + x(k)*v(k)
  end do
  r = r2**0.5d0


!  ########################################
!  ### from xyz to ecc and inc(radian)  ###
!  ########################################
  call mco_x2el (mstar+mpl,x(1),x(2),x(3),v(1),v(2),v(3),peri,ecc,inc,long,node,lmean)

! wave timescale
  tauwave = (mstar/mpl)*(mstar/(siggas*rdisk**2))*hgas**4/Omega
! eccentricity damping timescale
  tauecc = tauwave/0.78d0*(1d0 - 0.14d0*(ecc/hgas)**2 &
  + 6d-2*(ecc/hgas)**3 + 0.18d0*(ecc/hgas)*(inc/hgas)**2)
! inclination damping timescale
  tauinc = tauwave/0.544d0*(1d0 - 0.3d0*(inc/hgas)**2 &
  + 0.24d0*(inc/hgas)**3 + 0.14d0*(ecc/hgas)**2*(inc/hgas))
  pecc = ( 1d0 + (ecc/(2.25d0*hgas))**1.2d0 + &
  (ecc/(2.84d0*hgas))**6  )/(1.0d0 -(ecc/(2.02d0*hgas))**4 )
!  taumig = 2d0*tauwave/(2.7d0 - 1.1d0*index_p)/hgas**2*(pecc + pecc/abs(pecc) &
!  *(0.07d0*(inc/hgas) + 0.085d0*(inc/hgas)**4 -0.08d0*(ecc/hgas)*(inc/hgas)**2 ))
! calculate the type I migration prefactor based on two-component disk structures 
  call type1_prefactor (mpl,mstar,rdisk,siggas_vis,siggas_irr,siggas,hgas, rtrans,f1)

! only conside type I migration for all mass range 
  if (opt(12) == 1 ) then 
    f_tot = f1 
  end if 

! both conside type I & II 
  if (opt(12) == 2 .or. opt(12) == 3) then      
!  type II based on Kanagawa2018
    if (opt(12) == 2)   f2 = -1.d0 
!  type II  becomes infinite slow
    if (opt(12) == 3) f2 = -1d-10
! gap opening mass 
    m_gap = gap_opening(hgas,alpha_t,mstar) ! in earth mass 
    m_planet = mpl/K2/EARTH_MASS ! planet mass in earth mass !!!! unsolved, think about stellar-mass dependence 
    if (m_planet <m_gap) then 
      f_tot = f1
    else
    !fs  = 1.d0/(1.d0  + (m_planet/m_gap)**4)
    !f_tot  = f1*fs + f2*(1-fs)
    !f_tot = f_tot/(1.d0 + (m_planet/m_gap)**2)
      delta_m = m_gap/5.d0
      fs = exp(-(m_planet-m_gap)/delta_m)
      f_tot  = f1*fs + f2*(1.d0-fs)/(m_planet/m_gap)**2
    end if
  end if

!  print*, opt(12),alpha_t, m_gap, mpl,mpl/K2,mpl/K2/EARTH_MASS, mstar,K2

! migration  timescale f_tot>0, outward migration; f_tot<0, inward migration 
  taumig = 2d0*tauwave/f_tot/hgas**2
!  print*,'f1',t/365.25d0,f1,taumig/365.25
!  print*,'check',mpl/mstar,mstar/(siggas*rdisk**2),hgas,6.28/Omega/365.25
  
!  implement the acceralation terms due to disk torques 
  do  k = 1,3
       if (opt(12) /= 0) then  ! include migration 
         asemi(k) = v(k)/taumig
       else
         asemi(k) = 0d0
       end if 
       aecc(k) = -2*rv*x(k)/r**2/tauecc
       acc(k) = asemi(k) + aecc(k)
  end do
  ainc = -v(3)/tauinc
  acc(3) = acc(3) + ainc
  

  return
end subroutine type12



!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!%%%  Gas drag on planetesimal %%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
! Author: Beibei Liu
! Calculates the gas drag acceleration on a body (mainly planetesimal).
! The drag force describled in Mandell 2007


subroutine gasdrag (t,mstar,num, mpl,rphypl, x,v,acc)
!  mpl          = planet mass (in solar masses * K2)
!  rhypl          = planet radius 
!  mstar          = star mass (in solar masses * K2)
!  x      = coordinates (x,y,z) with respect to the central body [AU]
!  v      = velocities (vx,vy,vz) with respect to the central body [AU/day]
!  num      = current number of bodies
!  t          = current epoch [days]

  use physical_constant, only : PI, AU, MSUN, K2
  use mercury_globals, only : alpha_t, kap, mdot_gas0, L_s, Rdisk0, alpha, t0_dep,tau_dep
  implicit none

  
  ! Input
  integer, intent(in) :: num
  real(double_precision), intent(in) :: t, mpl,rphypl, mstar
  real(double_precision), intent(in) :: x(3), v(3)
  real(double_precision), intent(out) :: acc(3)
  real(double_precision) :: hgas, rdisk, siggas, Omega, vk
  real(double_precision) :: tmp, vg_x, vg_y, cosf, sinf
  real(double_precision) :: dv_x,dv_y,dv_z, dv, vk_x,vk_y,vg_over_vk
  real(double_precision) :: Cgas = 0.5d0 ! gas drag coefficient
  real(double_precision) :: rhopl  ! volumn density of the planet/planetesimal (in g/cm^3) 
  real(double_precision) :: rhogas, rhogas_mid ! volumn density of the disk gas (in g/cm^3) 
  real(double_precision) :: Rpl ! diamater of the planet/planetesimal  (! in AU)
  real(double_precision) :: etagas ! headwind factor \eta
  real(double_precision) :: rtrans, siggas_vis, siggas_irr ! headwind factor \eta
  integer :: k
!  planet's location in the disk r_p = (r,z) 
  rdisk = (x(1)**2 + x(2)**2 )**0.5d0

  !!!need to get a modt_gas(t)
!  call gasdisk to calculate the surface density and apsect ratio, headwind prefactor at planet's location  
  call gasdisk (t,rdisk,x(3),mstar,alpha,kap,hgas,etagas,siggas,siggas_vis,siggas_irr,rhogas,rtrans)
!  planet's radius 
  rhopl = 3d0*mpl/4d0/pi/rphypl**3/( AU * AU * AU * K2 / MSUN)
  !Rpl = (3.d0*mpl/K2*MSUN/PI/rhopl/4.d0)**(1d0/3d0)/AU ! (in AU)
  Rpl = rphypl ! (in AU)
  !print*, 'plt R',Rpl
!  gas drag prefactor
  tmp = - 3.d0/8.d0*Cgas*(rhogas/rhopl)/Rpl

!  angular velocity and Kepler velocity at the body's location  
  Omega = (mstar/rdisk**3)**0.5d0 ! in day^{-1}
  vk = (mstar/rdisk)**0.5d0 ! in day^{-1}
!  Keplerian velocity at planet's location (x,y)
  cosf  = x(1)/rdisk
  sinf  = x(2)/rdisk
  vk_x = - vk*sinf
  vk_y =   vk*cosf
! the vg(z,r) with respect to vk(r) Takeuchi-Lin 2002 eq(7)
  vg_over_vk = (x(1)**2 + x(2)**2 + x(3)**2)**0.5d0/rdisk &
  *(1.d0 + 0.5d0*hgas**2*(0-2.d0 -(2d0*1-1.0d0) &
  !*(1.d0 + 0.5d0*hgas**2*(index_p+index_q-2.d0 -(2d0*index_q-1.0d0) &
  /2.0d0*(x(3)/hgas/rdisk)**2))
! gas velocity at planet's location (x,y)
  vg_x = vg_over_vk*vk_x
  vg_y = vg_over_vk*vk_y
! relative velocity between gas and planet at planet's location (x,y,z)
  dv_x = v(1) - vg_x
  dv_y = v(2) - vg_y
  dv_z = v(3) 
  dv = (dv_x**2 + dv_y**2 + dv_z**2)**0.5d0
  acc(1) = tmp*dv_x*dv
  acc(2) = tmp*dv_y*dv
  acc(3) = tmp*dv_z*dv

end subroutine gasdrag





function fsmooth(r, rtrans)
  implicit none
  real(double_precision), intent(in) :: rtrans, r
  real(double_precision) :: fsmooth 
  fsmooth =  1./(1. + (r/rtrans)**4)
end function fsmooth


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!%%% two-component disk structure %%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! Author: Beibei Liu
! Calculates gas disk profile with inner viscously heated + outer stellar irradiation
! based on Liu 2019's population synthesis model 

subroutine gasdisk (t,r,z,mstar,alpha,kap,hgas,etagas,siggas,siggas_vis,siggas_irr,rhogas,rtrans)
!  t                     = current epoch [days]
!  r                     = disk radius [AU]
!  z                     = vertical distance [AU]
!  t0_dep                = onset time of disk disposal [Myr]
!  tau_dep                =  disk disposal timescale [Myr]
!  mdot_gas              = gas accretion rate at t [Msun/yr]
!  siggas                = gas surface desity at t  [K*MSUN/AU^2]
!  siggas_vis/irr        = gas surface desity at viscous/irradiation zone  [K*MSUN/AU^2]
!  rhogas_mid            = gas volumn desity in the midplane at r [g/cm^3]
!  rhogas                = gas volumn desity at (r,z) [g/cm^3]
!  etagas                = headwind prefactor 
!  hgas                  = gas disk aspect ratio 
!  rtrans                = transition radius between two regions  

  use physical_constant, only : PI, TWOPI, AU, MSUN, K2
  use mercury_globals, only : opt, L_s,mdot_gas0, t0_dep,tau_dep 

  implicit none
  real(double_precision), intent(in) :: t, r, z, mstar,alpha,kap
  real(double_precision), intent(out) :: etagas, siggas, rhogas,hgas
  real(double_precision), intent(out) :: rtrans, siggas_vis, siggas_irr
  real(double_precision) :: p, q, siggas0_vis, hgas0_vis, hgas_vis, T0_vis, T_vis
  real(double_precision) :: siggas0_irr,  hgas0_irr, hgas_irr, T0_irr, T_irr
  real(double_precision) :: fs, Temp,  rhogas_mid, etagas0
  real(double_precision) :: Omega,v_0,ts_0,mdot_gas
  real(double_precision) :: M_s, t0, td

!!! this will change in future !!!
  M_s = 1.0
!###gas viscous timescale ###
  p = - 15./14
  q = 2./7
  Omega = (mstar/r**3)**0.5d0 ! in day^{-1}
  hgas0_irr = 2.45e-2*(L_s/1.0)**(1./7)*(M_s/1.0)**(-4./7)
  hgas_irr = hgas0_irr*r**q
  v_0 =alpha*hgas_irr**2*Omega
  ts_0 = 1./(3d0*(2+p)**2)/v_0 ! day
  !########################################################
  !####### consider gas accretion rate time-dependent #####
  !########################################################
  if ( opt(11) ==1 ) then ! time evolution of mdot_gas
    !mdot_gas = mdot_gas0*(t/ts_0 +1)**(-(2.5+p)/(2+p))
    t0 = t0_dep*365.25d0*1d+6 ! day
    td = tau_dep*365.25d0*1d+6 ! day 
    if (t<t0)  then 
      mdot_gas =  mdot_gas0
    else        
      mdot_gas = mdot_gas0*exp(-(t-t0)/td)
    end if 
  else 
    mdot_gas =  mdot_gas0
  end if 


  !### viscous heated region ####
  p = - 0.375 ! power law index for surface density 
  q = -1./16 ! power law index for gas aspect ratio 
  ! surface density at 1AU and as a function of r  in g/cm^2
  siggas0_vis = 740*(mdot_gas/1d-8)**(1./2)*(M_s/1.d0)**(1./8)*(alpha/1d-3)**(-3./4)*(kap/1d-2)**(-1./4.)
  siggas_vis = siggas0_vis*r**p
  ! apsect ratio at 1AU and as a function of r
  hgas0_vis = 4.5e-2*(mdot_gas/1d-8)**(1./4)*(M_s/1.0)**(-5./16)*(alpha/1d-3)**(-1./8)*(kap/1d-2)**(1./8)
  hgas_vis = hgas0_vis*r**q
  ! temperature  at 1AU and as a function of r  in K
  T0_vis = 500.*(mdot_gas/1d-8)**(1./2)*(M_s/1.0)**(3./8)*(alpha/1d-3)**(-1./4)*(kap/1d-2)**(1./4)
  T_vis = T0_vis*r**(2*q-1)


!### stellar irradition region ####
  p = - 15./14
  q = 2./7
  ! surface density at 1AU and as a function of r  in g/cm^2
  siggas0_irr = 2500*(mdot_gas/1d-8)*(L_s/1.0)**(-2./7)*(M_s/1.0)**(9./14)*(alpha/1d-3)**(-1)
  siggas_irr = siggas0_irr*r**p
  ! apsect ratio at 1AU and as a function of r
  hgas0_irr = 2.45e-2*(L_s/1.0)**(1./7)*(M_s/1.0)**(-4./7)
  hgas_irr = hgas0_irr*r**q
  ! temperature  at 1AU and as a function of r  in K
  T0_irr = 150.*(L_s/1.0)**(2./7)*(M_s/1.0)**(-1./7)
  T_irr = T0_irr*r**(2*q-1)
  ! transition radius for two regions
  rtrans = (500./150.)**(56./39)*(mdot_gas/1d-8)**(28./39)*(M_s/1.0)**(29./39)& 
  *(alpha/1d-3)**(-14./39)*(L_s/1.0)**(-16./39)*(kap/1d-2)**(14./39)

!  combine of two regions
  if (opt(10) == 1) then ! include the inner viscously heated region
    fs= fsmooth(r, rtrans)
  else ! only irradiation
    fs = 0.0
  end if


  siggas = fs*siggas_vis + (1.0 - fs)*siggas_irr
  Temp = fs*T_vis + (1.0 - fs)*T_irr
  hgas = fs*hgas_vis + (1.0 - fs)*hgas_irr
  ! midplane volume density  in g/cm^3
  rhogas_mid = siggas/TWOPI**0.5d0/hgas/AU
  ! gas volume density rho(r,z)
  rhogas = rhogas_mid*exp(-z**2/(hgas*r)**2) ! in g/cm^3
  ! gas headwind prefactor at 1AU
  etagas0 = (2.d0 - p - q)/2.0d0
  ! gas headwind prefactor
  etagas = etagas0*hgas**2

  ! gas surface density at r  [in solar mass*K2]
  siggas = siggas/MSUN*AU**2*K2 ! in solar mass*K2
  siggas_vis = siggas_vis/MSUN*AU**2*K2 ! in solar mass*K2
  siggas_irr = siggas_irr/MSUN*AU**2*K2 ! in solar mass*K2

end subroutine gasdisk





!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
!%%% type 1 migration prefactor   %%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
! Author: Beibei Liu
! Calculates type I migration prefactor two-component disk structure  (inner viscously heated + outer stellar irradiation)
! based on Paardekoer2011 

subroutine type1_prefactor (mpl,mstar,r,siggas_vis,siggas_irr,siggas,hgas, rtrans,f1)
!  r                         =distance in the disk 
!  siggas_vis                = gas surface desity in viscously heated region   [solar mass*K2]
!  siggas_irr                = gas surface desity in stellar irradiation  [solar mass*K2]
!  siggas                    = gas surface desity  [solar mass*K2]
!  hgas                      = gas disk aspect ratio 
!  hgas_vis                  = gas disk aspect ratio in viscously heated region 
!  hgas_irr                  = gas disk aspect ratio in stellar irradiation region 
!  rtrans                    = transition radius  

  use physical_constant, only : PI, TWOPI, AU, MSUN, K2
  use mercury_globals, only : alpha_t
      
  implicit none
  real(double_precision), intent(in) :: r, rtrans, siggas_vis, siggas_irr, siggas
  real(double_precision), intent(in) :: hgas,mpl,mstar
  real(double_precision), intent(out) :: f1

  character(len=5) :: opt_corotation_damping 
  real(double_precision) :: p,pp,ppp,fp,pps,gpp,ppps,akppp
  real(double_precision) :: flb_vis,fhb_vis,fhe_vis,fcb_vis, fce_vis
  real(double_precision) :: flb_irr,fhb_irr,fhe_irr,fcb_irr, fce_irr
  real(double_precision) :: f1_vis,f1_irr, fs,f_ecc,f_inc


  opt_corotation_damping = 'False' 

  ! calculate the type 1 migration prefactor 
  ! saturation factor 
  p = 2./3.*(mpl/mstar)**0.75/((2.0*PI*alpha_t)**0.5*hgas**1.75)
  pp = p/(8.0/(45.0*PI))**0.5
  ppp = p/(28.0/(45.0*PI))**0.5
  fp = 1./(1. + (p/1.3)**2)
  pps = 1./(1. + pp**4)
  gpp = pps*(16.*pp**1.5/25.) + (1. - pps)*(1.-9.*pp**(-2.67)/25.)
  ppps = 1./(1. + ppp**4)
  akppp = ppps*(16.*ppp**1.5/25.) + (1. - ppps)*(1.-9.*ppp**(-2.67)/25.)
  !ftot
  !factors of each torques in viscous region: -0.375, -1.125       
  flb_vis = -4.375
  fhb_vis = 1.2375
  fhe_vis = 5.5018
  fcb_vis = 0.7875
  fce_vis = 1.17
  !factors of each torques in irr region: -3/7, -15/14             
  flb_irr = -3.08
  fhb_irr = 0.47
  fhe_irr = 0
  fcb_irr = 0.3
  fce_irr = 0.0
  if (opt_corotation_damping == 'True')   then 
  ! consider the corotation torque attenuation due to high eccentricity and inclination Bitsch&Kley2010,Fendyke&Nelson2014 
    f_ecc = 1 !exp(-ecc)  need to calculate ecc and inc!!! future 
    f_inc = 1 !
    f1_vis = flb_vis + (fhb_vis + fhe_vis*fp)*fp*gpp + (1 - akppp)*(fcb_vis + fce_vis)*f_ecc
    f1_irr = flb_irr + (fhb_irr + fhe_irr*fp)*fp*gpp + (1 - akppp)*(fcb_irr + fce_irr)*f_ecc
  else
    f1_vis = flb_vis + (fhb_vis + fhe_vis*fp)*fp*gpp + (1 - akppp)*(fcb_vis + fce_vis)
    f1_irr = flb_irr + (fhb_irr + fhe_irr*fp)*fp*gpp + (1 - akppp)*(fcb_irr + fce_irr)
  end if        
  fs = fsmooth(r,rtrans)
  f1 = (f1_vis* siggas_vis*fs + f1_irr*siggas_irr*(1-fs))/siggas 
  !!! ensure it is not zero
  if (abs(f1)<1d-10) Then 
    f1 = 1d-10
  end if
end subroutine type1_prefactor 



!###############################################
! #### calculate the gap opening mass  ####
!###############################################
function gap_opening(hgas,alpha,mstar)  ! in Earth mass  
!  hgas       =  gas disk aspect ratio 
!  alpha      =  turbulent viscos alpha
  use physical_constant, only : K2
  implicit none
  real(double_precision),intent(in) :: hgas, alpha, mstar
  real(double_precision) :: gap_opening
  gap_opening = 30d0*(hgas/5d-2)**2.5*(alpha/1d-3)**0.5*(mstar/K2)

end function gap_opening  ! in Earth mass

end module user_module

