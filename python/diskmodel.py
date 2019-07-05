# -*- coding: utf-8 -*- 
#!/usr/bin/env python
# disk model from Gauard&Lin2007, Liu2015
# type I migration torque from Paardekooper2011
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
import matplotlib
from  astrounit import *




### scale functions of disk accretion rate/luminosity and disk radius on the stellar mass ###
def mdot_gas_0(mdot_0, M_s):
    res = mdot_0*(M_s/1.0)**1.8 # [Msun/yr] Hartmann 1998 
    return res # [Msun/yr]

def Lum(M_s):
    res = (M_s/1.0)**2 # any reference
    return res # [Lsun]

def Rdisk_0(R0, M_s):
    res = R0*(M_s/1.0)**(0.0) # need to find some reference 
    return res # [AU]

# smooth function to combined viscous and irradiation disks
def fs(r,r_trans,opt_vis):
    if opt_vis == True:
        res = 1./(1. + (r/r_trans)**4)
    else: 
        res = 0.0
    return res 


# viscous disk: surface density, temperature, scale height
def siggas_vis(r,mdot_gas,M_s,alpha,kap):
    sig_vis0 = 740*(mdot_gas/1e-8)**(1./2)*(M_s/1.0)**(1./8)*(alpha/1e-3)**(-3./4)*(kap/0.01)**(-1./4.) # [g/cm^2]
    res = sig_vis0*r**-0.375
    return res # [g/cm^2]
def hgas_vis(r,mdot_gas,M_s,alpha,kap):
#    hoverr0 = (3./64./pi**2/gamma**0.5*(Rg*gamma/mu)**3*1e-2/sigma*(flux0*Msun/year)**2*(G*Msun)**(-5./2.)*AU**(-0.5)/alpha0)**(1./8.)
    hgas_vis0 = 4.5e-2*(mdot_gas/1e-8)**(1./4)*(M_s/1.0)**(-5./16)*(alpha/1e-3)**(-1./8)*(kap/0.01)**(1./8)
    res = hgas_vis0*r**(-1./16)
    return res
def T_vis(r,mdot_gas,M_s,alpha,kap):
    T_vis0 = 500.*(mdot_gas/1e-8)**(1./2)*(M_s/1.0)**(3./8)*(alpha/1e-3)**(-1./4)*(kap/0.01)**(1./4)
    res = T_vis0*r**-1.125
    return res #[K]
    
# irradition disk: surface density, temperature, scale height
def siggas_irr(r,mdot_gas,L_s,M_s,alpha):
    sig_irr0 = 2500*(mdot_gas/1e-8)*(L_s/1.0)**(-2./7)*(M_s/1.0)**(9./14)*(alpha/1e-3)**(-1) # [g/cm^2]
    res = sig_irr0*r**(-15./14)
    #sig_irr0 = 1340*(mdot_gas/1e-8)*(L_s/1.0)**(-1./4)*(M_s/1.0)**(1./2)*(alpha/1e-3)**(-1) # [g/cm^2]
    #res = sig_irr0*r**(-1)
    return res # [g/cm^2]
def hgas_irr(r,L_s,M_s):
    hgas_irr0 = 2.45e-2*(L_s/1.0)**(1./7)*(M_s/1.0)**(-4./7)
    #hgas_irr0 = 2.6e-2*(L_s/1.0)**(1./7)*(M_s/1.0)**(-4./7)
    #hgas_irr0 = 2.2e-2*(L_s/1.0)**(1./7)*(M_s/1.0)**(-4./7)
    res = hgas_irr0*r**(2./7)
    #hgas_irr0 = 3.35e-2*(L_s/1.0)**(1./8)*(M_s/1.0)**(-1./2)
    #res = hgas_irr0*r**(1./4)
    return res
def T_irr(r,L_s,M_s):
    T_irr0 = 150.*(L_s/1.0)**(2./7)*(M_s/1.0)**(-1./7)
    #T_irr0 = 170.*(L_s/1.0)**(2./7)*(M_s/1.0)**(-1./7)
    #T_irr0 = 125.*(L_s/1.0)**(2./7)*(M_s/1.0)**(-1./7)
    res = T_irr0*r**(-3./7) 
    #T_irr0 = 280.*(L_s/1.0)**(1./4)
    #res = T_irr0*r**(-1./2) 
    return res #[K]

def rtrans(mdot_gas,L_s,M_s,alpha,kap,opt_vis):
    if opt_vis == True:
        res = (500./150.)**(56./39)*(mdot_gas/1e-8)**(28./39)*(M_s/1.0)**(29./39)*(alpha/1e-3)**(-14./39)*(L_s/1.0)**(-16./39)*(kap/0.01)**(14./39)   # transit location 
     #   res = (500./170.)**(56./39)*(mdot_gas/1e-8)**(28./39)*(M_s/1.0)**(29./39)*(alpha/1e-3)**(-14./39)*(L_s/1.0)**(-16./39)*(kap/0.01)**(14./39)   # transit location 
        #res = (500./125.)**(56./39)*(mdot_gas/1e-8)**(28./39)*(M_s/1.0)**(29./39)*(alpha/1e-3)**(-14./39)*(L_s/1.0)**(-16./39)*(kap/0.01)**(14./39)   # transit location 
        #res = (500./280.)**(8./5)*(mdot_gas/1e-8)**(4./5)*(M_s/1.0)**(3./5)*(alpha/1e-3)**(-2./5)*(L_s/1.0)**(-2./5)*(kap/0.01)**(2./5)   # transit location 
    else:
        res = 0.0
    return res #[AU]


def rsnow(mdot_gas,L_s,M_s,alpha,kap,opt_vis):
    rsnow_vis = (50./17)**(8./9)*(mdot_gas/1e-8)**(4./9)*(M_s/1.0)**(1./3)*(alpha/1e-3)**(-2./9)*(kap/0.01)**(2./9)
    rsnow_irr = (15./17)**(7./3)*(M_s/1.0)**(-1./3)*(L_s/1.0)**(2./3)
    if opt_vis == True:
        res = max(rsnow_vis,rsnow_irr)
    else:
        res = rsnow_irr
    return res #[AU]



# surface density 
def siggas(r,mdot_gas,L_s,M_s,alpha,kap, opt_vis):
    r_trans = rtrans(mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    fsmooth = fs(r,r_trans,opt_vis)
    if opt_vis == True:
        res = fsmooth*siggas_vis(r,mdot_gas,M_s,alpha,kap) + (1-fsmooth)*siggas_irr(r,mdot_gas,L_s,M_s,alpha) 
    else:
        res = siggas_irr(r,mdot_gas,L_s,M_s,alpha)
    return res # [g/cm^2]

# apsect ratio 
def hgas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis):
    r_trans = rtrans(mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    fsmooth = fs(r,r_trans,opt_vis)
    if opt_vis == True:
        res = fsmooth*hgas_vis(r,mdot_gas,M_s,alpha,kap) + (1-fsmooth)*hgas_irr(r,L_s,M_s) 
    else:
        res = hgas_irr(r,L_s,M_s)   
    return res

#  headwind pressure prefactor: from Bitsch2018 Eq(9) [some errors]  
def eta_gas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis):
    r_trans = rtrans(mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    #eta_irr = -0.5*hgas_irr(1,L_s,M_s)**2*r**(2*2./7)*(2*2./7-15./14-2)
    #eta_vis = -0.5*hgas_vis(1,mdot_gas,M_s,alpha,kap)**2*r**(-2*.1/16)*(-2*1./16-3./8-2)
    eta_irr = 0.5*hgas_irr(r,L_s,M_s)**2*(-2./7+15./14+2)
    eta_vis = 0.5*hgas_vis(r,mdot_gas,M_s,alpha,kap)**2*(1./16+3./8+2)
    if opt_vis == True:
        fsmooth = fs(r,r_trans,opt_vis)
        res = fsmooth*eta_vis + (1-fsmooth)*eta_irr
    else:
        res = eta_irr
    return res

# temperature
def Tgas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis):
    r_trans = rtrans(mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    if opt_vis == True:
        fsmooth = fs(r,r_trans,opt_vis)
        res = fsmooth*T_vis(r,mdot_gas,M_s,alpha,kap) + (1-fsmooth)*T_irr(r,L_s,M_s)   
    else:
        res = T_irr(r,L_s,M_s)   
    return res  #[K]



# saturation pressure for water 
def Psat(T):
    ## cgs unit from Lichtenegger1991
    #A = 1.14e+13 
    #B = 6062
    ## SI unit from ?
    A = 6.034e+11 
    B = 5938.0
    res = A*np.exp(-B/T) 
    return res 


# water vapor pressure
def P_vapor(r,mdot_gas,L_s,M_s,alpha,kap, opt_vis,ceta):
    # in cgs unit ##
    '''
    T = Tgas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    sig = siggas(r,mdot_gas,L_s,M_s,alpha,kap, opt_vis) 
    h = hgas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    rhogas = sig/((2*np.pi)**0.5*h*(r*AU)) 
    Z_H2O = 1e-2 # metallicity of water 
    N_a = 6.022e+23 # avogadro constant [mol^-1]
    m_H2O = 18. # mass of water per mol [g/mol^{-1}]
    m_H2O = m_H2O/N_a # mean molecule weight of water [g]
    kB = 1.38e-16 # Boltzmann constant [ergK^-1]
    rho_vapor = rhogas*Z_H2O/m_H2O
    res = rho_vapor*kB*T
    '''

    T = Tgas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    sig = 10*siggas(r,mdot_gas,L_s,M_s,alpha,kap, opt_vis) # in SI unit 
    h = hgas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    rhogas = sig/((2*np.pi)**0.5*h*(r*AU/100.)) # AU/100 is in SI unit
    #Z_H2O = 1e-2 # metallicity of water 
    Z_H2O = ceta # metallicity of water 
    N_a = 6.022e+23 # avogadro constant [mol^-1]
    m_H2O = 1.8e-2 # mass of water per mol [kg/mol^{-1}]
    m_H2O = m_H2O/N_a # mean molecule weight of water [kg]
    kB = 1.38e-23 # Boltzmann constant [JK^-1]
    rho_vapor = rhogas*Z_H2O/m_H2O
    res = rho_vapor*kB*T
    return res 

def rsnow_new(r,mdot_gas,L_s,M_s,alpha,kap, opt_vis,ceta):
    T = Tgas(r,mdot_gas,L_s,M_s,alpha,kap,opt_vis)
    P1 = Psat(T)
    P2 = P_vapor(r,mdot_gas,L_s,M_s,alpha,kap, opt_vis,ceta)
    nlist = np.argmin( abs(P1/P2-1)) # when these two equals 
    res = r[nlist]
    return res #[AU]







# migration torque from Paardekooper2011
def torque(q,r,mdot_gas,L_s,M_s,alpha_v,alpha_d,kap,opt_vis):
    r_trans = rtrans(mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
    p = 2./3.*q**0.75/((2.0*np.pi*alpha_d)**0.5*hgas(r,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)**1.75)
    pp = p/(8.0/(45.0*np.pi))**0.5
    ppp = p/(28.0/(45.0*np.pi))**0.5
    fp = 1./(1. + (p/1.3)**2)
    pps = 1./(1. + pp**4)
    gpp = pps*(16.*pp**1.5/25.) + (1. - pps)*(1.-9.*pp**(-2.67)/25.)
    ppps = 1./(1. + ppp**4)
    akppp = ppps*(16.*ppp**1.5/25.) + (1. - ppps)*(1.-9.*ppp**(-2.67)/25.)
    # ftot
    #    factors of each torques in viscous region: -0.375, -1.125       
    flb_vis = -4.375
    fhb_vis = 1.2375
    fhe_vis = 5.5018
    fcb_vis = 0.7875
    fce_vis = 1.17
    #    factors of each torques in irr region: -3/7, -15/14             
    flb_irr = -3.08
    fhb_irr = 0.47
    fhe_irr = 0
    fcb_irr = 0.3
    fce_irr = 0.0
    ftot_vis = flb_vis + (fhb_vis + fhe_vis*fp)*fp*gpp + (1 - akppp)*(fcb_vis + fce_vis)
    ftot_irr = flb_irr + (fhb_irr + fhe_irr*fp)*fp*gpp + (1 - akppp)*(fcb_irr + fce_irr)
    # get torque
    factor = (r*AU)*(G*M_s*Msun)*(q/hgas(r,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis))**2
    torq_vis = ftot_vis*factor*siggas_vis(r,mdot_gas,M_s,alpha_v,kap)
    torq_irr = ftot_irr*factor*siggas_irr(r,mdot_gas,L_s,M_s,alpha_v)
    if opt_vis == True:
        fsmooth = fs(r,r_trans,opt_vis)
        torq = torq_vis*fsmooth + torq_irr*(1-fsmooth) # [cm g s] unit
    else:
        torq = torq_irr # [cm g s] unit
    # here choose only inward migration torque torque   
    #torq = ftot_irr*factor*siggas(r,mdot_gas,L_s,M_s,alpha_v,kap, opt_vis)
    f = torq/factor/siggas(r,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis) # normalized torque
    return torq, f

def tmig_1(mp,a,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis=True):
    """Computes the type I migration time for a given planetary mass at a given
    semi-major axis.
    
    mp: mass given in earth masses
    a : semi-major axis in AU
    mdot_gas: the gas accretion rate
    L_s: stellar luminosity
    M_s: stellar mass
    alpha_v: viscous alpha
    kap: disk opacity coefficient
    opt_vis: boolean statement for the inclusion of a modulation factor"""
    
    sigma_g = siggas(a,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)*(AU**2/Msun)
    hg      = hgas(a,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
    
    t = (M_s/(mp/msuntome))*(M_s/(sigma_g*a**2))*hg**2*(a**3/(4*np.pi**2*M_s))**(0.5)
    
    return t #[yr]

def tpeb(mp,a,mdot_gas,mdot_peb,L_s,M_s,alpha_v,kap,st,opt_vis=True):
    """Computes the pebble accretion time for a given planetary mass at a given
    semi-major axis. If the planets satisfies Racc>Hpeb we have 2D accretion
    while we have 3D if Racc<Hpeb.
    
    mp: mass given in earth masses
    a : semi-major axis in AU
    mdot_gas: the gas accretion rate
    mdot_peb: the pebble accretion rate
    L_s: stellar luminosity
    M_s: stellar mass
    alpha_v: viscous alpha
    kap: disk opacity coefficient
    st: stokes number of the pebbles in the disk
    opt_vis: boolean statement for the inclusion of a modulation factor"""
    
    hg = hgas(a,mdot_gas,L_s,M_s,alpha_v,kap,True)
    hpeb = np.sqrt(alpha_v/(alpha_v+st))*hg
    
    Mp = mp/msuntome
    
    Rh      = a*(Mp/(3*M_s))**(1/3)
    Racc    = Rh*st**(1/3)
    acc_2d  = Racc>hpeb*a
    acc_3d  = Racc<hpeb*a
    
    #We compute the pebble accretion efficiency
    eta = eta_gas(a,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
    
    eps_2d = 3e-3*(mp/1e-2)**(2/3)*(st/1e-2)**(-1/3)*(eta[acc_2d]/3e-3)**(-1)
    eps_3d = 4e-3*(mp/1e-2)*M_s**(-1)*(hpeb[acc_3d]/3e-3)**(-1)*(eta[acc_3d]/3e-3)**(-1)
    
    t         = np.zeros(len(a))
    t[acc_2d] = mp/(eps_2d*mdot_peb)
    t[acc_3d] = mp/(eps_3d*mdot_peb)
    
    return t#[yr]

def tpeb_alt(mp,a,mdot_gas,mdot_peb,L_s,M_s,alpha_v,kap,st,opt_vis=True):
    """Computes the pebble accretion time for a given planetary mass at a given
    semi-major axis. If the planets satisfies Racc>Hpeb we have 2D accretion
    while we have 3D if Racc<Hpeb.
    
    mp: mass given in earth masses
    a : semi-major axis in AU
    mdot_gas: the gas accretion rate
    mdot_peb: the pebble accretion rate
    L_s: stellar luminosity
    M_s: stellar mass
    alpha_v: viscous alpha
    kap: disk opacity coefficient
    st: stokes number of the pebbles in the disk
    opt_vis: boolean statement for the inclusion of a modulation factor"""
    
    hg = hgas(a,mdot_gas,L_s,M_s,alpha_v,kap,True)
    hpeb = np.sqrt(alpha_v/(alpha_v+st))*hg
    
    Mp = mp/msuntome
    
    Rh      = a*(Mp/(3*M_s))**(1/3)
    Racc    = Rh*st**(1/3)
    acc_2d  = Racc>hpeb*a
    acc_3d  = Racc<hpeb*a
    
    #We compute the pebble accretion efficiency
    eta = eta_gas(a,mdot_gas,L_s,M_s,alpha_v,kap,opt_vis)
    
    fset = np.exp(-0.07*(eta/2.5e-3)**2*(mp/0.01)**(-2/3)*(st/0.1)**(2/3))
    
    t         = np.zeros(len(a))
    t[acc_2d] = 9e-4*(mp/0.05)**(1/3)*(st/0.1)**(1/3)*(eta[acc_2d]/2.5e-3)*(mdot_peb/100)**(-1)
    t[acc_3d] = 9e-4*(hpeb[acc_3d]/4.2e-3)*(eta[acc_3d]/2.5e-3)*(mdot_peb/100)**(-1)
    
    return t#[yr]