#! /usr/bin/env python

import numpy as np
import shtns # necessary? imported within pyspharm
import pyspharm
from scipy import io
from scipy import interpolate
from matplotlib import pyplot as plt

# Code to solve the elastic sea level equation following 
# Kendall et al., 2005 and Austermann et al., 2015

# J. Austermann 2015
# Translated to Python and modified by A. Wickert, 2015

######################
# Parameters & Input #
######################

# Eventually, make these into variables that are defined at runtime
# e.g., func(maxdeg, rho_ice=916.7, rho_water=1000., rho_sed=2300., g=9.81)

# Specify maximum degree to which spherical transformations should be done
maxdeg = 128 # N

# parameters
rho_ice = 916.7
rho_water = 1000.
rho_sed = 2300.
g = 9.80616

# Using other library
#sh = shtns.sht(maxdeg, maxdeg)

# And with bindings
# grid, time step info
nlons = maxdeg * 2  # number of longitudes
ntrunc = maxdeg - 1
nlats = maxdeg # for gaussian grid.
a = rsphere = 6.37122e6 # earth radius

# setup up spherical harmonic instance, set lats/lons of grid
sh = pyspharm.Spharmt(nlons,nlats,ntrunc,rsphere,gridtype='gaussian')
# Check if needed
elons  = sh.lons * 180 / np.pi
lats   = sh.lats * 180 / np.pi
colats = 90 - lats
lon_out, lat_out = np.meshgrid(elons, colats)


# Variables needed later
M_e = 5.9742E24
# for calc_rot
omega = 7.292E-5
k_hydro = 0.934
CminA = 2.6E35
sqrt_32_15 = (32/15.)**0.5


# Precompute legendre polynomials
# --- this seems to be done already in the spherical harmonic library!
# see Schaeffer, 2015

# --------------------------------
# ICE
# --------------------------------

ice = io.loadmat('WAIS')
# ice_Ant
# ice_EAIS
# lat_WAIS
# lon_WAIS

# interpolate ice masks on common grid
# Consider RectBivariateSpline
f = interpolate.interp2d(ice['lon_WAIS'][0,:], \
                              ice['lat_WAIS'][:,0], ice['ice_Ant'])
ice_0 = f(elons, lats[::-1])[::-1]
f = interpolate.interp2d(ice['lon_WAIS'][0,:], \
                         ice['lat_WAIS'][:,0], ice['ice_EAIS'])
ice_j = f(elons, lats[::-1])[::-1]

del_ice = ice_j - ice_0

#ice_0[20:90,20:90] += 4E3
#ice_0[:,:] = 10E3


# --------------------------------
# DYNAMIC TOPOGRAPHY
# --------------------------------

del_DT = np.zeros(del_ice.shape)


# --------------------------------
# SEDIMENT
# --------------------------------

del_sed = np.zeros(del_ice.shape)


# --------------------------------
# TOPOGRAPHY
# --------------------------------

# load preloaded etopo (including ice) as topo_orig, lon_topo, lat_topo
topo = io.loadmat('gebco_08_15am')
topo_nlats, topo_nlons = topo['map_data'].shape
topo_colats = np.linspace(180./topo_nlats/2., 180 - 180./topo_nlats/2., topo_nlats)
topo_elons = np.linspace(360./topo_nlons/2., 360 - 360./topo_nlons/2., topo_nlons)

# interpolate topography grid onto Gauss Legendre Grid
f = interpolate.interp2d(topo_elons, topo_colats, topo['map_data'])
topo0 = f(elons, colats)


## Set up love number input

def love_lm(num, group='l'):
  """
  no dc shift also for gravity love number
  """
  h = np.hstack(( 0, num.squeeze() ))
  h_lm = [];
  if group == 'l':
    # This is the standard for Jerry's code
    for l in range(maxdeg+1):
      h_lm += ([h[l]] * l)
  elif group == 'm':
    # This works with the spherical harmonic library I am using here
    for m in range(maxdeg):
      h_lm += (list(h[:m+1]))
  h_lm = np.array(h_lm)
  return h_lm

# prepare love numbers in suitable format and calculate T_lm and E_lm 
# to calculate the fluid case, switch h_el to h_fl, k_el to k_fl and same
# for tidal love numbers
love = io.loadmat('SavedLN/LN_l90_VM2')
h_lm = love_lm(love['h_el'])
k_lm = love_lm(love['k_el'])
h_lm_tide = love_lm(love['h_el_tide'])
k_lm_tide = love_lm(love['k_el_tide'])

E_lm = 1 + k_lm - h_lm
# Need to have proper size for the in/out grid size; set higher-order values
# to 0
#E_lm = np.hstack(( E_lm, np.zeros(len(oc0_lm) - len(E_lm)) )) * 0

def get_tlm(maxdeg, group='l'):
  T_lm = []
  T = np.zeros(maxdeg+1)
  const = 4*np.pi*a**3/M_e
  if group == 'l':
    # This is the standard for Jerry's code
    for l in range(maxdeg+1):
      T[l] = const / (2 * l + 1.)
      T_lm += ([T[l]] * l) # IMPORTANT
  elif group == 'm':
    # This works with the spherical harmonic library I am using here, I think
    for m in range(maxdeg):
      T[m] = const / (2 * m + 1.)
      T_lm += (list(T[:m+1])) # IMPORTANT
  T_lm = np.array(T_lm)
  return T_lm

T_lm = get_tlm(maxdeg)
# Need to have proper size for the in/out grid size; set higher-order values
# to 0
#T_lm = np.hstack(( T_lm, np.zeros(len(oc0_lm) - len(T_lm)) )) * 0

E_lm_T = 1 + k_lm_tide - h_lm_tide
# Need to have proper size for the in/out grid size; set higher-order values
# to 0
#E_lm_T = np.hstack(( E_lm_T, np.zeros(len(oc0_lm) - len(E_lm_T)) )) * 0

# can switch this in if you want to exclude rotational effects
# E_lm_T = zeros(size(E_lm_T))

## Solve sea level equation (after Kendall 2005, Dalca 2013 & Austermann et al. 2015)

# Change this to make a better convergence criterion, eventually
k_max = 10   # maximum number of iterations
epsilon = 10**-4 # convergence criterion

# 0 = before
# j = after

# set up present-day topo and ocean function 
topo_0 = topo0.copy(); # already includes ice and dynamic topography

def sign_01(ingrid):
  out = -0.5*np.sign(ingrid)+0.5
  out = 0.5*np.sign(out-0.6)+0.5
  return out

oc_0 = sign_01(topo_0)

# set up topography and ocean function after the ice change
topo_j = topo_0 + del_ice; # del_ice is negative -> subtract ice that is melted
oc_j = sign_01(topo_j)

# INCLUDES IMAGINARY PARTS -- JUST TAKE REAL PART?

# calculate change in sediments and decompose into spherical harmonics
Sed_lm = sh.grdtospec(del_sed, norm='unity')

# expand ocean function into spherical harmonics
oc0_lm = sh.grdtospec(oc_0, norm='unity')

def calc_rot(L_lm, _k, _k_tide, group='l'):

  # extract degree 2:
  if group == 'l':
    L20 = L_lm[3]
    L21 = L_lm[4] 
    L22 = L_lm[5]
    print L22
  elif group == 'm':
    print "WARN: hard-coded max l,m; will break!"
    L20 = L_lm[2]
    L21 = L_lm[256+2]
    L22 = L_lm[512]
    print L22
  k_L = _k[1] # This has 256 values
  k_T = _k_tide[1] # This has 256 values

  I13 = sqrt_32_15 * np.pi * a**4 * np.real(L21)
  I23 = - sqrt_32_15 * np.pi * a**4 * np.imag(L21)

  II = I13 + 1j * I23

  # equation from Mitrovica and Wahr 2005

  m0 = II/CminA * (1 + k_L) / ( 1 - (k_T/k_hydro) )

  m1 = np.real(m0)
  m2 = np.imag(m0)
  m3 = 0

  m = np.hstack((m1, m2, m3))

  # calculate the perturbation to the rotational potential from Milne 1998

  La00 = a**2 * omega**2/3 * (np.sum(m**2) + 2*m3)
  La20 = a**2 * omega**2/(6*5**.5) * (m1**2 + m2**2 - 2*m3**2 - 4*m3)
  La21 = a**2 * omega**2/30**.5 * (m1*(1+m3) - 1j*m2*(1+m3))
  La22 = a**2 * omega**2/5**.5 * 24**.5 * ( (m2**2-m1**2) + 1j*2*m1*m2 )
  La2m1 = -1 * np.conj(La21)
  La2m2 = 1 * np.conj(La22)

  La_lm = np.zeros(L_lm.shape, dtype=np.complex128)

  La_lm[:6] += np.hstack((La00, 0, 0, La20, La21, La22))

  return La_lm
  
k = 0
chi = epsilon * 2

# Iterations k
# Iteration criterion epsilon
while (k < k_max) and (chi >= epsilon):

  #print chi

  # expand ocean function into spherical harmonics
  ocj_lm = sh.grdtospec(oc_j, norm='unity')

  # CHECK ICE MODEL 
  # check ice model for floating ice
  check1 = sign_01(-topo_j + ice_j)
  check2 = sign_01(+topo_j - ice_j) * \
    (sign_01(-ice_j*rho_ice - (topo_j - ice_j)*rho_water))
  
  ice_j_corr = check1 * ice_j + check2 * ice_j
  del_ice_corrected = ice_j_corr - ice_0
  
  deli_lm = sh.grdtospec(del_ice_corrected, norm='unity')
  
  # calculate topography correction
  TO = topo_0 * (oc_j-oc_0)
  # expand TO function into spherical harmonics
  TO_lm = sh.grdtospec(TO, norm='unity')
  
  # set up initial guess for sea level change
  if k == 0:
    # initial guess of sea level change is just to distribute the
    # ice over the oceans uniformally
    # (remember, harmonic (0,0) is the mean)
    delS_lm = ocj_lm/ocj_lm[0]*(-rho_ice/rho_water*deli_lm[0] + \
        TO_lm[0])
    # convert into spherical harmonics
    #delS_init = sh.spectogrd(delS_lm, norm='unity')
      
  # calculate loading term
  L_lm = rho_ice*deli_lm + rho_water*delS_lm + rho_sed*Sed_lm

  # calculate contribution from rotation
  La_lm = calc_rot(L_lm, love['k_el'], love['k_el_tide'])

  # calculate sea level perturbation
  # add ice and sea level and multiply with love numbers
  # DT doesn't load!
  # (note: this is curly L from Kendall et al. (2005), and has nothing to 
  #  do with the mathematical curl)
  # Also, truncate those arrays converted from higher-resolution to the
  # appropriate resolution for the harmonic operations
  # (don't know why they are going in/out at a resolution related to the 
  # spatial resolution)
  delSLcurl_lm_fl = E_lm * T_lm * (rho_ice*deli_lm + rho_water*delS_lm \
                                   + rho_sed*Sed_lm) + 1./g*E_lm_T * La_lm

  # convert to spherical harmonics and subtract terms that are part
  # of the topography to get the 'pure' sea level change
  delSLcurl_fl = sh.spectogrd( delSLcurl_lm_fl, norm='unity')
  delSLcurl = delSLcurl_fl - del_ice_corrected - del_DT - del_sed


  # compute and decompose RO
  RO = delSLcurl * oc_j
  RO_lm = sh.grdtospec(RO, norm='unity')

  # calculate eustatic sea level perturbation (delta Phi / g)
  delPhi_g = 1/ocj_lm[0] * (- rho_ice/rho_water*deli_lm[0] \
      - RO_lm[0] + TO_lm[0])

  # calculate overall perturbation of sea level over oceans
  # (spatially varying field and constant offset)
  delSL = delSLcurl + delPhi_g


  # update topography and ocean function
  topo_j = - delSL + topo_0
  oc_j = sign_01(topo_j)


  # calculate change in ocean height and decompose
  delS_new = delSL * oc_j -  topo_0 * (oc_j-oc_0)
  delS_lm_new = sh.grdtospec(delS_new, norm='unity')


  # calculate convergence criterion chi
  chi = np.abs( (np.sum(np.abs(delS_lm_new)) - np.sum(np.abs(delS_lm))) \
        / np.sum(np.abs(delS_lm)) )
        
  # update sea surface height
  delS_lm = delS_lm_new
        
  k += 1

if chi < epsilon:
  print 'Converged after iteration', k, 'Chi was', chi
else:
  print 'Did not yet converge.'
  print 'Finished iteration', k, '; Chi was', chi

# calculate the scaling to normalize the fingerprint (it's normalized to be
# one on average, when averaged over the final ocean basin). 
# calculate change in sea level over final ocean basin
del_scaling = (delSL + del_ice_corrected) * oc_j
# get the average of that when spreading the water over the whole globe
sca = sh.grdtospec(del_scaling, norm='unity')
# get the average of that when spreading the water only over the oceans.
scaling_fact = sca[0] / ocj_lm[0]


## Plot results

# We only want the sea level change cause by melted ice, so subtract
# del_ice
SL_change_rot = delSL + del_ice_corrected
# normalize it by the scaling factor
plotSL = SL_change_rot/scaling_fact

"""
# construct identical colormap to Mitrovica 2009 paper
MyColorMap = [238   44  37
    211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244;211 238 244
    173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235;173 224 235
    163 201 235
    111 147 201
    96  103 175
    74  102 176
    68  87  165
    58  84  163
    53  69  154
    44  47  137
    38  35  103
    19  15  54
    0   0   0
    0   0   0]
"""

# plot
plt.figure
plt.imshow(np.real(plotSL))
plt.colorbar()
plt.show()

"""
m_proj('hammer-aitoff','clongitude',0)
m_pcolor([lon_out(:,end/2+1:end)-360 lon_out(:,1:end/2)],lat_out,...
    [plotSL(:,end/2+1:end) plotSL(:,1:end/2)])
m_coast('color',[0 0 0])
m_grid('box','fancy','xticklabels',[],'yticklabels',[])
shading flat
colorbar
colormap(MyColorMap/255)
caxis( [-0.05, 1.6] )
title('Sea level fingerprint of West Antarctic Ice Sheet collapse')
"""
