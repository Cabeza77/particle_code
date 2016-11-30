#!/usr/bin/env python
# -*- coding: utf-8 -*-


##### IMPORTING MODULES ############################################################################################################
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Slider, Button
from random import choice
from string import ascii_letters

import astropy.constants as aconst
import astropy.units as u
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import shutil
import subprocess
import sys
import os


##### PHYSICAL PARAMETERS ##########################################################################################################
phys_mass = 1. * aconst.M_sun
phys_dist = 1. * u.au
phys_time = (np.sqrt( phys_dist**3. / (aconst.G * phys_mass) )).to(u.yr)


##### CHECK FOR ARGUMENTS ##########################################################################################################
try:
    fargo_dir = sys.argv[1]
except:
    sys.exit('You have to give the FARGO data diretory as first argument!')
try:
    dust_dir = sys.argv[2]
except:
    sys.exit('You have to give the dust data diretory as second argument!')
try:
    i_dust = int( sys.argv[3] )
except:
    sys.exit('You have to give the identifier of the dust particle as third argument!')


##### LOAD DIMS.DAT ################################################################################################################
# This is an output file which contains for example the number of grid points
try:
    dims = np.loadtxt(fargo_dir+'/dims.dat')
except:
    sys.exit('Could not load dims.dat')


##### LOAD USED_RAD.DAT ############################################################################################################
# This file contains the positions of the radial grid cell interfaces
try:
    used_rad = np.loadtxt(fargo_dir+'/used_rad.dat') * phys_dist
except:
    sys.exit('Could not load used_rad.dat!')


##### CREATE GRID ##################################################################################################################
N_t       = int( dims[5] )                            # Number of outputs
N_R       = int( dims[6] )                            # Number of radial gridpoint
N_theta   = int( dims[7] )                            # Number of angular grid points
R_cen     = 0.5 * (used_rad[:-1] + used_rad[1:] )     # Calculate radial cell centers
theta_cen = np.linspace( 0., 2.*const.pi, N_theta )   # Theta grid


##### LOAD PLANET0.DAT #############################################################################################################
# The planet data file
try:
    planet0 = np.loadtxt(fargo_dir+'/planet0.dat')
except:
    sys.exit('Could not load planet0.dat!')
dummy_i = -1
for i in range(planet0.shape[0]):
    if(planet0[i, 0] == dummy_i):
        planet0[i-1:planet0.shape[0]-1, :] = planet0[i:,:]
    dummy_i = planet0[i, 0]
# Cartesian coordinates of the planet and mass in stellar masses
time    = np.zeros(N_t+1) * u.yr
time[:] = np.nan
for it in range(N_t+1):
    try:
        time[it] = planet0[it, 7] * phys_time
    except:
        pass
        
        
##### CHECK DUST FILES #############################################################################################################
# Here we check which dust files exist and how many dust particles we have.
# ¡¡¡ All files need to have the same amount of dust particles !!!
# Boolean array of time steps where we have dust files
dust_exists = np.zeros( N_t+1, dtype=bool)
for i in range(N_t+1):
    if(os.path.isfile(dust_dir+'/dust'+repr(i)+'.dat')):
        dust_exists[i] = True
i_last = np.argmin(dust_exists) - 1

# Get number of dust particles and initialize dust array
N_dust = 0
if(np.any(dust_exists)):
    # Number of dust particles
    N_dust = (np.loadtxt(dust_dir+'/dust'+repr(np.where(dust_exists)[0][0])+'.dat', comments='#')).shape[0]
# Dust array
dust       = np.zeros( (N_t+1, 10) )
dust[:, :] = np.nan

if(i_dust < 1 or i_dust > N_dust):
    sys.exit('The particle identifier is out of range.')

# Read the dust properties
print('Loading files...')
for it in range(N_t+1):
    if(dust_exists[it]):
        dummy = np.genfromtxt(dust_dir+'/dust'+repr(it)+'.dat', comments='#', skip_header=i_dust, skip_footer=N_dust-i_dust)
        if(dummy[1]>0.):
            try:
                dust[it, 0] = dummy[ 1]   # R
            except:
                pass
            try:
                dust[it, 1] = dummy[ 2]   # phi
            except:
                pass
            try:
                dust[it, 2] = dummy[ 3]   # vR
            except:
                pass
            try:
                dust[it, 3] = dummy[ 4]   # vTheta
            except:
                pass
            try:
                dust[it, 4] = dummy[ 9]   # a
            except:
                pass
            try:
                dust[it, 5] = dummy[10]   # m
            except:
                pass
            try:
                dust[it, 6] = dummy[11]   # St
            except:
                pass
            try:
                dust[it, 7] = dummy[15]   # Sigma gas
            except:
                pass
            try:
                dust[it, 8] = dummy[13]   # T
            except:
                pass
            try:
                dust[it, 9] = dummy[14]   # xi
            except:
                pass
    printline = "\r      {:4.1f} %".format(100.0*it/N_t)
    print(printline),
# Modify phi to account for periodicity
abs_dphi      = np.abs( np.diff( dust[:, 1] ) )
abs_dphi_mean = abs_dphi[ np.isnan( abs_dphi ) == False ].mean()
abs_dphi_std  = abs_dphi[ np.isnan( abs_dphi ) == False ].std()
phi_mask      = np.hstack( [ abs_dphi > abs_dphi_mean + 3.*abs_dphi_std, [False] ] )
masked_phi    = np.ma.MaskedArray( dust[:, 1], phi_mask )
print('\n\r   ...done.')


##### PLOTTING #####################################################################################################################
fig = plt.figure( figsize=(16., 8.) )

gs0 = gs.GridSpec(1, 3, width_ratios=[1,1,2])
gs00 = gs.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs0[0])
gs01 = gs.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs0[1])
gs02 = gs.GridSpecFromSubplotSpec(1, 1, subplot_spec=gs0[2])

ax1  = fig.add_subplot(gs00[0, :])
plt.setp(ax1.get_xticklabels(), visible=False)

ax2  = fig.add_subplot(gs00[1, :], sharex=ax1)
plt.setp(ax2.get_xticklabels(), visible=False)

ax3  = fig.add_subplot(gs00[2, :], sharex=ax1)
plt.setp(ax3.get_xticklabels(), visible=False)

ax4  = fig.add_subplot(gs00[3, :], sharex=ax1)
plt.setp(ax4.get_xticklabels(), visible=False)

ax5  = fig.add_subplot(gs00[4, :], sharex=ax1)

ax6  = fig.add_subplot(gs01[0, :], sharex=ax1)
plt.setp(ax6.get_xticklabels(), visible=False)

ax7  = fig.add_subplot(gs01[1, :], sharex=ax1)
plt.setp(ax7.get_xticklabels(), visible=False)

ax8  = fig.add_subplot(gs01[2, :], sharex=ax1)
plt.setp(ax8.get_xticklabels(), visible=False)

ax9  = fig.add_subplot(gs01[3, :], sharex=ax1)
plt.setp(ax9.get_xticklabels(), visible=False)

ax10 = fig.add_subplot(gs01[4, :], sharex=ax1)

ax11 = fig.add_subplot(gs02[:, :], projection='polar')

linewidth = 1

ax1.plot(time.to(u.kyr).value, dust[:, 0], lw=linewidth, color='blue')
ax1.plot(time[0].to(u.kyr).value, dust[0, 0], 'o', markersize=8, color='#00FF00')
ax1.plot(time[i_last].to(u.kyr).value, dust[i_last, 0], 'o', markersize=8, color='#FF0000')
ax1.set_ylabel( r'$R\ \mathrm{[AU]}$' )
ax1.grid(b=True)

ax2.plot(time.to(u.kyr).value, (dust[:, 2]*u.cm/u.s).to(u.au/u.kyr).value, lw=linewidth, color='blue')
ax2.set_ylabel( r'$v_\mathrm{R}\ \mathrm{[AU/kyr]}$' )
range_max = np.amax( np.abs( (dust[:, 2][np.isnan( dust[:, 2] ) == False]*u.cm/u.s).to(u.au/u.kyr).value ) ) * 1.5
ax2.set_ylim(-10., 10.)
ax2.grid(b=True)

ax3.semilogy(time.to(u.kyr).value, dust[:, 6], lw=linewidth, color='purple')
ax3.set_ylabel( r'$\mathrm{St}$' )
ax3.grid(b=True)

ax4.semilogy(time.to(u.kyr).value, dust[:, 7], lw=linewidth, color='purple')
ax4.set_ylabel( r'$\Sigma_\mathrm{gas}\ \mathrm{[g/cm^2]}$' )
ax4.grid(b=True)

ax5.semilogy(time.to(u.kyr).value, dust[:, 4], lw=linewidth, color='green')
ax5.set_ylabel( r'$a\ \mathrm{[cm]}$' )
ax5.set_xlabel( r'$t\ \mathrm{[kyr]}$' )
ax5.grid(b=True)

#ax6.plot(time.to(u.kyr).value, dust[:, 1], lw=linewidth, color='blue')
ax6.plot(time.to(u.kyr).value, masked_phi, lw=linewidth, color='blue')
ax6.plot(time[0].to(u.kyr).value, dust[0, 1], 'o', markersize=8, color='#00FF00')
ax6.plot(time[i_last].to(u.kyr).value, dust[i_last, 1], 'o', markersize=8, color='#FF0000')
ax6.set_ylabel( r'$\varphi\ \mathrm{[rad]}$' )
ax6.set_ylim(0., 2.*const.pi)
ax6.grid(b=True)

ax7.plot(time.to(u.kyr).value, dust[:, 3]/np.sqrt( aconst.G * phys_mass / (dust[:, 0] * phys_dist) ).cgs.value, lw=linewidth, color='blue')
ax7.set_ylabel( r'$v_\varphi\ [v_\mathrm{K}]$' )
range_max = np.amax( np.abs( 1. - dust[:, 3][np.isnan( dust[:, 3] ) == False]/np.sqrt( aconst.G * phys_mass / (dust[:, 0][np.isnan( dust[:, 3] ) == False] * phys_dist) ).cgs.value ) ) * 1.5
ax7.set_ylim(1.-range_max, 1.+range_max)
ax7.grid(b=True)

ax8.plot(time.to(u.kyr).value, dust[:, 8], lw=linewidth, color='red')
ax8.set_ylabel( r'$T\ \mathrm{[K]}$' )
ax8.grid(b=True)

ax9.plot(time.to(u.kyr).value, dust[:, 9], lw=linewidth, color='red')
ax9.set_ylabel( r'$f_\mathrm{cryst}$' )
ax9.set_ylim(-0.1, 1.1)
ax9.grid(b=True)

ax10.semilogy(time.to(u.kyr).value, dust[:, 5], lw=linewidth, color='green')
ax10.set_ylabel( r'$m\ \mathrm{[g]}$' )
ax10.set_xlabel( r'$t\ \mathrm{[kyr]}$' )
ax10.grid(b=True)

ax11.plot(dust[:, 1], dust[:, 0])
ax11.plot(dust[0, 1], dust[0, 0], 'o', markersize=8, color='#00FF00')
ax11.plot(dust[i_last, 1], dust[i_last, 0], 'o', markersize=8, color='#FF0000')
ax11.set_rlim( -used_rad[0].to(u.au).value, used_rad[-1].to(u.au).value )
ax11.grid(b=True)

plt.tight_layout()

plt.draw()
plt.show()








 
