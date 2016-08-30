#!/usr/bin/env python
# -*- coding: utf-8 -*-


##### IMPORTING MODULES ############################################################################################################
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.widgets import Slider, Button
from random import choice
from string import ascii_letters

import astropy.constants as aconst
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants as const
import shutil
import subprocess
import sys
import os


##### CHECK FOR ARGUMENTS ##########################################################################################################
try:
    fargo_dir = sys.argv[1]
except:
    sys.exit('You have to give the FARGO data diretory as first argument!')
try:
    dust_dir = sys.argv[2]
except:
    sys.exit('You have to give the dust data diretory as second argument!')


##### LOAD DIMS.DAT ################################################################################################################
# This is an output file which contains for example the number of grid points
try:
    dims = np.loadtxt(fargo_dir+'/dims.dat')
except:
    sys.exit('Could not load dims.dat')


##### LOAD USED_RAD.DAT ############################################################################################################
# This file contains the positions of the radial grid cell interfaces
try:
    used_rad = np.loadtxt(fargo_dir+'/used_rad.dat') * u.AU
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
# Cartesian coordinates of the planet
planet_x = planet0[:,1] * u.AU
planet_y = planet0[:,2] * u.AU
# Trannsform to polar coordinates
planet_R     = np.sqrt( planet_x**2 + planet_y**2 )
planet_theta = np.arctan2( planet_y, planet_x )
        
        
##### LOAD DUST PARTICLES ##########################################################################################################
# Boolean array of time steps where there is dust
dust_exists = np.zeros( N_t+1, dtype=bool)
for i in range(N_t+1):
    if(os.path.isfile(dust_dir+'/dust'+repr(i)+'.dat')):
        dust_exists[i] = True
i_image = 0
# Load particle positions
if(np.any(dust_exists)):
    # Number of dust particles
    N_d = (np.loadtxt(dust_dir+'/dust'+repr(np.where(dust_exists)[0][0])+'.dat', comments='#')).shape[0]
    # Dust array
    dust = np.zeros( (N_d, 3) )
    dust[:, 2] = 1.e-100
if(dust_exists[i_image]):
    dummy = np.loadtxt(dust_dir+'/dust'+repr(i_image)+'.dat', comments='#')
    for j in range(N_d):
        dust[j, 0] = dummy[j,  1] # Radial position
        dust[j, 1] = dummy[j,  2] # Azimuthal position
        dust[j, 2] = dummy[j, 11] # Stokes number
    dust[ dust[:, 0] < R_cen[1].to(u.AU).value, 0 ] = 1.e6

##### PLOTTING #####################################################################################################################
# Create colormap
St_range = 0.1
color_dict = {'red':   ((0.00,          0.00,         0.00         ),
                        (0.50-St_range, 0.00,         0.00         ),
                        (0.50,          1.00,         1.00         ),
                        (0.50+St_range, 1.00,         1.00         ),
                        (1.00,          0.33,         0.33         )),
              'green': ((0.00,          0.00,         0.00         ),
                        (0.50-St_range, 0.00,         0.00         ),
                        (0.50,          0.75,         0.75         ),
                        (0.50+St_range, 0.00,         0.00         ),
                        (1.00,          0.00,         0.00         )),
              'blue':  ((0.00,          1.00,         1.00         ),
                        (0.50-St_range, 0.33,         0.33         ),
                        (0.50,          0.00,         0.00         ),
                        (0.50+St_range, 0.00,         0.00         ),
                        (1.00,          0.00,         0.00         )),
             }
cmap_stokes = LinearSegmentedColormap('Stokes number', color_dict)

# Particle display properties
particle_size  = 3.
particle_alpha = 1.

# Create plot
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))           # Use polar coordinate system
plt.subplots_adjust(bottom=0.25)                                      # Create space at the bottom for the slider

if(np.any(dust_exists)):
    scatter = ax.scatter(dust[:,1], dust[:,0], c=np.log10(dust[:,2]), s=particle_size, cmap=cmap_stokes, linewidths=0., alpha=particle_alpha, edgecolor='')
    St_max = np.ceil( np.log10( np.max( dust[:, 2] ) ) )
    St_min = St_max - 6.
    St_min = -2
    St_max = 2
    scatter.set_clim( [St_min, St_max] )
    cbar_ticks = np.arange( St_min, St_max+0.1, 1. )
    cbar = fig.colorbar(scatter)                                             # Show colorbar
    cbar.ax.set_ylabel('log Stokes Number')
    cbar.set_ticks(cbar_ticks)
    cbar.solids.set_edgecolor("face")

planet, = ax.plot(planet_theta[i_image], planet_R[i_image].to(u.AU), color='cyan', marker='o', markersize=8, markeredgecolor='black') # Star at planet position

ax.set_rmin(-R_cen[0].to(u.AU).value)                                 # Creating inner hole. Otherwise R=R_min would be in the center
ax.set_rmax(R_cen[-1].to(u.AU).value)                                 # Needed somehow to define outer limit.

ax.tick_params(axis='x', labelbottom='off')                           # Turns off the theta axis labelling
ax.tick_params(axis='y', colors='black', size=16)                     # Change radial labelling to white for better visibility


plt.grid(b=True)                                                      # Disable grid. Looks ugla with grid

# Create slider for the time
ax_time     = plt.axes([0.25, 0.1, 0.5, 0.03], axisbg='lightgoldenrodyellow')
slider_time = Slider(ax_time, 'timestep', 0.0, N_t, valinit=0,valfmt='%i')
ax._widgets = [slider_time] # Avoids garbage collection

# Create movie button
ax_button    = plt.axes([0.25, 0.04, 0.17, 0.04])
button_movie = Button(ax_button, 'Create movie', hovercolor='0.975')
ax._widgets += [button_movie] # Avoids garbage collection

plt.show()                  # Show plot


##### FUNCTIONS ####################################################################################################################
def update(val):
    i = int(np.floor(slider_time.val))  # Slider position
    
    if(dust_exists[i]):
        dummy = np.loadtxt(dust_dir+'/dust'+repr(i)+'.dat', comments='#')
        for j in range(N_d):
            dust[j, 0] = dummy[j,  1] # Radial position
            dust[j, 1] = dummy[j,  2] # Azimuthal position
            dust[j, 2] = dummy[j, 11] # Stokes number
        dust[ dust[:, 0] < R_cen[1].to(u.AU).value, 0 ] = 1.e6
    else:
        dust[:, 0] = 1.e6
            
    scatter.set_offsets( dust[:, 1::-1] )
    scatter.set_array( np.log10(dust[:,2]) )
    scatter.set_clim( [St_min, St_max] )
    
    # Update planet location    
    planet.set_data( (planet_theta[i], planet_R[i].to(u.AU)) )
    
    # Rescale axes. Don't know why this has to be done again.
    ax.set_rmin(-R_cen[0].to(u.AU).value)
    ax.set_rmax(R_cen[-1].to(u.AU).value)
    
    plt.draw() # Draw updates
    
def create_movie(event):
    dir_name = 'temporary_movie_files_'+''.join(choice(ascii_letters) for x in range(5))
    img_format = '.png'
    # Create the folder
    os.mkdir(dir_name)
    # Get current time step
    i_start = int(np.floor(slider_time.val))
    for j, i in enumerate(np.arange(i_start, N_t+1)):
        movie_fig, movie_ax = plt.subplots(subplot_kw=dict(projection='polar'))
        
        if(dust_exists[i]):
            dummy = np.loadtxt(dust_dir+'/dust'+repr(i)+'.dat', comments='#')
            for k in range(N_d):
                dust[k, 0] = dummy[k,  1] # Radial position
                dust[k, 1] = dummy[k,  2] # Azimuthal position
                dust[k, 2] = dummy[k, 11] # Stokes number
            dust[ dust[:, 0] < R_cen[1].to(u.AU).value, 0 ] = 1.e6
        else:
            dust[:, 0] = 1.e6
            
        movie_scatter = movie_ax.scatter(dust[:,1], dust[:,0], c=np.log10(dust[:,2]), s=particle_size, cmap=cmap_stokes, linewidths=0., alpha=particle_alpha, edgecolor='')
        movie_scatter.set_clim( [St_min, St_max] )
        movie_cbar_ticks = np.arange( St_min, St_max+0.1, 1. )
        movie_cbar = movie_fig.colorbar(scatter)                                             # Show colorbar
        movie_cbar.ax.set_ylabel('log Stokes Number')
        movie_cbar.set_ticks(movie_cbar_ticks)
        movie_cbar.solids.set_edgecolor("face")
            
        movie_planet, = movie_ax.plot(planet_theta[i], planet_R[i].to(u.AU), color='cyan', marker='o', markersize=8, markeredgecolor='black')
        movie_ax.set_rmin(-R_cen[0].to(u.AU).value)
        movie_ax.set_rmax(R_cen[-1].to(u.AU).value)
        movie_ax.tick_params(axis='x', labelbottom='off')
        movie_ax.tick_params(axis='y', colors='black')
        plt.grid(b=True)
        img_name = 'movie_plot_%05i%s'%(i, img_format)
        plt.savefig(dir_name+'/'+img_name)
        print 'Saving image: '+dir_name+'/'+img_name
        plt.close(movie_fig)
    
    # Create movie
    movie_name = 'movie.mp4'
    dummy_name = movie_name
    dummy_i = 0
    while os.path.isfile(dummy_name):
        dummy_i += 1
        dummy_name = movie_name.replace('.', '_%03i.'%dummy_i)
    movie_name = dummy_name
    subprocess.call(['ffmpeg','-start_number',repr(i_start),'-i',dir_name+os.sep+'movie_plot_%05d'+img_format,'-c:v','libx264','-crf','0','-pix_fmt','yuv420p',movie_name]);
    print 'Saving movie '+movie_name
    # Delete folder
    shutil.rmtree(dir_name)
    
slider_time.on_changed(update)          # This is executed when changing the slider
button_movie.on_clicked(create_movie)   # This is done by click on the movie button
