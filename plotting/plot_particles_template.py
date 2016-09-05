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


##### COSTUM COLORBARS #############################################################################################################
# Stokes number: blue - dark blue - yellow - red - dark red
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
# Particle radius: black - blue - purple - orange - green
color_dict = {'red':   ((0.00,          0.00,         0.00         ),
                        (0.50,          0.67,         0.67         ),
                        (0.75,          1.00,         1.00         ),
                        (1.00,          0.00,         0.00         )),
              'green': ((0.00,          0.00,         0.00         ),
                        (0.50,          0.00,         0.00         ),
                        (1.00,          1.00,         1.00         )),
              'blue':  ((0.00,          0.00,         0.00         ),
                        (0.25,          1.00,         1.00         ), 
                        (0.50,          0.67,         0.67         ),
                        (1.00,          0.00,         0.00         )),
             }
cmap_radius = LinearSegmentedColormap('Radius', color_dict)


##### OPTIONS ######################################################################################################################
column  = 9                           # Column in dust file you want plot
                                      #  0: particle identifier
                                      #  1: R
                                      #  2: theta
                                      #  3: R-velocity
                                      #  4: theta-velocity
                                      #  5: X
                                      #  6: Y
                                      #  7: X-velocity
                                      #  8: Y-velocity
                                      #  9: Radius [cm]
                                      # 10: Mass [g]
                                      # 11: Stokes number
                                      # 12: Stopping time [s]
                                      # 13: Temperature [K]
                                      # 14: Crystallinity fraction
val_min = 10.**(-0.5)                 # Minimum and maximum value you want
val_max = 10.**(0.5)                  #     to plot.
log     = True                        # Shall the scale be loagrithmic?

col_min = 1.e-4                       # Minimum and maximum value
col_max = 1.e2                        #     of colorbar
cmap    = cmap_radius                 # What colorbar do you want?
clabel  = r'Log Radius [cm]'          # Label of the colorbar. LaTeX allowed.

psize   = 3.                          # Particle point size in plot
palpha  = 1.                          # Particle alpha value in plot

spec_R  = True                        # Do want to specify your own R boundaries?
R_min   = 3.                          # Your chosen R boundaries, if spec_R = True.
R_max   = 100.                        #     (in AU)

col_pla = '#00FFFF'                   # Color of the planet
plsize  = 8                           # Size of the planet
show_La = True                        # Show the Lagrangian points?
col_La  = '#FF0000'                   # Color of the Lagrangian points
Lasize  = 6                           # Size of the Lagrangian points


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
    used_rad = np.loadtxt(fargo_dir+'/used_rad.dat')
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
planet_x = planet0[:,1]
planet_y = planet0[:,2]
# Trannsform to polar coordinates
planet_R     = np.sqrt( planet_x**2 + planet_y**2 )
planet_theta = np.arctan2( planet_y, planet_x )
        
        
##### CHECK DUST FILES #############################################################################################################
# Here we check which dust files exist and how many dust particles we have.
# ¡¡¡ All files need to have the same amount of dust particles !!!

# Boolean array of time steps where we have dust files
dust_exists = np.zeros( N_t+1, dtype=bool)
for i in range(N_t+1):
    if(os.path.isfile(dust_dir+'/dust'+repr(i)+'.dat')):
        dust_exists[i] = True

# Get number of dust particles and initialize dust array
N_dust = 0
if(np.any(dust_exists)):
    # Number of dust particles
    N_dust = (np.loadtxt(dust_dir+'/dust'+repr(np.where(dust_exists)[0][0])+'.dat', comments='#')).shape[0]
# Dust array
dust   = np.zeros( (N_dust, 3) )
    
    
##### FUNCTION TO LOAD DUST FILE ###################################################################################################
def load_dust_file(i_frame, i_col):
    if(not dust_exists[i_frame]):
        dust[:, 0] = 1.e100
        return
    dummy = np.loadtxt(dust_dir+'/dust'+repr(i_frame)+'.dat', comments='#')
    for j in range(N_dust):
        dust[j, 0] = dummy[j, 1]     # Radial position
        dust[j, 1] = dummy[j, 2]     # Azimuthal position
        dust[j, 2] = dummy[j, i_col] # Chosen column
    # If dust is out of boundaries, set R to huge value.
    dust[ dust[:, 0] < R_min, 0 ] = 1.e100
    # Chose column value in desired range
    dust[ dust[:, 2] < val_min, 0 ] = 1.e100
    dust[ dust[:, 2] > val_max, 0 ] = 1.e100
    # Convert to log scale if wanted
    if(log):
        dust[:, 2] = np.log10( dust[:, 2] )
    
    
##### SET OPTIONS ##################################################################################################################
# R boundaries
if(not spec_R):
    R_min = used_rad[0]
    R_max = used_rad[-1]
# Colorbar limits
if(log):
    clim = [ np.floor( np.log10( col_min ) ), np.ceil( np.log10( col_max ) ) ]
else:
    clim = [ col_min, col_max ]
    

##### PLOTTING #####################################################################################################################
# Load first frame
frame = 0
load_dust_file(frame, column)


# Create plot
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))           # Use polar coordinate system
plt.subplots_adjust(bottom=0.25)                                      # Create space at the bottom for the slider


scatter = ax.scatter(dust[:,1], dust[:,0], c=dust[:,2], s=psize, cmap=cmap, linewidths=0., alpha=palpha, edgecolor='')
scatter.set_clim( clim )
cbar = fig.colorbar(scatter)
cbar.ax.set_ylabel(clabel)
cbar.solids.set_edgecolor("face")

planet, = ax.plot(planet_theta[frame], planet_R[frame], color=col_pla, marker='o', markersize=plsize, markeredgecolor='black') # Star at planet position

if(show_La):
    L4, = ax.plot(planet_theta[frame]+const.pi/3., planet_R[frame], color=col_La, marker='o', markersize=Lasize, markeredgecolor='black')
    L5, = ax.plot(planet_theta[frame]-const.pi/3., planet_R[frame], color=col_La, marker='o', markersize=Lasize, markeredgecolor='black')

ax.set_rmin(-R_min)                                 # Creating inner hole. Otherwise R=R_min would be in the center
ax.set_rmax(R_max)                                  # Needed somehow to define outer limit.

ax.tick_params(axis='x', labelbottom='off')                           # Turns off the theta axis labelling
ax.tick_params(axis='y', colors='black', size=16)                     # Change radial labelling to white for better visibility


plt.grid(b=True)                                                      # Enable grid

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
    
    load_dust_file(i, column)
            
    scatter.set_offsets( dust[:, 1::-1] )
    scatter.set_array( dust[:,2] )
    
    # Update planet location    
    planet.set_data( (planet_theta[i], planet_R[i]) )
    
    if(show_La):
        L4.set_data( (planet_theta[i]+const.pi/3., planet_R[i]) )
        L5.set_data( (planet_theta[i]-const.pi/3., planet_R[i]) )
    
    # Rescale axes. Don't know why this has to be done again.
    ax.set_rmin(-R_min)
    ax.set_rmax(R_max)
    
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
        
        load_dust_file(i, column)
            
        movie_scatter = movie_ax.scatter(dust[:,1], dust[:,0], c=dust[:,2], s=psize, cmap=cmap, linewidths=0., alpha=palpha, edgecolor='')
        movie_scatter.set_clim( clim )
        movie_cbar = movie_fig.colorbar(movie_scatter)
        movie_cbar.ax.set_ylabel(clabel)
        movie_cbar.solids.set_edgecolor("face")
            
        movie_planet, = movie_ax.plot(planet_theta[i], planet_R[i], color=col_pla, marker='o', markersize=plsize, markeredgecolor='black')
        
        if(show_La):
            L4, = ax.plot(planet_theta[i]+const.pi/3., planet_R[i], color=col_La, marker='o', markersize=Lasize, markeredgecolor='black')
            L5, = ax.plot(planet_theta[i]-const.pi/3., planet_R[i], color=col_La, marker='o', markersize=Lasize, markeredgecolor='black')
        
        movie_ax.set_rmin(-R_min)
        movie_ax.set_rmax(R_max)
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
