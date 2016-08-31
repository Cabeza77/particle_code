#!/usr/bin/env python
# -*- coding: utf-8 -*-


##### IMPORTING MODULES ############################################################################################################
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.widgets import Slider, Button
from random import choice
from scipy.interpolate import interp2d
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
    data_dir = sys.argv[1]
except:
    sys.exit('You have to give the data diretory as first argument!')


##### LOAD DIMS.DAT ################################################################################################################
# This is an output file which contains for example the number of grid points
try:
    dims = np.loadtxt(data_dir+'/dims.dat')
except:
    sys.exit('Could not load dims.dat')


##### LOAD USED_RAD.DAT ############################################################################################################
# This file contains the positions of the radial grid cell interfaces
try:
    used_rad = np.loadtxt(data_dir+'/used_rad.dat') * u.AU
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
    planet0 = np.loadtxt(data_dir+'/planet0.dat')
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


##### LOAD DENSITIES ###############################################################################################################
# Load the gas densities
i_image = 0
sigma = np.zeros( (N_R, N_theta) ) * aconst.M_sun / u.AU**2 # Initialize density array
try:
    sigma[:, :] = np.fromfile(data_dir+'/gasdens'+repr(i_image)+'.dat').reshape(N_R, N_theta) * aconst.M_sun / u.AU**2
except:
    print 'Could not load gasdens'+repr(i_image)+'.dat'

# Function to interpolate data
def intpol(theta, r, data):
    f = interp2d(theta_cen, R_cen, data, kind='linear')
    return f(theta, r)
    

##### PLOTTING #####################################################################################################################
# Create amber-teal colormap
color_dict = {'red':   ((0.00, 0.10, 0.10),
                        (0.33, 0.10, 0.10),
                        (0.67, 1.00, 1.00),
                        (1.00, 1.00, 1.00)),
              'green': ((0.00, 0.10, 0.10),
                        (0.33, 0.10, 0.10),
                        (0.67, 0.50, 0.50),
                        (1.00, 1.00, 1.00)),
              'blue':  ((0.00, 0.10, 0.10),
                        (0.33, 0.50, 0.50),
                        (0.67, 0.10, 0.10),
                        (1.00, 1.00, 1.00))
             }
amber_teal = LinearSegmentedColormap('OrangeTeal1', color_dict)
amber_teal.set_under('#191919')
amber_teal.set_over('#FFFFFF')

# Constrain the colorbar
sigma_min  = np.min( sigma[sigma>0.].cgs.value )                                                  # Minimum (non-empty) value of sigma
sigma_max  = np.max( sigma.cgs.value )                                                            # Maximum value of sigma
levels     = np.linspace( np.floor(np.log10(sigma_min)), np.ceil(np.log10(sigma_max)), 100 )      # Levels of colorbar
cbar_ticks = np.arange( np.floor(np.log10(sigma_min)), np.ceil(np.log10(sigma_max))+0.1, 0.5 )    # Ticks of colorbar

# Create plot
fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))           # Use polar coordinate system
plt.subplots_adjust(bottom=0.25)                                      # Create space at the bottom for the slider

plot = ax.contourf(theta_cen, R_cen.to(u.AU), np.log10(sigma[ :, :].cgs.value), cmap=amber_teal, levels=levels, extend='both' )       # Filled contour plot
collections_list = plot.collections[:]                                                                                                # Collect the collections....lol
planet, = ax.plot(planet_theta[i_image], planet_R[i_image].to(u.AU), color='cyan', marker='o', markersize=8, markeredgecolor='black') # Cross at planet position

ax.set_rmin(-R_cen[0].to(u.AU).value)                                 # Creating inner hole. Otherwise R=R_min would be in the center
ax.set_rmax(R_cen[-1].to(u.AU).value)                                 # Needed somehow to define outer limit.

ax.tick_params(axis='x', labelbottom='off')                           # Turns off the theta axis labelling
ax.tick_params(axis='y', colors='white')                              # Change radial labelling to white for better visibility

cbar = fig.colorbar(plot)                                             # Show colorbar
cbar.ax.set_ylabel(r'log $\Sigma$  [ g/cm$^2$ ]')
cbar.set_ticks(cbar_ticks)
plt.grid(b=True)

# Create slider for the time
ax_time     = plt.axes([0.25, 0.1, 0.5, 0.03], axisbg='lightgoldenrodyellow')
slider_time = Slider(ax_time, 'timestep', 0.0, N_t, valinit=0,valfmt='%i')
ax._widgets = [slider_time] # Avoids garbage collection

def format_coord(theta, r):
    value = intpol(theta, r, sigma[ :, :].cgs)
    return r'theta=%1.4f rad, R=%1.4f AU, Sigma=%1.4f g/cm2' % (theta, r, value)

ax.format_coord = format_coord

# Create movie button
ax_button    = plt.axes([0.25, 0.04, 0.17, 0.04])
button_movie = Button(ax_button, 'Create movie', hovercolor='0.975')
ax._widgets += [button_movie] # Avoids garbage collection

plt.show()                  # Show plot


##### FUNCTIONS ####################################################################################################################
def update(val):
    i = int(np.floor(slider_time.val))  # Slider position
    
    try:
        sigma[:, :] = np.fromfile(data_dir+'/gasdens'+repr(i)+'.dat').reshape(N_R, N_theta) * aconst.M_sun / u.AU**2
    except:
        print 'Could not load gasdens'+repr(i)+'.dat'
    
    # Remove old collections and update contour plots
    for row in collections_list:
        ax.collections.remove(row)
        collections_list.remove(row)
    dummy = ax.contourf(theta_cen, R_cen.to(u.AU), np.log10(sigma[ :, :].cgs.value), cmap=amber_teal, levels=levels, extend='both' )
    for row in dummy.collections:
        collections_list.append(row)
    
    def format_coord_i(theta, r):
        value = intpol(theta, r, sigma[ :, :].cgs)
        return r'theta=%1.4f rad, R=%1.4f AU, Sigma=%1.4f g/cm2' % (theta, r, value)
    
    ax.format_coord = format_coord_i
    
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
        try:
            sigma[:, :] = np.fromfile(data_dir+'/gasdens'+repr(i)+'.dat').reshape(N_R, N_theta) * aconst.M_sun / u.AU**2
        except:
            print 'Could not load gasdens'+repr(i)+'.dat'
        movie_fig, movie_ax = plt.subplots(subplot_kw=dict(projection='polar'))
        movie_plot = movie_ax.contourf(theta_cen, R_cen.to(u.AU), np.log10(sigma[:, :].cgs.value), cmap=amber_teal, levels=levels, extend='both' )
        movie_planet, = movie_ax.plot(planet_theta[i], planet_R[i].to(u.AU), color='cyan', marker='o', markersize=8, markeredgecolor='black')
        movie_ax.set_rmin(-R_cen[0].to(u.AU).value)
        movie_ax.set_rmax(R_cen[-1].to(u.AU).value)
        movie_ax.tick_params(axis='x', labelbottom='off')
        movie_ax.tick_params(axis='y', colors='white')
        movie_cbar = movie_fig.colorbar(movie_plot)
        movie_cbar.ax.set_ylabel(r'log $\Sigma$  [ g/cm$^2$ ]')
        movie_cbar.set_ticks(cbar_ticks)
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
