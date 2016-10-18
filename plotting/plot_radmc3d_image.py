#!/usr/bin/env python

##### IMPORTING MODULES ############################################################################################################
from astropy.convolution import convolve, Gaussian2DKernel
from matplotlib.colors import LinearSegmentedColormap
import astropy.constants as aconst
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import radmc3dPy
import scipy.constants as const
import sys



##### CHECKING FOR ARGUMENTS #######################################################################################################
# Print empty line
print ''
try:
    imageFile = sys.argv[1]
except:
    sys.exit('You have to give the image file as first argument.')
try:
    FWHM = float(sys.argv[2]) * u.au
except:
    sys.exit('You have to give the FWHM (in au) as second argument.')
    
    
    
##### LOADING DATA #################################################################################################################
print 'Loading data...'
try:
    image = radmc3dPy.image.readImage(imageFile)
except:
    sys.exit('Could not read \''+imageFile+'\'.')
print '   ...Done.'



##### PROCESSING DATA ##############################################################################################################
print ''
print 'Processing data...'
X = image.x * u.cm
Y = image.y * u.cm
intensity = np.maximum( image.image[:, :, 0], 1.e-102 ) * u.erg / u.cm**2 / u.s / u.Hz / u.steradian
# Convolve image only if FWHM is larger than 0.
if(FWHM.value > 0.):
    stddev = (image.sizepix_x * u.cm)**(-1) * FWHM / 2. # Standard deviation
    kernel = Gaussian2DKernel(stddev.cgs.value)
    intensityConv = convolve(intensity.cgs.value, kernel, boundary='extend') * u.erg / u.cm**2 / u.s / u.Hz / u.steradian
else:
    intensityConv = intensity
print '   ...Done.'



##### PLOTTING #####################################################################################################################
print ''
print 'Plotting...'
# Defining amber-teal colormap
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
amber_teal = LinearSegmentedColormap('AmberTeal1', color_dict)
amber_teal.set_under('#191919')
amber_teal.set_over('#FFFFFF')
cmap='hot'

# Get global limits
margin    = 6.
intMax    = np.max( intensity.cgs.value )
intMin    = 10.**( np.ceil( np.log10( intMax ) ) - margin )
levelsInt = np.linspace( np.floor( np.log10( intMin ) ), np.ceil( np.log10( intMax ) ), 100 )    # Levels of intensity colorbar
cbarTicks = np.arange( np.floor( np.log10( intMin ) ), np.ceil( np.log10( intMax ) ) + 0.1, 1. ) # Ticks of colorbar

# Plotting
fig = plt.figure()

ax1 = fig.add_subplot(1, 2, 1, aspect='equal')
ax2 = fig.add_subplot(1, 2, 2, aspect='equal')

tickcolor = '#FFFFFF'
# Change color of major ticks to white
plt.setp(ax1.get_xticklines(), color=tickcolor)
plt.setp(ax1.get_yticklines(), color=tickcolor)
# Change color of minor ticks to white
ax1.tick_params(axis='both', which='minor', colors=tickcolor)
ax2.tick_params(axis='both', which='minor', colors=tickcolor)

ax1.set_title('Intensity map')
ax1.set_xlabel('X [AU]')
ax1.set_ylabel('Y [AU]')
ax1.set_xlim(X[0].to(u.au).value, X[-1].to(u.au).value)
ax1.set_ylim(Y[0].to(u.au).value, Y[-1].to(u.au).value)
ax1.grid(b=False)

ax2.set_title('Intensity map with FWHM = {:4.2f} AU'.format( FWHM.to(u.au).value ) )
ax2.set_xlabel('X [AU]')
ax2.set_ylabel('Y [AU]')
ax2.set_xlim(X[0].to(u.au).value, X[-1].to(u.au).value)
ax2.set_ylim(Y[0].to(u.au).value, Y[-1].to(u.au).value)
ax2.grid(b=False)

plot1 = ax1.contourf(X.to(u.au), Y.to(u.au), np.log10(intensity.cgs.value), cmap=cmap, levels=levelsInt, extend='both')
plot2 = ax2.contourf(X.to(u.au), Y.to(u.au), np.log10(intensityConv.cgs.value), cmap=cmap, levels=levelsInt, extend='both')

# Remove lines between the color shades. Important for exporting as .pdf-file
for c in plot1.collections:
    c.set_edgecolor("face")
for c in plot2.collections:
    c.set_edgecolor("face")
    
# Show colorbar
cbar1 = fig.colorbar(plot1, ax=ax1)
cbar1.ax.set_ylabel(r'log Intensity  [ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$ ]')
cbar1.set_ticks(cbarTicks)
cbar1.solids.set_edgecolor("face")
cbar2 = fig.colorbar(plot2, ax=ax2)
cbar2.ax.set_ylabel(r'log Intensity  [ erg cm$^{-2}$ s$^{-1}$ Hz$^{-1}$ ster$^{-1}$ ]')
cbar2.set_ticks(cbarTicks)
cbar2.solids.set_edgecolor("face")

plt.show()

# Maximize window and out it to foreground
try:
    figManager = plt.get_current_fig_manager()
    figManager.window.showMaximized()
    figManager.window.raise_()
except:
    pass
print '   ...Done.'
