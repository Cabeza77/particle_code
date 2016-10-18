#!/usr/bin/env python

##### IMPORTING MODULES ############################################################################################################
from scipy import interp, interpolate
import astropy.constants as aconst
import astropy.units as u
import ConfigParser
import numpy as np
import os
import scipy.constants as const
import sys



##### CHECKING FOR ARGUMENTS #######################################################################################################
print ''
try:
    fargoDataDir = sys.argv[1]
except:
    sys.exit('You have to give the FARGO data directory as first argument.')
try:
    dustDataDir = sys.argv[2]
except:
    sys.exit('You have to give the dust data directory as second argument.')
try:
    outDir = sys.argv[3]
except:
    sys.exit('You have to give the output directory as third argument.')
try:
    configFile = sys.argv[4]
except:
    sys.exit('You have to give the configuration file as forth argument.')
# Check if output directory already exists and is empty.
if(os.path.isdir(outDir)):
    if(not os.listdir(outDir)==[]):
        cont = raw_input('Output directory \''+outDir+'\' is not empty. Do you still want to continue? [n/Y] \n')
        print ''
        if(cont != 'Yes' and cont != 'Y'):
            sys.exit('Program aborted.')
else:
    cont = raw_input('Output directory \''+outDir+'\' does not exist. Shall I create it? [n/Y] \n')
    print ''
    if(cont != 'Yes' and cont != 'Y'):
        sys.exit('Program aborted.')
    else:
        os.mkdir(outDir)
# Loading config file
try:
    config = ConfigParser.ConfigParser()
    config.read(configFile)
except:
    sys.exit('Could not load \''+configFile+'\'.')
    
    
    
##### LOADING DATA #################################################################################################################
print 'Loading data...'
# Reading data from config file
try:
    i_snapshot     = config.getint   ('Simulation', 'i_snapshot')
    mu             = config.getfloat ('Simulation', 'mu') * u.g / u.mole
    M_disk         = config.getfloat ('Simulation', 'M_disk') * aconst.M_sun
    d2g            = config.getfloat ('Simulation', 'd2g') * u.one
    D              = config.getfloat ('Simulation', 'D') * u.au
    alpha          = config.getfloat ('Simulation', 'alpha') * u.one
    R_star         = config.getfloat ('Simulation', 'R_star') * aconst.R_sun
    M_star         = config.getfloat ('Simulation', 'M_star') * aconst.M_sun
    T_star         = config.getfloat ('Simulation', 'T_star') * u.K
    
    N_R            = config.getint  ('Grids',      'N_R')
    N_R_refine     = config.getint  ('Grids',      'N_R_refine')
    i_R_refine     = config.getint  ('Grids',      'i_R_refine')
    R_min          = config.getfloat('Grids',      'R_min') * u.au
    R_max          = config.getfloat('Grids',      'R_max') * u.au
   
    N_theta_hem    = config.getint  ('Grids',      'N_theta_hem')
    N_theta_coarse = config.getint  ('Grids',      'N_theta_coarse')
    theta_coarse   = config.getfloat('Grids',      'theta_coarse') * u.one
  
    N_phi          = config.getint  ('Grids',      'N_phi')
   
    lambda_1       = config.getfloat('Grids',      'lambda_1') * u.um
    lambda_2       = config.getfloat('Grids',      'lambda_2') * u.um
    lambda_3       = config.getfloat('Grids',      'lambda_3') * u.um
    lambda_4       = config.getfloat('Grids',      'lambda_4') * u.um
    N_lambda_1     = config.getint  ('Grids',      'N_lambda_1')
    N_lambda_2     = config.getint  ('Grids',      'N_lambda_2')
    N_lambda_3     = config.getint  ('Grids',      'N_lambda_3')
   
    N_phot         = config.getint  ('RadMC3D',    'N_phot')
    N_phot_scat    = config.getint  ('RadMC3D',    'N_phot_scat')
    scat_mode_max  = config.getint  ('RadMC3D',    'scat_mode_max')
    N_spec         = config.getint  ('RadMC3D',    'N_spec')
    a_min          = config.getfloat('RadMC3D',    'a_min') * u.cm
    a_max          = config.getfloat('RadMC3D',    'a_max') * u.cm
    rho_mat        = config.getfloat('RadMC3D',    'rho_mat') * u.g / u.cm**3.
    
    M_fargo        = config.getfloat('FARGO',      'M_FARGO') * aconst.M_sun
    d_fargo        = config.getfloat('FARGO',      'd_FARGO') * u.au
except:
    sys.exit('\''+configFile+'\' has an invalid format.')
# Computing total dust disk mass
M_dust = M_disk * d2g
# Loading dims.dat
try:
    dims = np.loadtxt(fargoDataDir+'/dims.dat')
except:
    sys.exit('Could not load dims.dat')
# Loading used_rad.dat
try:
    used_rad = np.loadtxt(fargoDataDir+'/used_rad.dat') * d_fargo
except:
    sys.exit('Could not load used_rad.dat!')
# Loading used_azi.dat
try:
    used_azi = np.loadtxt(fargoDataDir+'/used_azi.dat') * u.rad
except:
    sys.exit('Could not load used_azi.dat!')
# Extract fundamental parameters
N_R_fargo          = int( dims[6] )                                       # Number of radial gridpoint
N_phi_fargo        = int( dims[7] )                                       # Number of angular grid points
R_cen_fargo        = 0.5 * (used_rad[:-1] + used_rad[1:] )                # Calculate radial cell centers
phi_cen_fargo      = np.zeros(N_phi_fargo) * u.rad
phi_cen_fargo[:-1] = 0.5 * (used_azi[:-1] + used_azi[1:] )                # Calculate azimuthal cell centers
phi_cen_fargo[-1]  = 0.5 * (2.*const.pi*u.rad + used_azi[-1] )
# Load the gas density
sigmaGasFargo = np.zeros( (N_R_fargo, N_phi_fargo) ) * u.g / u.cm**2      # Initialize density array
try:
    sigmaGasFargo[:, :] = np.fromfile(fargoDataDir+'/gasdens'+repr(i_snapshot)+'.dat').reshape(N_R_fargo, N_phi_fargo) * M_fargo / d_fargo**2.
except:
    sys.exit( 'Could not load gasdens'+repr(i_snapshot)+'.dat' )
# Load the temperature
T_fargo = mu / aconst.R * aconst.G * M_fargo / d_fargo
tempFargo = np.zeros( (N_R_fargo, N_phi_fargo) ) * u.K  # Initialize temperature array
try:
    tempFargo[:, :] = np.fromfile(fargoDataDir+'/gasTemperature'+repr(i_snapshot)+'.dat').reshape(N_R_fargo, N_phi_fargo) * T_fargo
except:
    sys.exit( 'Could not load gasTemperature'+repr(i_snapshot)+'.dat' )
# Loading dust
try:
    dust   = np.loadtxt(dustDataDir+'/dust'+repr(i_snapshot)+'.dat', comments='#')
    N_dust = dust.shape[0]
except:
    sys.exit( 'Could not load dust'+repr(i_snapshot)+'.dat' )
print '   ...Done.'



##### PROCSSING DATA ###############################################################################################################
print ''
print 'Processing data...'
# The grid cell interfaces
# Radial grid
if(N_R_refine == 0):
    RInt = np.logspace( np.log10(R_min.to(u.au).value), np.log10(R_max.to(u.au).value), N_R+1, 10. ) * u.au
else:
    N_R_fine = N_R_refine * 2**i_R_refine
    coarse = np.logspace( np.log10(R_min.to(u.au).value), np.log10(R_max.to(u.au).value), num=N_R-N_R_fine+2, base=10. )
    fine   = np.logspace( np.log10(R_min.to(u.au).value), np.log10(coarse[1]), num=N_R_fine, base=10., endpoint=False )
    RInt   = np.concatenate( (fine, coarse[1:]), axis=0 ) * u.au
# Theta grid.
N_theta     = 2 * N_theta_hem
coarse      = np.linspace( 0., theta_coarse.value, N_theta_coarse, endpoint=False )
fine        = np.linspace(theta_coarse.value, const.pi/2., N_theta_hem-N_theta_coarse+1, endpoint=True)
theta_north = np.concatenate( (coarse, fine), axis=0 )
theta_south = theta_north[::-1][1:]
theta_south = const.pi - theta_south
thetaInt    = np.concatenate( (theta_north, theta_south), axis=0 ) * u.rad
# Phi grid
phiInt = np.linspace(0., 2.*const.pi, N_phi+1) * u.rad
# The particles
aInt = np.logspace( np.log10(a_min.cgs.value), np.log10(a_max.cgs.value), N_spec+1, 10. ) * u.cm

# The cell centers
RCen     = 0.5 * (     RInt[:-1] +     RInt[1:] )
thetaCen = 0.5 * ( thetaInt[:-1] + thetaInt[1:] )
phiCen   = 0.5 * (   phiInt[:-1] +   phiInt[1:] )
aCen     = 0.5 * (     aInt[:-1] +     aInt[1:] )

# The wavelength grid
lam1 = np.logspace( np.log10(lambda_1.to(u.um).value), np.log10(lambda_2.to(u.um).value), num=N_lambda_1, base=10., endpoint=False )
lam2 = np.logspace( np.log10(lambda_2.to(u.um).value), np.log10(lambda_3.to(u.um).value), num=N_lambda_2, base=10., endpoint=False )
lam3 = np.logspace( np.log10(lambda_3.to(u.um).value), np.log10(lambda_4.to(u.um).value), num=N_lambda_3, base=10., endpoint=True  )
lamb = np.concatenate( (lam1, lam2, lam3), axis=0 ) * u.um

# Kernel function to spread the particles
kernel = lambda R: 1. / (const.pi*D**2.) * np.exp( - R**2. / D**2. )

# Creating a 2D grid in R, phi and X, Y space for vector based calcualtions
phi_grid,       R_grid       = np.meshgrid( phiCen.cgs.value, RCen.to(u.au).value )
phi_grid_fargo, R_grid_fargo = np.meshgrid( phi_cen_fargo.cgs.value, R_cen_fargo.to(u.au).value ) # Used for griddata function
X_grid = R_grid * u.au * np.cos( phi_grid * u.rad )
Y_grid = R_grid * u.au * np.sin( phi_grid * u.rad )

# Sum up the dust particles within the final RadMC3D grid including the kernel
sigmaDustRadmc3d = np.zeros( (N_spec, N_R, N_phi) ) * u.g / u.cm**2.
sigmaDustRadmc3d[:, :, :] = 1.e-100 * u.g / u.cm**2
for idust in range(N_dust):
    R_dust = dust[idust, 1]  * u.au
    if(R_dust <= 0.0): continue # Don't include particles outside of the grid
    X_dust = dust[idust, 5]  * u.au
    Y_dust = dust[idust, 6]  * u.au
    R      = np.sqrt( (X_grid-X_dust)**2. + (Y_grid-Y_dust)**2. + (0.1*u.au)**2. )
    m_dust = dust[idust, 10] * u.g
    a_dust = dust[idust, 9]  * u.cm
    ia     = np.argmax( aInt > a_dust ) - 1
    sigmaDustRadmc3d[ia, :, :] = sigmaDustRadmc3d[ia, :, :] + m_dust * kernel(R)

# Average azimuthal grid cell size
phiAverage = np.average( phiInt[1:].cgs.value - phiInt[:-1].cgs.value ) * u.rad
# Calculating the grid cell sizes of the 2D and the 3D grid.
#A_grid_radmc3d = np.zeros( (N_R, N_phi) ) * u.cm**2.
V_grid_radmc3d = np.zeros( (N_R, N_theta, N_phi) ) * u.cm**3.
for iR in range(N_R):
    #A_grid_radmc3d[iR, :] = 0.5 * ( RInt[iR+1]**2. - RInt[iR]**2. ) * phiAverage
    for itheta in range(N_theta):
        for iphi in range(N_phi):
            V_grid_radmc3d[iR, itheta, iphi] = 1./3. * ( RInt[iR+1]**3. - RInt[iR]**3. ) * ( np.cos(thetaInt[itheta]) - np.cos(thetaInt[itheta+1]) ) * ( phiInt[iphi+1].cgs.value - phiInt[iphi].cgs.value )

## Computing total mass of all dust    
#M_part_tot = np.sum( np.sum( sigmaDustRadmc3d[:, :, :], axis=0 ) * A_grid_radmc3d )
## Convert dust surface density to desired dust disk mass
#sigmaDustRadmc3d = sigmaDustRadmc3d * M_dust / M_part_tot

# Interpolating temperature to RadMC3D grid
tempRadmc3d = interpolate.griddata( (phi_grid_fargo.flatten(), R_grid_fargo.flatten()), tempFargo.cgs.value.flatten(), (phi_grid, R_grid), method="nearest" ) * u.K
# Interpolating surface density to RadMC3D grid
sigmaGasRadmc3d = interpolate.griddata( (phi_grid_fargo.flatten(), R_grid_fargo.flatten()), sigmaGasFargo.cgs.value.flatten(), (phi_grid, R_grid), method="nearest" ) * u.g / u.cm**2.
# Calculating the sound speed on the RadMC3D grid
csRadmc3d = np.sqrt( aconst.k_B * tempRadmc3d / (mu.cgs.value*aconst.m_p) )
# Calculating the pressure scale height on the RadMC3D grid
HRadmc3d = csRadmc3d / np.sqrt( aconst.G * M_fargo / (R_grid*u.au)**3. )

# Function for calculating the Stokes number (only Epstein regime)
get_St = lambda a, sigma: const.pi / 2. * a * rho_mat / sigma

# Scale height of dust species
HDustRadmc3d = np.zeros( (N_spec, N_R, N_phi) ) * u.au
for ia in range( N_spec ):
    St = get_St( aCen[ia], sigmaGasRadmc3d )
    HDustRadmc3d[ia, :, :] = HRadmc3d * np.minimum( 1., np.sqrt( alpha / ( np.minimum( St, 0.5 ) * ( 1. + St**2. ) ) ) )

# z coordinate
zRadmc3d = np.zeros( (N_R, N_theta) ) * u.au
for iR in range(N_R):
    zRadmc3d[iR, :] = RCen[iR] * ( const.pi/2. - thetaCen.cgs.value ) # theta starts at North pole

# The actual RadMC3D density
rho_dust_radmc3d = np.zeros( (N_spec, N_R, N_theta, N_phi) ) * u.g / u.cm**3.
for ia in range(N_spec):
    for itheta in range(N_theta):
        for iphi in range(N_phi):
            # Enforcing a minimum value in case you want to plot de density logarithmically
            rho_dust_radmc3d[ia, :, itheta, iphi] = np.maximum( sigmaDustRadmc3d[ia, :, iphi] / ( np.sqrt(2.*const.pi) * HDustRadmc3d[ia, :, iphi] ) * np.exp( - 0.5 * (zRadmc3d[:, itheta]/HDustRadmc3d[ia, :, iphi]).cgs.value**2. ), 1.e-100 * u.g/u.cm**3. )
            
# Computing total mass of all dust    
M_tot = np.sum( np.sum( rho_dust_radmc3d, axis=0 ) * V_grid_radmc3d )
# Convert dust density to desired dust disk mass
rho_dust_radmc3d = rho_dust_radmc3d * M_dust / M_tot

## Coordinates of the star
coord_star = np.array([0., 0., 0.]) * u.cm
print '   ...Done.'



##### WRITING FILES ################################################################################################################
print ''
print 'Writing files...'

# amr_grid.inp
print r'   * amr_grid.inp'
f = open( outDir+'amr_grid.inp', 'w' )
f.write( '1' )
f.write( os.linesep )
f.write( '0' )
f.write( os.linesep )
f.write( '100' )
f.write( os.linesep )
f.write( '0' )
f.write( os.linesep )
if(N_phi > 1):
    f.write( '1  1  1' )
else:
    f.write( '1  1  0' )
f.write( os.linesep )
f.write( '{:d}  {:d}  {:d}'.format(len(RCen), len(thetaCen), len(phiCen) ) )
f.write( os.linesep )
for R in RInt:
    f.write( '{:13.6e}'.format( R.cgs.value ) )
    f.write( os.linesep )
for theta in thetaInt:
    f.write( '{:13.6f}'.format( theta.cgs.value ) )
    f.write( os.linesep )
for phi in phiInt:
    f.write( '{:13.6f}'.format( phi.cgs.value ) )
    f.write( os.linesep )
f.close()

# dust_density.inp
print r'   * dust_density.inp'
f = open( outDir+'dust_density.inp', 'w' )
f.write( '1' )
f.write( os.linesep )
f.write( '{:d}'.format( len(RCen)*len(thetaCen)*len(phiCen) ) )
f.write( os.linesep )
f.write( '{:d}'.format( len(aCen) ) )
f.write( os.linesep )
for ia in range(N_spec):
    print '      - Species {:d} of {:d}'.format(ia+1, N_spec)
    for iphi in range(len(phiCen)):
        for itheta in range(len(thetaCen)):
            for iR in range(len(RCen)):
                f.write( '{:13.6e}'.format( rho_dust_radmc3d[ia, iR, itheta, iphi].cgs.value ) )
                f.write( os.linesep )
f.close()

# radmc3d.inp
print r'   * radmc3d.inp'
f = open( outDir+'radmc3d.inp', 'w' )
f.write( 'nphot                    = {:d}'.format( N_phot ) )
f.write( os.linesep )
f.write( 'nphot_scat               = {:d}'.format( N_phot_scat ) )
f.write( os.linesep )
f.write( 'scattering_mode_max      = {:d}'.format( scat_mode_max ) )
f.write( os.linesep )
f.write( 'istar_sphere             = 1' )
f.write( os.linesep )
f.close()

# stars.inp
print r'   * stars.inp'
f = open( outDir+'stars.inp', 'w' )
f.write( '2' )
f.write( os.linesep )
f.write( '1  {:d}'.format( len(lamb) ) )
f.write( os.linesep )
f.write( '{:13.6e}  {:13.6e}  {:13.6e}  {:13.6e}  {:13.6e}'.format( R_star.cgs.value, M_star.cgs.value, coord_star[0].cgs.value, coord_star[1].cgs.value, coord_star[2].cgs.value ) )
f.write( os.linesep )
for lam in lamb:
    f.write( '{:13.6e}'.format( lam.to(u.um).value ) )
    f.write( os.linesep )
f.write( os.linesep )
f.write( '{:13.6f}'.format( -T_star.cgs.value ) )
f.close()

# wavelength_micron.inp
print r'   * wavelength_micron.inp'
f = open( outDir+'wavelength_micron.inp', 'w' )
f.write( '{:d}'.format( len(lamb) ) )
f.write( os.linesep )
for lam in lamb:
    f.write( '{:13.6e}'.format( lam.to(u.um).value ) )
    f.write( os.linesep )
f.close()

print '   ...Done.'
