#!/usr/bin/env python

##### IMPORTING MODULES ############################################################################################################
from bhmie_herbert_kaiser_july2012 import bhmie
from scipy import interp
import astropy.constants as aconst
import astropy.units as u
import ConfigParser
import f90nml
import numpy as np
import os
import scipy.constants as const
import sys



##### CHECKING FOR ARGUMENTS #######################################################################################################
# Print empty line
print ''
try:
    outDir = sys.argv[1]
except:
    sys.exit('You have to give the output directory as first argument.')
try:
    configFile = sys.argv[2]
except:
    sys.exit('You have to give the configuration file as second argument.')
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
    N_spec         = config.getint  ('RadMC3D',  'N_spec')
    a_min          = config.getfloat('RadMC3D',  'a_min') * u.cm
    a_max          = config.getfloat('RadMC3D',  'a_max') * u.cm
    material       = config.get     ('RadMC3D',    'material')
    rho_mat        = config.getfloat('RadMC3D',    'rho_mat') * u.g / u.cm**3
    N_ang          = config.getint  ('RadMC3D',    'N_ang')
    
    lambda_1       = config.getfloat('Grids',      'lambda_1') * u.um
    lambda_2       = config.getfloat('Grids',      'lambda_2') * u.um
    lambda_3       = config.getfloat('Grids',      'lambda_3') * u.um
    lambda_4       = config.getfloat('Grids',      'lambda_4') * u.um
    N_lambda_1     = config.getint  ('Grids',      'N_lambda_1')
    N_lambda_2     = config.getint  ('Grids',      'N_lambda_2')
    N_lambda_3     = config.getint  ('Grids',      'N_lambda_3')
except:
    sys.exit('\''+configFile+'\' has an invalid format.')
if(N_ang > 1000):
    sys.exit('N_ang has to be smaller than 1000.')
# Loading the optical constants
scriptDir = os.path.dirname(__file__)
try:
    dataM = np.loadtxt(scriptDir+'/opacities/'+material+'.lnk')
except:
    sys.exit('Could not load \'opacities/'+material+'.lnk\'.')
print '   ...Done.'



##### PROCSSING DATA ###############################################################################################################
print ''
print 'Processing data...'

# The particle grid
aInt = np.logspace( np.log10(a_min.cgs.value), np.log10(a_max.cgs.value), N_spec+1, 10. ) * u.cm
aCen = 0.5 * ( aInt[:-1] + aInt[1:] )

# The wavelength grid
lam1 = np.logspace( np.log10(lambda_1.to(u.um).value), np.log10(lambda_2.to(u.um).value), num=N_lambda_1, base=10., endpoint=False )
lam2 = np.logspace( np.log10(lambda_2.to(u.um).value), np.log10(lambda_3.to(u.um).value), num=N_lambda_2, base=10., endpoint=False )
lam3 = np.logspace( np.log10(lambda_3.to(u.um).value), np.log10(lambda_4.to(u.um).value), num=N_lambda_3, base=10., endpoint=True  )
lamb = np.concatenate( (lam1, lam2, lam3), axis=0 ) * u.um

# Angle grid
angle = np.linspace(0., const.pi, 2*N_ang-1, endpoint=True) * u.rad
mu    = np.cos(angle) * u.one

# Interpolation of the refractive index onto the wavelength grid
m_re = interp( lamb.to(u.um), dataM[:, 0], dataM[:, 1] )
m_im = interp( lamb.to(u.um), dataM[:, 0], dataM[:, 2] )
m    = m_re + 1.j * m_im    # This is the complex refractive index m = Re(m) + i*Im(m)
print '   ...Done.'



##### CALCULATING OPACITIES ########################################################################################################
print ''
print 'Calculating opacities...'
# Calculating the opacities for the particle species using bhmie
kappa_abs = np.zeros( (len(aCen), len(lamb)) ) * u.cm**2 / u.g
kappa_sca = np.zeros( (len(aCen), len(lamb)) ) * u.cm**2 / u.g
g_sca     = np.zeros( (len(aCen), len(lamb)) ) * u.one
xLim      = 1.e1 # The maximum size parameter for which the mie code is executed. For larger we only assume the geometric cross section.
Z         = np.zeros( (6, len(aCen), len(lamb), 2*N_ang-1) ) * u.cm**2 / u.g
for ia in range(len(aCen)):
    print '   * for species {:d} of {:d} with a = {:4.2e}.'.format(ia+1, len(aCen), aCen[ia].cgs)
    mass    = 4./3. * const.pi * rho_mat * aCen[ia]**3. # Mass of the particle
    sigGeo  = const.pi * aCen[ia]**2.                   # Geometric cross section of the particle
    for ilam in range(len(lamb)):
        # Size parameter x = k * a
        x = (2.*const.pi / lamb[ilam] * aCen[ia]).cgs.value
        if( (x<1000.) and (aCen[ia]<100.*u.um)  ): # Ignore scattering for x>1000 and a>0.1mm.
            convFactor = ( lamb[ilam]/(2.*const.pi) )**2. / mass # Conversion factor to convert the bhmie output to RadMC3D
            # bhmie returns [ S1, S2, Qext, Qsca, Qback, gsca ]
            S1, S2, Qext, Qsca, Qback, gsca = bhmie(x, m[ilam], N_ang)
            Qabs = Qext - Qsca
            S11 = 0.5 * ( np.absolute(S1[:])**2. + np.absolute(S2[:])**2. )
            S12 = 0.5 * ( np.absolute(S2[:])**2. - np.absolute(S1[:])**2. )
            S33 = ( S1[:] * S2[:].conjugate() ).real
            S34 = ( S1[:] * S2[:].conjugate() ).imag
            Z[0, ia, ilam, :] = S11 * convFactor
            Z[1, ia, ilam, :] = S12 * convFactor
            Z[2, ia, ilam, :] = S11 * convFactor
            Z[3, ia, ilam, :] = S33 * convFactor
            Z[4, ia, ilam, :] = S34 * convFactor
            Z[5, ia, ilam, :] = S33 * convFactor
        else:
            Qabs = 1.
            Qsca = 0.
            gsca = 0.
            Z[:, ia, ilam, :] = 0. * u.cm**2 / u.g
        kappa_abs[ia, ilam] = Qabs * sigGeo / mass
        kappa_sca[ia, ilam] = Qsca * sigGeo / mass
        g_sca    [ia, ilam] = gsca
        # Checking the validity of the Zs here
        # The integral over S11 has to be kappa_sca
        if(kappa_sca[ia, ilam] > 0.):
            checksum = 0. * u.cm**2 / u.g
            for i in range(1, 2*N_ang-1):
                checksum += 0.5 * ( Z[0, ia, ilam, i-1] + Z[0, ia, ilam, i] ) * np.absolute( mu[i] - mu[i-1] )
            checksum = checksum * 2. * const.pi
            error = np.absolute( checksum / kappa_sca[ia, ilam] - 1. )
            scaleFactor = kappa_sca[ia, ilam] / checksum
            Z[:, ia, ilam, :] = Z[:, ia, ilam, :] * scaleFactor
print '   ...Done.'



##### WRITING FILES ################################################################################################################
print ''
print 'Writing files...'

# dustopac.inp
print r'   * dustopac.inp'
f = open( outDir+'dustopac.inp', 'w' )
f.write( '2' )
f.write( os.linesep )
f.write( '{:d}'.format( N_spec ) )
f.write( os.linesep )
f.write( '------------------------------' )
f.write( os.linesep )
for ia in range(N_spec):
    f.write( '10' )
    f.write( os.linesep )
    f.write( '0' )
    f.write( os.linesep )
    f.write( 'spec{:s}'.format( str(ia).zfill(3) ) )
    f.write( os.linesep )
    f.write( '------------------------------' )
    f.write( os.linesep )
f.close()

# dustkapscatmat_spec*.inp
for ia in range(N_spec):
    print r'   * dustkapscatmat_spec'+repr(ia).zfill(3)+'.inp'
    f = open( outDir+'dustkapscatmat_spec'+repr(ia).zfill(3)+'.inp', 'w' )
    f.write( '# Opacities created with:' )
    f.write( os.linesep )
    f.write( '# bhmie_herbert_kaiser_july2012.py' )
    f.write( os.linesep )
    f.write( '# ( http://scatterlib.googlecode.com/files/bhmie_herbert_kaiser_july2012.py ).' )
    f.write( os.linesep )
    f.write( '#' )
    f.write( os.linesep )
    f.write( '# Optical constants from {:s}.lnk.'.format( material ) )
    f.write( os.linesep )
    f.write( '# Grain size: {:13.6e} cm.'.format( aCen[ia].cgs.value ) )
    f.write( os.linesep )
    f.write( '# Bulk density: {:13.6f} g/cm3.'.format( rho_mat.cgs.value ) )
    f.write( os.linesep )
    f.write( '#' )
    f.write( os.linesep )
    f.write( '1' )
    f.write( os.linesep )
    f.write( '{:d}'.format( len(lamb) ) )
    f.write( os.linesep )
    f.write( '{:d}'.format( 2*N_ang-1 ) )
    f.write( os.linesep )
    f.write( os.linesep )
    for ilam in range(len(lamb)):
        f.write( '{:13.6e}  {:13.6e}  {:13.6e}  {:13.6e}'.format( lamb[ilam].to(u.um).value, kappa_abs[ia, ilam].cgs.value, kappa_sca[ia, ilam].cgs.value, g_sca[ia, ilam].cgs.value ) )
        f.write( os.linesep )       
    f.write( os.linesep )
    for i in range(2*N_ang-1):
        f.write( '{:13.6f}'.format( angle[i].to(u.deg).value ) )
        f.write( os.linesep )
    f.write( os.linesep )
    for ilam in range(len(lamb)):
        for i in range(2*N_ang-1):
            f.write( '{:13.6e}  {:13.6e}  {:13.6e}  {:13.6e}  {:13.6e}  {:13.6e}'.format( Z[0, ia, ilam, i].cgs.value, Z[1, ia, ilam, i].cgs.value, Z[2, ia, ilam, i].cgs.value, Z[3, ia, ilam, i].cgs.value, Z[4, ia, ilam, i].cgs.value, Z[5, ia, ilam, i].cgs.value ) )
            f.write( os.linesep )
        f.write( os.linesep )
    f.close()

print '   ...Done.'
