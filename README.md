# particle_code
Particle code for FARGO-ADSG

This code is for calculating particle orbits in FARGO-ADSG simulations. It might also work for other FARGO versions, but it is not yet tested for them.



COMPILING
--------------------

To compile the code use

    make

If you want to compile the parallel version of the code use

    make BUILD=parallel
    
For additional debugging flags use

    make DEBUGMODE=1



EXECUTION
--------------------

To execute use

    ./particle_code input_file.nml
    
where you have to give the input file as first argument.

If you compiled in parallel mode use

    ./particle_code input_file.nml setthreads N
    
to run the code with N threads. If you don't give the number of threads, N=1 will be used.



WARNINGS
--------------------

The planet has to be on a fixed circular orbit! Other inertial forces that arise from a migrating planet on an eccentric orbit are not included, yet.



THE INPUT FILE
--------------------

The particle code needs an input file, whose location has to be given as command line argument. The input file is a simple FORTRAN namelist. It's variable are:

    character :: fargo_datadir

This is the location of the FARGO data directory relative to the path from which ./particle_code was called. The FARGO data needs to contain the files: dims.dat, used_rad.dat, planet0.dat, all gasdens###.dat, gasvrad###.dat, gasvtheta###.dat.


    character :: output_dir

The particle code saves its data files in this directory. The directory has to exist. There will be no checks, whether the directory is empty or not. Any existing files will be overwritten without mercy.
Every output file has one line of header that explains the meaning of the columns.


    integer :: i_start

The index of the first FARGO frame with which you want to start the simulation.


    integer :: i_stop

The index of the last FARGO frame with which you want to stop the simulation. All indices in between i_start and i_stop have to exist.


    double precision :: phys_dist

The physical distance in AU that represents one FARGO distance unit. The need for the Stokes numer of the particles breaks the scale-freeness of FARGO.


    double precision :: phys_mass
    
The physical mass in solar masses that represents one FARGO distance unit. This is also the mass of the central star.


    double precision :: adx
    
The adiabatic index of the gas.


    double precision :: alpha
    
The alpha viscosity parameter. This is used to add a random walk to the particles depending on their Stokes number due to turbulent motion of the gas. If this value is set to zero, no random walk will be added.


    double precision :: eps
    
The dust-to-gas ratio. For particle growth and fragmentation the code assumes a constant gas-to-dust ratio to calculate the collision rates.


    double presision :: mu
    
The mean molecular weight of the gas in proton masses.


    double precision :: smoothing

The smoothing length in units of pressure scale heights at the distance of the planet. This is the minimal distance a particle can have to the planet for the calculation of the gravitational potential in order to avoid division by zero.


    integer :: do_growth

Whether the code should do coagulation or not. 1: do coagulation, 0: don't do coagulation.


    integer :: do_frag

Whether the code should do fragmentation or not. 1: do fragmentation, 0: don't do fragmentation.


    double precision :: v_frag
    
The fragmentation velocity in cm/s. If the relative particle velocity is larger than v_frag the particles fragment. If not, they grow.


    double precision :: use_frag_ice
    
Do you want to use a different fragmentation velocity for ice particles?
0: no
1: yes.

    double precision :: v_frag_ice
    
Fragmentation velocity of ice particles in cm/s


    double precision :: T_ice
    
Temperature in K below which the particle is considered as ice particle.


    integer :: do_randomwalk

Whether the particles should do a random walk due to turbulence.
0: No.
1: Yes.


    double precision :: Ea
    
The activation energy of crystalization in units of the Boltzmann constant.
Meaning: Ea = activation energy / Boltzmann constant


    double precision :: nu_vib

The vibrational frequency of the solid state bonds in Hz


    double precision :: zeta

The initial volume fraction of crystal.


    integer :: N_dust

The number of dust particles you want to simulate


    double precision :: a_min

The minimum initial size of the particles in cm if a_dist_log=1.


    double precision :: a_max

The maximum initial particle size in cm if a_dist_log=1 and the monodisperse initial particle size in cm if a_dist_log=0.


    integer :: a_dist_log

The initial particle size distribution.
0: Monodisperse. All particles have size a_max
1: Logarithmically uniform between a_min and a_max


    double precision :: a_mono

The monomer size in cm. When do_fragmentation=1 this is the minimum size a particle can have.


    double precision :: rho_b

The bulk density of the particles in g/cm^3.


    integer :: dens_dist

The initial spatial distribution of the particles in the disk.
0: Uniformly in R-theta space. This leads naturally to a number density distribution of n(R) ~ 1/R
1: Following the surface density distribution of the starting frame i_start. This can be useful, if you want to start with a frame that already developed a gap. In this case you avoid placing particles in the gap. ATTENTION: This method uses rejection sampling to create the initial particle distribution. This can be very slow.
2: Uniformly in R-theta space in a ring centered on the planet of size R_ring.


    double precision :: R_ring

Size in AU of the ring for the initial particle distribution of dens_dist=2.



THE OUTPUT FILES
--------------------

The Output files dust0.dat of the code are saved in the output_dir directory. The first line is a header and explains the different columns. The columns are:

    # Identifier of the particle:
    # Radial distance in AU
    # Theta angle in rad
    # Radial velcity in cm/s
    # Azimuthal velocity in cm/s
    # X-coordinate in AU
    # Y-coordinate in AU
    # X-velocity in cm/s
    # Y-velocity in cm/s
    # Particle radius in cm
    # Particle mass in g
    # Stokes number
    # Stopping time in s
    # Temperature of the dust particle in K
    # Crystallinity fraction of the dust particle
    # Gas surface density at particle position in g/cm2



THE PLOTTING FILES
--------------------

    plot_density.py fargo_datadir
    
Plotting the FARGO density map assuming that the mass unit is 1 M_sun and the distance unit is 1 au.


    plot_temperature.py fargo_datadir
    
Plotting the FARGO temperature map assuming that the mass unit is 1 M_sun and the distance unit is 1 au and a mean moelcular weight of 2.3.


    plot_pressure.py fargo_datadir
    
Plotting the FARGO pressure map assuming that the mass unit is 1 M_sun and the distance unit is 1 au.


    plot_vorticity.py fargo_datadir
    
Plotting the vorticity assuming that the mass unit is 1 M_sun and the distance unit is 1 au.


    plot_vortensity.py fargo_datadir
    
Plotting the vortensity assuming that the mass unit is 1 M_sun and the distance unit is 1 au.


    plot_particles_radius.py fargo_datadir output_dir

Plots the particle positions in the disk colorcoded with their radius. Also shown are the positions of the planet and of the L4 and L5 Lagrange points.


    plot_particles_stokes.py fargo_datadir output_dir

Plots the particle positions in the disk colorcoded with their Stokes number. Also shown are the positions of the planet and of the L4 and L5 Lagrange points.


    plot_particles_crystallinity.py fargo_datadir output_dir

Plots the particle positions in the disk colorcoded with their crystallinity fraction. Also shown are the positions of the planet and of the L4 and L5 Lagrange points.
