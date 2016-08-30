# particle_code
Particle code for FARGO-ADSG

This code is for calculating particle orbits in FARGO-ADSG simulations. It might also work for other FARGO versions, but it is not yet tested for them.



COMPILING
--------------------

To compile the code use

    make

If you want to compile the parallel version of the code use

    make BUILD=parallel



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


    double precision :: alpha
    
The alpha viscosity parameter. This is used to add a random walk to the particles depending on their Stokes number due to turbulent motion of the gas. If this value is set to zero, no random walk will be added.


    double presision :: mu
    
The mean molecular weight of the gas in proton masses.


    double precision :: aspect

The aspect ratio H/R at the location of the FARGO distance R=1. The overall aspect ratio will be aspect(R) = aspect(1) * R^flaring_index


    double precision :: flaring_index

The flaring index of the disk. The disk's pressure scale height is then calculated as H(R) = aspect * R^(flaring_index+1)


    double precision :: smoothing

The smoothing length in units of pressure scale heights at the distance of the planet. This is the minimal distance a particle can have to the planet for the calculation of the gravitational potential in order to avoid division by zero.


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


    double precision :: rho_b

The bulk density of the particles in g/cm^3.


    integer :: dens_dist

The initial spatial distribution of the particles in the disk.
0: Uniformly in R-theta space. This leads naturally to a number density distribution of n~1/R^2
1: Following the surface density distribution of the starting frame i_start. This can be useful, if you want to start with a frame that already developed a gap. In this case you avoid placing particles in the gap. ATTENTION: This method uses rejection sampling to create the initial particle distribution. This can be very slow.
2: Uniformly in R-theta space in a ring centered on the planet of size R_ring.


    double precision :: R_ring

Size in AU of the ring for the initial particle distribution of dens_dist=2.



THE OUTPUT FILES
--------------------

The Output files dust0.dat of the code are saved in the output_dir directory. The first line is a header and explains the different columns. The columns are:

    # Identifier of the particle:
    # Radial distance
    # Theta angle
    # Radial velcity
    # Azimuthal velocity
    # X-coordinate
    # Y-coordinate
    # X-velocity
    # Y-velocity
    # Particle radius
    # Particle mass
    # Stokes number
    # Stopping time



THE PLOTTING FILES
--------------------

    plot_particles_stokes.py fargo_datadir output_dir

Plots the particle positions in the disk colorcoded with their Stokes number.
