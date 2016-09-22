module variables

    implicit none
    
! Parallel computing
    integer :: setthreads
    
! I/O files & directories
    
    ! Inputs namelist
    character(len=128) :: input_file
    ! FARGO data directory
    character(len=128) :: fargo_datadir
    ! Output directory
    character(len=128) :: output_dir
    
! Grid data

    ! Number of radial gridpoints
    integer :: N_R
    ! Number of azimuthal grid points
    integer :: N_theta
    ! Number of Fargo snapshots
    integer :: N_t
    
    ! Radial grid
    double precision, allocatable :: R(:)
    ! Interfaces of radial grid
    double precision, allocatable :: R_int(:)
    ! Theta grid
    double precision, allocatable :: theta(:)
    ! Interfaces of theta grid
    double precision, allocatable :: theta_int(:)
    
    ! Gas surface density
    double precision, allocatable :: sigma_gas(:, :, :)
    ! Gas radial velocity
    double precision, allocatable :: vR_gas(:, :, :)
    ! Gas azimuthal velocity
    double precision, allocatable :: vTheta_gas(:, :, :)
    ! Gas temperature
    double precision, allocatable :: T_gas(:, :, :)
    ! Adiabatic index
    double precision              :: adx
    
    ! Activation energy for crystallization
    double precision              :: Ea
    ! Vbrational frequency used for crystallization
    double precision              :: nu_vib
    ! Initial fraction of grwoth center
    double precision              :: zeta
    
    ! Planet position in cartesian coordinates
    double precision, allocatable :: x_planet(:)
    double precision, allocatable :: y_planet(:)
    ! Planet velocity in cartesian coordinates
    double precision, allocatable :: vX_planet(:)
    double precision, allocatable :: vY_planet(:)
    ! Planet position in polar coordinates
    double precision, allocatable :: R_planet(:)
    double precision, allocatable :: theta_planet(:)
    ! Planet velocity in polar coordinates
    double precision, allocatable :: vR_planet(:)
    double precision, allocatable :: vTheta_planet(:)
    ! Planet mass
    double precision, allocatable :: m_planet(:)
    
    ! Position and velocity of dust particles in cartesian coordinates
    double precision, allocatable :: x_dust(:)
    double precision, allocatable :: y_dust(:)
    double precision, allocatable :: vX_dust(:)
    double precision, allocatable :: vY_dust(:)
    ! Position and velocity of dust particles in polar coordinates
    double precision, allocatable :: R_dust(:)
    double precision, allocatable :: theta_dust(:)
    double precision, allocatable :: vR_dust(:)
    double precision, allocatable :: vTheta_dust(:)
    ! Specific angular momentum of dust
    double precision, allocatable :: L_dust(:)
    ! Mass of dust particles
    double precision, allocatable :: m_dust(:)
    ! Radius of dust particles
    double precision, allocatable :: a_dust(:)
    ! Temperature of dust particles
    double precision, allocatable :: T_dust(:)
    ! Crystallinity fraction of dust particles
    double precision, allocatable :: fc_dust(:)
    ! Local gas surface density at dust location
    double precision, allocatable :: loc_sigma_gas_dust(:)
    ! Change rate of particle radius
    double precision, allocatable :: dadt(:)
    ! Minimum and maximum initial size of dust particles and monomer radius
    double precision              :: a_min, a_max, a_mono
    ! Initial mass distribution of dust particles
    integer                       :: a_dist_log
    !Fragmentation velocity in cm/s
    double precision              :: v_frag
    ! Bulk density of dust particles
    double precision              :: rho_b
    ! Initial particle distribution
    integer                       :: dens_dist
    ! Box for particle distribution
    double precision :: R_ring
    ! Stokes number of particles
    double precision, allocatable :: St(:)
    ! Stopping time of particles
    double precision, allocatable :: tstop(:)
    
    ! Time of the snapshots
    double precision, allocatable :: time(:)
    
    ! Current values of variables
    double precision, allocatable :: cur_sigma_gas(:, :)
    double precision, allocatable :: cur_vR_gas(:, :)
    double precision, allocatable :: cur_vTheta_gas(:, :)
    double precision, allocatable :: cur_T_gas(:, :)
    double precision              :: cur_R_planet
    double precision              :: cur_theta_planet
    double precision              :: cur_X_planet
    double precision              :: cur_Y_planet
    double precision              :: cur_m_planet
    
! Simulation data
    
    ! First and last snapshot of dust simulation
    integer :: i_start, i_stop
    ! Number of dust particles
    integer :: N_dust
    ! Physical distance of the of the Fargo distance unit
    double precision :: phys_dist
    ! Physical mass of the of the Fargo mass unit
    double precision :: phys_mass
    ! Physical time of the of the Fargo time unit
    double precision :: phys_time
    ! Shall we do grwoth?
    integer :: do_growth
    ! Shall we do fragmentation?
    integer :: do_frag
    ! Shall we use a different fragmentation velocity for icy material
    integer :: use_frag_ice
    ! Fragmentation velocity of icy material
    double precision :: v_frag_ice
    ! Temperature below which we assume material to be ice
    double precision :: T_ice
    ! Shall we do a random walk?
    integer :: do_randomwalk
    
! Disk parameters

    ! Alpha viscosity parameter
    double precision :: alpha
    ! Mean molecular weight of gas
    double precision :: mu
    ! Smoothing length
    double precision :: smoothing
    ! Dust-to-gas ratio
    double precision :: eps

end module variables
