program particle_code

    use io
    use memory
    use variables
    use constants
    use interpolation
    use aerodynamics
    use coagulation

    implicit none
    
    ! Command line
    integer :: num_args
    
    ! Iteratives
    integer :: i
    integer :: iframe
    
    ! Random number generator in nrecip.f90
    double precision :: ran2
    ! Seed for random generator. Has to be negative
    integer :: iseed
    ! Function to check for validity
    logical :: number_invalid
    
    ! Time keeping
    integer :: time_start, time_finish, time_rate
    
    ! Dummies
    double precision :: dummy, dummy2, dr_dum, dt_dum
    !$ character(32)    :: dumstr
    
    ! Chunksize for loops in multithreading
    !$ integer :: chunksize
    
    ! Local gas parameters
    double precision :: loc_vR_gas, loc_vTheta_gas
    ! Local gravitational acceleration parameter
    double precision :: loc_aR_grav, loc_aTheta_grav
    
    ! Needed for the force calculations
    double precision :: d_pldu
    double precision :: Fx, Fy, Fr, Ftheta
    double precision :: H, nu_visc, l_visc, ran_theta
    double precision :: Xcms, Ycms, Rcms, thetaCms ! Center of mass coordinates
    double precision :: Fxin, Fyin                 ! Inertial forces
    
    ! Major time variable
    double precision :: cur_time
    ! Time step
    double precision :: dt
    
    ! Is grid logarithmic?
    logical :: is_grid_log = .FALSE.
    
    call write_welcome()
    
    num_args = command_argument_count()
    if(num_args .LT. 1) then
        call write_stop_message('You have to give the input file as command line argument')
        stop
    end if
    call get_command_argument(1, input_file)
    
    ! Start timer
    call system_clock( time_start, time_rate )
    
    ! Set seed for random generator
    iseed = -time_start

! ##### READ INPUTS AND ALLOCATE MEMORY ############################################################################################
    
    write(*,*) '# Reading inputs...'
    
    ! Read the inputs
    call read_inputs()
    
    ! Allocate memory
    call allocate_memory()
    
    ! Read data from files
    call read_data()
    
    if( (R_int(N_R+1) - R_int(N_R)) .GT. 2.d0*(R_int(2) - R_int(1)) ) then
        write(*,*) "          - Logarithmic radial grid detected."
        is_grid_log = .TRUE.
    end if
    
    write(*,*) '     ...done.'
    write(*,*)

    call write_information()
    
    !!$ Defaults number of threads for parallel computing. 1 = no parallelization
    !$ setthreads=1
    !$ if(num_args .EQ. 3 ) then
    !$     call get_command_argument(2, dumstr)
    !$     if(trim(dumstr)=='setthreads') then
    !$         call get_command_argument(3, dumstr)
    !$         read(dumstr,'(I6)') setthreads
    !$     end if
    !$ end if
    !$ setthreads = max(setthreads, 1)
    !$ write(*,'(1X, A)') achar(27)//'[95mYou are using the parallel version of'
    !$ if(setthreads==1) then
    !$     write(*,'(1X, A, I0, A)') 'the code with ', setthreads, ' thread.'
    !$ else
    !$     write(*,'(1X, A, I0, A)') 'the code with ', setthreads, ' threads.'
    !$ end if
    !$ write(*,*) achar(27)//'[0m'
    !$ call OMP_set_num_threads(setthreads)
    !$ chunksize = max( 1, N_dust / setthreads / 5 )

! ##### SET INITIAL CONDITIONS #####################################################################################################    
    write(*,*) '# Initializing...'

    ! Set initial sizes and masses of dust particles
    if(a_dist_log==1) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dummy) SCHEDULE(GUIDED, chunksize)
        do i=1, N_dust
            dummy = ran2(iseed)
            a_dust(i) = a_max**(dummy) / a_min**(dummy-1.d0)
        end do
!$OMP END PARALLEL DO
    else
        a_dust(:) = a_max
    end if
    m_dust = 4.d0/3.d0 * pi * rho_b * a_dust**3.
    
    ! Read first gas density frame
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasdens', i_start, 'dat')), sigma_gas(1, :, :))
    ! Read first gas velocity frames
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasvrad', i_start, 'dat')), vR_gas(1, :, :))
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasvtheta', i_start, 'dat')), vTheta_gas(1, :, :))
    ! Read gas temperature
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasTemperature', i_start, 'dat')), T_gas(1, :, :))
    ! Read gravitational accelerations if needed
    if(use_sg==1) then
        call read_frame( trim(fargo_datadir)//trim(make_filename('aR', i_start, 'dat')), aR_grav(1, :, :))
        call read_frame( trim(fargo_datadir)//trim(make_filename('aTheta', i_start, 'dat')), aTheta_grav(1, :, :))
    else
          aR_grav(1, :, :) = 0.d0
      aTheta_grav(1, :, :) = 0.d0
    end if
    
    ! Set initial positions of dust particles
    if(dens_dist==1) then
        write(*,*) achar(27)//'[36m        I''m setting the initial particle distribution with'
        write(*,*) '        the acceptance/rejection method. This may take'
        write(*,*) '        several minutes.'//achar(27)//'[0m'
    end if
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dummy) SCHEDULE(GUIDED, chunksize)
    do i=1, N_dust
        if(dens_dist==2) then
            R_dust(i)     = R_planet(i_start) - R_ring/2.d0 + ran2(iseed) * R_ring * 2.d0
            theta_dust(i) = 2.d0 * pi * ran2(iseed)
        else if(dens_dist==1) then
            do
                R_dust(i)     = ran2(iseed) * ( R(N_R)-R(1) ) + R(1)
                theta_dust(i) = 2.d0 * pi * ran2(iseed)
                dummy         = interp2d(theta_dust(i), R_dust(i), &
                                    & theta, R, sigma_gas(1, :, :)/maxval(sigma_gas(1, :, :)), N_theta, N_R)
                if(ran2(iseed) .LE. dummy) exit
            end do
        else
            R_dust(i)     = (R_int(N_R+1)-R_int(1))*ran2(iseed) + R_int(1)
            theta_dust(i) = 2.d0 * pi * ran2(iseed)
        end if
    end do
!$OMP END PARALLEL DO
    
    ! In carteesian coordinates
    X_dust = R_dust * cos( theta_dust )
    Y_dust = R_dust * sin( theta_dust )

    ! Initialize crystallinity fraction
    fc_dust = 0.d0

    ! Initiate the temperature, calculate Stokes number and stopping time of particles
    ! and with them the initial particle velocities
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dummy, dummy2) SCHEDULE(GUIDED, chunksize)
    do i=1, N_dust
        T_dust(i)             = interp2d(theta_dust(i), R_dust(i), theta, R,      T_gas(1, :, :), N_theta, N_R)
        dummy                 = interp2d(theta_dust(i), R_dust(i), theta, R,     vR_gas(1, :, :), N_theta, N_R)
        dummy2                = interp2d(theta_dust(i), R_dust(i), theta, R, vtheta_gas(1, :, :), N_theta, N_R)
        dummy2                = sqrt( (vR_dust(i)-dummy)**2.d0 + (vTheta_dust(i)-dummy2)**2.d0 ) ! relative gas-dust velocity
        loc_sigma_gas_dust(i) = interp2d(theta_dust(i), R_dust(i), theta, R,  sigma_gas(1, :, :), N_theta, N_R)
        St(i)                 = stokes_number( a_dust(i), loc_sigma_gas_dust(i), R_dust(i), dummy2, T_dust(i) )
        tstop(i)              = St(i) * sqrt( R_dust(i)**3.d0 )
        
        if(St(i) .LT. 1.d0) then ! Fully coupled to gas
            vR_dust(i)     = interp2d(theta_dust(i), R_dust(i), theta, R, vR_gas(1, :, :),     N_theta, N_R)
            vTheta_dust(i) = interp2d(theta_dust(i), R_dust(i), theta, R, vTheta_gas(1, :, :), N_theta, N_R)
        else ! Keplerian orbit
            vR_dust(i)     = 0.d0
            vTheta_dust(i) = 1./sqrt( R_dust(i) )
        end if
    end do
!$OMP END PARALLEL DO

    ! Specific angular momentum of the dust
    L_dust = R_dust * vTheta_dust
    ! X- and Y-component of dust velocity
    vX_dust = vR_dust * cos( theta_dust ) - vTheta_dust * sin( theta_dust ) 
    vY_dust = vR_dust * sin( theta_dust ) + vTheta_dust * cos( theta_dust )
    
    write(*,*) '     ...done.'
    write(*,*)
    
    write(*,*) '# Writing initial conditions...'
    ! Write initial conditions
    call write_dust_output(i_start)
    write(*,*) '     ...done.'
    write(*,*)

! ##### CALCULATIONS ###############################################################################################################

    write(*,*) '# Calculating particle orbits...'
    
    ! Set current time and frame
    iframe   = i_start
    cur_time = time(i_start)
    
    ! Load second data points for interpolation
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasdens',        iframe+1, 'dat')),  sigma_gas(2, :, :))
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasvrad',        iframe+1, 'dat')),     vR_gas(2, :, :))
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasvtheta',      iframe+1, 'dat')), vTheta_gas(2, :, :))
    call read_frame( trim(fargo_datadir)//trim(make_filename('gasTemperature', iframe+1, 'dat')),      T_gas(2, :, :))
    if(use_sg==1) then
        call read_frame( trim(fargo_datadir)//trim(make_filename('gassgaccr', i_start, 'dat')),         aR_grav(2, :, :))
        call read_frame( trim(fargo_datadir)//trim(make_filename('gassgacctheta', i_start, 'dat')), aTheta_grav(2, :, :))
    else
          aR_grav(2, :, :) = 0.d0
      aTheta_grav(2, :, :) = 0.d0
    end if
    
    TIME_LOOP: do
    
        ! The current gas parameters
        cur_sigma_gas   = (   sigma_gas(2, :, :) -   sigma_gas(1, :, :) ) / ( time(iframe+1)-time(iframe) ) &
                 & * ( cur_time - time(iframe) ) +   sigma_gas(1, :, :)
        cur_vR_gas      = (      vR_gas(2, :, :) -      vR_gas(1, :, :) ) / ( time(iframe+1)-time(iframe) ) &
                 & * ( cur_time - time(iframe) ) +      vR_gas(1, :, :)
        cur_vTheta_gas  = (  vTheta_gas(2, :, :) -  vTheta_gas(1, :, :) ) / ( time(iframe+1)-time(iframe) ) &
                 & * ( cur_time - time(iframe) ) +  vTheta_gas(1, :, :)
        cur_T_gas       = (       T_gas(2, :, :) -       T_gas(1, :, :) ) / ( time(iframe+1)-time(iframe) ) &
                 & * ( cur_time - time(iframe) ) +       T_gas(1, :, :)
        cur_aR_grav     = (     aR_grav(2, :, :) -      aR_grav(1, :, :) ) / ( time(iframe+1)-time(iframe) ) &
                 & * ( cur_time - time(iframe) ) +     aR_grav(1, :, :)
        cur_aTheta_grav = ( aTheta_grav(2, :, :) - aTheta_grav(1, :, :) ) / ( time(iframe+1)-time(iframe) ) &
                 & * ( cur_time - time(iframe) ) + aTheta_grav(1, :, :)
        
        ! The current planet position in cartesion coordinates
        cur_X_planet = ( X_planet(iframe+1)-X_planet(iframe) ) / ( time(iframe+1)-time(iframe) ) &
                & * ( cur_time-time(iframe) ) + X_planet(iframe)
        cur_Y_planet = ( Y_planet(iframe+1)-Y_planet(iframe) ) / ( time(iframe+1)-time(iframe) ) &
                & * ( cur_time-time(iframe) ) + Y_planet(iframe)
        ! The current position of the planet in polar coordinates
        cur_R_planet     = sqrt( cur_X_planet**2.d0 + cur_Y_planet**2.d0 )
        cur_theta_planet = atan2( cur_Y_planet, cur_X_planet )
        ! Catch for periodic theta
        if(cur_theta_planet .LT. 0.d0)      cur_theta_planet = cur_theta_planet + 2.d0 * pi
        if(cur_theta_planet .GT. 2.d0 * pi) cur_theta_planet = cur_theta_planet - 2.d0 * pi
        
        ! The current planet mass
        cur_m_planet = ( m_planet(iframe+1)-m_planet(iframe) ) / ( time(iframe+1)-time(iframe) ) &
                & * ( cur_time-time(iframe) ) + m_planet(iframe)
        
        ! The center of mass and inertial forces
        Xcms     = cur_X_planet * cur_m_planet / ( 1.d0 + cur_m_planet )
        Ycms     = cur_Y_planet * cur_m_planet / ( 1.d0 + cur_m_planet )
        Rcms     = sqrt( Xcms**2.d0 + Ycms**2.d0 )
        thetaCms = atan2( Ycms, Xcms )
        ! Catch for periodic theta
        if(thetaCms .LT. 0.d0)      thetaCms = thetaCms + 2.d0 * pi
        if(thetaCms .GT. 2.d0 * pi) thetaCms = thetaCms - 2.d0 * pi
        ! Inertial forces
        Fxin = - 1.d0 / cur_R_planet**3.d0 * Rcms * cos(thetaCms)
        Fyin = - 1.d0 / cur_R_planet**3.d0 * Rcms * sin(thetaCms)
        
        ! Set the timestep here
        ! Maximum allowed timestep until next snapshot
        dt = 1.d100
        ! Initialize growth rate
        dadt = 0.d0
        ! Set the maximum timestep, such that the particles cannot travel more than one cell.
        ! If Coagulation/fragmentation is included, limit the time step, such that the particle can change at maximum by it's
        ! own size
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, dr_dum, dt_dum) SCHEDULE(GUIDED, chunksize)
        do i=1, N_dust
            if(R_dust(i) .LE. 0.d0) goto 667
            ! Radial grid cell width
            if(is_grid_log) then
                dr_dum = R_dust(i) * ( log10(R_int(N_R+1)) - log10(R_int(1)) ) / N_R
            else
                dr_dum = ( R_int(N_R+1) - R_int(1) ) / N_R
            end if
            ! Azimuthal grid cell width
            dt_dum = R_dust(i) * 2.d0 * pi / N_theta
            ! Limit growth/fragmentation timestep
            if(do_growth==1 .OR. do_frag==1) then
                ! We need the stopping time
                loc_sigma_gas_dust(i) = interp2d(theta_dust(i), R_dust(i), theta, R, cur_sigma_gas,   N_theta, N_R)
                loc_vR_gas            = interp2d(theta_dust(i), R_dust(i), theta, R, cur_vR_gas,      N_theta, N_R)
                loc_vTheta_gas        = interp2d(theta_dust(i), R_dust(i), theta, R, cur_vTheta_gas,  N_theta, N_R)
                T_dust(i)             = interp2d(theta_dust(i), R_dust(i), theta, R, cur_T_gas,       N_theta, N_R)
                dummy                 = sqrt( (vR_dust(i)-loc_vR_gas)**2.d0 + (vTheta_dust(i)-loc_vTheta_gas)**2.d0 )
                St(i)                 = stokes_number( a_dust(i), loc_sigma_gas_dust(i), R_dust(i), dummy, T_dust(i) )
                dadt(i)               = coagfrag_rate( R_dust(i), loc_sigma_gas_dust(i), St(i), T_dust(i) )
                ! If particle fragments and already has minimum size, ignore to speed up programm
                if(dadt(i) .LT. 0.d0 .AND. a_dust(i)==a_mono) dadt(i) = 0.d0
            end if
  ! Make sure that only one thread is writing on dt  
  !$OMP ATOMIC
            dt = min(dt, abs(dr_dum/vR_dust(i)), abs(dt_dum/vTheta_dust(i)), abs(a_dust(i)/dadt(i)) )
667         continue
        end do
!$OMP END PARALLEL DO
        ! Add margin just to play safe...
        dt = 0.5d0 * dt
        ! Timestep shall not exceed next snapshot time
        dt = min(dt, time(iframe+1)-cur_time)

        ! Update particle values here
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(i, H, nu_visc, l_visc, ran_theta, loc_vR_gas, loc_vTheta_gas, d_pldu, Fx, Fy, Fr, Ftheta, dummy) &
!$OMP SCHEDULE(GUIDED, chunksize)
        do i=1, N_dust
            ! Continue if particle is out of grid
            if(R_dust(i) .LE. 0.d0) goto 666
            
            ! Do particle growth first
            if(do_growth==1 .OR. do_frag==1) then
                a_dust(i) = max( a_dust(i) + dt*dadt(i), a_mono )
                m_dust(i) = 4.d0/3.d0 * pi * rho_b * a_dust(i)**3.d0
            end if
         
            ! First part of leap-frogging to update R and theta
            theta_dust(i) = theta_dust(i) + L_dust(i) / R_dust(i)**2.d0 * dt/2.d0 ! Azimuthal coordinate
                R_dust(i) =     R_dust(i) + vR_dust(i)                  * dt/2.d0 ! Radial coordinate
            ! Catch for periodicity in theta
            if(theta_dust(i) .LT. 0.d0)      theta_dust(i) = theta_dust(i) + 2.d0 * pi
            if(theta_dust(i) .GT. 2.d0 * pi) theta_dust(i) = theta_dust(i) - 2.d0 * pi
            ! Convert new particle position into cartesian coordinates
            X_dust(i) = R_dust(i) * cos( theta_dust(i) )
            Y_dust(i) = R_dust(i) * sin( theta_dust(i) )
            
            ! New dust temperature
            T_dust(i) = interp2d(theta_dust(i), R_dust(i), theta, R, cur_T_gas, N_theta, N_R)
            
            ! Scale height at position of the particle
            H = sqrt( adx*T_dust(i)*R_dust(i)**3.d0 )
            
            ! Add turbulence as random walk
            if(do_randomwalk == 1) then
                ! Viscosity at the position of the particle
                ! nu = alpha * cs * H = alpha * H**2 * omega
                nu_visc = alpha * H**2.d0 / sqrt( R_dust(i)**3.d0 )
                ! Displacement length of the particle
                l_visc = sqrt( dt*nu_visc / (1.d0+St(i)**2.d0) )
                ! Displace in random direction
                ran_theta = ran2(iseed) * 2.d0 * pi
                ! Do the displacement
                X_dust(i) = X_dust(i) + l_visc * cos(ran_theta)
                Y_dust(i) = Y_dust(i) + l_visc * sin(ran_theta)
                ! Calculate the new polar coordinates
                R_dust(i)     = sqrt( X_dust(i)**2.d0 + Y_dust(i)**2.d0 )
                theta_dust(i) = atan2(Y_dust(i), X_dust(i))
                ! Catch for periodicity in theta
                if(theta_dust(i) .LT. 0.d0)      theta_dust(i) = theta_dust(i) + 2.d0 * pi
                if(theta_dust(i) .GT. 2.d0 * pi) theta_dust(i) = theta_dust(i) - 2.d0 * pi
                ! New dust temperature
                T_dust(i) = interp2d(theta_dust(i), R_dust(i), theta, R, cur_T_gas, N_theta, N_R)
                ! New scale height at the position of the particle
                H = sqrt( adx*T_dust(i)*R_dust(i)**3.d0 )
            end if
            
            ! If the particle left the boundaries we change the radial coordinate to
            ! -1.d0: Inner boundary
            !  0.d0: Outer boundary
            if(R_dust(i) .LT. R_int(1)) then
                R_dust(i) = -1.d0
                X_dust(i) = 0.d0
                Y_dust(i) = 0.d0
                goto 666
            end if
            if(R_dust(i) .GT. R_int(N_R+1)) then
                R_dust(i) =  0.d0
                X_dust(i) =  0.d0
                Y_dust(i) =  0.d0
                goto 666
            end if
            
            ! Local gas parameters at particle position
            loc_sigma_gas_dust(i) = interp2d(theta_dust(i), R_dust(i), theta, R, cur_sigma_gas,   N_theta, N_R)
            loc_vR_gas            = interp2d(theta_dust(i), R_dust(i), theta, R, cur_vR_gas,      N_theta, N_R)
            loc_vTheta_gas        = interp2d(theta_dust(i), R_dust(i), theta, R, cur_vTheta_gas,  N_theta, N_R)
            loc_aR_grav           = interp2d(theta_dust(i), R_dust(i), theta, R, cur_aR_grav,     N_theta, N_R)
            loc_aTheta_grav       = interp2d(theta_dust(i), R_dust(i), theta, R, cur_aTheta_grav, N_theta, N_R)
            
            ! Calculate particles Stokes number and stopping time
            dummy    = sqrt( (vR_dust(i)-loc_vR_gas)**2.d0 + (vTheta_dust(i)-loc_vTheta_gas)**2.d0 ) ! relative gas-dust velocity
            St(i)    = stokes_number( a_dust(i), loc_sigma_gas_dust(i), R_dust(i), dummy, T_dust(i) )
            tstop(i) = St(i) * sqrt( R_dust(i)**3.d0 )
            
            ! Calculate the planetary forces and update the velocities
            ! Planet-dust distance with smooting parameter
            d_pldu = sqrt( (X_dust(i)-cur_X_planet)**2.d0 + (Y_dust(i)-cur_Y_planet)**2.d0 + (smoothing * H)**2.d0 )
            ! The X- and Y-components of the Force
            Fx = - cos( theta_dust(i) ) / R_dust(i)**2.d0 &                     ! from star
                 - (X_dust(i)-cur_X_planet)/d_pldu * cur_m_planet/d_pldu**2.d0  ! from planet
            Fy = - sin( theta_dust(i) ) / R_dust(i)**2.d0 &                     ! from star
                 - (Y_dust(i)-cur_Y_planet)/d_pldu * cur_m_planet/d_pldu**2.d0  ! from planet
         
            ! Adding inertial forces
            Fx = Fx + Fxin
            Fy = Fy + Fyin
            ! Convert forces to polar coordiantes
            Fr     =  Fx * cos( theta_dust(i) ) + Fy * sin( theta_dust(i) )
            Ftheta = -Fx * sin( theta_dust(i) ) + Fy * cos( theta_dust(i) )
            ! Adding self gravity acceleration of the disk
            Fr     = Fr     + loc_aR_grav
            Ftheta = Ftheta + loc_aTheta_grav * R_dust(i)
               
            ! Updating the angular momentum
            L_dust(i) = L_dust(i) + dt * ( Ftheta + loc_vTheta_gas*R_dust(i)/tstop(i) )
            L_dust(i) = L_dust(i) / ( 1.d0 + dt/tstop(i) )
            ! Updating vR_dust
            vR_dust(i) = vR_dust(i) + dt * ( L_dust(i)**2.d0/R_dust(i)**3.d0 + Fr + loc_vR_gas/tstop(i) )
            vR_dust(i) = vR_dust(i) / ( 1.d0 + dt / tstop(i) )
            ! Updating the vTheta_dust
            vTheta_dust(i) = L_dust(i) / R_dust(i)

            ! Converting to cartesian coordinates
            vX_dust(i) = vR_dust(i) * cos( theta_dust(i) ) - vTheta_dust(i) * sin( theta_dust(i) )
            vY_dust(i) = vR_dust(i) * sin( theta_dust(i) ) + vTheta_dust(i) * cos( theta_dust(i) )
            
            ! Second part of leap-frogging to update R and theta
            theta_dust(i) = theta_dust(i) + L_dust(i) / R_dust(i)**2.d0 * dt/2.d0 ! Azimuthal coordinate
                R_dust(i) =     R_dust(i) + vR_dust(i)                  * dt/2.d0 ! Radial coordinate
            ! Catch for periodicity in theta
            if(theta_dust(i) .LT. 0.d0)      theta_dust(i) = theta_dust(i) + 2.d0 * pi
            if(theta_dust(i) .GT. 2.d0 * pi) theta_dust(i) = theta_dust(i) - 2.d0 * pi
            ! Convert new particle position into cartesian coordinates
            X_dust(i) = R_dust(i) * cos( theta_dust(i) )
            Y_dust(i) = R_dust(i) * sin( theta_dust(i) )
            
            ! If the particle left the boundaries we change the radial coordinate to
            ! -1.d0: Inner boundary
            !  0.d0: Outer boundary
            if(R_dust(i) .LT. R_int(1)) then
                R_dust(i) = -1.d0
                X_dust(i) = 0.d0
                Y_dust(i) = 0.d0
                goto 666
            end if
            if(R_dust(i) .GT. R_int(N_R+1)) then
                R_dust(i) =  0.d0
                X_dust(i) =  0.d0
                Y_dust(i) =  0.d0
                goto 666
            end if
            
            ! After all calculate the new crystallinity fraction
            dummy = T_dust(i) * mu / R_gas * G * phys_mass / phys_dist
            fc_dust(i) = min(1.d0, (fc_dust(i)**(1.d0/3.d0) + dt * 2.d0*zeta**(1.d0/3.d0)*nu_vib*exp(-Ea/dummy))**3.d0)            
            
666       continue

          ! Check for validity
          if( number_invalid(R_dust(i)) .OR. R_dust(i) .LT. -1.d0 ) then
  !$OMP CRITICAL
              call write_invalid_message('R_dust is invalid', R_dust(i))
              stop
  !$OMP END CRITICAL
          end if
          if( number_invalid(theta_dust(i)) .OR. theta_dust(i) .LT. 0.d0 .OR. theta_dust(i) .GT. 2.d0*pi ) then
  !$OMP CRITICAL
              call write_invalid_message('theta_dust is invalid', theta_dust(i))
              stop
  !$OMP END CRITICAL
          end if
          if( number_invalid(vR_dust(i)) ) then
  !$OMP CRITICAL
              call write_invalid_message('vR_dust is invalid', vR_dust(i))
              stop
  !$OMP END CRITICAL
          end if
          if( number_invalid(vTheta_dust(i)) ) then
  !$OMP CRITICAL
              call write_invalid_message('vTheta_dust is invalid', vTheta_dust(i))
              stop
  !$OMP END CRITICAL
          end if
          if( number_invalid(m_dust(i)) .OR. m_dust(i) .LE. 0.d0 ) then
  !$OMP CRITICAL
              call write_invalid_message('m_dust is invalid', m_dust(i))
              stop
  !$OMP END CRITICAL
          end if
          if( number_invalid(St(i)) .OR. St(i) .LT. 0.d0 ) then
  !$OMP CRITICAL
              call write_invalid_message('St is invalid', St(i))
              stop
  !$OMP END CRITICAL
          end if
          
          if( number_invalid(T_dust(i)) .OR. St(i) .LT. 0.d0 ) then
  !$OMP CRITICAL
              call write_invalid_message('T_dust is invalid', T_dust(i))
              stop
  !$OMP END CRITICAL
          end if

        end do
!$OMP END PARALLEL DO

        ! Update time
        cur_time = cur_time + dt
        
        ! When we arrive at a new snapshot, we need to write the output file and load a new frame
        if(cur_time == time(iframe+1)) then
            ! Write output
            call write_dust_output(iframe+1)
            ! Exit loop when we're at the end of the simulation
            if(iframe+1 == i_stop) exit TIME_LOOP
            ! If not, load new frame
              sigma_gas(1, :, :) =   sigma_gas(2, :, :)
                 vR_gas(1, :, :) =      vR_gas(2, :, :)
             vTheta_gas(1, :, :) =  vTheta_gas(2, :, :)
                  T_gas(1, :, :) =       T_gas(2, :, :)
                aR_grav(1, :, :) =     aR_grav(2, :, :)
            aTheta_grav(1, :, :) = aTheta_grav(2, :, :)
            call read_frame( trim(fargo_datadir)//trim(make_filename('gasdens',        iframe+1, 'dat')),  sigma_gas(2, :, :))
            call read_frame( trim(fargo_datadir)//trim(make_filename('gasvrad',        iframe+1, 'dat')),     vR_gas(2, :, :))
            call read_frame( trim(fargo_datadir)//trim(make_filename('gasvtheta',      iframe+1, 'dat')), vTheta_gas(2, :, :))
            call read_frame( trim(fargo_datadir)//trim(make_filename('gasTemperature', iframe+1, 'dat')),      T_gas(2, :, :))
            if(use_sg==1) then
                call read_frame( trim(fargo_datadir)//trim(make_filename('gassgaccr',     iframe+1, 'dat')),     aR_grav(2, :, :))
                call read_frame( trim(fargo_datadir)//trim(make_filename('gassgacctheta', iframe+1, 'dat')), aTheta_grav(2, :, :))
            else
                  aR_grav(2, :, :) = 0.d0
              aTheta_grav(2, :, :) = 0.d0
            end if
            iframe = iframe + 1
        end if
        
        ! If all particles are out of the simulation we stop
        if( maxval(R_dust) .LE. 0.d0 ) then
            write(*,*)
            call write_warning_message('All particles are out of the simulation domain')
            exit TIME_LOOP
        end if
    
    end do TIME_LOOP
    
    write(*,*) '     ...done.'
    write(*,*)
    
    call write_losses()

! ##### DEALLOCATE MEMORY ##########################################################################################################
    
    write(*,*) '# Deallocating memory...'
    ! Deallocate memory
    call deallocate_memory()
    write(*,*) '     ...done.'
    write(*,*)
    
    ! Stop timer
    call system_clock( time_finish, time_rate )
    
    call write_goodbye()
    
    write(*,'(1X, A, 1X, f9.1, 1X, A)') 'Elapsed time: ', 1.d0*(time_finish-time_start)/time_rate/3600.d0, 'hrs.'
    write(*,*)

end program particle_code
