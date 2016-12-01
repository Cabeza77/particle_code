module io

    implicit none
    
    public  :: make_filename, read_inputs, read_data, read_frame, write_dust_output, write_stop_message, write_warning_message, &
        & write_welcome, write_information, write_losses, write_invalid_message, read_dust_frame
    private :: check_inputs, read_dims_file, read_namelist, read_used_rad
    
    contains
    
        subroutine read_inputs()
        
          use variables
        
            implicit none
            
            ! Set standards
            i_start       = 0                       ! Starting snapshot of dust simulation
            i_stop        = 1000                    ! Last snapshot of dust simulation
            N_dust        = 10000                   ! Number of dust particles,
            a_min         = 1.d-4                   ! Minimum radius of dust particles in cm
            a_max         = 1.d1                    ! Maximum radius of dust particles in cm
            a_dist_log    = 1                       ! Initial size distribution of dust particles
                                                    !   1: Logarithmic between a_min and a_max
                                                    !   0: All particles have radius a_max
            a_mono        = 1.d-4                   ! Monomer size. This is the minimum radius a particle can reach, when
                                                    ! fragmentation is included.
            rho_b         = 3.3d0                   ! Bulk density of dust particles in g/cm**3
            dens_dist     = 1                       ! Initial particle distribution
                                                    !   2: Uniformly in ring of size R_ring centered on planet
                                                    !   1: Like the gas density of the starting frame
                                                    !   0: Uniformly distributed in R-theta-space (intrinsic 1/R dependency)
            do_growth     = 0                       ! Include particle growth?        1: Yes. 0: No.
            do_frag       = 0                       ! Include fragmentation?          1: Yes. 0: No.
            do_randomwalk = 0                       ! Include random particle walk?   1: Yes. 0: No.
            v_frag        = 100.d0                  ! Fragmentation velocity in cm/s
            use_frag_ice  = 0                       ! Shall we use a different velocity for ice?  1: Yes. 0: No.
            v_frag_ice    = 1000.d0                 ! Fragmentation velocity of ice material in cm/s
            T_ice         = 150.d0                  ! If colder that T_ice, we use v_frag_ice.
            R_ring        = 20.d0                   ! Size in AU of the ring for dens_dist=2 
            phys_dist     = 1.d0                    ! Physical distance in AU of the Fargo distance unit
            phys_mass     = 1.d0                    ! Physical mass in solar masses of the Fargo mass unit
            adx           = 1.4d0                   ! Adiabatic index
            alpha         = 1.d-3                   ! Turbulent alpha parameter for random walk. (alpha=0.d0 => No random walk)
            eps           = 1.d-2                   ! Dust-to-gas ratio
            mu            = 2.3                     ! Mean molecular weight of gas in proton masses
            smoothing     = 0.6d0                   ! Smooting length in units of scale height at planet position
            use_sg        = 0                       ! Using self-gravity of the disk?
            
            Ea            = 5330.d0                 ! Activation energy for crystallization (Ea/kb) in units of K
            nu_vib        = 2.2d13                  ! Vibrational frequency in s
            zeta          = 1.d-5                   ! Initial fraction of growth centers
            
            ! Read the namelist
            call read_namelist()
            
            ! Read dims.dat
            call read_dims_file()
            
            ! Check the inputs
            call check_inputs()
        
        end subroutine read_inputs
        
! ####################

        subroutine read_data()
        
            use constants
            use variables
        
            implicit none
            
            integer :: i
            
            ! Read used_rad.dat
            call read_used_rad()
            
            ! Create theta grid
            do i=1, N_theta+1
                theta_int(i) = ( i - 1 ) * 2.d0 * pi / N_theta
            end do
            theta = 0.5d0 * ( theta_int(2:) + theta_int(:N_theta) )
            
            ! Read planet0.dat
            call read_planet0()
        
        end subroutine read_data
        
! ####################
        
        subroutine read_namelist()
        
            use constants
            use variables
        
            implicit none
            
            integer :: ierror
            
            namelist / inputs / &
                & fargo_datadir, &
                & output_dir, &
                & i_start, i_stop, &
                & N_dust, &
                & a_min, a_max, a_dist_log, a_mono, do_growth, do_frag, v_frag, do_randomwalk, &
                & use_frag_ice, v_frag_ice, T_ice, &
                & rho_b, &
                & dens_dist, R_ring, &
                & phys_dist, phys_mass, &
                & alpha, mu, eps, adx, &
                & smoothing, use_sg, &
                & Ea, nu_vib, zeta
                
            open(unit=100, file=trim(input_file), delim='apostrophe', status='old', action='read', iostat=ierror)
            if(ierror==0) then
                read(100,nml=inputs)
            else
                call write_stop_message("Could not read '"//trim(input_file)//"'")
                stop
            end if
            close(100)
            
            phys_dist = phys_dist * AU
            phys_mass = phys_mass * M_sun
            phys_time = sqrt( phys_dist**3.d0 / (G*phys_mass) )
            
            R_ring    = R_ring * AU / phys_dist
            
            ! Convert T_ice in Fargo units
            T_ice = T_ice * R_gas/(mu*G) * phys_dist/phys_mass
        
        end subroutine read_namelist
        
! ####################

        subroutine read_dims_file()
        
            use variables
        
            implicit none
            
            integer :: ierror
            double precision :: dummy
            
            open(unit=100, file=trim(fargo_datadir)//'dims.dat', delim='apostrophe', status='old', action='read', iostat=ierror)
            if(ierror==0) then
                read(100,*) dummy, dummy, dummy, dummy, dummy, N_t, N_R, N_theta
            else
                call write_stop_message("Could not read '"//trim(fargo_datadir)//'dims.dat'//"'")
                stop
            end if
            close(100)
        
        end subroutine read_dims_file
        
! ####################

        subroutine read_used_rad()
        
            use variables
        
            implicit none
            
            integer :: ierror, i
            
            open(unit=100, file=trim(fargo_datadir)//'used_rad.dat', delim='apostrophe', status='old', action='read', iostat=ierror)
            if(ierror==0) then
                do i=1, N_R+1
                    read(100,*) R_int(i)
                end do
            else
                call write_stop_message("Could not read '"//trim(fargo_datadir)//'used_rad.dat'//"'")
                stop
            end if
            close(100)
            
            R = 0.5d0 * ( R_int(2:) + R_int(:N_R) )
        
        end subroutine read_used_rad
        
! ####################

        subroutine read_planet0()
        
            use variables
        
            implicit none
            
            integer :: ierror, dummyi
            double precision :: dummy, dummyx, dummyy, dummyvx, dummyvy, dummym, dummyt
            
            open(unit=100, file=trim(fargo_datadir)//'planet0.dat', delim='apostrophe', status='old', action='read', iostat=ierror)
            if(ierror==0) then
                do
                    read(100,*,iostat=ierror) dummyi, dummyx, dummyy, dummyvx, dummyvy, dummym, dummy, dummyt
                    x_planet(dummyi)  = dummyx
                    y_planet(dummyi)  = dummyy
                    vX_planet(dummyi) = dummyvx
                    vY_planet(dummyi) = dummyvy
                    vY_planet(dummyi) = dummyvy
                    m_planet(dummyi)  = dummym
                    time(dummyi)      = dummyt
                    ! Calculate polar coordiantes
                    R_planet(dummyi)      = sqrt( dummyx**2.d0 + dummyy**2.d0 )
                    theta_planet(dummyi)  = atan2( dummyy, dummyx )
                    if(ierror<0) then ! End of file => exit loop
                        exit
                    end if
                end do
            else
                call write_stop_message("Could not read '"//trim(fargo_datadir)//'planet0.dat'//"'")
                stop
            end if
            close(100)
        
        end subroutine read_planet0

! ####################

        subroutine read_dust_frame()
        
            use variables
            use constants
            
            implicit none
        
            integer :: ierror
            integer :: i
            
            integer :: idummy
            double precision :: dummy
            
            open(100, file=trim(output_dir)//trim(make_filename('dust', i_restart, 'dat')), delim='apostrophe', status='old', &
                & action='read', iostat=ierror)
            
                ! Skip header line
                read(100, *)
            
                do i=1, N_dust
                    read(100, *) idummy, R_dust(i), theta_dust(i), vR_dust(i), vTheta_dust(i), X_dust(i), Y_dust(i), vX_dust(i), &
                        & vY_dust(i), a_dust(i), m_dust(i), St(i), tstop(i), T_dust(i), fc_dust(i), loc_sigma_gas_dust(i)
                    
                    ! Convert to code units
                    R_dust(i)             = R_dust(i) * AU / phys_dist
                    vR_dust(i)            = vR_dust(i) * phys_time / phys_dist
                    vTheta_dust(i)        = vTheta_dust(i) * phys_time / phys_dist
                    X_dust(i)             = X_dust(i) * AU / phys_dist
                    Y_dust(i)             = Y_dust(i) * AU / phys_dist
                    vX_dust(i)            = vX_dust(i) * phys_time / phys_dist
                    vY_dust(i)            = vY_dust(i) * phys_time / phys_dist
                    tstop(i)              = tstop(i) / phys_time
                    T_dust(i)             = T_dust(i) * R_gas/mu * phys_dist/phys_mass / G
                    loc_sigma_gas_dust(i) = loc_sigma_gas_dust(i) * phys_dist**2.d0 / phys_mass
                end do
            
            close(100)
        
        end subroutine read_dust_frame

! ####################

        subroutine read_frame(filename, frame)
        
            use variables
        
            implicit none
            
            character(*),     intent(in)  :: filename
            double precision, intent(out) :: frame(N_theta, N_R)
            
            integer :: it, ir
            
            open(100, file=trim(filename), form='unformatted', access='direct', recl=8)
            
            do ir=1, N_R
                do it=1, N_theta
                    read(100, rec=(ir-1)*N_theta+it) frame(it, ir)
                end do
            end do
            
            close(100)
        
        end subroutine read_frame
        
! ####################

        character(64) function make_filename(base, num, ext)
        
            implicit none
            
            character(*) :: base
            character(*) :: ext
            integer      :: num
            
            character(8) :: i2str
            
            write(i2str,'(I0)') num
            
            make_filename = trim(base)//trim(i2str)//'.'//ext
        
        end function make_filename

! ####################

        subroutine check_inputs()
        
            use variables
        
            implicit none
            
            if(i_start .GE. i_stop) then
                call write_stop_message("i_start >= i_stop")
                stop
            end if
            
            if(i_stop .GT. N_t) then
                call write_stop_message("i_stop is larger than the maximum snapshot")
                stop
            end if
            
            if(i_start .LT. 0 .OR. i_stop .LT. 0) then
                call write_stop_message("i_start and/or i_stop is less than zero")
                stop
            end if
            
            if(N_dust < 1) then
                call write_stop_message("No dust particles: N_dust < 1")
                stop
            end if
            
            if(dens_dist .LT. 0 .OR. dens_dist .GT. 2) then
                call write_stop_message("Unknown dens_dist")
                stop
            end if
            
            if(a_min >= a_max .AND. a_dist_log==1) then
                call write_stop_message("a_min >= a_max")
                stop
            end if
            
            if(a_dist_log .LT. 0 .OR. a_dist_log .GT. 1) then
                call write_stop_message("Unknown a_dist_log")
                stop
            end if
            
            if(alpha .LT. 0.d0) then
                call write_stop_message("You are using a negative turbulent alpha")
                stop
            end if
        
            if(mu .LT. 0.d0) then
                call write_stop_message("Your mean molecular weight is negative")
                stop
            end if
            
            if(phys_mass .LT. 0.d0) then
                call write_stop_message("Your stellar mass is negative")
                stop
            end if
            
            if(phys_dist .LT. 0.d0) then
                call write_stop_message("Your physical distance is negative")
                stop
            end if
            
            if(R_ring .LT. 0.d0) then
                call write_stop_message("Your R_ring is negative")
                stop
            end if
            
            if(rho_b .LT. 0.d0) then
                call write_stop_message("Your particle density is negative")
                stop
            end if
            
            if(adx .LT. 1.d0) then
                call write_stop_message("Your adiabatic index is < 1.d0")
                stop
            end if
            
            if(smoothing .LT. 0.d0) then
                call write_stop_message("Your smoothing length is negative")
                stop
            end if
            
            if( (do_growth==1 .OR. do_frag==1) .AND. alpha==0.d0 ) then
                call write_stop_message("Growth/Fragmentation needs alpha>0")
                stop
            end if
            
            if( a_mono .GT. a_max ) then
                call write_stop_message("Your a_mono is larger than a_max")
                stop
            end if
            
            if( a_mono .GT. a_min .AND. a_dist_log==1 ) then
                call write_stop_message("Your a_mono is larger than a_min")
                stop
            end if
            
            if( do_randomwalk==1 .AND. alpha==0.d0 ) then
                call write_stop_message("Random walk needs alpha>0")
                stop
            end if
        
        end subroutine check_inputs

! ####################

        subroutine write_stop_message(msg)
            implicit none
            character(*) :: msg
            write(*,*) achar(27)//"[31mERROR:"//achar(27)//"[0m "//trim(msg)//"."
            write(*,*) achar(27)//"[31mProgram stopped!"//achar(27)//"[0m"
            write(*,*)
        end subroutine write_stop_message
        
! ####################

        subroutine write_invalid_message(msg, val)
            implicit none
            character(*)     :: msg
            double precision :: val
            write(*,'(1X, A, e10.3e3, A)') achar(27)//"[31mOoop! I made a boo-boo!"//achar(27)//"[0m "//trim(msg)//": ",val,"."
            write(*,*) achar(27)//"[31mProgram stopped!"//achar(27)//"[0m"
            write(*,*)
        end subroutine write_invalid_message
        
! ####################

        subroutine write_warning_message(msg)
            implicit none
            character(*) :: msg
            write(*,*) achar(27)//"[32mWARNING:"//achar(27)//"[0m "//trim(msg)//"."
            write(*,*)
        end subroutine write_warning_message
        
! ####################

        subroutine write_dust_output(i)
            use constants
            use variables
            implicit none
            integer, intent(in) :: i
            integer :: j
            integer :: ierror
            write(*,*) achar(27)//"[34m     ---> Writing file: '" &
                & //trim(output_dir)//trim(make_filename('dust', i, 'dat'))//"'."//achar(27)//"[0m"
            open(unit=100, file=trim(output_dir)//trim(make_filename('dust', i, 'dat')), status='replace', &
                & action='write', iostat=ierror)
            if(ierror.NE.0) then
                call write_stop_message("Could not write '"//trim(make_filename('dust', i, 'dat'))//"'")
                stop
            end if
            write(100,'(A1, 1X, A10, 15(1X, A19))') &
                & '#', 'Identifier', &
                & 'R', 'theta', &
                & 'v_R', 'v_theta', &
                & 'X', 'Y', &
                & 'v_X', 'v_Y', &
                & 'Particle radius', 'Mass', &
                & 'Stokes number', 'Stopping time', &
                & 'Temperature', 'Crystallinity', &
                & 'Gas density'
            do j=1, N_dust
                write(100, '(2X, I10, 4f20.10, 4f20.10, 2e20.10E3, 2e20.10E3, 2e20.10E3, e20.10E3)') &
                    & j, &
                    & R_dust(j)*phys_dist/AU, theta_dust(j), &
                    & vR_dust(j)*phys_dist/phys_time, vTheta_dust(j)*phys_dist/phys_time, &
                    & x_dust(j)*phys_dist/AU, y_dust(j)*phys_dist/AU, &
                    & vX_dust(j)*phys_dist/phys_time, vY_dust(j)*phys_dist/phys_time, &
                    & a_dust(j), m_dust(j), &
                    & St(j), tstop(j)*phys_time, &
                    & T_dust(j)*mu/R_gas*G*phys_mass/phys_dist, fc_dust(j), &
                    & loc_sigma_gas_dust(j)*phys_mass/phys_dist**2.d0
            end do
            close(100)
        end subroutine write_dust_output
        
! ####################

        subroutine write_welcome()
            implicit none
            write(*,*)
            write(*,*) ''//achar(27)//'[90m#######################################'
            write(*,*) '#                                     #'
            write(*,*) '#    '//achar(27)//'[31m####    #   #   #####   #####    '//achar(27)//'[90m#'
            write(*,*) '#    '//achar(27)//'[31m#   #   #   #   #         #      '//achar(27)//'[90m#'
            write(*,*) '#    '//achar(27)//'[31m#   #   #   #   #####     #      '//achar(27)//'[90m#'
            write(*,*) '#    '//achar(27)//'[31m#   #   #   #       #     #      '//achar(27)//'[90m#'
            write(*,*) '#    '//achar(27)//'[31m####    #####   #####     #      '//achar(27)//'[90m#'
            write(*,*) '#           for '//achar(27)//'[33mFARGO-ADSG            '//achar(27)//'[90m#'
            write(*,*) '#                v0.1                 #'
            write(*,*) '#                2016                 #'
            write(*,*) '#                                     #'
            write(*,*) '#######################################'
            write(*,*)
            write(*,*) '               '//achar(27)//'[0mwritten by'
            write(*,*) '   Sebastian Stammler, Sareh Ataiee &'
            write(*,*) '               Andras Zsom'
            write(*,*)
            write(*,*) '            '//achar(27)//'[36mreport bugs to:'
            write(*,*) '      stammler@uni-heidelberg.de'//achar(27)//'[0m'
            write(*,*)
        end subroutine write_welcome


        subroutine write_goodbye()
            implicit none
            write(*,*) achar(27)//'[31m#######################################'
            write(*,*) '     !!!   PROGRAM FINISHED   !!!      '
            write(*,*) '#######################################'//achar(27)//'[0m'
            write(*,*)
        end subroutine write_goodbye

! ####################
        
        subroutine write_information()
            use variables
            use constants
            implicit none
            write(*,*) achar(27)//'[32mWELCOME TO THE SHOW!'
            write(*,*) '     I detected the following inputs:'
            write(*,'(1X, A, I0)') '        # radial grid points:     ', N_R
            write(*,'(1X, A, I0)') '        # azimuthal grid points:  ', N_theta
            write(*,'(1X, A, I0)') '        # first snaphot:          ', i_start
            write(*,'(1X, A, I0)') '        # last snaphot:           ', i_stop
            if(do_restart==1) then
                write(*,*)
                write(*,*) '        ¡ This is a restarted simulation !'
                write(*,'(1X, A, I0)') '        # restart snaphot:        ', i_restart
            end if
            write(*,*)
            write(*,'(1X, A, I0)') '        # dust particles:         ', N_dust
            if(a_dist_log==1) then
                write(*,*) '        # initial particle sizes logarithmically spaced between a_min and a_max'
                write(*,'(1X, A, e9.3e2, A)') '             a_min = ', a_min, ' cm'
                write(*,'(1X, A, e9.3e2, A)') '             a_max = ', a_max, ' cm'
            else
                write(*,'(1X, A, e9.3e2, A)') '        # initial particle size:  ', a_max, ' cm'
            end if
            write(*,'(1X, A, e9.3e2, A)') '        # particle bulk density:  ', rho_b, ' g/cm3'
            write(*,*)
            if(do_growth==1) then
                write(*,*) '        # You want to do particle growth'
            end if
            if(do_frag==1) then
                write(*,*) '        # You want to do fragmentation'
                write(*,'(1X, A, f6.1, A)') '             v_frag     = ', v_frag, '    cm/s'
                if(use_frag_ice==1) then
                    write(*,'(1X, A, f6.1, A)') '             v_frag_ice = ', v_frag_ice, '    cm/s'
                    write(*,'(1X, A, f6.1, A)') '             T_ice      = ', T_ice*mu/R_gas*G*phys_mass/phys_dist, '    K'
                end if
                write(*,'(1X, A, e9.3e2, A)') '             a_mono     = ', a_mono, ' cm'
            end if
            if(do_frag==1 .OR. do_growth==1) then
                write(*,*)
            end if
            if(do_randomwalk==1) then
                write(*,*) '        # Your particles are doing a random walk'
                write(*,*)
            end if
            write(*,'(1X, A, f3.1)') '        # adiabatic index:        ', adx
            write(*,'(1X, A, f4.1, A)') '        # smoothing length:       ', smoothing*100.d0, ' % of local pressure scale height'
            write(*,'(1X, A, e9.3e2)') '        # alpha viscosity:        ', alpha
            write(*,*)
            write(*,'(1X, A, e9.3e2, A)') '        # activation energy:      ', Ea, ' K'
            write(*,'(1X, A, e9.3e2, A)') '        # vibrational frequency:  ', nu_vib, ' Hz'
            write(*,'(1X, A, e9.3e2)') '        # growth centers:         ', zeta
            write(*,*)
            write(*,'(1X, A, f4.1, A)') '        # one FARGO distance unit corresponds to ', phys_dist/AU, ' AU'
            write(*,'(1X, A, f4.1, A)') '        # the central star has a mass of ', phys_mass/m_sun, ' M_sun'
            write(*,*) '     If anything of this is wrong, you are screwed...'
            write(*,*) achar(27)//'[0m'
        end subroutine write_information
        
! ####################
        
        subroutine write_losses()
            use variables
            implicit none
            integer :: i, N_lost_in=0, N_lost_out=0
            do i=1, N_dust
                if(R_dust(i)==0.d0)   N_lost_out = N_lost_out+1
                if(R_dust(i).LT.0.d0) N_lost_in = N_lost_in+1
            end do
            write(*,*) achar(27)//'[32mHouston, we lost'
            if(N_lost_in == 1) then
                write(*,'(6X, I0, A)') N_lost_in, ' particle  through the inner and'
            else
                write(*,'(6X, I0, A)') N_lost_in, ' particles through the inner and'
            end if
            if(N_lost_out == 1) then
                write(*,'(6X, I0, A)') N_lost_out, ' particle  through the outer boundary!'
            else
                write(*,'(6X, I0, A)') N_lost_out, ' particles through the outer boundary!'
            end if
            write(*,*) achar(27)//'[0m'
        end subroutine write_losses

end module io
