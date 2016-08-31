module aerodynamics

    implicit none
    
    public :: stokes_number
    
    contains
    
    double precision function stokes_number(a, sigma, R, vRel, T)
    ! Calculates the Stokes number according Paardekooper (2007)
        use constants
        use variables, only: phys_dist, phys_mass, rho_b, mu, adx
    
        implicit none
        
        double precision, intent(in) :: a       ! Particle radius
        double precision, intent(in) :: sigma   ! Surface density
        double precision, intent(in) :: R       ! Radial position
        double precision, intent(in) :: vRel    ! Relative velocity gas-dust
        double precision, intent(in) :: T       ! Local gas temperature
        
        double precision :: cs            ! Dound speed
        double precision :: fd            ! Epstein drag coefficient
        double precision :: kd            ! Stokes drag coefficien
        double precision :: Kn            ! Knudsen number
        double precision :: rho_g         ! Midplane gas mass density
        double precision :: vK            ! Keplerian velocity
        
        double precision :: M             ! Mach number
        double precision :: mfp           ! Mean free path
        double precision :: H             ! Pressure scale height
        double precision :: omegaK        ! Keplerian frequency
        double precision :: Re            ! Reynolds number

        H      = sqrt( adx*T*R**3.d0 )                !
        omegaK = 1.d0 / sqrt(R**3.d0)                 ! This is all
        vK     = omegaK * R                           ! in FARGO
                                                      ! units!
        cs     = H * omegaK                           !
        M      = vRel / cs                            !
        
        rho_g  = sigma / (sqrt(2.d0*pi)*H)            ! Midplane gas mass density in FARGO units
        rho_g  = rho_g * phys_mass / phys_dist**3.d0  ! Midplane gas mass density in physical units
        mfp    = mu*m_p / (rho_g*sig_H2)              ! Mean free path in physical units
        
        Kn     = 0.5d0 * mfp/a                        !
        Re     = 3.d0 * sqrt( pi/8.d0 ) * M/Kn        ! These are all
        fd     = sqrt( 1.d0 + 9.d0*pi*M**2d0/128.d0 ) ! dimensionless quantities
        if(Re .LE. 500.d0) then                       !
            kd = 1.d0 + 0.15d0*Re**(0.687d0)          !
        else if(Re .LE. 1500.d0) then                 !
            kd = 3.96d-6 * Re**(2.4d0)                !
        else                                          !
            kd = 0.11d0 * Re                          !
        end if                                        !
        
        stokes_number = sqrt(pi/8.d0) * (3.d0*Kn+1.d0)**2.d0 / (9.d0*Kn**2.d0*fd+3.d0*Kn*kd) &
            & * a/(R*phys_dist) * rho_b/rho_g * vK/cs
    
    end function stokes_number

end module aerodynamics
