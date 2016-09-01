module coagulation

    implicit none
    
    public  :: coagfrag_rate
    private :: coag_rate, frag_rate, v_rel

    contains
    
    double precision function coagfrag_rate(R, sigma, St, T)
    ! Change in radius in physical units due to growth
        use constants
        use variables, only: do_growth, do_frag
    
        implicit none
        
        double precision, intent(in) :: R
        double precision, intent(in) :: sigma
        double precision, intent(in) :: St
        double precision, intent(in) :: T
        
        coagfrag_rate = 0.d0
        
        if(do_growth==1) then
            coagfrag_rate = coagfrag_rate + coag_rate(R, sigma, St, T)
        end if
        if(do_frag==1) then
            coagfrag_rate = coagfrag_rate + frag_rate(R, sigma, St, T)
        end if
        
    end function coagfrag_rate

!####################

    double precision function coag_rate(R, sigma, St, T)
    ! Change in radius in physical units due to growth
        use constants
        use variables, only: adx, alpha, rho_b, phys_dist, phys_mass, eps
    
        implicit none
        
        double precision, intent(in) :: R
        double precision, intent(in) :: sigma
        double precision, intent(in) :: St
        double precision, intent(in) :: T

        double precision :: H
        double precision :: omegaK
        double precision :: cs
        double precision :: dv
        double precision :: rho_dust
        double precision :: rho_bulk
        
        H        = sqrt( adx*T*R**3.d0 )                    ! ---
        omegaK   = 1.d0 / sqrt(R**3.d0)                     !  |
        cs       = H * omegaK                               !  |
                                                            !  |  This is in FARGO units
        dv       = v_rel(cs, St)                            !  |
                                                            !  |
        rho_dust = eps * sigma / (sqrt(2.d0*pi)*H)          !  |
        rho_bulk = rho_b / phys_mass * phys_dist**3.d0      ! ---
        
        coag_rate = rho_dust / rho_bulk * dv * phys_dist  ! This is in physical units
    end function coag_rate

!####################

    double precision function frag_rate(R, sigma, St, T)
    ! Change in radius in physical units due to fragmentation
        use constants
        use variables, only: alpha, adx, rho_b, phys_dist, phys_mass, eps, v_frag
    
        implicit none
        
        double precision, intent(in) :: R
        double precision, intent(in) :: sigma
        double precision, intent(in) :: St
        double precision, intent(in) :: T

        double precision :: H
        double precision :: omegaK
        double precision :: cs
        double precision :: dv
        double precision :: rho_dust
        double precision :: rho_bulk
        double precision :: vf
        double precision :: f
        
        H         = sqrt( adx*T*R**3.d0 )                    ! ---
        omegaK    = 1.d0 / sqrt(R**3.d0)                     !  |
        cs        = H * omegaK                               !  |
                                                             !  |  This is in FARGO units
        dv        = v_rel(cs, St)                            !  |
                                                             !  |
        rho_dust  = eps * sigma / (sqrt(2.d0*pi)*H)          !  |
        rho_bulk  = rho_b / phys_mass * phys_dist**3.d0      ! ---
        
        vf        = v_frag / sqrt( G*phys_mass/phys_dist )   ! Converting to FARGO units
        
        if(dv .GT. vf) then ! Only fragment, when relative velocity large enough
            f     = min( 1.d0, log(dv/vf) / log(5.d0) )      ! Loses at maximum its own sizeper collision
        end if
        
        frag_rate = - f * rho_dust / rho_bulk * dv * phys_dist
    end function frag_rate

!####################

    double precision function v_rel(cs, St)
    ! Relative velocity due to turbulence in FARGO units
        use variables, only: alpha
    
        implicit none
        
        double precision, intent(in) :: cs
        double precision, intent(in) :: St
    
        v_rel = cs * sqrt( 3.d0*alpha*(St+1.d0/St) )
    end function v_rel

end module coagulation
