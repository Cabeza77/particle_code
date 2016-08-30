module constants

    implicit none
    
    ! All units in the CGS system
    
    double precision, parameter :: AU     = 1.4960d13
    double precision, parameter :: k_b    = 1.38064852d-16
    double precision, parameter :: m_p    = 1.672621898d-24
    double precision, parameter :: M_sun  = 1.98855d33
    double precision, parameter :: pi     = acos(-1.d0)
    double precision, parameter :: R_gas  = 8.3144598d7
    double precision, parameter :: sig_H2 = 2.d-15
    double precision, parameter :: year   = 365.25d0 * 24.d0 * 3600.d0

end module constants
