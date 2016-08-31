module memory

    implicit none
    
    public allocate_memory, deallocate_memory
    
    contains
    
    subroutine allocate_memory()
        use variables
        implicit none
        
        allocate( R(N_R) )
        allocate( R_int(N_R+1) )
        allocate( theta(N_theta) )
        allocate( theta_int(N_theta+1) )
        
        allocate( sigma_gas(2, N_theta, N_R) )
        allocate( vTheta_gas(2, N_theta, N_R) )
        allocate( vR_gas(2, N_theta, N_R) )
        allocate( T_gas(2, N_theta, N_R) )
        
        allocate( x_planet(0:N_t) )
        allocate( y_planet(0:N_t) )
        allocate( vX_planet(0:N_t) )
        allocate( vY_planet(0:N_t) )
        allocate( R_planet(0:N_t) )
        allocate( theta_planet(0:N_t) )
        allocate( vR_planet(0:N_t) )
        allocate( vTheta_planet(0:N_t) )
        allocate( m_planet(0:N_t) )
        
        allocate( x_dust(N_dust) )
        allocate( y_dust(N_dust) )
        allocate( vX_dust(N_dust) )
        allocate( vY_dust(N_dust) )
        allocate( R_dust(N_dust) )
        allocate( theta_dust(N_dust) )
        allocate( vR_dust(N_dust) )
        allocate( vTheta_dust(N_dust) )
        allocate( L_dust(N_dust) )
        allocate( m_dust(N_dust) )
        allocate( a_dust(N_dust) )
        allocate( T_dust(N_dust) )
        allocate( fc_dust(N_dust) )
        allocate( St(N_dust) )
        allocate( tstop(N_dust) )
        
        allocate( time(0:N_t) )
        
        allocate( cur_sigma_gas(N_theta, N_R) )
        allocate( cur_vR_gas(N_theta, N_R) )
        allocate( cur_vTheta_gas(N_theta, N_R) )
        allocate( cur_T_gas(N_theta, N_R) )
        
    end subroutine allocate_memory
    
    subroutine deallocate_memory()
        use variables
        implicit none
        
        deallocate( R )
        deallocate( R_int )
        deallocate( theta )
        deallocate( theta_int )
        
        deallocate( sigma_gas )
        deallocate( vTheta_gas )
        deallocate( vR_gas )
        deallocate( T_gas )
        
        deallocate( x_planet )
        deallocate( y_planet )
        deallocate( vX_planet )
        deallocate( vY_planet )
        deallocate( R_planet )
        deallocate( theta_planet )
        deallocate( vR_planet )
        deallocate( vTheta_planet )
        deallocate( m_planet )
        
        deallocate( x_dust )
        deallocate( y_dust )
        deallocate( vX_dust )
        deallocate( vY_dust )
        deallocate( R_dust )
        deallocate( theta_dust )
        deallocate( vR_dust )
        deallocate( vTheta_dust )
        deallocate( L_dust )
        deallocate( m_dust )
        deallocate( a_dust )
        deallocate( T_dust )
        deallocate( fc_dust )
        deallocate( St )
        deallocate( tstop )
        
        deallocate( time )
        
        deallocate( cur_sigma_gas )
        deallocate( cur_vR_gas )
        deallocate( cur_vTheta_gas )
        deallocate( cur_T_gas )
        
    end subroutine deallocate_memory

end module memory
