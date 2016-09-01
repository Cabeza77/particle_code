module interpolation

    implicit none
    
    public interp2d
    
    contains
    
    recursive double precision function interp2d(x, y, xArr, yArr, mat, Nx, Ny)
    
        implicit none
        
        double precision, intent(in) :: x, y
        double precision, intent(in) :: xArr(Nx), yArr(Ny)
        double precision, intent(in) :: mat(Nx, Ny)
        integer,          intent(in) :: Nx, Ny

        integer          :: ix, iy
        double precision :: val00, val01, val10, val11
        double precision :: x0, x1, y0, y1
        double precision :: dum1, dum2

        ix=0
        iy=0

        call hunt(xArr, Nx, x, ix) 
        call hunt(yArr, Ny, y, iy)
        
        ! Catch out-of-range errors
        ! This is equivalent to extrapolating with the last valid slope
        if(ix <  1 ) ix = 1
        if(ix >= Nx) ix = Nx-1
        if(iy <  1 ) iy = 1
        if(iy >= Ny) iy = Ny-1
        
        val00 = mat(ix,   iy)
        val10 = mat(ix+1, iy)
        val01 = mat(ix,   iy+1)
        val11 = mat(ix+1, iy+1)
        
        x0 = xArr(ix)
        x1 = xArr(ix+1)
        y0 = yArr(iy)
        y1 = yArr(iy+1)
        
        dum1 = (val01-val00)/(y1-y0)*(y-y0)+val00
        dum2 = (val11-val10)/(y1-y0)*(y-y0)+val10
        
        interp2d = (dum2-dum1)/(x1-x0)*(x-x0)+dum1
    
    end function interp2d

end module interpolation
