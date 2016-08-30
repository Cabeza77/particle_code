!     --------------------------------------------------------------
!              NUMERICAL RECIPES RANDOM NUMBER GENERATOR
!     --------------------------------------------------------------
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      doubleprecision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
        & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
        & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      iy=0
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..


!     --------------------------------------------------------------
!              NUMERICAL RECIPES HUNT FUNCTION
!     --------------------------------------------------------------
      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      DOUBLEPRECISION x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1) return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
!  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..



! =============================================================================

! _____________________________________________________________________________
!                   FUNCTION: TEST VALIDITY OF NUMBER
!
!     0 = Number is okay
!     1 = Number is INF
!     2 = Number is NAN
! _____________________________________________________________________________
logical function number_invalid(a)
    implicit none
    doubleprecision a,b,c
    logical div,sub

    b=a*2.d0
    b=b/2.d0
    c=a-1.d100

    div = (b.eq.a)
    sub = (c.lt.a)

    if(div.and.sub) then
        number_invalid = .false.
    elseif(div) then
        number_invalid = .true.
    else
        number_invalid = .true.
    endif

    return

end function number_invalid
