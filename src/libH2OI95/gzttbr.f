      REAL(8) FUNCTION GZTTBR(delta, ier)

c     This gets the zeta function for Tbarr as a function of delta.
c     The result is used in routines VISCOS and THCOND to calculate
c     the respectie critical enhancement terms mubar2 and labar2.

c     Tbarr (the pertinent reference temperature factor) is 1.5. As
c     Tbar is defined as T/Tcr, Tbarr corresponds to a temperature of
c     1.5*Tcr, or 970.644K. One could calculate the requisite zeta
c     value from the equation-of-state model by calling routine EVAI95
c     for this temperature and the desired delta. However, that
c     seems unnecessarily inefficient for what is required. Here the
c     zeta function was evaluated for 970.644K and delta values
c     corresponding to pressures between 1.e-4 to 1000 MPa, and
c     polynomial curves were fitted to the results over two delta
c     ranges. These polynomials are not part of the original
c     viscosity model.

c     This function is called by:

c       VISCOS.f

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      REAL(8) delta

      INTEGER ier

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER i, j

      REAL(8) cfz(0:4,1:2)

      REAL(8) gx

c-----------------------------------------------------------------------

c     Coefficients for calculating gzttbr.

      DATA cfz(0:4,1) /       0.152966739132486d0,  0.127637911731621d0,
     $  0.098454718356579d0, -0.304726588463916d0,
     $  0.114242784534772d0  /

      DATA cfz(0:4,2) /       0.436894921256486d0, -0.256891980869883d0,
     $ -0.010426094984877d0,  0.032776263097026d0,
     $ -0.005480033061907d0  /

c-----------------------------------------------------------------------

      GZTTBR = 0.0d0

c     Set error flag.

      ier = 0

      j = 0

      IF (delta .LE. 2.0d0) THEN
        j = 1
      ELSEIF (delta .LE. 2.869d0) THEN
        j = 2
      ENDIF

      IF (j .LE. 0) THEN

c       Here delta is out of range.

        ier = 1
        GOTO 999
      ENDIF

      gx = cfz(0,j) + cfz(1,j)*delta

      DO i = 2,4
        gx = gx + cfz(i,j)*delta**i
      ENDDO

      GZTTBR = gx

  999 CONTINUE

      END