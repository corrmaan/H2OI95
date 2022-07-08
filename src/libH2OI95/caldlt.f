      SUBROUTINE CALDLT(nttyo, noutpt, qerr, qprnt1, qprnt2, qprnt3,
     $ qwrphi, pcr, rcnstw, rhocr, tcr, iter, tempk, tau, delta,
     $ rho, press, px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $ ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx)

c     Calculate properties including pressure as a function
c     of specified temperature and delta (reduced density). This
c     requires a direct evaluation of the equation of state.
c     No iteration is required, as when temperature and pressure
c     are specified.

c     This routine is indifferent to whether or not the value
c     of delta lies in the stable region. Values of delta that
c     lie between the delta of vapor and the delta of liquid on
c     the saturation curve may correspond to metastable or unstable
c     regions. All other code input options produce results only
c     in the stable field.

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qerr, qprnt1, qprnt2, qprnt3, qwrphi

      INTEGER nttyo, noutpt

      INTEGER iter

      REAL(8) pcr, rcnstw, rhocr, tcr

      REAL(8) delta, rho, tempk, tau

      REAL(8) press

      REAL(8) px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $ ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      LOGICAL qsilnx

      INTEGER ILNOBL

      INTEGER j2

      CHARACTER(LEN=8) unam8

      DATA qsilnx / .FALSE. /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      unam8 = 'H2O'
      j2 = ILNOBL(unam8)

      IF (qprnt1) THEN

        WRITE (noutpt,1150) unam8(1:j2), delta, rho

 1150   FORMAT(/6x,a,/9x,'delta  = ',e16.9,10x,'rho   = ',e16.9)

      ENDIF

c     The routine called here evaluates the equation-of-state
c     as a function of temperature and reduced density (delta).

c     Calling sequence substitutions:

c       qsilnx for qsilent

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, delta, rho, px, ax, ux, sx, hx,
     $ cvx, cpx, wx, mux, dtx, bsx, gx, vx, pdx, ptx, adx,
     $ gdx, avx, ktx, ztx, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilnx)

c     Calculations are complete for the last line read
c     from the input file. Write the results.

      WRITE (nttyo,1160)
      WRITE (noutpt,1160)

 1160 FORMAT(/8x,'Temp(K)',10x,'press(MPa)')

      WRITE (nttyo,1170) tempk,px
      WRITE (noutpt,1170) tempk,px

 1170 FORMAT(6x,f9.4,7x,e16.9)

      END