      SUBROUTINE CALSCP(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qwrphi, pcr, rcnstw, rhocr, tcr, press, ptripl,
     $ p1atm, ttripl, ts1atm, iter, betamx, bettol, btxtol, tempk,
     $ tau, rhosv, rhosl, deltsv, deltsl, pxv, uxv, sxv, hxv,
     $ cvxv, cpxv, wxv, muxv, dtxv, bsxv, axv, gxv, vxv, pdxv,
     $ ptxv, adxv, gdxv, avxv, ktxv, ztxv, pxl, uxl, sxl, hxl,
     $ cvxl, cpxl, wxl, muxl, dtxl, bsxl, axl, gxl, vxl, pdxl,
     $ ptxl, adxl, gdxl, avxl, ktxl, ztxl)

c     This routine calculates the saturation properties as a function
c     of specified pressure. This is done by iterating on pressure to
c     obtain the desired temperature. This routine calls CALSCT to
c     calculate the saturation pressure and densities (of vapor and
c     liquid). There is no easy way to calculate the requisite
c     derivative (dTsat/dPsat) in the refinement process.  Therefore,
c     the secant method is used instead of the Newton- Raphson method.

      IMPLICIT NONE

c     CALLING SEQUENCE VARIABLES.

      LOGICAL qerr, qfail, qprnt1, qprnt2, qprnt3, qwrphi

      INTEGER nttyo, noutpt

      INTEGER iter

      REAL(8) press

      REAL(8) betamx, bettol, btxtol

      REAL(8) ptripl, pcr, p1atm, rcnstw, rhocr, ttripl, ts1atm, tcr

      REAL(8) tempk, tau, rhosv, rhosl, deltsv, deltsl

      REAL(8) pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv, dtxv,
     $ bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv, ktxv, ztxv

      REAL(8) pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl, dtxl,
     $ bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl, ktxl, ztxl

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     LOCAL VARIABLES.

      LOGICAL qsilent

      INTEGER itermx, kter

      CHARACTER(LEN=8) unam8

      REAL(8) tempk0, psat, psat0, tsat, tsatp

      REAL(8) betmx0

      REAL(8) arelax, dtdp

      REAL(8) pdiff

      DATA qsilent / .TRUE. /

      DATA itermx / 35 /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      qerr = .FALSE.
      qfail = .FALSE.
      btxtol = bettol

      WRITE (nttyo,1050) press
      WRITE (noutpt,1050) press

 1050 FORMAT(/3x,'CALSCP: press  = ',e16.9,' MPa',/)

c     Calculate approximate saturation temperature to use
c     as a starting estimate.

      CALL APXSCP(noutpt, nttyo, qprnt1, rhocr, rhosv, rhosl,
     $ ttripl, ts1atm, tcr, ptripl, p1atm, pcr, press, tsat,
     $ tsatp, qerr)

      IF (qerr) GO TO 410

      tempk = tsat
      tau = tcr/tempk

      arelax = 1.0d0
      iter = 0

c     Using the current temperature estimate, calculate
c     the corresponding saturation properties. This will
c     include a calculated pressure to compare with the
c     specified pressure.

c     Calling sequence substitutions:

c       kter for iter
c       psat for press

      CALL CALSCT(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qsilent, qwrphi, pcr, rcnstw, rhocr, tcr, kter,
     $ betamx, bettol, btxtol, tempk, tau, psat, rhosv, rhosl,
     $ deltsv, deltsl, pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv,
     $ dtxv, bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv,
     $ ktxv, ztxv, pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl,
     $ dtxl, bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl,
     $ ktxl, ztxl)

      IF (qerr) GO TO 410

      pdiff = psat - press
      betamx = DABS(pdiff/press)

      WRITE (noutpt,1210) iter, tempk, betamx
      WRITE (nttyo,1210) iter, tempk, betamx

 1210 FORMAT(6x,'CALSCP: iter= ',i3,', tempk= ',f8.4,', betamx= ',
     $ 1pe12.5)

      IF (betamx .LE. btxtol) THEN

        IF (qprnt1) THEN

          WRITE (nttyo,1230) tempk, betamx, btxtol
          WRITE (noutpt,1230) tempk, betamx, btxtol

 1230     FORMAT(/3x,'CALSCP: ITERATION CONVERGED TO tempk= ',f8.4,
     $    ' (K).',/3x,'MAX NORM betamx = ',1pe10.3,
     $    ', TOLERANCE btxtol = ',e10.3)

        ENDIF

        GO TO 420
      ENDIF

      psat0 = psat
      tempk0 = tempk
      betmx0 = betamx

c     Make a small change in the temperature value to
c     generate a slope estimate.

      tempk = tempk0 + 2.0d0
      IF (tempk .GT. tcr) tempk = tempk0 - 2.0d0
      tau = tcr/tempk

c     The following is a return point for iteration.

  120 iter = iter + 1

c     Calculate the saturation properties for the
c     current temperature value.

c     Calling sequence substitutions:

c       kter for iter
c       psat for press

      CALL CALSCT(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qsilent, qwrphi, pcr, rcnstw, rhocr, tcr, kter,
     $ betamx, bettol, btxtol, tempk, tau, psat, rhosv, rhosl,
     $ deltsv, deltsl, pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv,
     $ dtxv, bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv,
     $ ktxv, ztxv, pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl,
     $ dtxl, bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl,
     $ ktxl, ztxl)

      IF (qerr) GO TO 410

      pdiff = psat - press
      betamx = DABS(pdiff/press)

      WRITE (noutpt,1210) iter, tempk, betamx
      WRITE (nttyo,1210) iter, tempk, betamx

      IF (betamx .LE. btxtol) THEN

        IF (qprnt1) THEN
          WRITE (nttyo,1230) tempk, betamx, btxtol
          WRITE (noutpt,1230) tempk, betamx, btxtol
        ENDIF

        GO TO 420
      ENDIF

      IF (iter .GE. itermx) THEN

c       Have done the maximum number of iterations.

        WRITE (nttyo,1240)
        WRITE (noutpt,1240)

 1240   FORMAT(/3x,'CALSCP: HAVE DONE THE MAXIMUM NUMBER OF',
     $  ' ITERATIONS.')

        GO TO 410
      ENDIF

c     Estimate the slope (dTsat/dPsat).

      dtdp = (tempk - tempk0)/(psat - psat0)

c     Choose the closest of the two active points as the
c     new base point.

      IF (betamx .LE. betmx0) THEN
        psat0 = psat
        tempk0 = tempk
        betmx0 = betamx
      ENDIF

c     Calculate a corrected temperature value.

      tempk = tempk0 + arelax*dtdp*(press - psat0)
      tau = tcr/tempk

      GO TO 120

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  410 CONTINUE

c     Iteration failed.

      qfail = .TRUE.

      WRITE (nttyo,1270)
      WRITE (noutpt,1270)

 1270 FORMAT(/3x,'CALSCP: ITERATION FAILED.',/)

      GO TO 999

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  420 CONTINUE

c     Calculations are complete for the last line read
c     from the input file.

c     Write the results.

      WRITE (nttyo,2400)
      WRITE (noutpt,2400)

 2400 FORMAT(/3x,'CALSCP key results:',//8x,'Temp(K)',10x,'press(MPa)')

      WRITE (nttyo,2410) tempk,psat
      WRITE (noutpt,2410) tempk,psat

 2410 FORMAT(6x,f9.4,7x,e16.9)

      WRITE (noutpt,1280) tempk, psat, tau, deltsv, rhosv,
     $ deltsl, rhosl

 1280 FORMAT(/6x,'Temp(K) = ',f9.4
     $ //9x,'psat   = ',e16.9,' MPa',6x,'tau    = ',e16.9,
     $ /9x,'deltsv = ',e16.9,10x,'rhosv  = ',e16.9,
     $ /9x,'deltsl = ',e16.9,10x,'rhosl  = ',e16.9)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 999  CONTINUE

      END