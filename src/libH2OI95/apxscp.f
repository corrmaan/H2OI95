      SUBROUTINE APXSCP(noutpt, nttyo, qprnt1, rhocr, rhosv, rhosl,
     $ ttripl, ts1atm, tcr, ptripl, p1atm, pcr, press, tsat,
     $ tsatp, qerr)

c     Evaluate the temperature (tsat) as a function of pressure along
c     the vapor-liquid equilibrium curve, using equation 2.5 of Wagner
c     and Pruss (2002). Also calculate the densities of the liquid and
c     vapor phases using equations 2.6 and 2.7 from the same source.
c     Results may not be fully consistent with IAPWS-95, but may serve
c     as close approximations or starting values for refinement.
c     This routine is similar to APXSCT, but evaluates the inverse
c     problem (tsat as a function of press instead of psat as a function
c     of tempk). Newton-Raphson iteration is used. This routine calls
c     APXSCT to calclate the pressure for a given temperature value and
c     the derivative of pressure with respect to temperature.

c       tsat = saturation temperature (K)
c       rhosv = density of vapor (kg/m3)
c       rhosl = density of liquid (kg/m3)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qprnt1, qerr

      INTEGER nttyo, noutpt

      REAL(8) rhosv, rhosl, rhocr, press, tcr, pcr

      REAL(8) ptripl, p1atm, ttripl, ts1atm

      REAL(8) tsat, tsatp

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      LOGICAL qsilent

      INTEGER iter, itermx

c     Note that betamx and bettol used here are local to this routine.

      REAL(8) pdiff, betamx, bettol, dtdp, psat, psatt, tempk, tempk0

      REAL(8) arelax

      DATA qsilent / .TRUE. /

      DATA itermx / 35 /

      DATA bettol / 1.0d-9 /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (qprnt1) THEN

        WRITE(nttyo,1100) press
        WRITE(noutpt,1100) press

 1100   FORMAT(/3x,'APXSCP: Input press= ',1pe12.5,' (MPa)',/)

      ENDIF

c     Check to see that the pressure is in the allowed range.

      qerr = .FALSE.
      IF (press.LT.ptripl .OR. press.GT.pcr) THEN
        qerr = .TRUE.
        GO TO 999
      ENDIF

c     Choose a starting temperature value. The exact value
c     is not that important.

      IF (press .GE. p1atm) THEN

c       Interpolate between 100C, 1.01325 bar and tcr, pcr.

        dtdp = (tcr - ts1atm)/(pcr - p1atm)
        tempk = ts1atm + dtdp*(press - p1atm)
      ELSE

c       Interpolate between the triple point and 100C, 1.013 bar

        dtdp = (ts1atm - ttripl)/(p1atm - ptripl)
        tempk = ttripl + dtdp*(press - ptripl)
      ENDIF

      iter = 0
      arelax = 0.25d0

  120 CALL APXSCT(noutpt, nttyo, qprnt1, qsilent, rhocr, rhosv,
     $ rhosl, tempk, tcr, pcr, psat, psatt, qerr)

      pdiff = psat - press
      betamx = DABS(pdiff/press)

      IF (qprnt1) THEN

        WRITE(nttyo,1110) iter, tempk, betamx
        WRITE(noutpt,1110) iter, tempk, betamx

 1110   FORMAT(6x,'APXSCP: iter= ',i3,', tempk= ',f8.4,', betamx= ',
     $  e12.5)

      ENDIF

      IF (betamx .LE. bettol) THEN

c       Have converged.

        tsat = tempk

        IF (qprnt1) THEN

          WRITE (nttyo,1130) tempk, betamx, bettol
          WRITE (noutpt,1130) tempk, betamx, bettol

 1130     FORMAT(/3x,'APXSCP: ITERATION CONVERGED TO tempk= ',f8.4,
     $    ' (K).',/3x,'MAX NORM betamx = ',1pe10.3,
     $    ', TOLERANCE bettol = ',e10.3)

        ENDIF

        GO TO 999

      ENDIF

      IF (iter .GE. itermx) THEN

c       Have done the maximum number of iterations.

        qerr = .TRUE.

        WRITE (nttyo,1150)
        WRITE (noutpt,1150)

 1150   FORMAT(/3x,'HAVE DONE THE MAXIMUM NUMBER OF ITERATIONS.')

        GO TO 999
      ENDIF

c     Estimate a new value for tempk.

      iter = iter + 1
      tempk0 = tempk
      tsatp = 1.0d0/psatt
      IF (iter .GE. 2) arelax = 1.0d0

c     Make a Newton-Raphson correction.

      tempk = tempk0 + arelax*tsatp*(press - psat)

      GO TO 120

  999 CONTINUE

      END