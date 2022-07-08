      REAL(8) FUNCTION TKMELT(press, noutpt, nttyo, ier)

c     This function gets the melting temperature of ice as a function
c     of pressure. This inverts the model for the melting pressure
c     as a function of temperature. That model is described in
c     IAPWS R14-08(2011), The International Association for the
c     Properties of Water and Steam, 2011, Revised Release on the
c     Pressure along the Melting and Sublimation Curves of Ordinary
c     Water Substance. This document may be found at:
c     http://www.iapws.org/relguide/MeltSub.html

c     The above model for the melting pressure is treated in
c     FUNCTION PMELT. The actual pressure must be specified in
c     order to resolve the problem of overlap in the temperature
c     ranges of the individual phases. See discussion in FUNCTION
c     PMELT.

c     Inversion of the model for the melting pressure is done
c     here using the secant method. This is chosen instead of the
c     Newton-Raphson method to avoid potential problems with slope
c     discontinuites at boundaries between the ice phases for
c     pressures above 208.566 MPa, which is the equilibrium pressure
c     for Ice Ih-Ice III-liquid. The corresponding equlibrium
c     temperature is 251.165K. Putative melting temperatures should
c     not be less than this for pressures above 208.566 MPa, nor
c     more than this for pressures less than this.

c     This function is called by:

c       H2OI95.f

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      REAL(8) press

      INTEGER ier

      INTEGER noutpt, nttyo

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER i, ier1, j2, j3, n, nrlx

      INTEGER ILNOBL

      REAL(8) PMELT

      REAL(8) t0, t1, p0, p1, tt, dp, dt, res, tlim0, tlim1, tx

      REAL(8) delt, rlx

      CHARACTER(LEN=24) ux24a, ux24b

      LOGICAL qwr

      DATA qwr / .FALSE. /

c-----------------------------------------------------------------------

      TKMELT = 0.0d0

c     Set error flag.

      ier = 0

c     Find the melting temperature that corresponds to the actual
c     pressure. The variables t0 and t1 are initial values for
c     secant method iteration.

      IF (press .LT. 208.566d0) THEN

c       In the Ice Ih field.

        tlim0 = 251.165d0
        tlim1 = 273.16d0

        t0 = tlim0
        t1 = tlim1

      ELSEIF (press.GE.208.566d0 .AND. press.LT.350.1d0) THEN

c       In the Ice III field.

        tlim0 = 251.165d0
        tlim1 = 256.164d0

        t0 = tlim0
        t1 = tlim1

      ELSEIF (press.GE.350.1d0 .AND. press.LT.632.4d0) THEN

c       In the Ice V field.

        tlim0 = 256.164d0
        tlim1 = 273.31d0

        t0 = tlim0
        t1 = tlim1

      ELSEIF (press.GE.632.4d0 .AND. press.LT.2216.0d0) THEN

c       In the Ice VI field.

        tlim0 = 273.31d0
        tlim1 = 355.0d0

        tx = tlim1 - tlim0
        t0 = tlim0 + 0.3d0*tx
        t1 = tlim0 + 0.7d0*tx

      ELSEIF (press.GE.2216.0d0 .AND. press.LE.10000.0d0) THEN

c       In the Ice VII field.
c       Note: the upper limit here of 10000 MPa is an arbitrary
c       cutoff suggested by Figure 1 from IAPWS R14, but this
c       is not part of the IAPWS R14 standard.

        tlim0 = 355.0d0
        tlim1 = 1000.0d0

        tx = tlim1 - tlim0
        t0 = tlim0 + 0.3d0*tx
        t1 = tlim0 + 0.7d0*tx

      ELSEIF (press.GT.20000.0d0) THEN

        ier = 1
        GOTO 990

      ENDIF

c     Calling sequence substitutions:
c       t1 for tempk
c       ier1 for ier

      p1 = PMELT(t1, press, noutpt, nttyo, ier1)

      IF (ier1 .GT. 0) THEN

        WRITE(ux24a,'(f16.9)') t1
        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(ux24b,'(f16.9)') press
        CALL LEJUST(ux24b)
        j3 = ILNOBL(ux24b)

        WRITE(nttyo,1010) ux24a(1:j2), ux24b(1:j3)
        WRITE(noutpt,1010) ux24a(1:j2), ux24b(1:j3)

 1010   FORMAT(/6x,'ERROR (TKMELT): PMELT failed for: t1 = ',a,'K,',
     $  ' press = ',a,' MPa')

        ier = 1
        GOTO 999

      ENDIF

      IF (qwr) THEN
        WRITE(noutpt,'( )')
      ENDIF

      n = 0

c     In general, t0 is the temperature to be tested nad
c     t1 is the other temperature used in refining t0.
c     Here test the starting t1 in case it is the solution.

      res = p1 - press

      IF (ABS(res) .LE. 1.0d-6) THEN
        TKMELT = t1

        IF (qwr) THEN

          WRITE(ux24a,'(f16.9)') t1
          CALL LEJUST(ux24a)
          j2 = ILNOBL(ux24a)

          WRITE(ux24b,'(f16.9)') res
          CALL LEJUST(ux24b)
          j3 = ILNOBL(ux24b)

          WRITE(noutpt,1020) n, ux24a(1:j2), ux24b(1:j3)

 1020     FORMAT(6x,'TKMELT: iter= ',i2,', tmelt= ',a,'K, residual= ',a,
     $    ' MPa')

        ENDIF

        GOTO 999
      ENDIF

  110 CONTINUE

c     Iterative loop, correcting t0.

c     Calling sequence substitutions:
c       t0 for tempk
c       ier1 for ier

      p0 = PMELT(t0, press, noutpt, nttyo, ier1)

      IF (ier1 .GT. 0) THEN

        WRITE(ux24a,'(f16.9)') t0
        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(ux24b,'(f16.9)') press
        CALL LEJUST(ux24b)
        j3 = ILNOBL(ux24b)

        WRITE(noutpt,1020) n, ux24a(1:j2), ux24b(1:j3)

        ier = 1
        GOTO 999

      ENDIF

      dp = p1 - p0
      dt = t1 - t0

      res = p0 - press

      IF (qwr) THEN

        WRITE(ux24a,'(f16.9)') t0
        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(ux24b,'(f16.9)') res
        CALL LEJUST(ux24b)
        j3 = ILNOBL(ux24b)

        WRITE(noutpt,1020) n, ux24a(1:j2), ux24b(1:j3)

      ENDIF

      IF (ABS(res) .LE. 1.0d-6) THEN
        TKMELT = t0
        GOTO 999
      ENDIF

      n = n + 1

      IF (n .GE. 100) GOTO 980

      delt = - (dt/dp)*res

      rlx = 1.0d0
      nrlx = 0

  210 tt = t0 + rlx*delt

c     Check to see that tt is within the specified temperature
c     limits. If not, apply under-relaxation.

      IF (tt.LT.tlim0 .OR. tt.GT.tlim1) THEN
        rlx = 0.5d0*rlx
        nrlx = nrlx + 1
        IF (nrlx .GT. 10) GOTO 980
        GOTO 210
      ENDIF

c     Determine whether tt is closer to t0 or t1.
c     Discard the more distant point.

      IF (abs(tt - t0) .LE. abs(tt - t1)) THEN

c       Make make t0 the new t1.

        t1 = t0
        p1 = p0

      ENDIF

c     Make tt the new t0.

      t0 = tt

      GOTO 110

  980 WRITE(ux24a,'(f16.9)') press
      CALL LEJUST(ux24a)
      j2 = ILNOBL(ux24a)

      WRITE(nttyo,1060) ux24a(1:j2)
      WRITE(noutpt,1060) ux24a(1:j2)

 1060 FORMAT(/6x,'ERROR (TKMELT): Iteration failed to find the melting',
     $ /9x,'temperature for ',a,' MPa')

      GOTO 999

  990 WRITE(ux24a,'(f16.9)') press

      CALL LEJUST(ux24a)
      j2 = ILNOBL(ux24a)

      WRITE(nttyo,1070) ux24a(1:j2)
      WRITE(noutpt,1070) ux24a(1:j2)

 1070 FORMAT(/6x,'ERROR (TKMELT): Cannot identify the correct',
     $ ' ice phase for ',/9x,a,' MPa')

  999 CONTINUE

      END