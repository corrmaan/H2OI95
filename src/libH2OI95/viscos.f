      REAL(8) FUNCTION VISCOS(delta, tau, pdx, ztx, tempk, press,
     $ rho, rhocr, rcnstw, tcr, pcr, ptripl, tmltx, mubar,
     $ nttyo, noutpt)

c     This function calculates the viscosity of water using the
c     "IAPWS Formulation 2008 for the Viscosity of Ordinary Water
c     Substance" (IAPWS R12-08). The corresponding journal reference
c     is: Huber M.L., Perkins R.A., Laesecke A., Friend D.G.,
c     Sengers J.V., Assael M.J., Metaxa I.N., Vogel E., Mares R.,
c     and Miyagawa K. (2009) New International Formulation for the
c     Viscosity of H2O. J. Phys. Chem. Ref. Data 38, 101-125.
c     Equations referenced below are from IAPWS R12-08.

c     This function is called by:

c       H2OI95.f

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      REAL(8) delta, tau, pdx, ztx, tempk, press, rho, rhocr, tcr, pcr,
     $ ptripl, rcnstw, tmltx, mubar

      INTEGER nttyo, noutpt

      INTEGER ILNOBL

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER i, j, j2, j3, j4, j5

      INTEGER ier

      CHARACTER(LEN=24) ux24a, ux24b, ux24c, ux24d, ux24e

      REAL(8) muref, sx0, sx1, sx2, dltm1, tbar, taum1, xx, zdx

      REAL(8) delchb, xicl, www, lcap, psid, ycap

      REAL(8) cfhi(0:3), cfhij(0:5,0:6)

      REAL(8) xmu, qcem1, qdem1, qc, qd, nuc, gamma, xicl0, gam0, tbarr

      REAL(8) qcxicl, qcfac, qdxicl, exx, xxx, yxx, yxx1, yxx2, yxx3,
     $ yxx4, ycf, yy1, yy2

      REAL(8) mubar0, mubar1, mubar2

      REAL(8) GZTTBR, TKMELT

      REAL(8) zettbr

      REAL(8) rhor, taur, tkr
      REAL(8) pxr, uxr, sxr, hxr, cvxr, cpxr, wxr, muxr, dtxr, bsxr,
     $ axr, gxr, vxr, pdxr, ptxr, adxr, gdxr, avxr, ktxr, ztxr

      LOGICAL qrchk

      LOGICAL qsilent, qwrphi, qprnt1, qprnt2, qprnt3

      DATA qsilent, qwrphi, qprnt1, qprnt2, qprnt3 / .FALSE., .FALSE.,
     $ .FALSE., .FALSE., .FALSE. /

c     The following are local hard-wired switches used in code
c     verification of viscosity calculations.

      LOGICAL qnomu2, qcbzet, qzetap, qwr1, qwr2, qwr3

      DATA qnomu2, qcbzet, qzetap, qwr1, qwr2, qwr3 / .FALSE., .FALSE.,
     $ .FALSE., .FALSE., .FALSE., .FALSE. /

c-----------------------------------------------------------------------

c     Here muref is the reference viscosity (Pa-s).

      DATA muref /1.00d-06 /

c     Coefficients for calculating Hi.

      DATA cfhi(0:3) / 1.67752d+00, 2.20462d+00, 0.6366564d+00,
     $ -0.241605d+00 /

c     Coefficients for calculating Hij.

      DATA cfhij(0:5,0) / 5.20094d-01, 8.50895d-02, -1.08374d+00,
     $ -2.89555d-01, 0.0d+00, 0.0d+00 /

      DATA cfhij(0:5,1) / 2.22531d-01, 9.99115d-01, 1.88797d+00,
     $ 1.26613d+00, 0.0d+00, 1.20573d-01 /

      DATA cfhij(0:5,2) / -2.81378d-01, -9.06851d-01, -7.72479d-01,
     $ -4.89837d-01, -2.57040d-01, 0.0d+00 /

      DATA cfhij(0:5,3) / 1.61913d-01, 2.57399d-01, 0.0d+00,
     $ 0.0d+00, 0.0d+00, 0.0d+00 /

      DATA cfhij(0:5,4) / -3.25372d-02, 0.0d+00, 0.0d+00,
     $ 6.98452d-02, 0.0d+00, 0.0d+00 /

      DATA cfhij(0:5,5) / 0.0d+00, 0.0d+00, 0.0d+00,
     $ 0.0d+00, 8.72102d-03, 0.0d+00 /

      DATA cfhij(0:5,6) / 0.0d+00, 0.0d+00, 0.0d+00,
     $ -4.356730d-03, 0.0d+00, -5.93264d-04 /

c     Coefficients for calculating the critical enhancement term
c     (mubar2).

      DATA xmu, qcem1, qdem1, nuc, gamma, xicl0, gam0, tbarr /
     $ 0.068d+00, 1.9d+00, 1.1d+00, 0.630d+00, 1.239d+00, 0.13d+00,
     $ 0.06d+00, 1.5d00 /

c-----------------------------------------------------------------------

c     Check range of validity.

      qrchk = .FALSE.

      IF (press.GT.0.0d0 .AND. press.LT.ptripl) THEN

        IF (tempk.GE.273.16d0 .AND. tempk.LE.1173.15d0) THEN
          qrchk = .TRUE.
        ENDIF

      ELSEIF (tmltx .GT. 0.0d0) THEN

        IF (press.GE.ptripl .AND. press.LE.300.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.1173.15d0) THEN
            qrchk = .TRUE.
          ENDIF

        ELSEIF (press.GT.300.0d0 .AND. press.LE.350.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.873.15d0) THEN
            qrchk = .TRUE.
          ENDIF

        ELSEIF (press.GT.350.0d0 .AND. press.LE.500.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.433.15d0) THEN
            qrchk = .TRUE.
          ENDIF

        ELSEIF (press.GT.500.0d0 .AND. press.LE.1000.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.373.15d0) THEN
            qrchk = .TRUE.
          ENDIF

        ENDIF

      ENDIF

      IF (.NOT.qrchk) THEN

c       Failed check on valid temperature-pressure range.
c       Issue an appropriate warning.

        IF (tmltx .LE. 0.0d0) THEN

          WRITE(noutpt,1020)
          WRITE(nttyo,1020)

 1020     FORMAT(/6x,'WARNING (VISCOS): Could not calculate the',
     $    /9x,'melting temperature. Hence could not check the',
     $    /9x,'range of validity of the viscosity model.')

        ELSEIF (.NOT.qrchk) THEN

          WRITE(noutpt,1050)
          WRITE(nttyo,1050)

 1050     FORMAT(/6x,'WARNING (VISCOS): Outside the range of validity',
     $    /9x,' of the viscosity model')
        ENDIF

      ENDIF

c     Note: Although this viscosity model must be used with the
c     IAPWS-95 equation-of-state model, the language employed is
c     slightly different. "Tbar" (tbar) is equivalent to 1/tau. Also,
c     "rhobar" is equivalent to delta. Insofar as possible, we will
c     use tau and delta in evaluating the appropriate equations.
c     "Tbar" is only necessary in one instance, in the calculation
c     of mubar0.

c     Here the scaled viscosity (mubar) depends on three factors,
c     mubar0, mubar1, and mubar2. Note that
c     mubar = mubar1*mubar2*mubar3, not mubar1 + mubar2 + mubar3.

      tbar = 1.0d0/tau

c     Calculate the viscosity in the dilute gas limit (mubar0).

      mubar0 = 1.0d0

      sx0 = cfhi(0) + cfhi(1)*tau
      DO i = 2,3
        sx0 = sx0 + cfhi(i)*(tau**i)
      ENDDO

      mubar0 = 100*SQRT(tbar)/sx0

c     Calculate the contribution to viscosity due to finite density
c     (mubar1).

      mubar1 = 1.0d0

      dltm1 = delta - 1.0d0
      taum1 = tau - 1.0d0

      sx1 = 0.0d0
      DO i = 0,5

        sx2 = cfhij(i,0) + cfhij(i,1)*dltm1
        DO j = 2,6
          sx2 = sx2 + cfhij(i,j)*(dltm1**j)
        ENDDO
        sx1 = sx1 + (taum1**i)*sx2
      ENDDO
      xx = delta*sx1
      mubar1 = EXP(xx)

c     Calculate the critical enhancement (mubar2), which is only
c     signficant in a very small region of density and temperature
c     around the critical point.

      mubar2 = 1.0d0

c     If not using mubar2, leave its value at unity.

      IF (qnomu2) GO TO 800

c     Get zeta at the reference temperature (1.5*Tcr, or 970.644K).

      zettbr = 0.0d0
      ztxr = 0.0d0

      IF (qcbzet .OR. qzetap) THEN

c       Estimate zeta at the reference temperature using fitted
c       polynomials). This is experimental.

        zettbr = GZTTBR(delta, ier)

        IF (ier .GT. 0) THEN

c         Error.

          WRITE(nttyo,1210)

 1210     FORMAT(/6x,'WARNING (VISCOS): Could not estimate the',
     $    ' value of zeta for the reference',/9x,'temperature',
     $    ' of 1.5*Tcr using a fitted polynomial.')

          zettbr = 0.0d0

        ENDIF

      ENDIF

      IF (qcbzet .OR. .NOT.qzetap .OR. zettbr.EQ.0.0d0) THEN

c       Estimate the reference zeta value by directly evaluating
c       the EOS model. This is more accurate.

        tkr = tbarr*tcr
        taur = 1.0d0/tbarr
        rhor = rho

c       Calling sequence substitutions:
c         tkr for tempk
c         taur for tau
c         rhor for rho
c         ztxr for ztx
c         Plus many others

        CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $  press, tkr, taur, delta, rhor, pxr, axr, uxr, sxr, hxr,
     $  cvxr, cpxr, wxr, muxr, dtxr, bsxr, gxr, vxr, pdxr, ptxr, adxr,
     $  gdxr, avxr, ktxr, ztxr, qprnt1, qprnt2, qprnt3,
     $  qwrphi, qsilent)

      ENDIF

      IF (qwr1) THEN

        IF (qcbzet) THEN

c         Compare the reference zeta value obtained by the two methods.

          zdx = zettbr - ztxr

c         Write zettbr and ztxr on the standard output.

          WRITE(ux24a,'(f16.9)') zettbr
          CALL LEJUST(ux24a)
          j2 = ILNOBL(ux24a)

          WRITE(ux24b,'(f16.9)') ztxr
          CALL LEJUST(ux24b)
          j3 = ILNOBL(ux24b)

          WRITE(ux24c,'(f16.9)') zdx
          CALL LEJUST(ux24c)
          j4 = ILNOBL(ux24c)

          WRITE(nttyo,1220) ux24a(1:j2), ux24b(1:j3), ux24c(1:j4)

 1220     FORMAT(/6x,'VISCOS: ztxr   = ',a,' (Poly)',/
     $    14x,'ztxr   = ',a,' (EOS)',/14x,'diff   = ',a)

        ENDIF

      ENDIF

      IF (qzetap .AND. zettbr.GT.0.0d0) THEN

        ztxr = zettbr
        ux24e = '(used, Poly)'

      ELSE

        ux24e = '(used, EOS)'

      ENDIF

      IF (qwr1) THEN

c       Print the zeta value for the refrence temperature.

        WRITE(ux24b,'(f16.9)') ztxr
        CALL LEJUST(ux24b)
        j3 = ILNOBL(ux24b)
        j4 = ILNOBL(ux24e)

        WRITE(nttyo,1250) ux24b(1:j3), ux24e(1:j4)

 1250   FORMAT(/6x,'VISCOS: ztxr   = ',a,1x,a)

      ENDIF

c     Get delchb (DELTA chibar, equation 21).

      delchb = delta*( ztx - (ztxr*tbarr*tau) )
      IF (delchb .LT. 0.0d0) delchb = 0.0d0

c     Get xicl (equation 20).

      xxx = delchb/gam0
      exx = nuc/gamma
      xicl = xicl0*( xxx**exx )

      IF (qwr2) THEN

c       Write xicl on the standard output.

        WRITE(ux24a,'(f16.9)') xicl

        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(nttyo,1260) ux24a(1:j2)

 1260   FORMAT(/6x,'VISCOS: xicl = ',a)

      ENDIF

c     Get psid (equation 17).

      qd = 1.0d0/qdem1
      qdxicl = qd*xicl

      xx = ( 1.0d0 + qdxicl**2 )**(-0.5d0)
      psid = ACOS(xx)

c     Get www (w in equation 19).

      qc = 1.0d0/qcem1
      qcxicl = qc*xicl
      qcfac = ABS( (qcxicl - 1.0d0)/(qcxicl + 1.0d0) )

      www = SQRT(qcfac) * TAN(0.5d0*psid)

c     Get lcap (L(w) in equation 18).

      lcap = 0.0d0
      IF (qcxicl .GT. 1.0d0) THEN
        lcap = LOG( (1.0d0 + www)/(1.0d0 - www) )
      ELSE
        lcap = 2.0d0*ATAN( ABS(www) )
      ENDIF

c     Get ycap (Y is equations 15 and 16).

      ycap = 0.0d0

      IF (xicl .GE. 0.0d0) THEN

        IF (xicl .LE. 0.3817016416d0) THEN
          yxx  = 1.0d0 - qcxicl + qcxicl**2
     $           - ( 765.0d0/504.0d0 )*qdxicl**2
          ycap  = 0.2d0*qcxicl*(qdxicl**5)*yxx
        ELSE
          yxx1 = ( 1.0d0/12.0d0 )*SIN(3*psid)
          yxx2 = -( 1.0d0/(4.0d0*qcxicl) )*SIN(2.0d0*psid)
          yxx3 = ( 1.0d0/qcxicl**2 )*( 1.0d0 - 1.25d0*qcxicl**2 )
     $           *SIN(psid)
          ycf = -( 1.0d0/qcxicl**3 )
          yy1 = ( 1.0d0 - 1.5d0*qcxicl**2 )
          yy2 = ( ABS(qcxicl**2 - 1.0d0) )**1.5d0
          yxx4 = ycf*( yy1*psid - yy2*lcap )

          ycap = yxx1 + yxx2 + yxx3 + yxx4
        ENDIF

        mubar2 = EXP(xmu*ycap)

      ELSE

c       Error.

        WRITE(nttyo,1270)
        WRITE(noutpt,1270)

 1270   FORMAT(/6x,'ERROR (VISCOS): YCAP could not be calculated',
     $  ' because',/9x,'xicl was out of range. Consequently, the',
     $  ' critical ehnancement',/9x,'term will be set to unity.')

      ENDIF

  800 CONTINUE

      mubar = mubar0*mubar1*mubar2
      VISCOS = mubar*muref

      IF (qwr3) THEN

c       Write mubar0, mubar1, mubar2, and mubar on
c       the standard output.

        WRITE(ux24a,'(f16.9)') mubar0
        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(ux24b,'(f16.9)') mubar1
        CALL LEJUST(ux24b)
        j3 = ILNOBL(ux24b)

        WRITE(ux24c,'(f16.9)') mubar2
        CALL LEJUST(ux24c)
        j4 = ILNOBL(ux24c)

        WRITE(ux24d,'(f16.9)') mubar
        CALL LEJUST(ux24d)
        j5 = ILNOBL(ux24d)

        WRITE(nttyo,1280) ux24a(1:j2), ux24b(1:j3),
     $  ux24c(1:j4), ux24d(1:j5)

 1280   FORMAT(/6x,'VISCOS: mubar0 = ',a,/14x,'mubar1 = ',a,
     $  /14x,'mubar2 = ',a,/14x,'mubar  = ',a)

      ENDIF

      END