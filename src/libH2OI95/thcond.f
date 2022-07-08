      REAL(8) FUNCTION THCOND(delta, tau, pdx, ztx, tempk, press,
     $ rho, rhocr, rcnstw, tcr, pcr, ptripl, cpx, cvx, tmltx,
     $ mubar, nttyo, noutpt)

c     This function calculates the thermal conductivity of water
c     using the "IAPWS Formulation 2011 for the Thermal Conductivity
c     of Ordinary Water Substance" (IAPWS R15-11). The corresponding
c     journal reference is: Huber M.L., Perkins R.A., Friend D.G.,
c     Sengers J.V., Assael M.J., Metaxa I.N., Miyagawa K., Hellmann R.,
c     and Vogel E. (2012) New International Formulation for the
c     Thermal Conductivity of H2O. J. Phys. Chem. Ref. Data 41,
c     033102, 1-23.  Equations referenced below are from
c     IAPWS R15-11.

c     Note that the dimensionless viscosity (mubar) is required
c     as an input. FUNCTION VISCOS must be called first to provide
c     this.

c     This function is called by:

c       H2OI95.f

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      REAL(8) delta, tau, pdx, ztx, tempk, press, rho, rhocr, tcr, pcr,
     $ ptripl, rcnstw, cpx, cvx, tmltx, mubar

      INTEGER nttyo, noutpt

      INTEGER ILNOBL

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER i, j, k, j2, j3, j4, j5

      INTEGER ier

      CHARACTER(LEN=24) ux24a, ux24b, ux24c, ux24d, ux24e

      REAL(8) laref, sx0, sx1, sx2, dltm1, tbar, taum1, xx, zdx

      REAL(8) delchb, xicl, yx, lcap, zcap, cpxbar, kappa, kapem1

      REAL(8) cflk(0:4), cflij(0:4,0:5)

      REAL(8) lacap, qdem1, nuc, gamma, xicl0, gam0, tbarr

      REAL(8) qd, exx, xxx, zxx, zxx1, zxx2, zxx3, pi

      REAL(8) labar, labar0, labar1, labar2

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
c     verification of thermal conductivity calculations.

      LOGICAL qnola2, qcbzet, qzetap, qwr1, qwr2, qwr3

      DATA qnola2, qcbzet, qzetap, qwr1, qwr2, qwr3 / .FALSE., .FALSE.,
     $ .FALSE., .FALSE., .FALSE., .FALSE. /

c-----------------------------------------------------------------------

c     Here laref is the reference thermal conductivity (W/K-m).

      DATA laref /1.00d-03 /

c     Coefficients (L(k)) for calculating labar1.

      DATA cflk(0:4) / 2.443221d-03, 1.323095d-02, 6.770357d-03,
     $ -3.454586d-03, 4.096266d-04 /

c     Coefficients (L(ij)) for calculating labar2.

      DATA cflij(0:4,0) / 1.60397357d0, 2.33771842d0, 2.19650529d0,
     $ -1.21051378d0, -2.7203370d0 /

      DATA cflij(0:4,1) / -0.646013523d0, -2.78843778d0, -4.54580785d0,
     $ 1.60812989d0, 4.57586331d0 /

      DATA cflij(0:4,2) / 0.111443906d0, 1.53616167d0, 3.55777244d0,
     $ -0.621178141d0, -3.18369245d0 /

      DATA cflij(0:4,3) / 0.102997357d0, -0.463045512d0, -1.40944978d0,
     $ 0.0716373224d0, 1.1168348d0 /

      DATA cflij(0:4,4) / -0.0504123634d0, 0.0832827019d0,
     $ 0.275418278d0, 0.0d0, -0.19268305d0 /

      DATA cflij(0:4,5) / 0.00609859258d0, -0.00719201245d0,
     $ -0.0205938816d0, 0.0d0, 0.012913842d0 /

c     Coefficients for calculating the critical enhancement term
c     (labar2).

      DATA lacap, qdem1, nuc, gamma, xicl0, gam0, tbarr /
     $ 177.8514d0, 0.40d0, 0.630d0, 1.239d0, 0.13d0, 0.06d0, 1.5d0 /

c-----------------------------------------------------------------------

c     Check range of validity.

      qrchk = .FALSE.

      IF (press.GT.0.0d0 .AND. press.LT.ptripl) THEN

        IF (tempk.GE.273.16d0 .AND. tempk.LE.1173.15d0) THEN
          qrchk = .TRUE.
        ENDIF

      ELSEIF (tmltx .GT. 0.0d0) THEN

        IF (press.GE.ptripl .AND. press.LE.100.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.1173.15d0) THEN
            qrchk = .TRUE.
          ENDIF

        ELSEIF (press.GT.100.0d0 .AND. press.LE.250.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.874.0d0) THEN
            qrchk = .TRUE.
          ENDIF

        ELSEIF (press.GT.250.0d0 .AND. press.LE.687.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.573.0d0) THEN
            qrchk = .TRUE.
          ENDIF

        ELSEIF (press.GT.687.0d0 .AND. press.LE.785.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.403.0d0) THEN
            qrchk = .TRUE.
          ENDIF

        ELSEIF (press.GT.785.0d0 .AND. press.LE.1000.0d0) THEN

          IF (tempk.GE.tmltx .AND. tempk.LE.348.0d0) THEN
            qrchk = .TRUE.
          ENDIF

        ENDIF

      ENDIF

      IF (.NOT.qrchk) THEN

c       Failed check on valid temperature-pressure range.
c       Issue an appropriate warning.

        IF (tmltx .LE. 0) THEN

          WRITE(noutpt,1020)
          WRITE(nttyo,1020)

 1020     FORMAT(/6x,'WARNING (THCOND): Could not calculate the',
     $    /9x,'melting temperature. Hence could not check the range',
     $    /9x,'of validity of the thermal conductivity model.')

        ELSEIF (.NOT.qrchk) THEN

          WRITE(noutpt,1050)
          WRITE(nttyo,1050)

 1050     FORMAT(/6x,'WARNING (THCOND): Outside the range of validity',
     $    /9x,'of the thermal conductivity model')

        ENDIF

      ENDIF

c     Note: Although this thermal conductivity model must be used with
c     the IAPWS-95 equation-of-state model, the language employed is
c     slightly different. "Tbar" (tbar) is equivalent to 1/tau. Also,
c     "rhobar" is equivalent to delta. Insofar as possible, we will
c     use tau and delta in evaluating the appropriate equations.
c     "Tbar" is only necessary in one instance, in the calculation
c     of labar0.

c     Here the scaled thermal conductivity (labar) depends on three
c     parts, labar0, labar1, and labar2. Note that
c     labar = labar0*labar1 + labar2. This is not analogous to
c     the representation of the viscosity, as the critical
c     enhancement is a term, not a factor.

      tbar = 1.0d0/tau

c     Calculate the viscosity in the dilute gas limit (labar0).

      labar0 = 1.0d0

      sx0 = cflk(0) + cflk(1)*tau
      DO k = 2,4
        sx0 = sx0 + cflk(k)*(tau**k)
      ENDDO

      labar0 = SQRT(tbar)/sx0

c     Calculate the contribution to viscosity due to finite density
c     (labar1).

      labar1 = 1.0d0

      dltm1 = delta - 1.0d0
      taum1 = tau - 1.0d0

      sx1 = 0.0d0
      DO i = 0,4

        sx2 = cflij(i,0) + cflij(i,1)*dltm1
        DO j = 2,5
          sx2 = sx2 + cflij(i,j)*(dltm1**j)
        ENDDO
        sx1 = sx1 + (taum1**i)*sx2
      ENDDO
      xx = delta*sx1
      labar1 = EXP(xx)

c     Calculate the critical enhancement (labar2), which is only
c     signficant in a very small region of density and temperature
c     around the critical point.

      labar2 = 0.0d0

c     If not using labar2, go to 800.

      IF (qnola2) GO TO 800

c     Get zeta at the reference temperature (1.5*Tcr, or 970.644K).

      zettbr = 0.0d0
      ztxr = 0.0d0

      IF (qcbzet .OR. qzetap) THEN

c       Estimate zeta at the reference temperature using fitted
c       polynomials.

        zettbr = GZTTBR(delta, ier)

        IF (ier .GT. 0) THEN

          WRITE(nttyo,1210)

 1210     FORMAT(/6x,'WARNING (THCOND): Could not estimate the',
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

 1220     FORMAT(/6x,'THCOND: zettbr =  ',a,' (Poly)',/
     $    14x,'  ztxr =  ',a,' (EOS)',/14x,'  diff = ',a)

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

1250    FORMAT(/6x,'THCOND: ztxr   = ',a,1x,a)

      ENDIF

c     Get delchb (DELTA chibar, equation 23).

      delchb = delta*( ztx - (ztxr*tbarr*tau) )
      IF (delchb .LT. 0.0d0) delchb = 0.0d0

c     Get xicl (equation 22).

      xxx = delchb/gam0
      exx = nuc/gamma
      xicl = xicl0*( xxx**exx )

      IF (qwr2) THEN

c       Write xicl on the standard output.

        WRITE(ux24a,'(f16.9)') xicl
        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(nttyo,1260) ux24a(1:j2)

 1260   FORMAT(/6x,'THCOND: xicl = ',a)

      ENDIF

c     Get yx (equation 20).

      qd = 1.0d0/qdem1
      yx = qd*xicl

c     Get zcap (Z in equation 19).

      IF (yx .LT. 1.2d-07) THEN
        zcap = 0.0d0
        labar2 = 0.0d0
      ELSE

        kappa = cpx/cvx
        kapem1 = 1.0d0/kappa
        pi = 4*ATAN(1.0d0)
        cpxbar = cpx/rcnstw

        zxx = 2.0d0/(pi*yx)
        zxx1 = (1.0d0 - kapem1)*ATAN(yx) + kapem1*yx
        zxx3 = yx**(-1) + ( yx**2/(3*delta**2) )
        zxx2 = 1.0d0 - EXP( -1.0d0/zxx3 )
        zcap = zxx*(zxx1 - zxx2)

        labar2 = lacap*delta*cpxbar*tbar*zcap/mubar

      ENDIF

  800 CONTINUE

      labar = labar0*labar1 + labar2

      THCOND = labar*laref

      IF (qwr3) THEN

c       Write labar0, labar1, labar2, and labar on
c       the standard output.

        WRITE(ux24a,'(f16.9)') labar0
        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(ux24b,'(f16.9)') labar1
        CALL LEJUST(ux24b)
        j3 = ILNOBL(ux24b)

        IF (labar2.LT.1.0d+06) THEN
          WRITE(ux24c,'(f16.9)') labar2
        ELSE
          WRITE(ux24c,'(1pe16.9)') labar2
        ENDIF
        CALL LEJUST(ux24c)
        j4 = ILNOBL(ux24c)

        IF (labar.LT.1.0d+06) THEN
          WRITE(ux24d,'(f16.9)') labar
        ELSE
          WRITE(ux24d,'(1pe16.9)') labar
        ENDIF
        CALL LEJUST(ux24d)
        j5 = ILNOBL(ux24d)

        WRITE(nttyo,1280) ux24a(1:j2), ux24b(1:j3),
     $  ux24c(1:j4), ux24d(1:j5)

 1280   FORMAT(/6x,'THCOND: labar0 = ',a,/14x,'labar1 = ',a,
     $  /14x,'labar2 = ',a,/14x,'labar  = ',a)

      ENDIF

      END