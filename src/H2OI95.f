      program H2OI95

c     Version 1.1, Build 118

c     2020/05/18

c     This version information must be incorporated into string
c     string variables named uverst and ubuild, which are used to
c     write this information on the output files.

c     Version 1.0 calculates the thermodynamic properties of water
c     using the equation of state model of Wagner and Pruss (2002).
c     It calculates the standard thermochemical properties from
c     that model in conjunction with the final CODATA
c     recommendations (Cox et al., 1989).

c     Version 1.1 incorporates routines to calculate additional
c     properties of water, including viscosity, thermal conductivity,
c     surface tension, dielectric constant, Born functions, and
c     Debye-Huckel constants.

c     This is a stand-alone program written as a research code.
c     It is not optimized to support codes such as SUPCRT92
c     or reactive transport codes. In particular, it has
c     a lot of I/O that would be unncessary in that context.

      IMPLICIT NONE

      LOGICAL qdelta, qerr, qfail, qheadr, qrho, qpress, qprnt1,
     $ qprnt2, qprnt3, qsilent, qpsatc, qtsatc, qwrphi, qrhog

      LOGICAL qhdr1t, qhdr1p, q1phsc

      INTEGER nttyi, nttyo, ninpt, noutpt, nmcsvf, nccsvf, nxcsvf

      INTEGER ILNOBL

      INTEGER i, ier, icount, iter, j, j2, j3, j4, j5

      CHARACTER(LEN=80) uline, ulinex, ux80

      CHARACTER(LEN=24) udescr

      CHARACTER(LEN=24) ux24a, ux24b, ux24c

      CHARACTER(LEN=16) uverst, ubuild

      CHARACTER(LEN=16) ushdes

      CHARACTER(LEN=8) ux8, unam8

      REAL(8) delta, tau, rho, tempk, press, psat, psatt, tsat,
     $ tsatp, datum1, datum2, datum3, rhog

      REAL(8) rhosv, rhosl, deltsv, deltsl

      REAL(8) px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $ ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx, vsx, thcx,
     $ stx, epsx, deptp, deppt, mubar, agme, agm10, aphi,
     $ adhv, adhh, ahrt, bgm, bdhv, bdhh, bhrt, pmltx, tmltx

      REAL(8) pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv, dtxv,
     $ bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv, ktxv, ztxv,
     $ vsxv, thcxv, stxv, epsxv, deptpv, depptv, mubarv, agmev,
     $ agm10v, aphiv, adhvv, adhhv, ahrtv, bgmv, bdhvv, bdhhv,
     $ bhrtv, pmltxv, tmltxv

      REAL(8) pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl, dtxl,
     $ bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl, ktxl, ztxl,
     $ vsxl, thcxl, stxl, epsxl, deptpl, depptl, mubarl, agmel,
     $ agm10l, aphil, adhvl, adhhl, ahrtl, bgml, bdhvl, bdhhl,
     $ bhrtl, pmltxl, tmltxl

      REAL(8) uxcu, sxcu, hxcu, cvxcu, cpxcu, dtxcu, axcu, gxcu, vxcu

      REAL(8) uxcuv, sxcuv, hxcuv, cvxcuv, cpxcuv, dtxcuv, axcuv,
     $ gxcuv, vxcuv

      REAL(8) uxcul, sxcul, hxcul, cvxcul, cpxcul, dtxcul, axcul,
     $ gxcul, vxcul

      REAL(8) tcr, rhocr, pcr, rcnstw, xmcapw, betamx

      REAL(8) ttripl, ts1atm, ptripl, p1atm, bettol, btxtol,
     $ rhotol, rhosvx, rhoslx, rxx, psctol

      REAL(8) VISCOS, THCOND, SURFTE, PMELT, TKMELT

c     Values for triple point, normal boiling point, and critical
c     point parameters and the specific gas constant are from
c     Wagner and Pruss (2002), p. 395.

      DATA ttripl / 273.16d0 /, ts1atm / 373.124d0 /
      DATA ptripl / 0.000611655d0 /, p1atm / 1.01325d-1 /

      DATA bettol / 1.0d-9 /

      DATA rhotol / 1.0d-8 /

      DATA psctol / 1.0d-5 /

c     Note: the value of rcnstw (kJ kg^-1 K^-1) is specific to the
c     IAPWS-95 model and should not be replaced by any other value
c     in evaluating the model. The constant mW (xmcapw) is the mass
c     of one mole of water (kg mol^-1). The product of these two
c     constants is the gas constant R (kJ mol^-1 K^-1). Note
c     that 1 kJ = 1 m^3 kPa.

      DATA tcr / 647.096d+00 /, rhocr / 322.0d+00 /,
     $ pcr / 22.064d0 /, rcnstw / 0.46151805d+00 /,
     $ xmcapw / 0.018015268d0 /

      DATA uverst / '1.1' /
      DATA ubuild / '118' /

      DATA unam8 / 'H2O' /

      DATA nttyi /5/, nttyo /6/

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Open the main output file.

      CALL OPENOU(noutpt,nttyo,'output','formatted',80,noutpt)

c     Open the input file.

      CALL OPENIN(noutpt,nttyo,'input','formatted',ninpt)

c     Open the main .csv output file.

      CALL OPENOU(noutpt,nttyo,'mtab.csv','formatted',400,nmcsvf)

c     Open the thermochemical .csv output file.

      CALL OPENOU(noutpt,nttyo,'ctab.csv','formatted',400,nccsvf)

c     Open the extra properties .csv output file.

      CALL OPENOU(noutpt,nttyo,'xtab.csv','formatted',400,nxcsvf)

c     Write the code name and version information on the output files.

      j2 = ILNOBL(uverst)
      j3 = ILNOBL(ubuild)
      WRITE (nttyo,1000) uverst(1:j2),ubuild(1:j3)
      WRITE (noutpt,1000) uverst(1:j2),ubuild(1:j3)

 1000 FORMAT(1x,'H2OI95, Version ',a,', Build ',a,/)

      WRITE (nmcsvf,1005) uverst(1:j2),ubuild(1:j3)
      WRITE (nccsvf,1005) uverst(1:j2),ubuild(1:j3)
      WRITE (nxcsvf,1005) uverst(1:j2),ubuild(1:j3)

 1005 FORMAT('H2OI95,Version ',a,',Build ',a,/)

c     Write header for the main .csv file (mtab.csv).

      WRITE (nmcsvf,1010)

 1010 FORMAT('Substance,Type,Temp(K),Press(MPa),rho(kg/m3),u(kJ/kg),',
     $ 'h(kJ/kg),s(kJ/kg/K),a(kJ/kg),g(kJ/kg),v(m3/kg),',
     $ 'cv(kJ/kg/K),cp(kJ/kg/K),w(m/s),mu(K/MPa),dt(kJ/kg/MPa,',
     $ 'bs(K/MPa),iter,betamx,btxtol,rhog(kg/m3)',/)

c     Write header for the thermochemical .csv file (ctab.csv).

      WRITE (nccsvf,1020)

 1020 FORMAT('Substance,Type,Temp(K),Press(MPa),rho(kg/m3),u(kJ/mol),',
     $ 'h(kJ/mol),s(J/mol/K),a(kJ/mol),g(kJ/mol),v(m3/mol),',
     $ 'cv(J/mol/K),cp(J/mol/K),w(m/s),mu(K/MPa),dt(kJ/mol/MPa,',
     $ 'bs(K/MPa)',/)

c     Write header for the extra properties csv file (xtab.csv).

      WRITE (nxcsvf,1025)

 1025 FORMAT('Substance,Type,Temp(K),Press(MPa),rho(kg/m3),av(/K),',
     $ 'kt(/MPa),zt(),vs(Pa-s),thc(W/K-m),st(mN/m),eps(),',
     $ 'deptp((/K),deppt((/MPa),agme((kg/mol)^0.5),',
     $ 'agm10((kg/mol)^0.5),aphi((kg/mol)^0.5),',
     $ 'AV(cm^3-kg^0.5-mol^-1.5),',
     $ 'AH/RT((kg/mol)^0.5),bgm((kg/mol)^0.5 A^-1),',
     $ 'BV(cm^3 kg^0.5 mol^-1.5 A^-1),',
     $ 'BH/RT(kg^0.5 mol^-0.5 A^-1),','Pmelt(MPa),','Tmelt(K)',/)

      qheadr = .FALSE.
      qrho = .FALSE.
      qdelta = .FALSE.
      qpress = .FALSE.
      qrhog = .FALSE.
      qprnt1 = .FALSE.
      qprnt2 = .FALSE.
      qprnt3 = .FALSE.
      qpsatc = .FALSE.
      qtsatc = .FALSE.
      qwrphi = .FALSE.
      qhdr1t = .FALSE.
      qhdr1p = .FALSE.

      icount = 0
      iter = 0
      betamx = 0.0d0

c     Below is a return point for reading a new line from the
c     input file.

  120 READ(ninpt,'(a80)',END=300) uline
      icount = icount + 1
      qerr = .FALSE.
      qfail = .FALSE.
      j2 = ILNOBL(uline)

c     Here bettol is the normal convergence tolerance.
c     Under certain conditions, it may be necessary to
c     loosen the tolerance actually applied, which is
c     btxtol.

      btxtol = bettol

c     Make a copy (ulinex) of the current line. This copy will
c     be examined and parsed as needed.

      ulinex = uline

c     Convert the copy to lower case.

      CALL LOCASE(ulinex)

c     Left-justify the copy.

      CALL LEJUST(ulinex)
      j3 = ILNOBL(ulinex)
      j3 = MIN(j3,72)

      IF (.NOT.qheadr) THEN

c       The header line has not been previously encountered.
c       If the current line is not a header line, then write
c       it to the output file now.

        i = INDEX(ulinex,'tempk')
        j = INDEX(ulinex,'press')
        IF (i.NE.1 .AND. j.NE.1) THEN

          WRITE (noutpt,1030) uline(1:j2)

 1030     FORMAT(1x,a)

        ENDIF
      ENDIF

c     Ignore empty lines.

      IF (j3 .LE. 0) GO TO 120

c     Ignore comment lines.

      IF (ulinex(1:1) .EQ. '#') GO TO 120
      IF (ulinex(1:1) .EQ. '*') GO TO 120

      IF (.NOT. qheadr) THEN

c       Process print option strings.

        IF (ulinex(1:16) .EQ. 'showphi') THEN
          qwrphi = .TRUE.
          GO TO 120
        ENDIF

        IF (ulinex(1:16) .EQ. 'showdetails1') THEN
          qprnt1 = .TRUE.
          GO TO 120
        ENDIF

        IF (ulinex(1:16) .EQ. 'showdetails2') THEN
          qprnt1 = .TRUE.
          qprnt2 = .TRUE.
          GO TO 120
        ENDIF

        IF (ulinex(1:16) .EQ. 'showdetails3') THEN
          qprnt1 = .TRUE.
          qprnt2 = .TRUE.
          qprnt3 = .TRUE.
          GO TO 120
        ENDIF

c       Check if the current line is the data header.
c       This must be one of the following:

c         'tempk   rho'
c         'tempk   delta'
c         'tempk   press   rhog'
c         'tempk   psat'
c         'press   tsat'

c       where rhog denotes optional starting values
c       of density (rho) for Newton-Raphson iteration.
c       Here press, psat, rho, and rhog may have units
c       appended (e.g., press(MPa, rho(kg/m3)).

        IF (ulinex(1:5).EQ.'tempk') THEN
          qhdr1t = .TRUE.
        ELSEIF (ulinex(1:5).EQ.'press') THEN
          qhdr1p = .TRUE.
        ELSE
          GO TO 120
        ENDIF

c       The current line is the data header.

        ux80 = uline
        CALL LEJUST(ux80)
        j5 = ILNOBL(ux80)

        WRITE (nttyo,1035) ux80(1:j5)
        WRITE (noutpt,1035) ux80(1:j5)

 1035   FORMAT(1x,'The header line contains "',a,'"',/)

c       Remove optional units.

        j4 = INDEX(ulinex,'press(mpa)')
        IF (j4 .GT. 0) THEN
          i = j4 + 5
          j = i + 4
          ulinex(i:j) = ' '
        ENDIF

        j4 = INDEX(ulinex,'psat(mpa)')
        IF (j4 .GT. 0) THEN
          i = j4 + 4
          j = i + 4
          ulinex(i:j) = ' '
        ENDIF

        j4 = INDEX(ulinex,'rho(kg/m3)')
        IF (j4 .GT. 0) THEN
          i = j4 + 3
          j = i + 6
          ulinex(i:j) = ' '
        ENDIF

        j4 = INDEX(ulinex,'rhog(kg/m3)')
        IF (j4 .GT. 0) THEN
          i = j4 + 4
          j = i + 6
          ulinex(i:j) = ' '
        ENDIF

        IF (qhdr1t) THEN

          ulinex(1:5) = ' '
          CALL LEJUST(ulinex)

          IF (ulinex(1:3) .EQ. 'rho') THEN
            qrho = .TRUE.
            qheadr = .TRUE.
            GO TO 120
          ELSEIF (ulinex(1:5) .EQ. 'delta') THEN
            qdelta = .TRUE.
            qheadr = .TRUE.
            GO TO 120
          ELSEIF (ulinex(1:5) .EQ. 'press') THEN
            qpress = .TRUE.
            qheadr = .TRUE.

            ulinex(1:5) = ' '
            CALL LEJUST(ulinex)
            IF (ulinex(1:4) .EQ. 'rhog') qrhog = .TRUE.

            GO TO 120
          ELSEIF (ulinex(1:4) .EQ. 'psat') THEN
            qpsatc = .TRUE.
            qheadr = .TRUE.
            GO TO 120
          ELSE
            CALL LEJUST(uline)
            j2 = ILNOBL(uline)

            WRITE (nttyo,1040) uline(1:j2)
            WRITE (noutpt,1040) uline(1:j2)

 1040       FORMAT(/6x,'ERROR (H2OI95): The data header on the input',
     $      ' file is not valid.',/9x,'It must be of the general form',
     $      ' "tempk   rho", "tempk   delta",',
     $      /9x,'"tempk   press   (rhog)", "tempk   psat",',
     $      ' or "press   tsat".',/9x,'Units may be included.',
     $      ' The last line read was',/9x,'"',a,'".')

            STOP
          ENDIF

        ELSEIF (qhdr1p) THEN

          ulinex(1:5) = ' '
          CALL LEJUST(ulinex)

          IF (ulinex(1:4) .EQ. 'tsat') THEN
            qtsatc = .TRUE.
            qheadr = .TRUE.
            GO TO 120
          ELSE
            CALL LEJUST(uline)
            j2 = ILNOBL(uline)
            WRITE (nttyo,1040) uline(1:j2)
            WRITE (noutpt,1040) uline(1:j2)
            STOP
          ENDIF

        ENDIF

      ENDIF

c     At this point should have found the data header.

      ux24a = ' '
      ux24b = ' '
      ux24c = ' '

c     Get the first datum.

      CALL LEJUST(uline)
      j2 = INDEX(uline,' ')
      ux24a(1:j2) = uline(1:j2)
      READ(ux24a,'(d16.9)') datum1

c     Get the second datum.

      uline(1:j2) = ' '
      CALL LEJUST(uline)
      j2 = INDEX(uline,' ')
      ux24b(1:j2) = uline(1:j2)
      READ(ux24b,'(d16.9)') datum2

c     Get the third datum, if any.

      uline(1:j2) = ' '
      CALL LEJUST(uline)
      j2 = INDEX(uline,' ')
      IF (j2 .GT. 0) ux24c(1:j2) = uline(1:j2)
      READ(ux24c,'(d16.9)') datum3

      IF (qhdr1t) THEN
        tempk = datum1

        IF (tempk.LT.273.16d0 .OR. tempk.GT.1273.00d0) THEN

          WRITE (nttyo,1050) tempk
          WRITE (noutpt,1050) tempk

 1050     FORMAT(/6x,'WARNING (H2OI95): Temp = ',f9.4,' K. This is',
     $    ' outside',/9x,'the nominal range of the IAPWS-95 model.',
     $    ' Continuing.')

        ENDIF

      ELSEIF (qhdr1p) THEN
        press = datum1

        IF (press.LT.1.0d-12 .OR. press.GT.1000.0d0) THEN

          WRITE (nttyo,1060) press
          WRITE (noutpt,1060) press

 1060     FORMAT(/6x,'WARNING (H2OI95): Press = ',1pe12.5,' MPa.',
     $    ' This is outside',/9x,'the nominal range of the IAPWS-95',
     $    ' model. Continuing.')

        ENDIF

      ENDIF

c     Initialize q1phsc. This will be set to .TRUE. if a point not
c     specified as on the saturation curve happens to be on or very
c     close to the saturation curve. The surface tension will then
c     be calculated and reported.

      q1phsc = .FALSE.

c     Initialize the fluid type descriptor.

      udescr = 'unknown'

      IF (qhdr1t) THEN

c       The first datum is temperature. There are a number of possible
c       options for this case.

        tau = tcr/tempk

        IF (qdelta) THEN
          delta = datum2
          rho = rhocr*delta
        ELSEIF (qrho) THEN
          rho = datum2
          delta = rho/rhocr
        ENDIF

        IF (qrho .OR. qdelta) THEN

c         Calculate properties including pressure as a function
c         of specified temperature and delta. This requires only
c         a direct evaluation of the equation of state.

          CALL CALDLT(nttyo, noutpt, qerr, qprnt1, qprnt2, qprnt3,
     $    qwrphi, pcr, rcnstw, rhocr, tcr, iter, tempk, tau, delta,
     $    rho, press, px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $    ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx)

          press = px

c         Find a description of the H2O fluid type.

c         Check for unstable fluid.

          IF (pdx .LT. 0.0d0) THEN

            WRITE (nttyo,1070)
            WRITE (noutpt,1070)

 1070       FORMAT(/6x,'WARNING (H2OI95): PRESSURE DERIVATIVE WITH',
     $      /9x,'RESPECT TO DELTA IS LESS THAN ZERO')

            udescr = 'unstable'
            GO TO 210
          ENDIF

          IF (px .LE. 0.0d0) THEN

            WRITE (nttyo,1100)
            WRITE (noutpt,1100)

 1100       FORMAT(/6x,'WARNING (H2OI95): PRESSURE IS NOT POSITIVE')

            udescr = 'unrealistic'
            GO TO 210
          ENDIF

          CALL FDESCR(nttyo, noutpt, qprnt1, qprnt2, qprnt3,
     $    qwrphi, pcr, rcnstw, rhocr, tcr, bettol, btxtol, psat,
     $    rho, rhosv, rhosl, rhotol, tempk, tau, press, udescr)

c         Is the current point on or very close to the saturation
c         curve?

          IF (psat .GT. 0.0d0) THEN
            IF (tempk.GE.ttripl .AND. tempk.LE.tcr) THEN
              IF (ABS(press - psat) .LE. psctol) q1phsc = .TRUE.
            ENDIF
          ENDIF

c         If in the P-T quadrant where the saturation curve exists,
c         the description found by FDESCR pertains to the stable
c         form at the specified temperature and pressure. Check
c         for metastability.

          IF (tempk.LE.tcr .AND. press.LE.pcr) THEN

            rhosvx = rhosv + rhotol
            rhoslx = rhosl - rhotol

            IF (rho.GE.rhosvx .AND. rho.LE.rhoslx) THEN

c             Can be sure at this point that the fluid is metastable.
c             Discrimination of metastable vapor from metastable
c             liquid here is somewhat arbitrary.

              rxx = 0.5d0*(rhosv + rhosl)
              IF (rho .LE. rxx) THEN
                udescr = 'metastable vapor'
              ELSE
                udescr = 'metastable liquid'
              ENDIF

            ENDIF

          ENDIF

  210     CONTINUE

        ELSEIF (qpress) THEN

c         Calculate properties as a function of specified
c         temperature and pressure. This requires finding the
c         appropriate delta value by iteration.

          press = datum2
          rhog = 0.0d0
          if (qrhog) rhog = datum3

          IF (press.LT.1.0d-12 .OR. press.GT.1000.0d0) THEN
            WRITE (nttyo,1060) press
            WRITE (noutpt,1060) press
          ENDIF

          CALL CALPRE(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $    qprnt3, qwrphi, qrhog, rhog, pcr, rcnstw, rhocr, tcr, iter,
     $    betamx, bettol, rhotol, btxtol, udescr, press, tempk, tau,
     $    delta, rho, psat, px, ux, sx, hx, cvx, cpx, wx, mux, dtx,
     $    bsx, ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx)

          IF (psat .GT. 0.0d0) THEN
            IF (tempk.GE.ttripl .AND. tempk.LE.tcr) THEN
              IF (ABS(press - psat) .LE. psctol) q1phsc = .TRUE.
            ENDIF
          ENDIF

        ELSEIF (qpsatc) THEN

c         Calculate the saturation curve properties for the
c         specified temperature. In this option, datum2 is not used.
c         Note: here we allow the temperature slightly out of
c         range to accommodate 0C.

          IF (tempk.LT.273.15d0 .OR. tempk.GT.tcr) THEN

c           Error, specified temperature is not in the allowed range.

            WRITE (nttyo,1210) tempk
            WRITE (noutpt,1210) tempk

 1210       FORMAT(/6x,'ERROR (H2OI95): Temp = ',f9.4,' K. This is',
     $      /9x,'out of range for the saturation curve.')

            qerr = .TRUE.

          ELSE

c           Calculate the saturation curve properties.

            qsilent = .FALSE.

            CALL CALSCT(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $      qprnt3, qsilent, qwrphi, pcr, rcnstw, rhocr, tcr, iter,
     $      betamx, bettol, btxtol, tempk, tau, press, rhosv, rhosl,
     $      deltsv, deltsl, pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv,
     $      dtxv, bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv,
     $      ktxv, ztxv, pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl,
     $      dtxl, bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl,
     $      ktxl, ztxl)

          ENDIF

        ELSE

c         Have encountered an unrecognized input.

          WRITE (nttyo,1250)
          WRITE (noutpt,1250)

 1250     FORMAT(/6x,'ERROR (H2OI95): Do not recognize how to process',
     $    /9x,'the current input line.')

          GO TO 280

        ENDIF

      ELSEIF (qhdr1p) THEN

        IF (qtsatc) THEN

c         Calculate the saturation curve properties for the
c         specified pressure. In this option, datum2 is not used.

          IF (press.LT.ptripl .OR. press.GT.pcr) THEN

c           Error, specified pressure is not in the allowed range.

            WRITE (nttyo,1260) press
            WRITE (noutpt,1260) press

 1260       FORMAT(/6x,'ERROR (H2OI95): Press = ',1pe12.5,' MPa.',
     $      /9x,'this is out of range for the saturation curve.')

            qerr = .TRUE.

          ELSE

c           Calculate the saturation curve properties for the
c           specified temperature. In this option, datum2 is not used.

            CALL CALSCP(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $      qprnt3, qwrphi, pcr, rcnstw, rhocr, tcr, press, ptripl,
     $      p1atm, ttripl, ts1atm, iter, betamx, bettol, btxtol, tempk,
     $      tau, rhosv, rhosl, deltsv, deltsl, pxv, uxv, sxv, hxv,
     $      cvxv, cpxv, wxv, muxv, dtxv, bsxv, axv, gxv, vxv, pdxv,
     $      ptxv, adxv, gdxv, avxv, ktxv, ztxv, pxl, uxl, sxl, hxl,
     $      cvxl, cpxl, wxl, muxl, dtxl, bsxl, axl, gxl, vxl, pdxl,
     $      ptxl, adxl, gdxl, avxl, ktxl, ztxl)

          ENDIF

        ELSE

c         Have encountered an unrecognized input.

          WRITE (nttyo,1250)
          WRITE (noutpt,1250)
          GO TO 280

        ENDIF

      ENDIF

c     Create short descriptions of fluid type to write on the
c     .csv output files.

      ushdes = udescr(1:16)
      j4 = ILNOBL(udescr)
      IF (udescr(1:j4) .EQ. 'compressed liquid') THEN
        ushdes = 'comp. liquid'
      ELSEIF (udescr(1:j4) .EQ. 'supercritical fluid') THEN
        ushdes = 'supercritical'
      ENDIF

      IF (qerr .OR. qfail) THEN

        WRITE (nttyo,1270)
        WRITE (noutpt,1270)

 1270   FORMAT(/6x,'WARNING (H2OI95): CALCULATION FAILED. NO',
     $  /9x,'RESULTS WERE OBTAINED.')

        j4 = ILNOBL(ushdes)
        j5 = ILNOBL(unam8)

        IF (qhdr1t) THEN

          IF (qrho .OR. qdelta) THEN

            WRITE (nmcsvf,1272) unam8(1:j5), ushdes(1:j4), tempk, rho
            WRITE (nccsvf,1272) unam8(1:j5), ushdes(1:j4), tempk, rho
            WRITE (nxcsvf,1272) unam8(1:j5), ushdes(1:j4), tempk, rho

 1272       FORMAT(a,',',a,',',f11.6,',FAILED,',1pe16.9,',')

          ELSEIF (qpress) THEN

            WRITE (nmcsvf,1273) unam8(1:j5), ushdes(1:j4), tempk,
     $      press, rhog

 1273       FORMAT(a,',',a,',',f11.6,',',1pe16.9,
     $      ',FAILED,,,,,,,,,,,,,,,,',e16.9,',')

            WRITE (nccsvf,1274) unam8(1:j5), ushdes(1:j4), tempk,
     $      press
            WRITE (nxcsvf,1274) unam8(1:j5), ushdes(1:j4), tempk,
     $      press

 1274       FORMAT(a,',',a,',',f11.6,',',1pe16.9,
     $      ',FAILED,')

          ELSEIF (qpsatc) THEN

            WRITE (nmcsvf,1275) unam8(1:j5), tempk
            WRITE (nccsvf,1275) unam8(1:j5), tempk
            WRITE (nxcsvf,1275) unam8(1:j5), tempk

 1275       FORMAT(a,',vapor,',f11.6,',FAILED,')

            WRITE (nmcsvf,1276) unam8(1:j5), tempk
            WRITE (nccsvf,1276) unam8(1:j5), tempk
            WRITE (nxcsvf,1276) unam8(1:j5), tempk

 1276       FORMAT(a,',liquid,',f11.6,',FAILED,')

          ENDIF

        ELSEIF (qhdr1p) THEN

          IF (qtsatc) THEN

            WRITE (nmcsvf,1277) unam8(1:j5), press
            WRITE (nccsvf,1277) unam8(1:j5), press
            WRITE (nxcsvf,1277) unam8(1:j5), press

 1277       FORMAT(a,',vapor,FAILED,',1pe16.9,',')

            WRITE (nmcsvf,1278) unam8(1:j5), press
            WRITE (nccsvf,1278) unam8(1:j5), press
            WRITE (nxcsvf,1278) unam8(1:j5), press

 1278       FORMAT(a,',liquid,FAILED,',1pe16.9,',')

          ENDIF

        ENDIF
        GO TO 280

      ENDIF

      IF (.NOT.qpsatc .AND. .NOT.qtsatc) THEN

c       Have a single phase.

c       Calculate thermodynamic functions on the standard
c       thermochemical scale where appropriate. Convert units
c       to those commonly given in steam tables.

        CALL CTHERM(nttyo, noutpt, unam8, tempk, press,
     $  px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $  ax, gx, vx, uxcu, sxcu, hxcu, cvxcu, cpxcu,
     $  dtxcu, axcu, gxcu, vxcu)

c       Calculate additional properties.

c       Melting pressure and melting temperature.

        pmltx = PMELT(tempk, press, noutpt, nttyo, ier)
        tmltx = TKMELT(press, noutpt, nttyo, ier)

c       Viscosity.

        vsx = VISCOS(delta, tau, pdx, ztx, tempk, press,
     $  rho, rhocr, rcnstw, tcr, pcr, ptripl, tmltx, mubar,
     $  nttyo, noutpt)

c       Thermal conductivity.

        thcx = THCOND(delta, tau, pdx, ztx, tempk, press,
     $  rho, rhocr, rcnstw, tcr, pcr, ptripl, cpx, cvx, tmltx,
     $  mubar, nttyo, noutpt)

c       Surface tension.

        stx = 0.0d0
        IF (q1phsc) THEN
          stx = SURFTE(tempk, tcr, q1phsc, qpsatc, qtsatc,
     $    nttyo)
        ENDIF

c       Dielectric constant.

        CALL DIELEC(rho, rhocr, tempk, tcr, rcnstw, xmcapw,
     $  avx, ktx, ptx, nttyo, epsx, deptp, deppt, agme,
     $  agm10, aphi, adhv, adhh, ahrt, bgm, bdhv, bdhh,
     $  bhrt)

c       Write results to the output file, the main .csv file,
c       and the thermochemical .csv file.

        CALL WROTAB(nttyo, noutpt, nmcsvf, nccsvf, nxcsvf, unam8,
     $  udescr, ushdes, tempk, delta, rho, px, ux, hx, sx, ax,
     $  gx, vx, cvx, cpx, wx, mux, dtx, bsx, pdx, avx, ktx,
     $  ztx, vsx, thcx, stx, epsx, deptp, deppt, agme, agm10,
     $  aphi, adhv, adhh, ahrt, bgm, bdhv, bdhh, bhrt, pmltx,
     $  tmltx, uxcu, hxcu, sxcu, axcu, gxcu, vxcu, cvxcu,
     $  cpxcu, dtxcu, iter, betamx, btxtol, qrhog, rhog)

      ELSE

c       Have vapor and liquid phases.

c       First process the vapor.

        udescr = 'vapor'
        ushdes = udescr(1:16)

c       Calculate thermodynamic functions on the standard
c       thermochemical scale where appropriate. Convert units
c       to those commonly given in steam tables.

c       Calling sequence substitutions:

c         pxv for px
c         uxv for ux
c         hxv for hx
c         sxv for sx
c         axv for ax
c         gxv for gx
c         vxv for vx
c         cvxv for cvx
c         cpxv for cpx
c         wxv for wx
c         muxv for mux
c         dtxv for dtx
c         bsxv for bsx
c         uxcuv for uxcu
c         sxcuv for sxcu
c         hxcuv for hxcu
c         cvxcuv for cvxcu
c         cpxcuv for cpxcu
c         dtxcuv for dtxcu
c         axcuv for axcu
c         gxcuv for gxcu
c         vxcuv for vxcu

c       Calculate thermodynamic functions on the standard
c       thermochemical scale where appropriate. Convert units
c       to those commonly given in steam tables.

        CALL CTHERM(nttyo, noutpt, unam8, tempk, press,
     $  pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv, dtxv, bsxv,
     $  axv, gxv, vxv, uxcuv, sxcuv, hxcuv, cvxcuv, cpxcuv,
     $  dtxcuv, axcuv, gxcuv, vxcuv)

c       Calculate additional properties.

c       Melting pressure and melting temperature mean nothing to
c       vapor. In the future, one might calculate the sublimation
c       pressure and sublimation temperature here.

        pmltxv = 0.0d0
        tmltxv = 0.0d0

c       Viscosity.

c       Calling sequence substitutions:

c         ztxv for ztx
c         rhosv for rho
c         deltsv for delta
c         tmltxv for tmltx
c         mubarv for mubar

        vsxv = VISCOS(deltsv, tau, pdx, ztxv, tempk, press,
     $  rhosv, rhocr, rcnstw, tcr, pcr, ptripl, tmltxv, mubarv,
     $  nttyo, noutpt)

c       Thermal conductivity.

c       Calling sequence substitutions:

c         ztxv for ztx
c         rhosv for rho
c         deltsv for delta
c         cpxv for cpx
c         cvxv for cvx
c         tmltxv for tmltx
c         mubarv for mubar

        thcxv = THCOND(deltsv, tau, pdx, ztxv, tempk, press,
     $  rhosv, rhocr, rcnstw, tcr, pcr, ptripl, cpxv, cvxv, tmltxv,
     $  mubarv, nttyo, noutpt)

c       Surface tension. The value to be calculated will be
c       associated with the liquid (technically it is associated
c       with both).

        stxv = 0.0d0

c       Dielectric constant.

c       Calling sequence substitutions:

c         rhosv for rho
c         avxv for avx
c         ktxv for ktx
c         ptxv for ptx
c         epsxv for epsx
c         deptpv for deptp
c         depptv for deppt
c         agmev for agme
c         agm10v for agm10
c         aphiv for aphi
c         adhvv for adhv
c         adhhv for adhh
c         ahrtv for ahrt
c         bgmv for bgm
c         bdhvv for bdhv
c         bdhhv for bdhh
c         bhrtv for bhrt

        CALL DIELEC(rhosv, rhocr, tempk, tcr, rcnstw, xmcapw,
     $  avxv, ktxv, ptxv, nttyo, epsxv, deptpv, depptv, agmev,
     $  agm10v, aphiv, adhvv, adhhv, ahrtv, bgmv, bdhvv, bdhhv,
     $  bhrtv)

c       Write results to the output file, the main .csv file,
c       and the thermochemical .csv file.

c       Calling sequence substitutions:

c         deltsv for delta
c         rhosv for rho
c         pxv for px
c         uxv for ux
c         hxv for hx
c         sxv for sx
c         axv for ax
c         gxv for gx
c         vxv for vx
c         cvxv for cvx
c         cpxv for cpx
c         wxv for wx
c         muxv for mux
c         dtxv for dtx
c         bsxv for bsx
c         pdxv for pdx
c         avxv for avx
c         ktxv for ktx
c         ztxv for ztx
c         vsxv for vsx
c         thcxv for thcx
c         stxv for stx
c         epsxv for epsx
c         deptpv for deptp
c         depptv for deppt
c         agmev for agme
c         agm10v for agm10
c         aphiv for aphi
c         adhvv for adhv
c         adhhv for adhh
c         ahrtv for ahrt
c         bgmv for bgm
c         bdhvv for bdhv
c         bdhhv for bdhh
c         bhrtv for bhrt
c         pmltxv for pmltx
c         tmltxv for tmltx
c         uxcuv for uxcu
c         sxcuv for sxcu
c         hxcuv for hxcu
c         cvxcuv for cvxcu
c         cpxcuv for cpxcu
c         dtxcuv for dtxcu
c         axcuv for axcu
c         gxcuv for gxcu
c         vxcuv for vxcu

        CALL WROTAB(nttyo, noutpt, nmcsvf, nccsvf, nxcsvf, unam8,
     $  udescr, ushdes, tempk, deltsv, rhosv, pxv, uxv, hxv, sxv, axv,
     $  gxv, vxv, cvxv, cpxv, wxv, muxv, dtxv, bsxv, pdxv, avxv, ktxv,
     $  ztxv, vsxv, thcxv, stxv, epsxv, deptpv, depptv, agmev, agm10v,
     $  aphiv, adhvv, adhhv, ahrtv, bgmv, bdhvv, bdhhv, bhrtv, pmltxv,
     $  tmltxv, uxcuv, hxcuv, sxcuv, axcuv, gxcuv, vxcuv, cvxcuv,
     $  cpxcuv, dtxcuv, iter, betamx, btxtol, qrhog, rhog)

c       Now process the liquid.

        udescr = 'liquid'
        ushdes = udescr(1:16)

c       Calculate thermodynamic functions on the standard
c       thermochemical scale where appropriate. Convert units
c       to those commonly given in steam tables.

c       Calling sequence substitutions:

c         pxl for px
c         uxl for ux
c         hxl for hx
c         sxl for sx
c         axl for ax
c         gxl for gx
c         vxl for vx
c         cvxl for cvx
c         cpxl for cpx
c         wxl for wx
c         muxl for mux
c         dtxl for dtx
c         bsxl for bsx
c         uxcul for uxcu
c         sxcul for sxcu
c         hxcul for hxcu
c         cvxcul for cvxcu
c         cpxcul for cpxcu
c         dtxcul for dtxcu
c         axcul for axcu
c         gxcul for gxcu
c         vxcul for vxcu
c         bettol for btxtol

c       Calculate thermodynamic functions on the standard
c       thermochemical scale where appropriate. Convert units
c       to those commonly given in steam tables.

        CALL CTHERM(nttyo, noutpt, unam8, tempk, press,
     $  pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl, dtxl, bsxl,
     $  axl, gxl, vxl, uxcul, sxcul, hxcul, cvxcul, cpxcul,
     $  dtxcul, axcul, gxcul, vxcul)

c       Write results to the output file, the main .csv file,
c       and the thermochemical .csv file.

c       Calculate additional properties.

c       Melting pressure and melting temperature.

        pmltxl = PMELT(tempk, press, noutpt, nttyo, ier)
        tmltxl = TKMELT(press, noutpt, nttyo, ier)

c       Viscosity.

c       Calling sequence substitutions:

c         ztxl for ztx
c         rhosl for rho
c         deltsl for delta
c         tmltxl for tmltx
c         mubarl for mubar

        vsxl = VISCOS(deltsl, tau, pdx, ztxl, tempk, press,
     $  rhosl, rhocr, rcnstw, tcr, pcr, ptripl, tmltxl, mubarl,
     $  nttyo, noutpt)

c       Thermal conductivity.

c       Calling sequence substitutions:

c         ztxl for ztx
c         rhosl for rho
c         deltsl for delta
c         cpxl for cpx
c         cvxl for cvx
c         tmltxl for tmltx
c         mubarl for mubar

        thcxl = THCOND(deltsl, tau, pdx, ztxl, tempk, press,
     $  rhosl, rhocr, rcnstw, tcr, pcr, ptripl, cpxl, cvxl, tmltxl,
     $  mubarl, nttyo, noutpt)

c       Surface tension.

        stxl = SURFTE(tempk, tcr, q1phsc, qpsatc, qtsatc,
     $  nttyo)

c       Dielectric constant.

c       Calling sequence substitutions:

c         rhosl for rho
c         avxl for avx
c         ktxl for ktx
c         ptxl for ptx
c         epsxl for epsx
c         deptpl for deptp
c         depptl for deppt
c         agmel for agme
c         agm10l for agm10
c         aphil for aphi
c         adhvl for adhv
c         adhhl for adhh
c         ahrtl for ahrt
c         bgml for bgm
c         bdhvl for bdhv
c         bdhhl for bdhh
c         bhrtl for bhrt

        CALL DIELEC(rhosl, rhocr, tempk, tcr, rcnstw, xmcapw,
     $  avxl, ktxl, ptxl, nttyo, epsxl, deptpl, depptl, agmel,
     $  agm10l, aphil, adhvl, adhhl, ahrtl, bgml, bdhvl, bdhhl,
     $  bhrtl)

c       Write results to the output file, the main .csv file,
c       and the thermochemical .csv file.

c       Calling sequence substitutions:

c         deltsl for delta
c         rhosl for rho
c         pxl for px
c         uxl for ux
c         hxl for hx
c         sxl for sx
c         axl for ax
c         gxl for gx
c         vxl for vx
c         cvxl for cvx
c         cpxl for cpx
c         wxl for wx
c         muxl for mux
c         dtxl for dtx
c         bsxl for bsx
c         pdxl for pdx
c         avxl for avx
c         ktxl for ktx
c         ztxl for ztx
c         vsxl for vsx
c         thcxl for thcx
c         stxl for stx
c         epsxl for epsx
c         deptpl for deptp
c         depptl for deppt
c         agmel for agme
c         agm10l for agm10
c         aphil for aphi
c         adhvl for adhv
c         adhhl for adhh
c         ahrtl for ahrt
c         bgml for bgm
c         bdhvl for bdhv
c         bdhhl for bdhh
c         bhrtl for bhrt
c         pmltxl for pmltx
c         tmltxl for tmltx
c         uxcul for uxcu
c         sxcul for sxcu
c         hxcul for hxcu
c         cvxcul for cvxcu
c         cpxcul for cpxcu
c         dtxcul for dtxcu
c         axcul for axcu
c         gxcul for gxcu
c         vxcul for vxcu

        CALL WROTAB(nttyo, noutpt, nmcsvf, nccsvf, nxcsvf, unam8,
     $  udescr, ushdes, tempk, deltsl, rhosl, pxl, uxl, hxl, sxl, axl,
     $  gxl, vxl, cvxl, cpxl, wxl, muxl, dtxl, bsxl, pdxl, avxl, ktxl,
     $  ztxl, vsxl, thcxl, stxl, epsxl, deptpl, depptl, agmel, agm10l,
     $  aphil, adhvl, adhhl, ahrtl, bgml, bdhvl, bdhhl, bhrtl, pmltxl,
     $  tmltxl, uxcul, hxcul, sxcul, axcul, gxcul, vxcul, cvxcul,
     $  cpxcul, dtxcul, iter, betamx, bettol, qrhog, rhog)

      ENDIF

  280 CONTINUE

c     Write a break line.

      WRITE (nttyo,1280)
      WRITE (noutpt,1280)

 1280 FORMAT(/1x,'----------------------------------------')

c     Go back to process another line from the input file.

      GO TO 120

  300 CONTINUE

c     Have now read all lines from the input file.

      IF (.NOT.qheadr) THEN

        WRITE(nttyo,1320)
        WRITE(noutpt,1320)

 1320   FORMAT(/6x,'ERROR (H2OI95): Did not find the data header',
     $  ' on the input file.',/9x,'It must be "tempk   rho",',
     $  ' "tempk   delta", "tempk   press",',/9x,
     $  '"tempk   psat", or "press   tsat".')

      ENDIF

      WRITE(ux8,'(i8)') icount
      CALL LEJUST(ux8)
      j2 = ILNOBL(ux8)

      WRITE(nttyo,1330) ux8(1:j2)
      WRITE(noutpt,1330) ux8(1:j2)

 1330 FORMAT(/3x,'Have read ',a,' lines from the input file.',
     $ /3x,'Done.',/)

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE CTHERM(nttyo, noutpt, unam8, tempk, press,
     $ px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $ ax, gx, vx, uxcu, sxcu, hxcu, cvxcu, cpxcu,
     $ dtxcu, axcu, gxcu, vxcu)

c     Convert entropy, internal energy, enthalpy, Helmholtz energy,
c     and Gibbs energy values from the equation of state scale to the
c     standard thermochemical scale. Convert these and other equation
c     of state results from kilogram units to moles. Results remain
c     in Joule, meter, and MPa units.

      IMPLICIT NONE

c     Calling sequence variables.

      INTEGER nttyo, noutpt

      CHARACTER(LEN=8) unam8

      REAL(8) tempk, press

      REAL(8) rho, px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $ ax, gx, vx

      REAL(8) rhocu, uxcu, sxcu, hxcu, cvxcu, cpxcu, dtxcu, axcu,
     $ gxcu, vxcu

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      REAL(8) acoks, acoku, acokh, acoka, acokg

c     The following data are in kilogram units.

      DATA acoks / 3.5156150d0 /, acoku / -15767.19391d0 /,
     $ acokh / -15970.89538d0 /, acoka / -11906.84446d0 /,
     $ acokg / -12110.54592d0 /

      REAL(8) xx

      REAL(8) htripl, mwH2O

c     Note: the htripl value was computed using this software.

      DATA htripl / 0.000611782d0 /

      DATA mwH2O / 18.01528d0 /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate the entropy, internal energy, enthalpy, Helmholtz
c     energy, and Gibbs energy on the standard thermochemical scale
c     using offset parameters. Results are still in kilogram units.

      sxcu = sx + acoks
      uxcu = ux + acoku
      hxcu = (hx - htripl) + acokh
      axcu = ux - tempk*sxcu + acoka
      gxcu = hx - tempk*sxcu + acokg

c     Calculate the entropy, internal energy, enthalpy, Helmholtz
c     energy, Gibbs energy, volume, heat capacity at constant volume,
c     heat capacity at constant pressure, and isothermal throttling
c     coefficient (dt) in molar units.

      xx = mwH2O/1000.0d0

      sxcu = xx*sxcu
      uxcu = xx*uxcu
      hxcu = xx*hxcu
      axcu = xx*axcu
      gxcu = xx*gxcu

      vxcu = xx*vx
      cvxcu = xx*cvx
      cpxcu = xx*cpx
      dtxcu = xx*dtx

c     Adjust entropy and heat capacity values from kJ units to J.

      xx = 1000.0d0

      sxcu = xx*sxcu
      cvxcu = xx*cvx
      cpxcu = xx*cpx

c     The following have no units changes at this point: wx, mux, bsx


c     At this point, all thermodynamic functions are in units commonly
c     tabulated in modern steam tables. Pressure is still in MPa and
c     length is still in meters.

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE CALPRE(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qwrphi, qrhog, rhog, pcr, rcnstw, rhocr, tcr, iter,
     $ betamx, bettol, rhotol, btxtol, udescr, press, tempk, tau,
     $ delta, rho, psat, px, ux, sx, hx, cvx, cpx, wx, mux, dtx,
     $ bsx, ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx)

c     Find the thermodynamic properties of water at given temperature
c     and pressure. The problem reduces to finding the value of reduced
c     density (delta) that is consistent with the desired pressure.
c     The Newton-Raphson method is employed. Under-relaxation techniques
c     are required to assure that delta or other requisite parameters
c     do not take on out-of-bounds values in the iteration process.
c     Small negative values of calculated pressure are okay. Zero or
c     negative values for calculated "pdx" (pressure derivative with
c     respect to delta) imply the unstable zone and must be avoided.

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qerr, qfail, qprnt1, qprnt2, qprnt3, qwrphi, qrhog

      CHARACTER(LEN=24) udescr

      INTEGER nttyo, noutpt

      INTEGER iter

      REAL(8) press, tempk, tau, delta, rho, rhog, psat

      REAL(8) pcr, rcnstw, rhocr, tcr

      REAL(8) betamx, bettol, btxtol, rhotol

      REAL(8) px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx,
     $ ax, gx, vx, pdx, ptx, adx, gdx, avx, ktx, ztx

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      LOGICAL qsilent

      INTEGER ILNOBL

      INTEGER icdlta, itermx, irelax, irlxmx, j2

      CHARACTER(LEN=8) ux8

      REAL(8) rhoidg, rhosup, rhocpa, rhosv, rhosl, rxx, psatt

      REAL(8) arelax, delta0, deltas, delts0, deltol, dltamx, dltsx,
     $ pdiff, px0, rho0, sbetmx

      DATA qsilent / .TRUE. /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      itermx = 60
      IF (qrhog) itermx = 80
      btxtol = bettol
      deltol = 1.0d-12
      icdlta = 0
      qerr = .FALSE.
      qfail = .FALSE.

      WRITE (nttyo,1050) tempk, tau, press
      WRITE (noutpt,1050) tempk, tau, press

 1050 FORMAT(/3x,'CALPRE: Temp(K) = ',f9.4,13x,'tau    = ',e16.9,
     $ /11x,'press   = ',e16.9,' MPa',/)

      IF (qrhog) THEN

        WRITE (nttyo,1060) rhog
        WRITE (noutpt,1060) rhog

 1060   FORMAT(11x,'rhog   = ',e16.9,' kg/m3 (starting value)',/)

      ENDIF

      iter = 0
      rho = 0.0d0

c     Obtain a description (udescr) of the H2O fluid.
c     Ignore the qerr error flag.

      CALL FDESCR(nttyo, noutpt, qprnt1, qprnt2, qprnt3,
     $ qwrphi, pcr, rcnstw, rhocr, tcr, bettol, btxtol, psat,
     $ rho, rhosv, rhosl, rhotol, tempk, tau, press, udescr)

      IF (udescr(1:8) .EQ. 'vapor   ') THEN

c       Vapor: assume ideal gas behavior.

c       Ideal gas.
        rho = 1000.0d0*press/(tempk*rcnstw)

c       Alternate: ideal gas correction to saturation
c       density (slightly different).
c       rho = (press/psat)*rhosv

      ELSEIF (udescr(1:8) .EQ. 'liquid  ') THEN

c       Liquid: use a liquid-like density.

c       The liquid density on the saturation curve.
        rho = rhosl

c       Something slightly greater than the saturation
c       density.
c       rho = 1.0002d0*rhosl

      ELSEIF (udescr(1:8) .EQ. 'compress') THEN

c       Estimate the density of compressed liquid.

c       Ideal gas correction to critical point density.
c       This produces one of two false solutions at 500K
c       and 25 MPa. Do not use this.
c       rho = (press/pcr)*(tcr/tempk)*rhocr

c       Twice the ideal gas correction to the critical
c       point density.
c       rho = 2.0d0*(press/pcr)*(tcr/tempk)*rhocr

c       The saturated liquid density is a minimum value.

        rho = 1.10d0*rhosl

c       Close to the upper limit for this field
c       (T near the triple point, 1000 MPa).
c       For higher pressure, a higher value might
c       be needed.
c       rho = 1250.0d0

        rho = DMAX1(rho, rhosl)
        rho = DMIN1(rho, 1400.0d0)

      ELSEIF (udescr(1:8).EQ.'supercri') THEN

c       Estimate the density of supercritical fluid.

c       Ideal gas.
c       rho = 1000.0d0*press/(tempk*rcnstw)

c       SUPCRT92 estimate, about 15% higher than ideal gas.
c       rho = 2500.0d0*press/tempk

c       Ideal gas correction to critical point density.
c       rho = (press/pcr)*(tcr/tempk)*rhocr

c       Twice the ideal gas correction to the critical
c       point density.
        rho = 2.0d0*(press/pcr)*(tcr/tempk)*rhocr

c       Close to the upper limit for this P, T field.
c       (T near the critical point, 1000 MPa).
c       Calculated value of 1068.7 kg/m3 is rounded up.
c       For higher pressure, a higher value might
c       be needed.
c       rho = 1100.0d0

        rho = DMIN1(rho, 1100.0d0)

      ELSEIF (udescr(1:8).EQ.'hot vapo') THEN

c       Estimate the density of hot vapor.

c       Ideal gas.
        rhoidg = 1000.0d0*press/(tempk*rcnstw)

c       SUPCRT92 estimate, about 15% higher than ideal gas.
        rhosup = 2500.0d0*press/tempk

c       Ideal gas correction to critical point density.
        rhocpa = (press/pcr)*(tcr/tempk)*rhocr

c       The upper limit for this field, the critical pressure
c       (rhocr), 22.064 MPa.
c       rho = rhocr

        IF (press .LE. 1.00d0) THEN
          rho = rhoidg
        ELSEIF (press .LE. 18.0d0) THEN
          rho = rhosup
        ELSE
          rho = rhocpa
        ENDIF

        rho = DMIN1(rho, rhocr)

      ELSE

c       The H2O fluid type could not be determined.

        WRITE (nttyo,1070)
        WRITE (noutpt,1070)

 1070   FORMAT(/6x,'WARNING (CALPRE): The H2O fluid type could not be',
     $  ' determined.',/9x,'A good starting estimate of density',
     $  ' could not be established.')

c       Will try four times the critical density.

        rho = 4*rhocr

      ENDIF

      IF (qrhog) THEN

c       Use the user-specified starting value for rho.

        rho = rhog
      ENDIF

      delta = rho/rhocr

      irelax = 0
      irlxmx = 20
      arelax = 0.25d0
      sbetmx = 0.0d0
      betamx = 0.0d0
      j2 = ILNOBL(udescr)

      WRITE (nttyo,1080) udescr(1:j2)
      WRITE (noutpt,1080) udescr(1:j2)

 1080 FORMAT(6x,'This appears to be ',a,'.',/)

c     Below is the return point for iterating on delta to
c     obtain a desired pressure.

  250 CONTINUE

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, delta, rho, px, ax, ux, sx, hx,
     $ cvx, cpx, wx, mux, dtx, bsx, gx, vx, pdx, ptx, adx,
     $ gdx, avx, ktx, ztx, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilent)

c     Test to see if the pressure has converged to the desired value.

      pdiff = px - press
      sbetmx = pdiff/press
      betamx = DABS(sbetmx)

      WRITE (noutpt,1250) iter, px, sbetmx
      WRITE (nttyo,1250) iter, px, sbetmx

 1250 FORMAT(6x,'CALPRE: iter= ',i3,', px= ',1pe16.9,', sbetmx= ',
     $ e12.5)

      IF (qprnt1) THEN

        WRITE (noutpt,1260) delta, rho
        WRITE (nttyo,1260) delta, rho

 1260   FORMAT(/9x,'delta  = ',e16.9,10x,'rho    = ',e16.9)

        WRITE (noutpt,1290) px, pdx
        WRITE (nttyo,1290) px, pdx

 1290   FORMAT(9x,'px     = ',e16.9,' MPa',6x,'pdx    = ',e16.9)

      ENDIF

      IF (betamx .LE. btxtol) THEN

c       Have converged.

        IF (qprnt1) THEN

          WRITE (nttyo,1670) betamx, btxtol
          WRITE (noutpt,1670) betamx, btxtol

 1670     FORMAT(/3x,'CALPRE: ITERATION CONVERGED.',
     $    /3x,'MAX NORM betamx =',1pe10.3,', TOLERANCE btxtol = ',
     $    e10.3,/)

        ENDIF

        GO TO 420
      ENDIF

      IF (iter .GE. itermx) THEN

c       Have done the maximum number of iterations.

        WRITE (nttyo,1370)
        WRITE (noutpt,1370)

 1370   FORMAT(/3x,'CALPRE: HAVE DONE THE MAXIMUM NUMBER OF',
     $  ' ITERATIONS.')

        GO TO 410
      ENDIF

c     Estimate a new value for delta (and rho).

      iter = iter + 1
      irelax = 0
      delta0 = delta
      rho0 = rho
      px0 = px
      deltas = -pdiff/pdx
      delts0 = deltas
      dltamx = DABS(deltas)

      IF (qprnt1) THEN

        WRITE (noutpt,1372) dltamx
        WRITE (nttyo,1372) dltamx

 1372   FORMAT(9x,'dltamx = ',e16.9,/)

      ENDIF

      IF (iter .GE. 5) THEN
        IF (dltamx .LE. deltol) icdlta = icdlta + 1
      ENDIF

      IF (icdlta .GE. 5) THEN

        IF (betamx .LE. 1.0d-6) THEN

c         Have pseudo-converged.

          WRITE (nttyo,1375) betamx, btxtol, dltamx, deltol
          WRITE (noutpt,1375) betamx, btxtol, dltamx, deltol

 1375     FORMAT(/3x,'CALPRE: ITERATION PSEUDO-CONVERGED.',
     $    /6x,'MAX NORM betamx = ',1pe10.3,', TOLERANCE btxtol = ',
     $    e10.3,/6x,'MAX NORM dltamx = ',e10.3,
     $    ', TOLERANCE deltol = ',1pe10.3)

          GO TO 420

        ENDIF

      ENDIF

      IF (icdlta .GE. 10) THEN

c       Iteration is failing to result in any improvement.

        WRITE (nttyo,1377) betamx, btxtol, dltamx, deltol
        WRITE (noutpt,1377) betamx, btxtol, dltamx, deltol

 1377   FORMAT(/3x,'CALPRE: ITERATION IS NOT LEADING TO IMPROVEMENT.',
     $  /6x,'MAX NORM betamx = ',1pe10.3,', TOLERANCE btxtol = ',
     $  e10.3,/6x,'MAX NORM dltamx = ',e10.3,
     $  ', TOLERANCE deltol = ',1pe10.3)

        GO TO 410

      ENDIF

      dltsx = 0.10d0*delta0
      IF (deltas .GT. dltsx) THEN

c       Clamp deltas to keep delta from large fractional increase.

        IF (qprnt1) THEN

          WRITE (nttyo,1277)
          WRITE (noutpt,1277)

 1277     FORMAT(/3x,'CLAMPING CHANGE IN DELTA TO AVOID LARGE',
     $    ' FRACTIONAL INCREASE.')

        ENDIF

        deltas = dltsx
      ENDIF

      dltsx = -0.40d0*delta0
      IF (deltas .LT. dltsx) THEN

c       Clamp deltas to keep delta from large fractional decrease.

        IF (qprnt1) THEN

          WRITE (nttyo,1278)
          WRITE (noutpt,1278)

 1278     FORMAT(/3x,'CLAMPING CHANGE IN DELTA TO AVOID LARGE',
     $    ' FRACTIONAL DECREASE.')

        ENDIF

        deltas = dltsx
      ENDIF

c     Newton-Raphson correction.

  220 delta = delta0 + deltas

      IF (delta .LT. 1.0d-15) then

c       Under-relax to keep delta from being too small. Only
c       positive values are physical.

        IF (irelax .GE. irlxmx) GO TO 405
          IF (qprnt1) THEN

          WRITE (nttyo,1330)
          WRITE (noutpt,1330)

 1330     FORMAT(/3x,'UNDER-RELAXING TO AVOID NEGATIVE DELTA.')

        ENDIF

        irelax = irelax + 1
        deltas = arelax*deltas
        GO TO 220
      ENDIF

      rho = rhocr*delta

      GO TO 250

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  405 CONTINUE

c     Have hit the limit for under-relaxation steps.

      WRITE (ux8,'(i8)') irlxmx
      j2 = ILNOBL(ux8)
      CALL LEJUST(ux8)

      WRITE (nttyo,1350) ux8(1:j2)
      WRITE (noutpt,1350) ux8(1:j2)

 1350 FORMAT(/3x,'CALPRE: HAVE HIT THE UNDER-RELAXATION LIMIT OF ',a,
     $ /3x,' TIMES PER ITERATION.')

c     Continue iteration.

      rho = rhocr*delta

      GO TO 250

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  410 CONTINUE

c     Iteration failed.

      qfail = .TRUE.

      WRITE (nttyo,1690)
      WRITE (noutpt,1690)

 1690 FORMAT(/3x,'CALPRE: ITERATION FAILED.',/)

      GO TO 999

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  420 CONTINUE

c     Calculations are complete for the last line read
c     from the input file. Write the results.

      WRITE (nttyo,2400)
      WRITE (noutpt,2400)

 2400 FORMAT(/8x,'Temp(K)',10x,'press(MPa)')

      WRITE (nttyo,2410) tempk,px
      WRITE (noutpt,2410) tempk,px

 2410 FORMAT(6x,f9.4,7x,e16.9)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE CALSCT(nttyo, noutpt, qerr, qfail, qprnt1, qprnt2,
     $ qprnt3, qsilent, qwrphi, pcr, rcnstw, rhocr, tcr, iter,
     $ betamx, bettol, btxtol, tempk, tau, press, rhosv, rhosl,
     $ deltsv, deltsl, pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv,
     $ dtxv, bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv,
     $ ktxv, ztxv, pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl,
     $ dtxl, bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl,
     $ ktxl, ztxl)

c     This routine calculates the saturation properties as a function
c     of specified temperature. This is done using Newton-Raphson
c     iteration to refine values of pressure, vapor density, and
c     liquid density, starting with results obtained using approximate
c     equations included by Wagner and Pruss (2002) in their description
c     of the IAPWS-95 model.

      IMPLICIT NONE

c     CALLING SEQUENCE VARIABLES.

      LOGICAL qerr, qfail, qprnt1, qprnt2, qprnt3, qsilent, qwrphi

      INTEGER nttyo, noutpt

      INTEGER iter

      REAL(8) pcr, rcnstw, rhocr, tempk, tau, tcr

      REAL(8) betamx, bettol, btxtol

      REAL(8) press, rhosv, rhosl, deltsv, deltsl

      REAL(8) pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv, dtxv,
     $ bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv, ktxv, ztxv

      REAL(8) pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl, dtxl,
     $ bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl, ktxl, ztxl

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     LOCAL VARIABLES.

      INTEGER k_par

      PARAMETER (k_par = 3)

      LOGICAL qrestl, qxiter, qsilnx

      INTEGER ILNOBL

      INTEGER i, itermx, irelax, j, j2, kdim, kmax, icutv,
     $ icutl, icutvq, icutlq, ncut

      CHARACTER(LEN=8) ux8, unam8

      REAL(8) delta, rho, rhosv0, rhosl0, psat, psat0, psatt

      REAL(8) delta0, dltsvq, dltslq, dltsv0, dltsl0, pxm, px0, pdiff

      REAL(8) aamatr(k_par,k_par), alpha(k_par), deltas(k_par),
     $ beta(k_par), betmx0

      REAL(8) xx, xxx

      REAL(8) arelax, dix, dltx

      REAL(8) bettl1, bettl2, bettl3

      DATA qsilnx / .TRUE. /

      DATA bettl1 / 1.0d-8 /
      DATA bettl2 / 1.0d-7 /
      DATA bettl3 / 1.0d-6 /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (.NOT.qsilent) THEN

        WRITE (nttyo,1010) tempk, tau
        WRITE (noutpt,1010) tempk, tau

 1010   FORMAT(/3x,'CALSCT: Temp(K) = ',f9.4,13x,'tau    = ',e16.9,/)

      ENDIF

      kmax = k_par
      kdim = kmax

      qerr = .FALSE.
      qfail = .FALSE.
      qrestl = .FALSE.
      qxiter = .FALSE.
      itermx = 35
      btxtol = bettol

      IF (tempk.LE.298.15d0) THEN
        qrestl = .TRUE.
        btxtol = bettl1
        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1020) tempk, bettl1
          WRITE (noutpt,1020) tempk, bettl1

 1020     FORMAT(6x,'Temp(K) = ',f10.5,', LIES BETWEEN THE TRIPLE',
     $    ' POINT and 298.15K.',/6x,'APPLYING LOOSER CONVERGENCE',
     $    ' TOLERANCE bettl1 = ',1pe10.3,'.',/)

        ENDIF
      ENDIF

      IF (tempk.GT.647.090d0 .AND. tempk.LT.tcr) THEN
        IF (.NOT. qsilent) THEN

          WRITE (nttyo,1040)
          WRITE (noutpt,1040)

 1040     FORMAT(6x,'AS PRESENTLY TUNED, SATURATION CURVE ITERATION',
     $    /6x,'WILL NOT CONVERGE NORMALLY FOR TEMPERATURES GREATER',
     $    /6x,'THAN 647.090 K AND LESS THAN THE CRITICAL TEMPERATURE',
     $    /6x,'OF 647.096 K. A MAXIMUM OF FIVE ITERATIONS WILL BE',
     $    ' DONE.',/)

        ENDIF
        qxiter = .TRUE.
      ENDIF

      IF (tempk.GT.647.082d0) THEN
        qrestl = .TRUE.
        btxtol = bettl2
        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1050) tempk, bettl2
          WRITE (noutpt,1050) tempk, bettl2

 1050     FORMAT(6x,'Temp(K) = ',f10.5,', VERY CLOSE TO THE CRITICAL',
     $    ' TEMPERATURE.',/6x,'APPLYING LOOSER CONVERGENCE TOLERANCE',
     $    ' bettl1 = ',1pe10.3,'.',/)

        ENDIF
      ENDIF

      arelax = 1.0d0

      rhosv = 0.0d0
      rhosl = 0.0d0
      press = 0.0d0

c     Calculate approximate saturation pressure and
c     corresponding densities of liquid and vapor.
c     These results are not those of the IAPWS-95 model
c     itself, but can serve as starting estimates.

      CALL APXSCT(noutpt, nttyo, qprnt1, qsilent, rhocr, rhosv,
     $ rhosl, tempk, tcr, pcr, psat, psatt, qerr)

      IF (qerr) THEN

c       Error, specified temperature is not in the allowed range.

        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1210) tempk
          WRITE (noutpt,1210) tempk

 1210     FORMAT(/6x,'ERROR (CALSCT): Temp = ',f9.4,'K. This is out',
     $    /9x,'of range for the saturation curve.')

        ENDIF

        qfail = .TRUE.
        GO TO 999
      ENDIF

      icutv = 0
      icutl = 0
      icutvq = 0
      icutlq = 0

      deltsv = rhosv/rhocr
      deltsl = rhosl/rhocr

c     Save the values from the approximation.

      dltsvq = deltsv
      dltslq = deltsl

  120 CONTINUE

      IF (qprnt1 .AND. .NOT.qsilent) THEN

        WRITE (nttyo,1240) tempk, psat, tau, deltsv, rhosv,
     $  deltsl, rhosl
        WRITE (noutpt,1240) tempk, psat, tau, deltsv, rhosv,
     $  deltsl, rhosl

 1240   FORMAT(/6x,'CALSCT: Starting values for Temp(K) = ',f9.4,
     $  //9x,'psat   = ',e16.9,' MPa',6x,'tau    = ',e16.9,
     $  /9x,'deltsv = ',e16.9,10x,'rhosv  = ',e16.9,
     $  /9x,'deltsl = ',e16.9,10x,'rhosl  = ',e16.9,/)

      ENDIF

      press = psat
      rhosv0 = rhosv
      rhosl0 = rhosl

      iter = 0
      irelax = 0
      betmx0 = 1.0d+100

      psat0 = psat
      dltsv0 = deltsv
      dltsl0 = deltsl

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Below is the return point to refine the saturation
c     curve properties.

  250 CONTINUE

c     First calculate the vapor properties by calling EVAI95
c     with the vapor delta.

      unam8 = 'Vapor'
      j2 = ILNOBL(unam8)

      IF (qprnt2) THEN

        WRITE (noutpt,1250) unam8(1:j2), deltsv, rhosv

 1250   FORMAT(6x,a,/9x,'delta  = ',e16.9,10x,'rho   = ',e16.9,/)

      ENDIF

c     Calling sequence substitutions:

c       deltsv for delta
c       rhosv sfor rho
c       pxv for px
c       axv for ax
c       uxv for ux
c       sxv for sx
c       hxv for hx
c       cvxv for cvx
c       cpxv for cpx
c       wxv for wx
c       muxv for mux
c       dtxv for dtx
c       bsxv for bsx
c       gxv for gx
c       vxv for vx
c       pdxv for pdx
c       ptxv for ptx
c       adxv for adx
c       gdxv for gdx
c       avxv for avx
c       ktxv for ktx
c       ztxv for ztx
c       qsilnx for qsilent

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, deltsv, rhosv, pxv, axv, uxv, sxv, hxv,
     $ cvxv, cpxv, wxv, muxv, dtxv, bsxv, gxv, vxv, pdxv, ptxv, adxv,
     $ gdxv, avxv, ktxv, ztxv, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilnx)

      IF (iter .LE. 0) THEN
        IF (.NOT.qsilent .AND.qprnt1) THEN

          WRITE (nttyo,1270) pxv, pdxv
          WRITE (noutpt,1270) pxv, pdxv

 1270     FORMAT(/9x,'pxv    = ',e16.9,' MPa',6x,'pdxv   = ',e16.9)

        ENDIF
      ENDIF

c     Now calculate the liquid properties by calling EVAI95
c     with the liquid delta.

      unam8 = 'Liquid'
      j2 = ILNOBL(unam8)

      IF (qprnt2) THEN
        WRITE (noutpt,1250) unam8(1:j2), deltsl, rhosl
      ENDIF

c     Calling sequence substitutions:

c       deltsl for delta
c       rhosl sfor rho
c       pxl for px
c       axl for ax
c       uxl for ux
c       sxl for sx
c       hxl for hx
c       cvxl for cvx
c       cpxl for cpx
c       wxl for wx
c       muxl for mux
c       dtxl for dtx
c       bsxl for bsx
c       gxl for gx
c       vxl for vx
c       pdxl for pdx
c       ptxl for ptx
c       adxl for adx
c       gdxl for gdx
c       avxl for avx
c       ktxl for ktx
c       ztxl for ztx
c       qsilnx for qsilent

      CALL EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, deltsl, rhosl, pxl, axl, uxl, sxl, hxl,
     $ cvxl, cpxl, wxl, muxl, dtxl, bsxl, gxl, vxl, pdxl, ptxl, adxl,
     $ gdxl, avxl, ktxl, ztxl, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilnx)

      IF (iter .LE. 0) THEN
        IF (.NOT.qsilent .AND.qprnt1) THEN

          WRITE (nttyo,1292) pxl, pdxl
          WRITE (noutpt,1292) pxl, pdxl

 1292     FORMAT(9x,'pxl    = ',e16.9,' MPa',6x,'pdxl   = ',e16.9)

        ENDIF
      ENDIF

      IF (qprnt2) THEN
        xx = (axv - axl)

        WRITE (noutpt,2420) axv, axl, xx
        WRITE (nttyo,2420) axv, axl, xx

 2420   FORMAT(/9x,'axv    = ',e16.9,' kJ/kg',4x,'axl   = ',e16.9,
     $  /9x,'adif   = ',e16.9,' kJ/kg')

        WRITE (noutpt,2430) adxv, adxl
        WRITE (nttyo,2430) adxv, adxl

 2430   FORMAT(9x,'adxv   = ',e16.9,' kJ/kg',4x,'adxl  = ',e16.9)

      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c       The pdx for vapor cannot be negative.
c       Under-relax to prevent this.

        IF (pdxv. lt. 0.0d0) THEN

          IF (iter .LE. 0) THEN
            IF (icutv .ge. 30) GO TO 410
            icutv = icutv + 1
            rhosv = 0.995d0*rhosv
            deltsv = rhosv/rhocr
            IF (qprnt1) THEN

              WRITE (nttyo,1140)
              WRITE (noutpt,1140)

 1140         FORMAT(/3x,'CALSCT: DECREASING INITIAL VAPOR DENSITY',
     $        ' TO AVOID NEGATIVE PDX',/)

            ENDIF
            GO TO 120
          ELSE
            IF (icutv .ge. 30) GO TO 410
            icutv = icutv + 1
            IF (qprnt1) THEN

              WRITE (nttyo,1150)
              WRITE (noutpt,1150)

 1150         FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID NEGATIVE',
     $        ' PDX FOR VAPOR',/)

            ENDIF
            arelax = 0.25d0
            GO TO 150
          ENDIF
        ENDIF
        icutv = 0

c       The pdx for liquid cannot be negative.
c       Under-relax to prevent this.

        IF (pdxl. lt. 0.0d0) THEN

          IF (iter .LE. 0) THEN
            IF (icutl .ge. 30) GO TO 410
            icutl = icutl + 1
            rhosl = 1.001d0*rhosl
            deltsl = rhosl/rhocr
            IF (qprnt1) THEN

              WRITE (nttyo,1160)
              WRITE (noutpt,1160)

 1160         FORMAT(/3x,'CALSCT: INCREASING INITIAL LIQUID DENSITY',
     $        ' TO AVOID NEGATIVE PDX',/)

            ENDIF
            GO TO 120
          ELSE
            IF (icutl .ge. 30) GO TO 410
            icutl = icutl + 1
            IF (qprnt1) THEN

              WRITE (nttyo,1170)
              WRITE (noutpt,1170)

 1170         FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID NEGATIVE',
     $        ' PDX FOR LIQUID',/)

            ENDIF
            arelax = 0.25d0
            GO TO 150
          ENDIF
        ENDIF
        icutl = 0

c       The revised delta for vapor cannot be less than a good
c       fraction of the value obtained from the intial approximation.

        IF (deltsv. lt. 0.90d0*dltsvq) THEN

          IF (icutvq .ge. 30) GO TO 410
          icutvq = icutvq + 1
          arelax = 0.25d0

          IF (qprnt1) THEN

            WRITE (nttyo,1180)
            WRITE (noutpt,1180)

 1180       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID VAPOR DELTA',
     $      /3x,'TOO MUCH LESS THAN THE APPROXIMATION VALUE.',/)

          ENDIF

          GO TO 150

        ENDIF
        icutvq = 0

c       The revised delta for liquid cannot be greater than a small
c       fraction of the value obtained from the intial approximation.

        IF (deltsl. gt. 1.10d0*dltslq) THEN

          IF (icutlq .ge. 30) GO TO 410
          icutlq = icutlq + 1
          IF (qprnt1) THEN

            WRITE (nttyo,1190)
            WRITE (noutpt,1190)

 1190       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID LIQUID DELTA',
     $      /3x,'TOO MUCH GREATER THAN THE APPROXIMATION VALUE.',/)

          ENDIF
          arelax = 0.25d0
          GO TO 150
        ENDIF
        icutlq = 0

        GO TO 170

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  150 CONTINUE

      IF (qprnt2) THEN
        IF (arelax .NE. 1.0d0) THEN
          WRITE (noutpt,1530) arelax
        ENDIF
      ENDIF

      xxx = 0.0d0

      DO j = 1,3
        deltas(j) = arelax*deltas(j)
        xx = DABS(deltas(j))
        xxx = DMAX1(xxx,xx)
      ENDDO

      psat = psat0 + deltas(3)
      deltsv = dltsv0 + deltas(1)
      deltsl = dltsl0 + deltas(2)
      rhosv = rhocr*deltsv
      rhosl = rhocr*deltsl
      GO TO 250

  170 CONTINUE

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Have obtained valid (outside the unstable zone) vapor and
c     liquid properties for the current iteration. Improve the
c     calculated saturation properties by solving three equations
c     in three unknowns. The equations are all in terms of pressure.
c     The unknowns to be found are psat, deltsv, and deltsl.

c     Calculate the Maxwell criterion pressure
c     (Gibbs energy equality expressed through the
c     Helmholtz energies and the pressure)

      dix = (1.0d0/deltsl) - (1.0d0/deltsv)
      dltx = deltsl - deltsv
      IF (.NOT.qsilent .AND. qprnt2) THEN

        WRITE (noutpt,1295) dix, dltx

 1295   FORMAT(9x,'dix    = ',e16.9,'    ',6x,'dltx   = ',e16.9)

      ENDIF

      IF (DABS(dix) .GT. 1.0d-15) THEN
        IF (DABS(axv - axl) .GT. 1.0d-15) THEN

c         Normal calculation, result in kPa.

          pxm = rhocr*( axv - axl )/dix

c         Convert from kPa to MPa.

          pxm = 0.001d0*pxm
        ELSE

c         There is no difference in the Helmholtz
c         energies of the vapor and the liquid.

          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1300)
            WRITE (noutpt,1300)

 1300       FORMAT(/3x,'CALSCT: THE TWO HELMHOLTZ ENERGIES ARE EQUAL.',
     $      /3x,'THE MAXWELL EQUATION CANNOT BE EVALUATED.')

          ENDIF

          GO TO 410
        ENDIF
      ELSE

c       Exception intended for the critical point.

        IF (DABS(tempk - tcr).LE.1.0d-10 .AND.
     $    DABS(deltsv - 1.0d0).LE.1.0d-10 .AND.
     $    DABS(deltsl - 1.0d0).LE.1.0d-10) THEN

c         Am at the critical point.

          pxm = pxv
        ELSE

c         Not at the critical point, but the vapor
c         and liquid densities have converged.

          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1310)
            WRITE (noutpt,1310)

 1310       FORMAT(/3x,'CALSCT: THE TWO DELTA VALUES ARE EQUAL.',
     $      /3x,'THE MAXWELL EQUATION CANNOT BE EVALUATED.')

          ENDIF

          GO TO 410
        ENDIF
      ENDIF

      IF (iter .LE. 0) THEN
        IF (.NOT.qsilent .AND. qprnt1) THEN

          WRITE (nttyo,1320) pxm
          WRITE (noutpt,1320) pxm

 1320     FORMAT(/6x,'Maxwell',/9x,'pxm    = ',e16.9,' MPa',/)

        ENDIF
      ENDIF

c     Calculate residual functions.

      alpha(1) = pxm - psat
      alpha(2) = pxv - psat
      alpha(3) = pxl - psat

      beta(1) = DABS(alpha(1)/psat)
      beta(2) = DABS(alpha(2)/psat)
      beta(3) = DABS(alpha(3)/psat)

      betamx = DMAX1(beta(1), beta(2), beta(3))

      IF (qprnt2) THEN

        WRITE (noutpt,1330) (alpha(i), i = 1,kdim)

 1330   FORMAT(/6x,'Residuals',
     $  /9x,'alpha(1) = ',e16.9,1x,'(Maxwell)',
     $  /9x,'alpha(2) = ',e16.9,1x,'(Vapor)',
     $  /9x,'alpha(3) = ',e16.9,1x,'(Liquid)')

        WRITE (noutpt,1332) (beta(i), i = 1,kdim)

 1332   FORMAT(/6x,'Relative Residuals',
     $  /10x,'beta(1) = ',e16.9,1x,'(Maxwell)',
     $  /10x,'beta(2) = ',e16.9,1x,'(Vapor)',
     $  /10x,'beta(3) = ',e16.9,1x,'(Liquid)')

      ENDIF

      IF (.NOT.qsilent) THEN

        WRITE (noutpt,1335) iter, psat, betamx
        WRITE (nttyo,1335) iter, psat, betamx

 1335   FORMAT(6x,'CALSCT: iter= ',i3,', psat= ',1pe16.9,', betamx= ',
     $  e12.5)

      ENDIF

c     Note: using a convergence tolerance below 1.0d-11
c     may lead to non-convergence due to the limitations
c     of 64-bit arithmetic.

      IF (betamx .LE. btxtol) THEN

c       Iteration has converged.

        press = psat

        IF (.NOT.qsilent .AND. qprnt1) THEN

          WRITE (nttyo,1360) press, betamx, btxtol
          WRITE (noutpt,1360) press, betamx, btxtol

 1360     FORMAT(/3x,'CALSCT: ITERATION CONVERGED TO press= ',1pe16.9,
     $    ' (MPa).',/3x,'MAX NORM betamx = ',e10.3,
     $    ', TOLERANCE btxtol = ',1pe10.3)

        ENDIF

        IF (deltsv .LE. 1.0d-12) THEN

          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1365)
            WRITE (noutpt,1365)

 1365       FORMAT(/3x,'CALSCT: THE SOLUTION IS A FALSE ONE. THE VAPOR',
     $      /3x,'DENSITY AND PRESSURE ARE NEAR-ZERO.',/)

          ENDIF

          GO TO 410
        ENDIF

        GO TO 420

      ELSEIF (iter .GE. itermx) THEN

c       Have done the maximum number of iterations.

        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1370)
          WRITE (noutpt,1370)

 1370     FORMAT(/3x,'CALSCT: HAVE DONE THE MAXIMUM NUMBER OF',
     $    ' ITERATIONS.')

        ENDIF

        GO TO 410

      ELSEIF (qxiter .AND. iter.GE.5) THEN

c       Have done the maximum number of iterations.
c       Report calculated results with warning.

        IF (.NOT.qsilent) THEN

          WRITE (nttyo,1375) betamx, btxtol
          WRITE (noutpt,1375) betamx, btxtol

 1375     FORMAT(/3x,'CALSCT: CALCULATED RESULTS ARE APPROXIMATE.',
     $    /3x,'MAX NORM betamx = ',1pe10.3,', TOLERANCE btxtol = ',
     $    e10.3)

        ENDIF

        GO TO 420

      ELSE

c       Make corrections and do another iteration.

        iter = iter + 1

c       The Jacobian matrix J here is aamatr(kdim,kdim).

        DO i = 1,3
          DO j = 1,3
            aamatr(i,j) = 0.30
          ENDDO
        ENDDO

        aamatr(1,1) = -(1.0d0/(dix*deltsv**2)) + (adxv/(axv - axl))
        aamatr(1,2) =  (1.0d0/(dix*deltsl**2)) - (adxl/(axv - axl))
        aamatr(1,1) = pxm*aamatr(1,1)
        aamatr(1,2) = pxm*aamatr(1,2)
        aamatr(1,3) = -1.0d0

        aamatr(2,1) = pdxv
        aamatr(2,2) = 0.0d0
        aamatr(2,3) = -1.0d0

        aamatr(3,1) = 0.0d0
        aamatr(3,2) = pdxl
        aamatr(3,3) = -1.0d0

        IF (qprnt3) THEN

          WRITE (noutpt,1410)

 1410     FORMAT(/3x,'Starting matrix')

          WRITE (noutpt,1420) ((aamatr(i,j), j = 1,kdim),
     $    alpha(i), i = 1,kdim)

 1420     FORMAT(/6x,'Jacobian and residual',
     $    /9x,1pe12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $    /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $    /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5)

        ENDIF

c       Since this matrix is only 3x3, solve using Gaussian
c       elimination with partial pivoting. There is no need
c       to use the "Jordan" part of the Gauss-Jordan algorithm
c       for a matrix of this size.

c       The simultaneous equations have the form:

c         aamatr(i,1)*deltas(1) + aamatr(i,2)*deltas(2)
c           aamatr(i,3)*deltas(3) = -alpha(i), i = 1,3

        CALL GAUSSE(noutpt, nttyo, qerr, qprnt3, kdim,
     $  kmax, aamatr, alpha, deltas)

        IF (qerr) THEN
          IF (.NOT.qsilent) THEN

            WRITE (nttyo,1430)
            WRITE (noutpt,1430)

 1430       FORMAT(/6x,'ERROR (CALSCT): Gaussian elimination failed.')

          ENDIF

          GO TO 410
        ENDIF

c       Reverse the sign of deltas to be consistent
c       with definition of correction.

        DO j = 1,kdim
          deltas(j) = -deltas(j)
        ENDDO

c       IF (qprnt2) THEN
c
c         WRITE (noutpt,1510)
c
c1510     FORMAT(/6x,'Corrections')
c
c         DO j = 1,kdim
c
c           WRITE (noutpt,1520) j, deltas(j)
c
c1520       FORMAT(13x,i2,3x,1pe16.9)
c
c         ENDDO
c       ENDIF

        ncut = 0
        arelax = 1.0d0

        IF (qprnt2) THEN
          IF (arelax .NE. 1.0d0) THEN

            WRITE (noutpt,1530) arelax

 1530       FORMAT(/9x,'arelax = ',e16.9)

          ENDIF
        ENDIF

c       Save current values.

        psat0 = psat
        dltsv0 = deltsv
        dltsl0 = deltsl

c       Make corrections.

  300   CONTINUE
        IF (qprnt2) THEN

          WRITE (noutpt,1510)

 1510     FORMAT(/6x,'Corrections')

          DO j = 1,kdim

            WRITE (noutpt,1520) j, deltas(j)

 1520       FORMAT(13x,i2,3x,1pe16.9)

          ENDDO
        ENDIF

        psat = psat0 + deltas(3)
        deltsv = dltsv0 + deltas(1)
        deltsl = dltsl0 + deltas(2)

c       The delta for liquid cannot be less than the delta for vapor.
c       Under-relax to prevent this.

        IF (deltsl. lt. deltsv) THEN
          IF (ncut .ge. 20) GO TO 410
          ncut = ncut + 1
          IF (qprnt1) THEN

            WRITE (nttyo,1540)
            WRITE (noutpt,1540)

 1540       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID VAPOR-LIQUID',
     $      ' DENSITY INVERSION',/)

          ENDIF
          arelax = 0.25d0
          GO TO 320
        ENDIF

c       Corrected delta values must be positive to avoid
c       a singularity in the equation-of-state model equations.

        IF (deltsv.LE.0.0d0 .OR. deltsl.LE.0.0d0) THEN
          IF (ncut .ge. 20) GO TO 410
          ncut = ncut + 1
          IF (qprnt1) THEN

            WRITE (nttyo,1580)
            WRITE (noutpt,1580)

 1580       FORMAT(/3x,'CALSCT: UNDER-RELAXING TO AVOID NEGATIVE',
     $      ' DELTA.',/)

          ENDIF
          arelax = 0.25d0
          GO TO 320
        ENDIF

        GO TO 330

  320   CONTINUE
        xxx = 0.0d0
        DO j = 1,3
          deltas(j) = arelax*deltas(j)
          xx = DABS(deltas(j))
          xxx = DMAX1(xxx,xx)
        ENDDO
        GO TO 300

  330   rhosv = rhocr*deltsv
        rhosl = rhocr*deltsl

c       Test for vapor-liquid density inversion.
c       That is, is the vapor more dense than the liquid?

        IF (qprnt2) THEN

          WRITE (noutpt,1600) psat, deltsv, deltsl

 1600     FORMAT(/6x,'Corrected variables',
     $    /9x,'psat   = ',e16.9,' MPa',
     $    /9x,'deltsv = ',e16.9,/9x,'deltsl = ',e16.9)

        ENDIF

        GO TO 250

      ENDIF

      GO TO 420

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  410 CONTINUE

c     Iteration failed.

      qfail = .TRUE.
      IF (.NOT.qsilent) THEN

        WRITE (nttyo,1690)
        WRITE (noutpt,1690)

 1690   FORMAT(/3x,'CALSCT: ITERATION FAILED.',/)

      ENDIF

      GO TO 999

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  420 CONTINUE

c     Calculations are complete for the last line read
c     from the input file. Write the results.

      IF (.NOT.qsilent) THEN

        WRITE (nttyo,2400)
        WRITE (noutpt,2400)

 2400   FORMAT(/3x,'CALSCT key results:',//8x,'Temp(K)',10x,
     $  'press(MPa)')

        WRITE (nttyo,2410) tempk,psat
        WRITE (noutpt,2410) tempk,psat

 2410   FORMAT(6x,f9.4,7x,e16.9)

        WRITE (noutpt,2415) tempk, psat, tau, deltsv, rhosv,
     $   deltsl, rhosl

 2415   FORMAT(/6x,'Temp(K) = ',f9.4
     $  //9x,'psat   = ',e16.9,' MPa',6x,'tau    = ',e16.9,
     $  /9x,'deltsv = ',e16.9,10x,'rhosv  = ',e16.9,
     $  /9x,'deltsl = ',e16.9,10x,'rhosl  = ',e16.9)

      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 999  CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE FDESCR(nttyo, noutpt, qprnt1, qprnt2, qprnt3,
     $ qwrphi, pcr, rcnstw, rhocr, tcr, bettol, btxtol, psat,
     $ rho, rhosv, rhosl, rhotol, tempk, tau, press, udescr)

c     Given the temperature and pressure, find the appropriate
c     description of the H2O fluid. A problem may occur if the
c     pressure is equal or nearly equal to the saturation pressure.
c     Here comparing the pressure with the saturation pressure
c     pressure may lead to the wrong description, as vapor and
c     liquid coexist at the saturation pressure. It then becomes
c     neccesary to compare the fluid density with the saturated
c     vapor and saturated liquid densities. If the density is
c     known (as following a CALDLT calculation), it will be used.
c     If it is not known (as in a CALPRE calculation), it will
c     not be used. In that case, the results obtained here will
c     determine the starting density estimate, thus in essence
c     choosing "vapor" or "liquid" for pressures close to the
c     saturation pressure.

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     CALLING SEQUENCE VARIABLES.

      LOGICAL qprnt1, qprnt2, qprnt3, qwrphi

      CHARACTER(LEN=24) udescr

      INTEGER nttyo, noutpt

      REAL(8) tempk, press

      REAL(8) pcr, rcnstw, rhocr, tau, tcr

      REAL(8) psat, rho, rhosv, rhosl, rhotol

      REAL(8) bettol, btxtol

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     LOCAL VARIABLES.

      LOGICAL qsilent, qerrx, qfailx

      INTEGER iter

      REAL(8) betamx, btxsav, deltsv, deltsl

      REAL(8) pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv, dtxv,
     $ bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv, ktxv, ztxv

      REAL(8) pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl, dtxl,
     $ bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl, ktxl, ztxl

      REAL(8) ptest, btest, rtestl, rtestv

      DATA qsilent / .TRUE. /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      udescr = 'unknown'
      qerrx = .FALSE.
      qfailx = .FALSE.
      psat = 0.0d0
      rhosv = 0.0d0
      rhosl = 0.0d0

      IF (tempk.LT.273.15d0) THEN

c       Note that the allowed temperature range has been extended a bit
c       on the low end to include 0C.

        WRITE (nttyo,1010)
        WRITE (noutpt,1010)

 1010   FORMAT(/6x,'WARNING (FDESCR): Temperature is too low to',
     $  ' assign',/9x,'a description to the H2O fluid.')

        GO TO 999
      ENDIF

      IF (tempk .LE. tcr) THEN

        btxsav = btxtol

c       Calculate the saturation curve properties.

c       Calling sequence substitutions:

c         psat for press
c         qerrx for qerr
c         qfailx for qfail

        CALL CALSCT(nttyo, noutpt, qerrx, qfailx, qprnt1, qprnt2,
     $  qprnt3, qsilent, qwrphi, pcr, rcnstw, rhocr, tcr, iter,
     $  betamx, bettol, btxtol, tempk, tau, psat, rhosv, rhosl,
     $  deltsv, deltsl, pxv, uxv, sxv, hxv, cvxv, cpxv, wxv, muxv,
     $  dtxv, bsxv, axv, gxv, vxv, pdxv, ptxv, adxv, gdxv, avxv,
     $  ktxv, ztxv, pxl, uxl, sxl, hxl, cvxl, cpxl, wxl, muxl,
     $  dtxl, bsxl, axl, gxl, vxl, pdxl, ptxl, adxl, gdxl, avxl,
     $  ktxl, ztxl)

        IF (qerrx) THEN

          WRITE (nttyo,1020)
          WRITE (noutpt,1020)

 1020     FORMAT(/6x,'WARNING (FDESCR): The call to CALSCT to',
     $    ' calculate',/9x,'the saturation pressure failed.')

          rhosv = 0.0d0
          rhosl = 0.0d0
        ENDIF

c       Reset the convergence tolerance.

        btxtol = btxsav

        IF (press .GT. pcr) THEN
          udescr = 'compressed liquid'

          IF (qerrx) THEN

            WRITE (nttyo,1030)
            WRITE (noutpt,1030)

 1030       FORMAT(/6x,'WARNING (FDESCR): Because the call to CALSCT',
     $      /9x,'failed, an arbitrary liquid-like density will be',
     $      /9x,'assigned as a starting value for compressed liquid.')

            rhosl = 1.05d0
          ENDIF

        ELSE

          IF (qerrx) THEN

            WRITE (nttyo,1040)
            WRITE (noutpt,1040)

 1040       FORMAT(/6x,'WARNING (FDESCR): Because the call to CALSCT',
     $      /9x,'failed, the vapor and liquid states cannot be',
     $      /9x,'distinguished from one another. Liquid is assigned',
     $      /9x,'arbitrarily.')

            udescr = 'liquid'
            rhosl = 1.05d0
          ELSE

            IF (press .GE. psat) THEN
              udescr = 'liquid'
            ELSE
              udescr = 'vapor'
            ENDIF

c           Use density (rho) if available and pressure is close
c           to psat.

            IF (rho .GT. 0.0d0) THEN

              ptest = (press - psat)/psat
              btest = 10*btxtol

              IF (DABS(ptest) .LE. btest) THEN

c               Here press is very close to psat.
c               Use rho to determine vapor or liquid.

                rtestl = (rho - rhosl)/rhosl
                rtestv = (rho - rhosv)/rhosv

                IF (DABS(rtestl) .LE. rhotol) THEN
                  udescr = 'liquid'
                ELSEIF (DABS(rtestv) .LE. rhotol) THEN
                  udescr = 'vapor'
                ELSE

                  WRITE (nttyo,1050)
                  WRITE (noutpt,1050)

 1050             FORMAT(/6x,'WARNING (FDESCR): Could not use density',
     $            ' rho to discriminate',/6x,'vapor from liquid.')

                ENDIF

              ENDIF

            ENDIF

          ENDIF

        ENDIF

      ELSE

        IF (press .GT. pcr) THEN
          udescr = 'supercritical fluid'
        ELSE
          udescr = 'hot vapor'
        ENDIF

      ENDIF

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE WROTAB(nttyo, noutpt, nmcsvf, nccsvf, nxcsvf, unam8,
     $ udescr, ushdes, tempk, delta, rho, px, ux, hx, sx, ax,
     $ gx, vx, cvx, cpx, wx, mux, dtx, bsx, pdx, avx, ktx,
     $ ztx, vsx, thcx, stx, epsx, deptp, deppt, agme, agm10,
     $ aphi, adhv, adhh, ahrt, bgm, bdhv, bdhh, bhrt, pmltx,
     $ tmltx, uxcu, hxcu, sxcu, axcu, gxcu, vxcu, cvxcu,
     $ cpxcu, dtxcu, iter, betamx, btxtol, qrhog, rhog)

c     This routine writes summary results to the output file, the main
c     .csv file, and the thermochemical .csv file.

      LOGICAL qrhog

      INTEGER nttyo, noutpt, nmcsvf, nccsvf, nxcsvf

      INTEGER iter

      CHARACTER(LEN=8) unam8

      CHARACTER(LEN=24) udescr

      CHARACTER(LEN=16) ushdes

      REAL(8) tempk, delta, rho, betamx, btxtol, rhog

      REAL(8)  px, ux, sx, hx, ax, gx, vx, cvx, cpx, wx, mux, dtx,
     $ bsx, pdx, avx, ktx, ztx, vsx, thcx, stx, epsx, deptp, deppt,
     $ agme, agm10, aphi, adhv, adhh, ahrt, bgm, bdhv, bdhh, bhrt,
     $ pmltx, tmltx

      REAL(8)  uxcu, sxcu, hxcu, axcu, gxcu, vxcu, cvxcu, cpxcu, dtxcu

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      INTEGER ILNOBL

      INTEGER j2, j3, j4

      REAL(8) asize

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      j2 = ILNOBL(unam8)
      j3 = ILNOBL(udescr)

      WRITE (nttyo,1010) unam8(1:j2), udescr(1:j3)
      WRITE (noutpt,1010) unam8(1:j2), udescr(1:j3)

 1010 FORMAT(/7x,a,'(',a,')')

      WRITE (nttyo,1020)
      WRITE (noutpt,1020)

 1020 FORMAT(/16x,'delta',11x,'rho(kg/m3)')

      WRITE (nttyo,1030) delta, rho
      WRITE (noutpt,1030) delta, rho

 1030 FORMAT(10x,e16.9,2x,e16.9)

c     Write the usual equation-of-state thermodynamic results.

      WRITE (noutpt,1040)

 1040 FORMAT(/9x,'Thermodynamic Results')

      WRITE (noutpt,1050)

 1050 FORMAT(/15x,'u(kJ/kg)',10x,'h(kJ/kg)',9x,'s(kJ/kg/K)')

      WRITE (noutpt,1060) ux, hx, sx

 1060 FORMAT(8x,3(2x,e16.9))

      WRITE (noutpt,1070)

 1070 FORMAT(/15x,'a(kJ/kg)',10x,'g(kJ/kg)',9x,'v(m3/kg)')
      WRITE (noutpt,1060) ax, gx, vx

      WRITE (noutpt,1080)

 1080 FORMAT(/13x,'cv(kJ/kg/K)',7x,'cp(kJ/kg/K)',9x,'w(m/s)')
      WRITE (noutpt,1060) cvx, cpx, wx

      WRITE (noutpt,1090)

 1090 FORMAT(/14x,'mu(K/MPa)',7x,'dt(kJ/kg/MPa)',7x,'bs(K/MPa)')

      WRITE (noutpt,1060) mux, dtx, bsx

      WRITE (noutpt,1100)

 1100 FORMAT(/15x,'av(/K)',12x,'kt(/MPa)')

      WRITE (noutpt,1060) avx, ktx

c     Write thermochemical results.

      WRITE (noutpt,1140)

 1140 FORMAT(/9x,'Thermochemical Results')

      WRITE (noutpt,1150)

 1150 FORMAT(/15x,'g(kJ/mol)',8x,'h(kJ/mol)',8x,'s(J/mol/K)')

      asize = DMAX1(DABS(gxcu), DABS(hxcu), DABS(sxcu))

      IF (asize .LE. 1.0d+06) THEN

        WRITE (noutpt,1160) gxcu, hxcu, sxcu

 1160   FORMAT(12x,f12.4,5x,f12.4,5x,f12.4)

      ELSE

        WRITE (noutpt,1170) gxcu, hxcu, sxcu

 1170   FORMAT(12x,1pe12.5,5x,e12.5,6x,e12.5)

      ENDIF

      WRITE (noutpt,1180)

 1180 FORMAT(/15x,'a(kJ/mol)',8x,'u(kJ/mol)',8x,'v(m3/mol)')

      asize = DMAX1(DABS(axcu), DABS(uxcu), DABS(vxcu))
      IF (asize .LE. 1.0d+06) THEN

        IF (DABS(vxcu) .GT. 0.10d0) THEN

          WRITE (noutpt,1160) axcu, uxcu, vxcu

        ELSE

          WRITE (noutpt,1165) axcu, uxcu, vxcu

 1165     FORMAT(12x,f12.4,5x,f12.4,5x,1pe12.5)

        ENDIF

      ELSE

        WRITE (noutpt,1170) axcu, uxcu, vxcu

      ENDIF

      WRITE (noutpt,1190)

 1190 FORMAT(/13x,'cv(J/mol/K)',7x,'cp(J/mol/K)',5x,'dt(kJ/mol/Mpa')

      asize = DMAX1(DABS(cvxcu), DABS(cpxcu))
      IF (asize .LE. 1.0d+06) THEN

        WRITE (noutpt,1200) cvxcu, cpxcu, dtxcu

 1200   FORMAT(12x,f12.4,5x,f12.4,6x,1pe12.5)
      ELSE

        WRITE (noutpt,1170) cvxcu, cpxcu, dtxcu

      ENDIF

      IF (px .LE. 0.0d0) THEN

        WRITE (nttyo,1202)
        WRITE (noutpt,1202)

 1202   FORMAT(/6x,'WARNING (WROTAB): PRESSURE IS NOT POSITIVE.',
     $  /9x,'HAVE A NON-PHYSICAL SOLUTION.')

      ENDIF

      IF (pdx .LT. 0.0d0) THEN

        WRITE (nttyo,1205)
        WRITE (noutpt,1205)

 1205   FORMAT(/6x,'WARNING (WROTAB): PRESSURE DERIVATIVE WITH RESPECT',
     $  /9x,'TO DENSITY IS LESS THAN ZERO. HAVE AN UNSTABLE SOLUTION.')

      ENDIF

      IF (wx .LE. 0.0d0) THEN

        WRITE (nttyo,1210)
        WRITE (noutpt,1210)

 1210   FORMAT(/6x,'WARNING (WROTAB): SPEED OF SOUND VALUE IS NOT',
     $  ' POSITIVE.',/9x,'MAY HAVE A NON-PHYSICAL SOLUTION.')

      ENDIF

      IF (cpx .LE. 0.0d0) THEN

        WRITE (nttyo,1215)
        WRITE (noutpt,1215)

 1215   FORMAT(/6x,'WARNING (WROTAB): HEAT CAPACITY (CP) VALUE IS NOT',
     $  ' POSITIVE.',/9x,'MAY HAVE A NON-PHYSICAL SOLUTION.')

      ENDIF

      IF (sxcu .LE. 0.0d0) THEN

        WRITE (nttyo,1220)
        WRITE (noutpt,1220)

 1220   FORMAT(/6x,'WARNING (WROTAB): THERMOCHEMICAL ENTROPY VALUE',
     $  ' IS NOT POSITIVE.',/9x,'MAY HAVE A NON-PHYSICAL SOLUTION.')

      ENDIF

      j4 = ILNOBL(ushdes)

c     Write results to the main .csv file (mtab.csv).

      IF (.NOT.qrhog) THEN

c       Note that qrhog (specified initial value for density
c       when iterating to find the desired pressure) is of
c       of limited interest).

        WRITE (nmcsvf,1520) unam8(1:j2), ushdes(1:j4), tempk, px,
     $  rho, ux, hx, sx, ax, gx, vx, cvx, cpx, wx, mux, dtx, bsx,
     $  iter, betamx, btxtol

 1520   FORMAT(a,',',a,',',f11.6,',',1pe16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',
     $  i3,',',e10.3,',',e10.3)

      ELSE

        WRITE (nmcsvf,1522) unam8(1:j2), ushdes(1:j4), tempk, px,
     $  rho, ux, hx, sx, ax, gx, vx, cvx, cpx, wx, mux, dtx, bsx,
     $  iter, betamx, btxtol, rhog

 1522   FORMAT(a,',',a,',',f11.6,',',1pe16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',
     $  i3,',',e10.3,',',e10.3,',',e16.9)

      ENDIF

c     Write results to the thermochemical .csv file (ctab.csv).

        WRITE (nccsvf,1540) unam8(1:j2), ushdes(1:j4), tempk, px, rho,
     $  uxcu, hxcu, sxcu, axcu, gxcu, vxcu, cvxcu, cpxcu, wx, mux,
     $  dtxcu, bsx

 1540   FORMAT(a,',',a,',',f11.6,',',1pe16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9)

c     Write results to the extra properties .csv file (xtab.csv).

c       adhh = AV, ahrt = AV/RT
c       bdhh = BV, bhrt = BV/RT

        WRITE (nxcsvf,1550) unam8(1:j2), ushdes(1:j4), tempk, px, rho,
     $  avx, ktx, ztx, vsx, thcx, stx, epsx, deptp, deppt, agme,
     $  agm10, aphi, adhv, ahrt, bgm, bdhv, bhrt, pmltx, tmltx

 1550   FORMAT(a,',',a,',',f11.6,',',1pe16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',
     $  e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9,',',e16.9)

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE EVAI95(noutpt, nttyo, tcr, rhocr, pcr, rcnstw,
     $ press, tempk, tau, delta, rho, px, ax, ux, sx, hx,
     $ cvx, cpx, wx, mux, dtx, bsx, gx, vx, pdx, ptx, adx,
     $ gdx, avx, ktx, ztx, qprnt1, qprnt2, qprnt3,
     $ qwrphi, qsilent)

c     Evaluate the basic equation of state, which is written as
c     a function of:

c       reduced density: delta = rho/rhoc
c       inverse reduced temperature: tau = tcr/tempk

c     The pressure px is an output. To obtain properties as
c     a function of temperature and pressure, it is necessary to
c     iterate on pressure (adjusting density).

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qprnt1, qprnt2, qprnt3, qwrphi, qsilent

      INTEGER nttyo, noutpt

      REAL(8) tcr, rhocr, pcr, rcnstw
      REAL(8) delta, tau, rho, tempk, press
      REAL(8) ax, px, ux, sx, hx, cvx, cpx, wx, mux, dtx, bsx, gx,
     $ vx, pdx, ptx, adx, gdx, avx, ktx, ztx

      REAL(8) cin0(1:8), cigam0(4:8)

c     Original Wagner and Pruss (2002) values for cin0(1) and cin0(2).

c     DATA cin0(1:8) / -8.32044648201d+00, 6.6832105268d+00,
c    $ 3.00632d+00, 0.012436d+00, 0.97315d+00, 1.27950d+00,
c    $ 0.96956d+00, 0.24873d+00 /

c     IAPWS (2016) revised values for cin0(1) and cin0(2).

      DATA cin0(1:8) / -8.3204464837497d+00, 6.6832105275932d+00,
     $ 3.00632d+00, 0.012436d+00, 0.97315d+00, 1.27950d+00,
     $ 0.96956d+00, 0.24873d+00 /

      DATA cigam0(4:8) / 1.28728967d+00, 3.53734222d+00,
     $ 7.74073708d+00, 9.24437796d+00, 27.5075105d+00 /

c     Note: cigam0(1:3) are not used in the EOS model.

c     Coefficients for the residual part.

c     First part, indices 1-51.

      INTEGER c1c(8:54), c1d(1:54)

      REAL(8) c1t(1:54), c1n(1:54)

      DATA c1c(8:51) / 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     $  2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     $  3, 3, 3, 3, 4, 6, 6, 6, 6 /

c     Note: c1c(1:7) are not used in the EOS model.

      DATA c1d(1:51) / 1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5,
     $ 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9,
     $ 9, 9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6 /

      DATA c1t(1:51) / -0.5d0, 0.875d0, 1.0d0, 0.5d0, 0.75d0, 0.375d0,
     $ 1.0d0, 4.0d0, 6.0d0, 12.0d0, 1.0d0, 5.0d0, 4.0d0, 2.0d0,
     $ 13.0d0, 9.0d0, 3.0d0, 4.0d0, 11.0d0, 4.0d0, 13.0d0, 1.0d0,
     $ 7.0d0, 1.0d0, 9.0d0, 10.0d0, 10.0d0, 3.0d0, 7.0d0, 10.0d0,
     $ 10.0d0, 6.0d0, 10.0d0, 10.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0,
     $ 8.0d0, 6.0d0, 9.0d0, 8.0d0, 16.0d0, 22.0d0, 23.0d0, 23.0d0,
     $ 10.0d0, 50.0d0, 44.0d0, 46.0d0, 50.0d0 /

      DATA c1n(1:51) /  0.12533547935523d-01,
     $  0.78957634722828d+01, -0.87803203303561d+01,
     $  0.31802509345418d+00, -0.26145533859358d+00,
     $ -0.78199751687981d-02,  0.88089493102134d-02,
     $ -0.66856572307965d+00,  0.20433810950965d+00,
     $ -0.66212605039687d-04, -0.19232721156002d+00,
     $ -0.25709043003438d+00,  0.16074868486251d+00,
     $ -0.40092828925807d-01,  0.39343422603254d-06,
     $ -0.75941377088144d-05,  0.56250979351888d-03,
     $ -0.15608652257135d-04,  0.11537996422951d-08,
     $  0.36582165144204d-06, -0.13251180074668d-11,
     $ -0.62639586912454d-09, -0.10793600908932d+00,
     $  0.17611491008752d-01,  0.22132295167546d+00,
     $ -0.40247669763528d+00,  0.58083399985759d+00,
     $  0.49969146990806d-02, -0.31358700712549d-01,
     $ -0.74315929710341d+00,  0.47807329915480d+00,
     $  0.20527940895948d-01, -0.13636435110343d+00,
     $  0.14180634400617d-01,  0.83326504880713d-02,
     $ -0.29052336009585d-01,  0.38615085574206d-01,
     $ -0.20393486513704d-01, -0.16554050063734d-02,
     $  0.19955571979541d-02,  0.15870308324157d-03,
     $ -0.16388568342530d-04,  0.43613615723811d-01,
     $  0.34994005463765d-01, -0.76788197844621d-01,
     $  0.22446277332006d-01, -0.62689710414685d-04,
     $ -0.55711118565645d-09, -0.19905718354408d+00,
     $  0.31777497330738d+00, -0.11841182425981d+00 /

c     Second part, indices 52-54.

      INTEGER c2d(52:54), c2t(52:54), c2alph(52:54), c2beta(52:54),
     $ c2eps(52:54)

      REAL(8) c2n(52:54), c2gamm(52:54)

c     Note: c2c(52:54) are not used in the EOS model.

      DATA c2d(52:54) / 3, 3, 3 /
      DATA c2t(52:54) / 0, 1, 4 /
      DATA c2n(52:54) / -0.31306260323435d+02, 0.31546140237781d+02,
     $ -0.25213154341695d+04 /

      DATA c2alph(52:54) / 20, 20, 20 /
      DATA c2beta(52:54)  / 150, 150, 250 /
      DATA c2gamm(52:54) / 1.21d0, 1.21d0, 1.25d0 /
      DATA c2eps(52:54)   / 1, 1, 1 /

c     Third part, indices 55-56.

      INTEGER c3ccap(55:56), c3dcap(55:56)

      REAL(8) c3a(55:56), c3b(55:56), c3bcap(55:56), c3n(55:56),
     $ c3acap(55:56), c3beta(55:56)

      DATA c3a(55:56) / 3.5d0, 3.50d0 /
      DATA c3b(55:56) / 0.85d0, 0.95d0 /
      DATA c3bcap(55:56) / 0.2d0, 0.2d0 /
      DATA c3n(55:56) / -0.14874640856724d0, 0.31806110878444d0 /
      DATA c3ccap(55:56) / 28, 32 /
      DATA c3dcap(55:56) / 700, 800 /
      DATA c3acap(55:56) / 0.32d0, 0.32d0 /
      DATA c3beta(55:56) / 0.3d0, 0.3d0 /

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      REAL(8) phi0, phi0d, phi0dd, phi0t, phi0tt, phi0dt

      REAL(8) phir, phird, phirdd, phirt, phirtt, phirdt

      REAL(8) delsq, tausq

      REAL(8) egam0(4:8), xegam0(4:8)

      REAL(8) deldi(1:54),tauti(1:54)
      REAL(8) deldm1(1:54), deldm2(1:54)
      REAL(8) delci(8:51), dmci(8:51), edci(8:51)
      REAL(8) tautm1(1:51), tautm2(1:51)

      REAL(8) exxi(52:54),dade(52:54),tbtg(52:54)

      REAL(8) TH(55:56), HK(55:56), HKK(55:56), BEX(55:56)
      REAL(8) WW(55:56), WWd(55:56), WWdd(55:56)
      REAL(8) HH(55:56), HHd(55:56), HHdd(55:56), HHt(55:56),
     $ HHtt(55:56), HHdt(55:56)
      REAL(8) PS(55:56), PSd(55:56), PSdd(55:56), PSt(55:56),
     $ PStt(55:56), PSdt(55:56)

      REAL(8) XC3B, XC3B2

      REAL(8) dm1, dm1sq, tm1, tm1sq

      REAL(8) epxc, rtx, x1, x2, x3, x4, xsum, xxt, ztxc

      INTEGER i

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      epxc = 100*EPSILON(delta)

c     For 64-bit arithmetic, 100*EPS is approximately 2.2204d-14.

c     Avoid a singularity for delta too close to zero.

      IF (delta .LT. epxc) THEN
        delta = epxc
        rho = rhocr*delta
      ENDIF

c     Avoid a singularity at the critical density (at any temperature).
c     If delta = rho/rhocr is unity, then delta - 1 is zero. This
c     cannot be used with (delta - 1)**n where n is negative, or in a
c     division (x/(d -1)).

      x1 = 1.0d0 - epxc
      x2 = 1.0d0 + epxc

      IF (delta.GT.x1 .AND. delta.LT.x2) THEN

c       The density is too close to the critical density.

        IF (qprnt1) THEN
          IF (.NOT.qsilent) THEN

            WRITE(nttyo,1110)
            WRITE(noutpt,1110)

 1110       FORMAT(/3x,'TOO CLOSE TO THE CRITICAL DENSITY,',
     $      ' ADJUSTING.',/)

          ENDIF
        ENDIF

        IF (delta .LT. 1.0d0) THEN
          delta = x1
        ELSE
          delta = x2
        ENDIF

        rho = rhocr*delta

      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Zero parameters to be calculated.

      px = 0.0d0
      ax = 0.0d0
      ux = 0.0d0
      sx = 0.0d0
      hx = 0.0d0
      cvx = 0.0d0
      cpx = 0.0d0
      wx = 0.0d0
      mux = 0.0d0
      dtx = 0.0d0
      bsx = 0.0d0
      gx = 0.0d0
      vx = 0.0d0
      pdx = 0.0d0
      ptx = 0.0d0
      adx = 0.0d0
      gdx = 0.0d0
      avx = 0.0d0
      ktx = 0.0d0
      ztx = 0.0d0

c     Calculate some common pieces.

      delsq = delta**2
      tausq = tau**2
      rtx = rcnstw*tempk

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate the ideal-gas parts.

c     Calculate common pieces first.

      DO i = 4,8
        egam0(i) = DEXP( -cigam0(i)*tau )
        xegam0(i) = 1.0d0 - egam0(i)
      ENDDO

c     Now calculate phi0 and its partial derivatives.

      phi0 = DLOG(delta) + cin0(1) + cin0(2)*tau
     $ + cin0(3)*DLOG(tau)
      DO i = 4,8
        phi0 = phi0 + cin0(i)*DLOG(xegam0(i))
      ENDDO

      phi0d = 1.0d0/delta

      phi0dd = -1.0d0/delsq

      phi0t = cin0(2) + (cin0(3)/tau)
      DO i = 4,8
        phi0t = phi0t + cin0(i)*cigam0(i)*( (xegam0(i)**(-1)) -1.0d0 )
      ENDDO

      phi0tt = -cin0(3)/tausq
      DO i = 4,8
        phi0tt = phi0tt
     $  - cin0(i)*(cigam0(i)**2)*egam0(i)*( xegam0(i)**(-2) )
      ENDDO

      phi0dt = 0.0d0

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate the residual parts.

c     Calculate common pieces.

      dm1 = delta - 1.0d0
      dm1sq = dm1**2
      tm1 = tau - 1.0d0
      tm1sq = tm1**2

      DO i = 1,51
        deldi(i) = delta**c1d(i)
        tauti(i) = tau**c1t(i)
        deldm1(i) = delta**(c1d(i) - 1.0d0)
        deldm2(i) = delta**(c1d(i) - 2.0d0)
        tautm1(i) = tau**(c1t(i) - 1.0d0)
        tautm2(i) = tau**(c1t(i) - 2.0d0)
      ENDDO

      DO i = 8,51
        delci(i) = delta**c1c(i)
        edci(i) = DEXP(-delci(i))
        dmci(i) = c1d(i) - c1c(i)*delci(i)
      ENDDO

      DO i = 52,54
        deldi(i) = delta**c2d(i)
        tauti(i) = tau**c2t(i)
        deldm1(i) = delta**(c2d(i) - 1.0d0)
        deldm2(i) = delta**(c2d(i) - 2.0d0)
        exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
     $  - c2beta(i)*(tau - c2gamm(i))**2)
        dade(i) = (c2d(i)/delta) -2.0d0*c2alph(i)*(delta - c2eps(i))
        tbtg(i) = (c2t(i)/tau) -2.0d0*c2beta(i)*(tau - c2gamm(i))
      ENDDO

      DO i = 55,56
        XC3B = 1.0d0/( 2.0d0*c3beta(i) )

        BEX(i) = XC3B - 1.0d0
        TH(i) = (1.0d0 - tau) + c3acap(i)*dm1sq**XC3B
        PS(i) = DEXP( -c3ccap(i)*dm1sq - c3dcap(i)*tm1**2 )

        WW(i) = TH(i)**2 + c3bcap(i)*dm1sq**c3a(i)
        HH(i) = WW(i)**c3b(i)
        HK(i) = WW(i)**(c3b(i) - 1.0d0)
        HKK(i) = WW(i)**(c3b(i) - 2.0d0)

        WWd(i) = dm1*( c3acap(i)*TH(i)*(2.0d0/c3beta(i))*dm1sq**BEX(i)
     $  + 2.0d0*c3bcap(i)*c3a(i)*dm1sq**(c3a(i) - 1.0d0) )

        x1 = 4.0d0*c3bcap(i)*c3a(i)*(c3a(i) - 1.0d0)
     $  *dm1sq**(c3a(i) - 2.0d0)

        XC3B2 = (1.0d0/c3beta(i))**2
        x2 = 2.0d0*(c3acap(i)**2)*XC3B2*( dm1sq**BEX(i) )**2

        x3 = c3acap(i)*TH(i)*( 4.0d0/c3beta(i) )*BEX(i)
     $  *dm1sq**(BEX(i) - 1.0d0)

        xsum = x1 + x2 + x3
        WWdd(i) = (WWd(i)/dm1) + dm1sq*xsum

        HHd(i) = c3b(i)*HK(i)*WWd(i)

        HHdd(i) = c3b(i)*( HK(i)*WWdd(i) + (c3b(i) - 1.0d0)
     $  *HKK(i)*( WWd(i)**2 ) )

        HHt(i) = -2.0d0*TH(i)*c3b(i)*HK(i)

        HHtt(i) = 2.0d0*c3b(i)*HK(i) + 4.0d0*(TH(i)**2)*c3b(i)
     $  *(c3b(i) - 1.0d0)*HKK(i)

        x1 = 2.0d0/c3beta(i)
        HHdt(i) = -c3acap(i)*c3b(i)*x1*HK(i)*dm1*(dm1sq**BEX(i))
     $  - 2.0d0*TH(i)*c3b(i)*(c3b(i) - 1.0d0)*HKK(i)*WWd(i)

        PSd(i) = -2.0d0*c3ccap(i)*dm1*PS(i)
        PSdd(i) = ( 2.0d0*c3ccap(i)*dm1sq - 1.0d0 )
     $  *2.0d0*c3ccap(i)*PS(i)
        PSt(i) = -2.0d0*c3dcap(i)*tm1*PS(i)
        PStt(i) = ( 2.0d0*c3dcap(i)*tm1sq -1 )*2.0d0*c3dcap(i)*PS(i)
        PSdt(i) = 4.0d0*c3ccap(i)*c3dcap(i)*dm1*tm1*PS(i)

      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Now calculate phir and its partial derivatives.

c     Calculate phir.

      phir = 0.0d0

      DO i =1,7

c       deldi(i) = delta**c1d(i)
c       tauti(i) = tau**c1t(i)

        phir = phir + c1n(i)*deldi(i)*tauti(i)
      ENDDO

      DO i = 8,51

c       delci(i) = delta**c1c(i)
c       deldi(i) = delta**c1d(i)
c       tauti(i) = tau**c1t(i)
c       dmci(i) = c1d(i) - c1c(i)*delci(i)
c       edci(i) = DEXP(-delci(i))

        phir = phir + c1n(i)*deldi(i)*tauti(i)*edci(i)
      ENDDO

      DO i = 52,54

c       deldi(i) = delta**c2d(i)
c       tauti(i) = tau**c2t(i)
c       exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
c    $  - c2beta(i)*(tau - c2gamm(i))**2)

        phir = phir + c2n(i)*deldi(i)*tauti(i)*exxi(i)
      ENDDO

      DO i = 55,56

c       XC3B = 1.0d0/( 2.0d0*c3beta(i) )
c       TH(i) = (1.0d0 - tau) + c3acap(i)*dm1sq**XC3B
c       WW(i) = TH(i)**2 + c3bcap(i)*dm1sq**c3a(i)
c       HH(i) = WW(i)**c3b(i)
c       PS(i) = DEXP( -c3ccap(i)*dm1sq - c3dcap(i)*tm1**2 )

        phir = phir + c3n(i)*HH(i)*delta*PS(i)
      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phird.

      phird = 0.0d0

      DO i =1,7

c       deldm1(i) = delta**(c1d(i) - 1.0d0)
c       tauti(i) = tau**c1t(i)

        phird = phird + c1n(i)*c1d(i)*deldm1(i)*tauti(i)
      ENDDO

      DO i = 8,51

c       delci(i) = delta**c1c(i)
c       edci(i) = DEXP(-delci(i))
c       deldm1(i) = delta**(c1d(i) - 1.0d0)
c       tauti(i) = tau**c1t(i)
c       dmci(i) = c1d(i) - c1c(i)*delci(i)

        phird = phird + c1n(i)*edci(i)*deldm1(i)*tauti(i)*dmci(i)
      ENDDO

      DO i = 52,54

c       deldi(i) = delta**c2d(i)
c       tauti(i) = tau**c2t(i)
c       exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
c    $  - c2beta(i)*(tau - c2gamm(i))**2)
c       dade(i) = (c2d(i)/delta) -2.0d0*c2alph(i)*(delta - c2eps(i))

        phird = phird + c2n(i)*deldi(i)*tauti(i)*exxi(i)*dade(i)
      ENDDO

      DO i = 55,56

c       XC3B = 1.0d0/( 2.0d0*c3beta(i) )
c       BEX(i) = XC3B - 1.0d0
c       TH(i) = (1.0d0 - tau) + c3acap(i)*dm1sq**XC3B
c       WW(i) = TH(i)**2 + c3bcap(i)*dm1sq**c3a(i)
c       HH(i) = WW(i)**c3b(i)
c       PS(i) = DEXP( -c3ccap(i)*dm1sq - c3dcap(i)*tm1**2 )
c       HK(i) = WW(i)**(c3b(i) - 1.0d0)
c       PSd(i) = -2.0d0*c3ccap(i)*dm1*PS(i)
c       WWd(i) = dm1*( c3acap(i)*TH(i)*(2.0d0/c3beta(i))*dm1sq**BEX(i)
c    $  + 2.0d0*c3bcap(i)*c3a(i)*( dm1sq**(c3a(i) - 1.0d0) ) )
c       HHd(i) = c3b(i)*HK(i)*WWd(i)

        phird = phird + c3n(i)*( HH(i)*(PS(i) + delta*PSd(i))
     $  + HHd(i)*delta*PS(i) )
      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirdd.

      phirdd = 0.0d0

      DO i =1,7
        phirdd = phirdd + c1n(i)*c1d(i)*(c1d(i) - 1.0d0)
     $  *deldm2(i)*tauti(i)
      ENDDO

      DO i = 8,51

c       delci(i) = delta**c1c(i)
c       edci(i) = DEXP(-delci(i))
c       dmci(i) = c1d(i) - c1c(i)*delci(i)

        xxt = dmci(i)*(dmci(i) - 1.0d0) - (c1c(i)**2)*delci(i)
        phirdd = phirdd + c1n(i)*edci(i)*deldm2(i)*tauti(i)*xxt
      ENDDO

      DO i = 52,54

c       exxi(i) = DEXP( -c2alph(i)*(delta - c2eps(i))**2
c    $  - c2beta(i)*(tau - c2gamm(i))**2)

        x1 = -2.0d0*c2alph(i)*deldi(i)
        x2 = 4.0d0*(c2alph(i)**2)*deldi(i)*(delta - c2eps(i))**2
        x3 = -4.0d0*c2d(i)*c2alph(i)*deldm1(i)*(delta - c2eps(i))
        x4 = c2d(i)*(c2d(i) - 1.0d0)*deldm2(i)
        xsum = x1 + x2 + x3 + x4
        phirdd = phirdd + c2n(i)*tauti(i)*exxi(i)*xsum
      ENDDO

      DO i = 55,56
        x1 = HH(i)*( 2.0d0*PSd(i) + delta*PSdd(i) )
        x2 = 2.0d0*HHd(i)*( PS(i) + delta*PSd(i) )
        x3 = HHdd(i)*delta*PS(i)
        xsum = x1 + x2 + x3
        phirdd = phirdd + c3n(i)*xsum
      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirt.

      phirt = 0.0d0

      DO i =1,7
        phirt = phirt + c1n(i)*c1t(i)*deldi(i)*tautm1(i)
      ENDDO

      DO i = 8,51
        phirt = phirt + c1n(i)*c1t(i)*deldi(i)*tautm1(i)*edci(i)
      ENDDO

      DO i = 52,54
        phirt = phirt + c2n(i)*deldi(i)*tauti(i)*exxi(i)*tbtg(i)
      ENDDO

      DO i = 55,56
        xsum = HHt(i)*PS(i) + HH(i)*PSt(i)
        phirt = phirt + c3n(i)*delta*xsum
      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirtt.

      phirtt = 0.0d0

      DO i =1,7
        phirtt = phirtt + c1n(i)*c1t(i)*(c1t(i) - 1.0d0)
     $  *deldi(i)*tautm2(i)
      ENDDO

      DO i = 8,51
        phirtt = phirtt + c1n(i)*c1t(i)*(c1t(i) - 1.0d0)*deldi(i)
     $  *tautm2(i)*edci(i)
      ENDDO

      DO i = 52,54
        x1 = tbtg(i)**2
        x2 = -( c2t(i)/tausq ) - 2.0d0*c2beta(i)
        xsum = x1 + x2
        phirtt = phirtt + c2n(i)*deldi(i)*tauti(i)*exxi(i)*xsum
      ENDDO

      DO i = 55,56
        xsum = HHtt(i)*PS(i) + 2.0d0*HHt(i)*PSt(i) + HH(i)*PStt(i)
        phirtt = phirtt + c3n(i)*delta*xsum
      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate phirdt.

      phirdt = 0.0d0

      DO i =1,7
        phirdt = phirdt + c1n(i)*c1d(i)*c1t(i)*deldm1(i)*tautm1(i)
      ENDDO

      DO i = 8,51
        phirdt = phirdt + c1n(i)*c1t(i)*deldm1(i)*tautm1(i)
     $  *dmci(i)*edci(i)
      ENDDO

      DO i = 52,54
        xsum = dade(i)*tbtg(i)
        phirdt = phirdt + c2n(i)*deldi(i)*tauti(i)*exxi(i)*xsum
      ENDDO

      DO i = 55,56
        x1 = HH(i)*( PSt(i) + delta*PSdt(i) )
        x2 = delta*HHd(i)*PSt(i)
        x3 = HHt(i)*( PS(i) + delta*PSd(i) )
        x4 = HHdt(i)*delta*PS(i)
        xsum = x1 + x2 + x3 + x4
        phirdt = phirdt + c3n(i)*xsum
      ENDDO

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Calculate thermodynamic functions.

c     Helmholtz energy. The value is in kJ/kg-K.

      ax = rtx*( phi0 + phir )

c     Pressure. The value here is in kPa.

      px = rho*rtx*( 1.0d0 + delta*phird )

c     Internal energy.

      ux = rtx*tau*( phi0t + phirt )

c     Entropy.

      sx = rcnstw*( tau*(phi0t + phirt) - phi0 - phir )

c     Enthalpy.

      hx = rtx*( 1.0d0 + tau*(phi0t + phirt) + delta*phird )

c     Gibbs energy.

      gx = rtx*( 1.0d0 + phi0 + phir + delta*phird )

c     Alternate formulas for the Gibbs energy.

c     gx = hx - tempk*sx
c     gx = ax + hx - ux

c     Volume.

      vx = 1.0d0/rho

c     Isochoric heat capacity.

      cvx = -rcnstw*tausq*( phi0tt + phirtt )

c     Isobaric heat capacity.

      x1 = ( 1.0d0 + delta*phird - delta*tau*phirdt )**2
      x2 = 1.0d0 + 2.0d0*delta*phird + delsq*phirdd
      IF (DABS(x2) .GT. 1.0d-15) THEN
        cpx = cvx + rcnstw*(x1/x2)
      ELSE
        cpx = 1.0d+100
      ENDIF

c     Speed of sound.

      x1 = ( 1.0d0 + delta*phird - delta*tau*phirdt )**2
      x2 = tausq*( phi0tt + phirtt )
      x3 = x1/x2
      xxt = rtx*( 1.0d0 + 2.0d0*delta*phird
     $ + delsq*phirdd - x3 )
      IF (xxt .GT. 0.0d0) THEN
        wx = sqrt(xxt)
      ELSE
        wx = 0.0d0
      ENDIF

c     Joule-Thomsen coefficient.

      x1 = delta*phird + delsq*phirdd + delta*tau*phirdt
      x2 = ( 1.0d0 + delta*phird - delta*tau*phirdt )**2
      x3 = ( phi0tt + phirtt )*( 1.0d0 + 2.0d0
     $ *delta*phird + delsq*phirdd )
      mux = ( - x1/( x2 - tausq*x3 ) )/(rcnstw*rho)

c     Isothermal throttling coefficient.

      x1 = 1.0d0 + delta*phird - delta*tau*phirdt
      x2 = 1.0d0 + 2.0d0*delta*phird + delsq*phirdd
      dtx = ( 1.0d0 - ( x1/x2 ) )/rho

c     Isentropic temperature-pressure coefficient.

      x1 = 1.0d0 + delta*phird - delta*tau*phirdt
      x2 = x1**2
      x3 = ( phi0tt + phirtt )*( 1.0d0 + 2.0d0*delta*phird
     $ + delsq*phirdd )
      bsx = ( x1/( x2 - tausq*x3 ) )/(rcnstw*rho)

c     Derivative of pressure with respect to delta
c     (needed to perform Newton-Raphson iteration
c     to matched desired pressure). Here the value
c     is in kPa. Recall that:
c     px = rho*rtx*( 1.0d0 + delta*phird )

      pdx = ( px/delta ) + delta*rhocr*rtx*( phird + delta*phirdd )

c     Derivative of pressure with respect to tau (needed to calculate
c     the thermal expansion coefficient). Here the value is again
c     in kPa.

      ptx = ( -px/tau ) + px*delta*phirdt/(1.0d0 + delta*phird)

c     Compressibility. Here the value is in /Kpa.

      ktx = 1.0d0/(delta*pdx)

c     Calculate ztx (needed to calculate viscosity). Note: pcr here
c     is in MPa, but pdx is still in /kPa. Hence there is the need
c     to include the factor of 1000 kPa/MPa. Note that ztx itself
c     is dimensionless, so pressure unit correction needs to be made
c     here.

      ztx = 1000.0d0*pcr/pdx

c     An alternative formula is: ztx = 1000.0d0*delta*pcr*ktx
c     This can be useful for getting zeta from a calculator that gives
c     the compressibility (but not zeta).

c     Thermal expansion coefficient (thermal expansivity).
c     This calculation is based on the Maxwell relation:
c     (del P/del T) at const V = alpha/kappa

      avx = ktx*ptx*( -tau/tempk )

c     Parts needed to calculate residuals and Jacobian elements
c     if refining saturation properties at fixed temperature.

c     Helmholtz energy.

c     ax = rtx*( phi0 + phir )

      adx = rtx*( phi0d + phird )

c     Gibbs energy.

c     gx = rtx*( 1.0d0 + phi0 + phir + delta*phird )

      gdx = rtx*( phi0d + 2.0d0*phird + delta*phirdd )

      IF (qwrphi) THEN
        IF (.NOT.qsilent) THEN

          WRITE(noutpt,1120) tempk, rho, delta, tau

 1120     FORMAT(/7x,'tempk = ',e16.9,3x,'rho   = ',e16.9,
     $    /7x,'delta = ',e16.9,3x,'tau   = ',e16.9,//)

          WRITE(noutpt,1130) phi0, phir, phi0d, phird, phi0dd, phirdd,
     $    phi0t, phirt, phi0tt, phirtt, phi0dt, phirdt

 1130     FORMAT(/12x,'phi functions',/
     $    /5x,'phi0   = ',e16.9,5x,'phir   = ',e16.9,
     $    /5x,'phi0d  = ',e16.9,5x,'phird  = ',e16.9,
     $    /5x,'phi0dd = ',e16.9,5x,'phirdd = ',e16.9,
     $    /5x,'phi0t  = ',e16.9,5x,'phirt  = ',e16.9,
     $    /5x,'phi0tt = ',e16.9,5x,'phirtt = ',e16.9,
     $    /5x,'phi0dt = ',e16.9,5x,'phirdt = ',e16.9,//)

        ENDIF
      ENDIF

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Unit conversions. Direct evaluation of the IAPWS-95
c     model equations gives some of the results in in odd
c     units. Conversions are therefore desirable for most
c     usage. Prior to conversion, the units are determined
c     by the units adopted by the model for temperature (K),
c     density (kg/m3), and the gas constant (kJ/kg/K).
c     See Table 3.

c       Note that density is input as kg/m3.
c       Volume is therefore in m3/kg.

c       Pressure comes out in kPA. Divide by 1000 to obtain
c       pressure in MPa. One could divide instead by 100 to
c       obtain pressure in bars. Here we will use MPa.
c       The pdx (the partial derivative of pressure with respect
c       to delta) and ktx (the isothermal compressibility)
c       must also be corrected to be consistent with
c       pressure in MPa.

c       Internal energy, enthalpy, and Helmholtz energy come
c       out in kJ/kg.

c       Entropy and the two heat capacity functions come out
c       in kJ/kg/K.

c       Gibbs energy is therefore in kJ/kg.

c       Speed of sound comes out in sqrt(kJ/kg). Multiply by
c       sqrt(1000) to obtain the speed of sound in m/s.

c       Joule-Thomson coefficient comes out in K-m3/kJ.
c       These units are equivalent to the usual K/MPa.

c       Isothermal throttling coefficient comes out in m3/kg.
c       Divide by 1000 to obtain the result in the usual
c       kJ/kg/MPa.

c       Isentropic temperature-pressure coefficient comes out
c       K-m3/kJ (the same units as for the Joule-Thomson coefficient).
c       These units are equivalent to the usual K/MPa.

      px = 0.001d0*px
      pdx = 0.001d0*pdx
      ptx = 0.001d0*ptx
      ktx = 1000.0d0*ktx

      wx = wx*sqrt(1000.0d0)
      dtx = 0.001d0*dtx

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE APXSCT(noutpt, nttyo, qprnt1, qsilent, rhocr, rhosv,
     $ rhosl, tempk, tcr, pcr, psat, psatt, qerr)

c     Evaluate the pressure (psat) as a function of temperature along
c     the vapor-liquid equilibrium curve, using equation 2.5 of Wagner
c     and Pruss (2002). Also calculate the densities of the liquid and
c     vapor phases using equations 2.6 and 2.7 from the same source.
c     Results are not fully consistent with IAPWS-95 to high precision,
c     but may serve as close approximations or starting values for
c     refinement.

c       psat = saturation pressure (MPa)
c       rhosv = density of vapor (kg/m3)
c       rhosl = density of liquid (kg/m3)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qprnt1, qsilent, qerr

      INTEGER nttyo, noutpt

      REAL(8) delta, rhosv, rhosl, rhocr, tempk, tcr, pcr

      REAL(8) psat, psatt

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

c     Model coefficients.

      REAL(8) c1a(1:6), c1b(1:6), c1c(1:6)
      REAL(8) e1b(1:6), e1c(1:6)

      DATA c1a(1:6) / -7.85951783d0, 1.84408259d0, -11.7866497d0,
     $ 22.6807411d0, -15.9618719d0, 1.80122502d0 /

      DATA c1b(1:6) / 1.99274064d0, 1.09965342d0, -0.510839303d0,
     $ -1.75493479d0, -45.5170352d0, -6.74694450d+05 /

      DATA c1c(1:6) / -2.03150240d0, -2.68302940d0, -5.38626492d0,
     $ -17.2991605d0, -44.7586581d0, -63.9201063d0 /

      LOGICAL qfirst

      SAVE qfirst, e1b, e1c

      DATA qfirst / .TRUE. /

      REAL(8) varth, x1, x2

      INTEGER i

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF (qprnt1 .AND. .NOT.qsilent) THEN

        WRITE(nttyo,1100) tempk
        WRITE(noutpt,1100) tempk

 1100   FORMAT(/3x,'APXSCT: Input tempk= ',f8.4' (K)',/)

      ENDIF

c     If this is the first call, calculate the exponents for the
c     density equations.

      IF (qfirst) THEN

        e1b(1) = 1.0d0/3.0d0
        e1b(2) = 2.0d0/3.0d0
        e1b(3) = 5.0d0/3.0d0
        e1b(4) = 16.0d0/3.0d0
        e1b(5) = 43.0d0/3.0d0
        e1b(6) = 110.0d0/3.0d0

        e1c(1) = 2.0d0/6.0d0
        e1c(2) = 4.0d0/6.0d0
        e1c(3) = 8.0d0/6.0d0
        e1c(4) = 18.0d0/6.0d0
        e1c(5) = 37.0d0/6.0d0
        e1c(6) = 71.0d0/6.0d0

        qfirst = .FALSE.

      ENDIF

c     Zero parameters to be calculated.

      psat = 0.0d0
      rhosv = 0.0d0
      rhosl = 0.0d0

c     Check to see that the temperature is in the allowed range.

      qerr = .FALSE.
      IF (tempk.LT.273.15d0 .OR. tempk.GT.tcr) THEN
        qerr = .TRUE.
        GO TO 999
      ENDIF

      varth = 1.0d0 - (tempk/tcr)

c     Saturation pressure.

      x1 = c1a(1)*varth + c1a(2)*varth**1.5d0 + c1a(3)*varth**3.0d0
     $ + c1a(4)*varth**3.5d0 + c1a(5)*varth**4.0d0
     $ + c1a(6)*varth**7.5d0

      x2 = ( tcr/tempk )*x1

      psat = pcr*DEXP(x2)

c     Derivative of saturation pressure with respect to temperature.

      x1 = c1a(1) + 1.5d0*c1a(2)*varth**0.5d0
     $ + 3.0d0*c1a(3)*varth**2.0d0 + 3.5d0*c1a(4)*varth**2.5d0
     $ + 4.0d0*c1a(5)*varth**3.0d0 + 7.5d0*c1a(6)*varth**6.5d0

      x2 = DLOG( psat/pcr ) + x1

      psatt = -( psat/tempk )*x2

c     Density of liquid.

      x1 = 1.0d0 + c1b(1)*varth**(1.0d0/3.0d0)
     $ + c1b(2)*varth**(2.0d0/3.0d0) + c1b(3)*varth**(5.0d0/3.0d0)
     $ + c1b(4)*varth**(16.0d0/3.0d0) + c1b(5)*varth**(43.0d0/3.0d0)
     $ + c1b(6)*varth**(110.0d0/3.0d0)

      rhosl = rhocr*x1

c     Density of vapor.

      x1 = c1c(1)*varth**(2.0d0/6.0d0) + c1c(2)*varth**(4.0d0/6.0d0)
     $ + c1c(3)*varth**(8.0d0/6.0d0) + c1c(4)*varth**(18.0d0/6.0d0)
     $ + c1c(5)*varth**(37.0d0/6.0d0) + c1c(6)*varth**(71.0d0/6.0d0)

      rhosv = rhocr*DEXP(x1)

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE GAUSSE(noutpt, nttyo, qerr, qprnt3, kdim,
     $ kmax, aamatr, alpha, deltas)

c     Solve the matrix equation aamatr*deltas = alpha using Gaussian
c     elimination with partial pivoting.

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

c     Calling sequence variables.

      LOGICAL qerr, qprnt3

      INTEGER nttyo, noutpt

      INTEGER kdim, kmax

      REAL(8) aamatr(kmax,kmax), alpha(kmax), deltas(kmax)

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c     Local variables.

      INTEGER i, j, k, kk, kopt

      REAL(8) aamx, abmx, amx, axx

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      qerr = .FALSE.

      DO k = 1,kdim - 1

c       Find the optimal pivot row in rows k through kdim
c       Unless row k happens to be the optimal row, exchange
c       it for the optimal row.

        kopt = k
        aamx = DABS(aamatr(k,k))

        DO kk = k + 1, kdim
          abmx = DABS(aamatr(kk,k))
          IF (abmx .GT. aamx) THEN
            kopt = kk
            aamx = abmx
          ENDIF
        ENDDO

        IF (aamx .LE. 1.0d-14) THEN

          WRITE(nttyo,1120)
          WRITE(noutpt,1120)

 1120     FORMAT(/6x,'ERROR (GAUSSE): Have a singular matrix',/)

          qerr = .FALSE.
          GO TO 999
        ENDIF

        IF (kopt .GT. k) THEN

c         Exchange rows k and kopt in the extended matrix
c         (don't forget the right-hand-side part).

          DO j =  k,kdim
            amx = aamatr(k,j)
            aamatr(k,j) = aamatr(kopt,j)
            aamatr(kopt,j) = amx
          ENDDO
          axx = alpha(k)
          alpha(k) = alpha(kopt)
          alpha(kopt) = axx

          IF (qprnt3) THEN

            WRITE (noutpt,1130)

 1130       FORMAT(/3x,'Pivot (row exchange)')

            WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $      alpha(i), i = 1,kdim)

 1140       FORMAT(/6x,'Jacobian and residual',
     $      /9x,1pe12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $      /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5,
     $      /9x,e12.5,3x,e12.5,3x,e12.5,' | ',e12.5)

          ENDIF

        ENDIF

        amx = aamatr(k,k)
        aamatr(k,k) = 1.0d0
        DO j =  k + 1,kdim
          aamatr(k,j) = aamatr(k,j)/amx
        ENDDO
        alpha(k) = alpha(k)/amx

        IF (qprnt3) THEN

          WRITE (noutpt,1150)

 1150     FORMAT(/3x,'Scale the pivot row for unit diagonal',
     $    ' entry')

          WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $    alpha(i), i = 1,kdim)

        ENDIF

        DO kk = k + 1,kdim

          axx = aamatr(kk,k)

          DO j = k,kdim
            aamatr(kk,j) = aamatr(kk,j) - axx*aamatr(k,j)
          ENDDO
          alpha(kk) = alpha(kk) - axx*alpha(k)

        ENDDO

        IF (qprnt3) THEN

          WRITE (noutpt,1160)

 1160     FORMAT(/3x,'Zero lower row entries in the pivot',
     $    ' column')

          WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $    alpha(i), i = 1,kdim)

        ENDIF

      ENDDO

c     Scale the last row for unit diagonal entry.

      amx = aamatr(k,k)
      aamatr(k,k) = 1.0d0
      DO j =  k + 1,kdim
        aamatr(k,j) = aamatr(k,j)/amx
      ENDDO
      alpha(k) = alpha(k)/amx

      IF (qprnt3) THEN
        WRITE (noutpt,1150)
        WRITE (noutpt,1140) ((aamatr(i,j), j = 1,kdim),
     $  alpha(i), i = 1,kdim)
      ENDIF

c     Find the solution vector.

      DO k = kdim,1,-1
        deltas(k) = alpha(k)
        DO kk = k + 1,kdim
          deltas(k) = deltas(k) - aamatr(k,kk)*deltas(kk)
        ENDDO
      ENDDO

c     deltas(3) = alpha(3)
c     deltas(2) = ( alpha(2) - aamatr(2,3)*deltas(3) )
c     deltas(1) = ( alpha(1) - aamatr(1,2)*deltas(2)
c    $ - aamatr(1,3)*deltas(3) )

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE OPENIN(noutpt,nttyo,ufiln,uform,ilu)

c     This subroutine opens an input type file. The file must already
c     exist. An unused logical unit number is obtained.

c     This subroutine is called by:

c       H2OI95

c-----------------------------------------------------------------------

c     Input:

c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c       ufiln  = the name of the file to open
c       uform  = the file format, 'formatted' or 'unformatted'

c     Output:

c       ilu     = the logical unit number of opened file

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

c     Note- the character length for the variables ufiln and uform
c     can not be exactly specified.

      INTEGER ilu,noutpt,nttyo

      CHARACTER(LEN=*) ufiln,uform

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j2,nerr

      INTEGER ILNOBL

      LOGICAL qex

      CHARACTER(LEN=8) uformo

c-----------------------------------------------------------------------

      j2 = ILNOBL(ufiln)

c     Check to make sure the file exists.

      INQUIRE(file=ufiln,exist=qex,formatted=uformo)
      IF (.NOT.qex) THEN
        IF (noutpt .GT. 0) WRITE (noutpt,1000) ufiln(1:j2)

        WRITE (nttyo,1000) ufiln(1:j2)

 1000   FORMAT(/6x,'ERROR (OPENIN): The file "',a,'"',
     $  /9x,'does not exist.')

        STOP
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     Check the file format.

      IF (uform.EQ.'formatted' .AND. uformo.EQ.'no') THEN

        IF (noutpt .GT. 0) WRITE (noutpt,1010) ufiln(1:j2)

        WRITE (nttyo,1010) ufiln(1:j2)

 1010   FORMAT(/6x,'ERROR (OPENIN): The file "',a,'"',
     $  /9x,'should be formatted, but it is not.')

        STOP
      ENDIF

      IF (uform.EQ.'unformatted' .AND .uformo.EQ.'yes') THEN

        IF (noutpt .GT. 0) WRITE (noutpt,1020) ufiln(1:j2)

        WRITE (nttyo,1020) ufiln(1:j2)

 1020   FORMAT(/6x,'ERROR (OPENIN): The file "',a,'"',
     $  /9x,'should be unformatted, but it is not.')

        STOP
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     Get the next available logical unit number.

      CALL GETLU(ilu,nerr)
      IF (nerr .NE. 0) THEN
        IF (noutpt .GT. 0) WRITE (noutpt,1050) ufiln(1:j2)

        WRITE (nttyo,1050) ufiln(1:j2)

 1050   FORMAT(/6x,'ERROR (OPENIN): No logical unit number',
     $  /9x,'is available for the file "',a,'".')

        STOP
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      OPEN(ilu,file=ufiln,form=uform,status='old',err=10)
      GO TO 999

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   10 CONTINUE
      IF (noutpt .GT. 0) WRITE (noutpt,1060) ufiln(1:j2)

      WRITE (nttyo,1060) ufiln(1:j2)

 1060 FORMAT(/6x,"ERROR (OPENIN): Can't open the file",
     $ /9x,'"',a,'".')

      STOP

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE OPENOU(noutpt,nttyo,ufiln,uform,nrecl,ilu)

c     This subroutine opens an output type file. If a file of the
c     same name already exists, it is first destroyed. An unused
c     logical unit number is obtained.

c     This subroutine is called by:

c       H2OI95

c-----------------------------------------------------------------------

c     Input:

c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c       ufiln  = the name of the file to open
c       uform  = the file format, 'formatted' or 'unformatted'
c       nrecl  = the record length (number of characters per line)

c     Output:

c       ilu     = the logical unit number of opened file

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      INTEGER noutpt,nttyo

      INTEGER ilu,nrecl

c     Note- the character length for the variables ufiln and uform
c     can not be exactly specified.

      CHARACTER(LEN=*) ufiln,uform

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j2,nerr

      INTEGER ILNOBL

      LOGICAL qex

      CHARACTER(LEN=8) ustat

c-----------------------------------------------------------------------

      j2 = ILNOBL(ufiln)

c     See if a file of the same name already exists. If so,
c     destroy it. This makes the logical unit number available.
c     If a file of the same name does not exist, get the next
c     available logical unit number.

      INQUIRE(file=ufiln,exist=qex)
      IF (qex) THEN
        ustat = 'old'
        CALL GETLU(ilu,nerr)
        IF (nerr .NE. 0) THEN

          IF (noutpt .GT. 0) WRITE (noutpt,1000) ustat,ufiln(1:j2)

          WRITE (nttyo,1000) ustat,ufiln(1:j2)

 1000     FORMAT(/6x,'ERROR (OPENOU): No logical unit number',
     $    /9x,'is available to open the ',a3,' file "',a,'".')

          STOP
        ENDIF
        OPEN(ilu,file=ufiln,status=ustat,err=10)
        CLOSE(ilu,status='delete',err=15)
      ELSE
        ustat = 'new'
        CALL GETLU(ilu,nerr)
        IF (nerr .NE. 0) THEN
          IF (noutpt .GT. 0) WRITE (noutpt,1000) ustat,ufiln(1:j2)
          WRITE (nttyo,1000) ustat,ufiln(1:j2)
          STOP
        ENDIF
      ENDIF

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     Open the new file.

      IF (nrecl .GT. 0) THEN

c       Use the specified record length.

        OPEN(ilu,file=ufiln,form=uform,status='new',recl=nrecl,err=10)
      ELSE

c       The record length is not specified. Open the file at the
c       default record length.

        OPEN(ilu,file=ufiln,form=uform,status='new',err=10)
      ENDIF
      GO TO 999

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   10 IF (noutpt .GT. 0) WRITE (noutpt,1010) ustat,ufiln(1:j2)

      WRITE (nttyo ,1010) ustat,ufiln(1:j2)

 1010 FORMAT(/6x,"ERROR (OPENOU): Can't open the ",a3,' copy',
     $ /9x,'of the file "',a,'".')

      STOP

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

   15 IF (noutpt .GT. 0) WRITE (noutpt,1020) ufiln(1:j2)

      WRITE (nttyo,1020) ufiln(1:j2)

 1020 FORMAT(/6x,"ERROR (OPENOU): Can't delete the old copy",
     $ /9x,'of the file "',a,'".')

      STOP

c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE GETLU(nlu,nerr)

c     This subroutine finds a currently unused unit number.

c     This subroutine is called by:

c       OPENIN.f
c       OPENOU.f

c-----------------------------------------------------------------------

c     Input:

c       None

c     Output:

c       nlu    = first currently unused unit number
c       nerr   = error flag:
c                  = 0   Okay
c                  = 1   Error

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      INTEGER nlu,nerr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER iumax,iumin

      LOGICAL qopen

c-----------------------------------------------------------------------

      data iumax,iumin /40,0/

c-----------------------------------------------------------------------

      nerr = 0

c     Loop through all valid file numbers, beginning with the largest.

      DO nlu = iumax,iumin,-1
        INQUIRE(unit=nlu,opened=qopen)
        IF (.NOT.qopen) GO TO 999
      ENDDO

      nerr = 1

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      INTEGER FUNCTION IFNOBL(ustr)

c     This subroutine finds the position of the first non-blank
c     character in the string ustr.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       IFNOBL = the position of the first non-blank character

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      CHARACTER(LEN=*) ustr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j,nchars

c-----------------------------------------------------------------------

c     Get the length of the string variable.

      nchars = len(ustr)

c     Find the first non-blank character.

      IFNOBL = 0
      DO j = 1,nchars
        IF (ustr(j:j) .NE. ' ') THEN
          IFNOBL = j
          GO TO 999
        ENDIF
      ENDDO

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      INTEGER FUNCTION ILNOBL(ustr)

c     This subroutine finds the position of the last non-blank character
c     in the string ustr.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       ILNOBL = the position of the first non-blank character

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      CHARACTER(LEN=*) ustr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j,nchars

c-----------------------------------------------------------------------

c     Get the length of the string variable.

      nchars = len(ustr)

c     Find the first non-blank character.

      ILNOBL = 0
      DO j = nchars,1,-1
        IF (ustr(j:j) .NE. ' ') THEN
          ILNOBL = j
          GO TO 999
        ENDIF
      ENDDO

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE LEJUST(ustr)

c     This subroutine left-justifies the non-blank portion of the string
c     ustr.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       ustr   = the output string variable

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      CHARACTER(LEN=*) ustr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j,jj,jbl,j1,nchars

      INTEGER IFNOBL

c-----------------------------------------------------------------------

c     Get the length of the string variable.

      nchars = len(ustr)

c     Get the position of the first non-blank character and the number
c     of blanks on the left-hand-side.

      j1 = IFNOBL(ustr)
      jbl = j1 - 1

      IF (jbl .GT. 0) THEN
        DO jj = j1,nchars
          j = jj - jbl
          ustr(j:j) = ustr(jj:jj)
        ENDDO
        DO j = nchars - jbl + 1,nchars
          ustr(j:j) = ' '
        ENDDO
      ENDIF

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE LOCASE(ustr)

c     This subroutine converts a string from upper case to lower case.

c     This subroutine is called by:

c       Any

c-----------------------------------------------------------------------

c     Input:

c       ustr   = the input string variable

c     Output:

c       ustr   = the output string variable

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      CHARACTER(LEN=*) ustr

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER idel,j,nchars

      CHARACTER(LEN=1) u1

c-----------------------------------------------------------------------

      idel = ichar('A') - ichar('a')
      IF (idel .NE. 0) THEN
        nchars = len(ustr)
        DO j = 1,nchars
          u1 = ustr(j:j)
          IF (u1.GE.'A' .AND. u1.LE.'Z')
     $    ustr(j:j) = char(ichar(u1) - idel)
        ENDDO
      ENDIF

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      REAL(8) FUNCTION PMELT(tempk, press, noutpt, nttyo, ier)

c     This function gets the melting pressure of ice as a function
c     of temperature. This model is described in IAPWS R14-08(2011),
c     The International Association for the Properties of Water and
c     Steam, 2011, Revised Release on the Pressure along the Melting
c     and Sublimation Curves of Ordinary Water Substance. This
c     document may be found at:
c     http://www.iapws.org/relguide/MeltSub.html

c     Five ice phases are covered here. The melting pressure is
c     not a single-valued function of temperature as there is some
c     overlap in the temperature ranges of the individual phases.
c     There is no overlap in the temperature ranges of Ices III,
c     V, VI, and VII, which together span the range 251.165-715K.
c     The melting pressure is continuous and monotonically increasing
c     over this range, albeit with discontinuities in slope at
c     the triple points where two ice phases and liquid are in
c     equilibrium. The problem comes in with Ice Ih, whose
c     temperature range completely overlaps that of Ice III and
c     partially overlaps that of Ice V. For a temperature In the
c     range for Ice Ih, there are two possible melting pressures.

c     The possible ambiguity here in the meaning of melting pressure
c     is not present if the temperature is greater than or equal to
c     the triple point temperature of 273.16K, or if the pressure is
c     greater than or equal to 208.566 MPa (the triple point pressure
c     for Ice Ih-Ice III-liquid). If neither of these conditions
c     are satisfied, then the Ice Ih-liquid curve will be used.
c     To deal with the pressure condition noted above, this function
c     assumes that an actual pressure is specified.

c     IAPWS R14-08 makes reference to the melting temperature for
c     pressure between the triple point pressure and 300 MPa, which
c     includes the region of ambiguity. However, it does not discuss
c     the ambiguity or make any suggestion of how to deal with it.
c     presents the overall model.

c     This function is called by:

c       H2OI95.f
c       TKMELT.f

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      REAL(8) tempk, press

      INTEGER ier

      INTEGER noutpt, nttyo

c-----------------------------------------------------------------------

c     Local variable declarations.

      LOGICAL qwr

      INTEGER i

      INTEGER ILNOBL

      REAL(8) cfa(1:3), cfb(1:3)

      REAL(8) theta, pimelt, px

      INTEGER j2, j3

      CHARACTER(LEN=24) ux24a, ux24b

      DATA qwr / .FALSE. /

c-----------------------------------------------------------------------

c     Coefficients for calculating the melting pressure of Ice Ih.

      DATA cfa(1:3) / 0.119539337d+07, 0.808183159d+05,
     $ 0.333826860d+04 /

      DATA cfb(1:3) / 0.300000d+01, 0.257500d+02, 0.103750d+03 /

c-----------------------------------------------------------------------

      PMELT = 0.0d0
      pimelt = 0.0d0

c     Set error flag.

      ier = 0

      IF (press .LE. 208.566d0) THEN

        IF (tempk.GE.251.165d0 .AND. tempk.LE.273.16d0) THEN

c         Ice Ih.

          theta = tempk/273.16d0

          pimelt = 1.0d0
          DO i = 1,3
            pimelt = pimelt + cfa(i)*( 1 - theta**cfb(i) )
          ENDDO

          PMELT = pimelt*611.647d-06
        ELSE

          ier = 1

        ENDIF

      ELSE

        IF (tempk.GE.251.165d0 .AND. tempk.LT.256.164d0) THEN

c         Ice III.

          theta = tempk/251.165d0

          pimelt = 1 - 0.299948d0*( 1 - theta**60 )

          PMELT = pimelt*208.566d0

        ELSEIF (tempk.GE.256.164d0 .AND. tempk.LT.273.31d0) THEN

c         Ice V.

          theta = tempk/256.164d0

          pimelt = 1 - 1.18721d0*( 1 - theta**8 )

          PMELT = pimelt*350.1d0

        ELSEIF (tempk.GE.273.31d0 .AND. tempk.LT.355.0d0) THEN

c         Ice VI.

          theta = tempk/273.31d0

          pimelt = 1 - 1.07476d0*( 1 - theta**4.6d0 )

          PMELT = pimelt*632.4d0

        ELSEIF (tempk.GE.355.0d0 .AND. tempk.LT.715.0d0) THEN

c         Ice VII.

          theta = tempk/355.0d0

          px =   0.173683d+01*( 1 - theta**(-1) )
     $         - 0.544606d-01*( 1 - theta**5 )
     $         + 0.806106d-07*( 1 - theta**22 )

          pimelt = DEXP(px)
          PMELT = pimelt*2216.0d0

        ELSEIF (tempk.GT.715.0d0 .AND. tempk.LE.2000.0d0) THEN

c         This is out of range. Ice VII, extrapolated.

          theta = tempk/355.0d0

          px =   0.173683d+01*( 1 - theta**(-1) )
     $         - 0.544606d-01*( 1 - theta**5 )
     $         + 0.806106d-07*( 1 - theta**22 )

          pimelt = DEXP(px)
          PMELT = pimelt*2216.0d0

        ELSE

          ier = 1

        ENDIF

      ENDIF

      IF (ier .GT. 0) THEN

        WRITE(ux24a,'(f16.9)') tempk

        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(nttyo,1010) ux24a(1:j2)
        WRITE(noutpt,1010) ux24a(1:j2)

 1010   FORMAT(/6x,'WARNING (PMELT): The temperature of ',a,'K is',
     $  ' out of range.')

      ENDIF

      IF (qwr) THEN

c       Write PMELT on the standard output.

        WRITE(ux24a,'(f16.9)') tempk
        CALL LEJUST(ux24a)
        j2 = ILNOBL(ux24a)

        WRITE(ux24b,'(f16.9)') PMELT
        CALL LEJUST(ux24b)
        j3 = ILNOBL(ux24b)

        WRITE(nttyo,1050) ux24a(1:j2), ux24b(1:j3)

 1050   FORMAT(/6x,'PMELT: tempk = ',a,'K, PMELT = ',a,' MPa')

      ENDIF

  999 CONTINUE

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      REAL(8) FUNCTION SURFTE(tempk, tcr, q1phsc, qpsatc, qtsatc,
     $ nttyo)

c     This function calculates the surface tension of water on
c     the saturation curve. The model equation is a simple function
c     of temperature. The model is described in "Revised Release
c     on Surface Tension of Ordinary Water Substance" (IAPWS
c     R1-76, 2014). . The corresponding journal reference is
c     Vargaftik N. B., Volkov B. N., and Voljak L. D. (1983)
c     International Tables of the Surface Tension of Water.
c     J. Phys. Chem. Ref. Data 12, 817-820. The equation
c     referenced below is from IAPWS R1-76.

c     This function is called by:

c       H2OI95.f

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      LOGICAL q1phsc, qpsatc, qtsatc

      REAL(8) tempk, tcr

      INTEGER nttyo

      INTEGER ILNOBL

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j2

      CHARACTER(LEN=24) ux24

      REAL(8) stx, xtau

      REAL(8) bcap, bsma, xmu

c     The following is a local hard-wired switch used in code
c     verification of surface tension calculations.

      LOGICAL qwr1

      DATA qwr1 / .FALSE. /

c-----------------------------------------------------------------------

c     Model constants.

      DATA bcap, bsma, xmu / 235.8d0, -0.625d0, 1.256d0 /

c-----------------------------------------------------------------------

      stx = 0.d0

      IF (q1phsc .OR. qpsatc .OR. qtsatc) THEN

        xtau = 1.0d0 - (tempk/tcr)

        stx = bcap*(xtau**xmu)*(1.0d0 + bsma*xtau)

      ENDIF

      IF (qwr1) THEN

c       Write stx on the standard output.

        WRITE(ux24,'(f16.9)') stx
        CALL LEJUST(ux24)
        j2 = ILNOBL(ux24)

        WRITE(nttyo,1250) ux24(1:j2)

 1250   FORMAT(/6x,'SURFTE: stx = ',a)

      ENDIF

      SURFTE = stx

      END

c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SUBROUTINE DIELEC(rho, rhocr, tempk, tcr, rcnstw, xmcapw,
     $ avx, ktx, ptx, nttyo, epsx, deptp, deppt, agme,
     $ agm10, aphi, adhv, adhh, ahrt, bgm, bdhv, bdhh,
     $ bhrt)

c     This subroutine calculates the dielectric constant of water
c     (epsx) as a function of temperature and density. The model is
c     described in "Release on the Static Dielectric Constant of
c     Ordinary Water Substance for Temperatures from 238 K to 873 K
c     and Pressures up to 1000 MPa" (IAPWS R8-97, 1997). The
c     corresponding journal reference is Fernandez D. P.,
c     Goodwin A. R. H., Lemmon E. W., Levelt Sengers J. M. H.,
c     and Williams R. C. (1997) A Formulation for the Static
c     Permittivity of Water and Steam at Temperatures from 238 K
c     to 873 K at Pressures up to 1200 MPa, Including Derivatives
c     and Debye-HÅckel Coefficients. J. Phys. Chem. Ref. Data 26,
c     1125-1166.  The equations used for the dielectric constant
c     here are taken directly from IAPWS R8-97.

c     This routine also calculates some closely related water
c     properties. Equations for The temperature and pressure
c     derivatives of the dielectric constant and the Debye-Huckel
c     "A" parameters are taken from Fernandez et al. (1997). These
c     equations are not sanctioned by the IAPWS. Equations for the
c     Debye-Huckel "B" parameters and some Born functions are taken
c     from Helgeson H. C. and Kirkham D. H. (1974) Theoretical
c     Prediction of the Thermodynamic Behavior of Aqueous Electrolytes
c     at High Pressures and Temperatures: II. Debye-Huckel Parameters
c     for Activity Coefficients and Relative Partial Molal Properties.
c     Am. J. Sci. 274, 1199-1251. The values of associated constants
c     are updated to values consistent with Fernandez et al. (1997).
c

c     This routine is called by:

c       H2OI95.f

c-----------------------------------------------------------------------

      IMPLICIT NONE

c-----------------------------------------------------------------------

c     Calling sequence variable declarations.

      INTEGER nttyo

      REAL(8) rho, rhocr, tempk, tcr, rcnstw, xmcapw, avx, ktx, ptx

      REAL(8) epsx, deptp, deppt, agme, agm10, aphi, adhv, adhh, ahrt,
     $ bgm, bdhv, bdhh, bhrt

c-----------------------------------------------------------------------

c     Local variable declarations.

      INTEGER j2

      CHARACTER(LEN=24) ux24

      REAL(8) delta, tau, rhom, rhocrm, gha, pi, eps0, tfx228,
     $ acap, bcap

      REAL(8) al10, cx, cxcf, ttx, ttxt, xx, xxx, yyy, xsum, xsqrt,
     $ rtxm, rtxmv, bdhvx

      REAL(8) dghar, dghat, acap1, acap2, bcap1, ccap, x44b

      REAL(8) drpt, drtp, dptr, deprt, deptr

      REAL(8) dbgtp, dbgpt, dleptp, dleppt

      INTEGER ILNOBL

      REAL(8) ncap(1:12), jh(1:11)

      INTEGER h

      INTEGER ih(1:11)

      REAL(8) alpmmp, xmudp, xkbolt, avogad, elemch, xctesu,
     $ echesu, rhogc, xkblte, bgmx, bdhhx

c     The following is a local hard-wired switch used in code
c     verification of dielectric constant calculations.

      LOGICAL qwr1

      DATA qwr1 / .FALSE. /

c-----------------------------------------------------------------------

c     Model constants.

      DATA ncap(1:12) /  0.978224486826d0,   -0.957771379375d0,
     $                   0.237511794148d0,    0.714692244396d0,
     $                  -0.298217036956d0,   -0.108863472196d0,
     $                   0.949327488264d-01, -0.980469816509d-02,
     $                   0.165167634970d-04,  0.937359795772d-04,
     $                  -0.123179218720d-09,  0.196096504426d-02 /

      DATA ih(1:11) / 1, 1, 1, 2, 3, 3, 4, 5, 6, 7, 10 /

      DATA jh(1:11) / 0.25, 1.0d0, 2.5d0, 1.5d0, 1.5d0, 2.5d0,
     $ 2.0d0, 2.0d0, 5.0d0, 0.5d0, 10.0d0 /

c     The following constant values (alpmmp through xctesu) are from
c     Table 3 of Fernandez et al. (1997). A subset of these constants
c     is used in IAPWS R8-92 (1997).

      DATA alpmmp, xmudp, xkbolt, avogad
     $ / 1.636d-40, 6.138D-30, 1.380658d-23, 6.0221367d+23 /

c     The elementary charge (elemch) is in Coulombs (C).

      DATA elemch / 1.6021773310d-19 /

c     The following factor (xctesu) converts from Coulombs (C) to esu.
c     In Table 3 of Fernandez et al. (1997) this is buried in the
c     expression for calculating the permittivity of free space, eps0.

      DATA xctesu  / 2997924580.0d0 /

c     Slightly different values for some of these constants were used
c     by Helgeson and Kirkham (1974b).

c     Relevant 2019 redefinitions (for information only):

c       The elementary charge e is 1.602176634 x 10^-19 coulomb (C)
c       The Boltzmann constant k is 1.380649 x 10^-23 J K^-1
c       The Avogadro constant NA is 6.02214076 x 10^23 mol^-1

c-----------------------------------------------------------------------

      epsx = 0.0d0
      rhom = rho/xmcapw
      rhocrm = rhocr/xmcapw

c     Get the permittivity of free space (eps0).

      pi = 4*ATAN(1.0d0)

c     Note: the constant below is xctesu/10.

      cx = (299792458d0)**2
      xx = 4.0d-07*pi*cx
      eps0 = 1.0d0/xx

      delta = rhom/rhocrm
      tau = tcr/tempk

c     Get the Harris and Alder g factor (gha).

      tfx228 = (tempk/228.0d0) - 1.0d0

      gha = 1.0d0
      DO h = 1,11
        gha = gha + ncap(h)*(delta**(ih(h)))*(tau**(jh(h)))
      ENDDO
      ttx = tfx228**(-1.2d0)
      gha = gha + ncap(12)*delta*ttx

c     Get A (acap) and B (bcap).

      acap = ( avogad*(xmudp**2)*rhom*gha )/( eps0*xkbolt*tempk )
      bcap = ( avogad*alpmmp*rhom )/( 3*eps0 )

c     Get the dielectric constant (epsx).

      xsum = 9.0d0 + 2*acap + 18*bcap + acap**2 + 10*acap*bcap
     $ + 9*bcap**2
      xsqrt = SQRT(xsum)
      epsx = ( 1.0d0 + acap + 5*bcap + xsqrt )/( 4.0d0 - 4*bcap )

      IF (qwr1) THEN

c       Write epsx on the standard output.

        WRITE(ux24,'(f16.9)') epsx
        CALL LEJUST(ux24)
        j2 = ILNOBL(ux24)

        WRITE(nttyo,1250) ux24(1:j2)

 1250   FORMAT(/6x,'DIELEC: epsx = ',a)

      ENDIF

c     Get the Debye-Huckel A(phi) (aphi) and A(gamma) parameters
c     (agme, agme10). Units are kg^0.5 mol^-0.5 for all three.

      xxx = 2*pi*avogad*rhom*xmcapw
      yyy = elemch**2/( 4*pi*epsx*eps0*xkbolt*tempk )
      agme = SQRT(xxx)*(yyy**(1.5d0))

      al10 = LOG(10.0d0)
      agm10 = agme/al10

      aphi = agme/3.0d0

      IF (qwr1) THEN

c       Write agme on the standard output.

        WRITE(ux24,'(f16.9)') agme
        CALL LEJUST(ux24)
        j2 = ILNOBL(ux24)

        WRITE(nttyo,1260) ux24(1:j2)

 1260   FORMAT(/6x,'DIELEC: agme = ',a)

      ENDIF

c     Get the partial derivative of the Harris and Alder g factor
c     with respect to density at constant temperature.

      dghar = 0.0d0
      DO h = 1,11
        dghar = dghar
     $  + (ncap(h)*ih(h)/rhocrm)*(delta**(ih(h) - 1))*(tau**(jh(h)))
      ENDDO
      ttx = tfx228**(-1.2d0)
      dghar = dghar + (ncap(12)/rhocrm)*ttx

c     Get the partial derivative of the Harris and Alder g factor
c     with respect to temperature at constant density.

      dghat = 0.0d0
      DO h = 1,11
        dghat = dghat + ncap(h)*(delta**(ih(h)))
     $  *jh(h)*(tau**(jh(h)- 1))*(-tau/tempk)
      ENDDO
      ttxt = ( -1.2d0*tfx228**(-2.2d0) )/228.0d0
      dghat = dghat + ncap(12)*delta*ttxt

c     Get A1 (acap1), B1 (bcap1), C (ccap), and A2 (acap2).
c     Note: A1 is the partial of A with respect to rhom.
c     A2 is the partial of A with respect to temperature.
c     B1 is the partial of B with respect to rhom.
c     There is not B2 analogous to A2 (it would be zero).

      acap1 = (acap/rhom) + (acap/gha)*dghar
      acap2 = -(acap/tempk) + (acap/gha)*dghat

      bcap1 = bcap/rhom

      ccap = 9 + 2*acap + 18*bcap + acap**2 + 10*acap*bcap
     $ + 9*(bcap**2)

      cxcf = 0.5d0*(ccap**(-0.5d0))
      x44b = 4.0d0 - 4*bcap

c     Get the partial of eps with respect to rhom at
c     constant T (deprt).

      xx = 2*acap1 + 18*bcap1 + 2*acap*acap1
     $ + 10*(acap1*bcap + acap*bcap1) + 18*bcap*bcap1
      xxx = acap1 + 5*bcap1 + cxcf*xx
      deprt = (4*bcap1*epsx/x44b) + (xxx/x44b)

c     Get the partial of eps with respect to T at
c     constant rhom (deptr).

      xx = 2 + 2*acap + 10*bcap
      xxx= acap2 + cxcf*acap2*xx
      deptr = xxx/x44b

c     Get the partial of rhom with respect to P at constant T.
c     Use the compressibility (ktx) obtained from the EOS model.

      drpt = rhom*ktx

c     Get the partial of eps with respect to P at constant
c     T (deppt).

      deppt = deprt*drpt

c     Will need the partial of P with respect to T at constant
c     rhom (dptr). This is related to "ptx" as returned from
c     the evaluation of the EOS model. Here "ptx" is the
c     partial of P with respect to tau (tcr/tempk).

      dptr = -tau*ptx/tempk

c     Get the partial of eps with respect to T at constant
c     P (deptp).

      deptp = deptr - deprt*dptr*drpt

c     Need to use "molar" R below, where R = rcnstw*xmcapw
c     has units of kJ mol^-1 K^-1. Multiply by 1000 to get
c     R in J mol^-1 K^-1, for which the same numerical value
c     applies to R in cm^3 MPa mol^-1 K^-1.

      rtxm = rcnstw*xmcapw*tempk
      rtxmv = 1000*rtxm

c     Get the Debye-Huckel AV constant, adhv.
c     Units are cm^3 kg^1/2 mol^-3/2 (Fernandez et al., 1997).

      adhv = 2*aphi*rtxmv*( (3*deppt/epsx) - (drpt/rhom) )

c     Get the partial of rhom with respect to T at constant P.
c     Use the thermal expansivity (avx) obtained from the EOS model.

      drtp = -rhom*avx

c     Get the Debye-Huckel AH constant, adhh.
c     This constant is often tabulated as AH/RT. That is
c     expressed below as ahrt. Units for AH/RT are
c     are kg^1/2 mol^-1/2 (Fernandez et al., 1997).
c     Units for AH are kJ kg^1/2 mol^-3/2.

      xx = 1.0d0 + (tempk*deptp/epsx) - (tempk*drtp/(3*rhom))
      ahrt = -6*aphi*xx

      adhh = rtxm*ahrt

c     Get the Debye-Huckel B(gamma) parameter (bgm).
c     Units are kg^0.5 mol^-0.5 cm^-1 (Helgeson and
c     Kirkham, 1974b). The values tabulated by Helgeson
c     and Kirkham in their Table 2 are consistent with units
c     of kg^0.5 mol^-0.5 Angstrom^-1 (see below).

c     The tabulated values are described by these authors
c     as "B(gamma) in (kg^0.5 mol^-0.5 cm^-1) x 10^-8".
c     Technically these values correspond to
c     B(gamma) x 10^+8 for the units specified ("... cm-1")
c     in the table's title. There is a similar issue
c     with power of 10 multipliers in regard to BV in their
c     Table 12 and BH in their Table 4.

      echesu = xctesu*elemch

c     Get Boltzmann's constant in erg K^-1. The constant
c     represented by xkbolt is in J K^-1.

      xkblte = 1.0d+07*xkbolt
      rhogc = 0.001d0*rho

      xxx = 8*pi*avogad*rhogc*(echesu**2)
      yyy = 1000*epsx*xkblte*tempk

      bgmx = SQRT(xxx/yyy)

c     Get B(gamma) in units of kg^1/2 mol^-1/2 Angstrom^-1.
c     This is more traditional (e.g., Helgeson and Kirkham, 1974b)
c     for usage in calculating activity coefficients.
c     1 Angstrom(A) = 10^-8 cm. The multiplier below is
c     10^-8 cm Angstrom^-1.

      bgm = 1.0d-08*bgmx

c     Get the Debye-Huckel BV constant, bdhvx.
c     Units are cm^2 kg^1/2 mol^-3/2, equivalent to
c     cm^3 kg^1/2 mol^-3/2 cm^-1.

      dleppt = deppt/epsx
      dbgpt = (bgmx/2)*( ktx - dleppt )
      bdhvx = 2*al10*rtxmv*dbgpt

c     bdhv is BV in units of cm^3 kg^1/2 mol^-3/2 Angstrom^-1.

      bdhv = 1.0d-08*bdhvx

c     Get the Debye-Huckel BH constant, bdhhx.
c     Units are kJ kg^1/2 mol^-3/2 cm^-1,

      dleptp = (1.0d0/epsx)*deptp
      dbgtp = -(bgmx/2)*(dleptp + (1.0d0/tempk) + avx )
      bdhhx = 2*al10*rtxm*tempk*dbgtp

c     bdhh is BH in units of kJ kg^1/2 mol^-3/2 Angstrom^-1,

      bdhh = 1.0d-08*bdhhx

c     This parameter is often tabulated as BH/RT. That is
c     expressed below as bhrt. Units below are
c     kg^1/2 mol^-1/2 Angstrom^-1.

      bhrt = bdhh/rtxm

      END
