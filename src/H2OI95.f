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