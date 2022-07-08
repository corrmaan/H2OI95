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