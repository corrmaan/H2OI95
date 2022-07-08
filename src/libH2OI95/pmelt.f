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