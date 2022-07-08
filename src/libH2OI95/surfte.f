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