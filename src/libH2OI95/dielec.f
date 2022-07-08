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
c     and Debye-H�ckel Coefficients. J. Phys. Chem. Ref. Data 26,
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