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