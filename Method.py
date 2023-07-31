# LMTD / e-NUT METHOD FOR HX PROJECT
# # Carlos Eduardo Bibow Corrêa

import math
import CoolProp.CoolProp as CP
import ReNu
import HXDim
import Input
from Input import mhtotal, mctotal, Tcin, Thin, pcin, phin
from HXDim import Atot_c, Atot_h, d_c, d_h, L_c, L_h, t_c, t_h, k, n_c, w, Afree_c, Afree_h



fluid = 'Water'
air = 'Air'

Thout = Thin  # ASSUMING INLET AND OUTLET TEMPERATURES ARE THE SAME, AS ESTIMATED VALUE
Tcout = Tcin  # ASSUMING INLET AND OUTLET TEMPERATURES ARE THE SAME, AS ESTIMATED VALUE
error = 0.5  # INITIAL VALUE
est = error
tolerance = 1E-5  # ERROR TOLERANCE

while error > tolerance:

    mhin = Input.mhtotal / HXDim.n_h
    mcin = Input.mctotal / HXDim.n_c

    Thmed = (Thout + Thin) / 2
    Tcmed = (Tcout + Tcin) / 2

    # PROPRIEDADES
    Hot_fluid = 'Water'
    Cold_fluid = 'Air'

    cp1h = CP.PropsSI('C', 'T', Thmed + 273, 'P', phin * 1e5, Hot_fluid)
    mu1h = CP.PropsSI('V', 'T', Thmed + 273, 'P', phin * 1e5, Hot_fluid)
    rho1h = CP.PropsSI('D', 'T', Thmed + 273, 'P', phin * 1e5, Hot_fluid)
    k1h = CP.PropsSI('L', 'T', Thmed + 273, 'P', phin * 1e5, Hot_fluid)

    cp1c = CP.PropsSI('C', 'T', Tcmed + 273, 'P', pcin * 1e5, Cold_fluid)
    mu1c = CP.PropsSI('V', 'T', Tcmed + 273, 'P', pcin * 1e5, Cold_fluid)
    rho1c = CP.PropsSI('D', 'T', Tcmed + 273, 'P', pcin * 1e5, Cold_fluid)
    k1c = CP.PropsSI('L', 'T', Tcmed + 273, 'P', pcin * 1e5, Cold_fluid)

    Prh = (mu1h * cp1h) / k1h
    Prc = (mu1c * cp1c) / k1c

    Cc = cp1c * mctotal
    Ch = cp1h * mhtotal

    if Cc < Ch:
        Cmin = Cc
        Cmax = Ch
    else:
        Cmin = Ch
        Cmax = Cc

    Cr = Cmin / Cmax

    Thout = Thin - est * (Cmin * (Thin - Tcin)) / Ch
    Tcout = Tcin + est * (Cmin * (Thin - Tcin)) / Cc

    REc = (d_c * mctotal) / (mu1c * Afree_c)  # REYNOLDS FOR COLD CHANNEL
    REh = (d_h * mhtotal) / (mu1h * Afree_h)  # REYNOLDS FOR HOT CHANNEL

    REL = 2300  # MAXIMUM REYNOLDS VALUE FOR LAMINAR FLOW
    RET = 10000  # MINIMUM REYNOLDS VALUE FOR TURBULENT FLOW

    lthh = (L_h/d_h)/(REc*Prc)
    lthc = (L_c/d_c)/(REh*Prh)
    fh = (1 / 4) * (1.8 * math.log(REh) - 1.5) ** (-2)
    fc = (1 / 4) * (1.8 * math.log(REc) - 1.5) ** (-2)

    # Para o ramal quente
    if REh <= REL:
        NUh = ReNu.calc_nusselt_laminar(Prh, REh, d_h, L_h)

    elif (REh > REL) and (REh < RET):
        NUh = (((ReNu.calc_nusselt_laminar(Prh, REh, d_h, L_h))**6) + ((ReNu.calc_nusselt_turbulent(fh, REh, Prh))**6))**(1/6)

    else:
        NUh = ReNu.calc_nusselt_turbulent(fh, REh, Prh)

    # Para o ramal frio
    if REh <= REL:
        NUc = ReNu.calc_nusselt_laminar(Prc, REc, d_c, L_c)

    elif (REh > REL) and (REh < RET):
        NUc = (((ReNu.calc_nusselt_laminar(Prc, REc, d_c, L_c))**6) + ((ReNu.calc_nusselt_turbulent(fc, REc, Prc))**6))**(1/6)

    else:
        NUc = ReNu.calc_nusselt_turbulent(fc, REc, Prc)

    hh = (k1h * NUh)/d_h
    hc = (k1c * NUc)/d_c
    # ---------------------------- APROXIMAÇÃO DE ALETAS
    Ah = d_h
    Ac = d_c

    Palh = 2 * (L_h + ((t_h - Ah)/2))
    Atralh = L_h * ((t_h - Ah)/2)
    Palc = 2 * (L_c + ((t_c - Ac)/2))
    Atralc = L_c * ((t_c - Ac)/2)

    Mhot = (Palh * hh)/(k * Atralh)
    Mcold = (Palc * hc)/(k * Atralc)

    Etah = math.tanh((Mhot * Ah/2)/(Mhot * Ah/2))
    Etac = math.tanh((Mcold * Ac/2)/(Mcold * Ac/2))

    Aletahvc = 0.5 * Ah * L_h
    Aletacvc = 0.5 * Ac * L_c
    Awethvc = d_h * math.pi * L_h/2
    Awetcvc = d_h * math.pi * L_c/2
    Awall = L_c * L_h * ((2*n_c) + 1)

    etah = 1 - ((2 * Aletahvc/Awethvc) * (1 - Etah))
    etac = 1 - ((2 * Aletacvc/Awetcvc) * (1 - Etac))

    Rhconv = 1 / (Atot_h * hh * etah)
    Rcconv = 1 / (Atot_c * hc * etac)
    Rcond = w/(k * Awall)

    Rtot = Rhconv + Rcconv + Rcond
    UA = 1/Rtot

    NUT = UA/Cmin
    Cr = Cmin/Cmax

    epslon = 1 - math.exp((1/Cr)*(NUT**0.22)*(math.exp(-Cr*(NUT**0.78))-1))

    error = abs(epslon - est)

    est = epslon

q = epslon * Cmin * (Thin - Tcin)
qc = mctotal * cp1c * (Tcin - Tcout)
qh = mctotal * cp1h * (Thin - Thout)
    
print("est, cp1h, q, qc, qh", est, cp1h, q, qc, qh)
