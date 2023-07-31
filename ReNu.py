# REYNOLDS CORRELATIONS AND NUSSELT VALUE

def calc_nusselt_laminar(Prh, Reh, d, L):
    NUlam = ((4.354**3) + (0.6**3) + ((1.953 * ((Prh*Reh*(d/L))**(1/3)) - 0.6)** 3) + 0.924 * ((Prh)**(1/3)) * ((Reh*(d/L))**(1/2)))**(1/3)
    return NUlam

def calc_nusselt_turbulent(fh, Reh, Prh):
    NUturb = ((fh / 8) * (Reh - 1000) * Prh) / (1 + 12.7 * ((fh / 8)** 0.5) * ((Prh) ** (2 / 3)) - 1)  # NUSSELT TURBULENTO
    return NUturb

def calc_nusselt_transition(calc):
    NUtran = (((calc_nusselt_laminar)**6) + ((calc_nusselt_turbulent)**6))**(1/6)