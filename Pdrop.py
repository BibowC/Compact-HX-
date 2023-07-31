 # MODELO DE QUEDA DE PRESSÃO

Gh = mhin / Atot/2
Gc = mctotal / Atot/2
rhomh = ((1/rho1h)+(1/rho2h))/2
rhomc = ((1/rho1c)+(1/rho2c))/2
if (REh > 2300):
 Keh = 0.75
 Kch = 0.47
elif (REh<2300):
 Keh = 0.69
 Kch = 0.99

if (REc > 2300):
 Kec = 0.75
 Kcc = 0.47
elif (REc<2300):
 Kec = 0.69
 Kcc = 0.99

#Ramal quente
dPCoreh = (Gh**2 / 2*rho1h) * ((4*fh * L/d * rhomh + ((rho1h/rho2h) - 1) +
(1 - Poro**2) - (rho1h/rho2h)*(1 - Poro**2) + Kch + Keh*(rho1h/rho2h)))
dPCorec = (Gc**2 / 2*rho1c) * ((4*fc * L/d * rhomc + ((rho1c/rho2c) - 1) +
(1 - Poro**2) - (rho1c/rho2c)*(1 - Poro**2) + Kcc + Kec*(rho1c/rho2c)))