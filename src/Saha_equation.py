import numpy as np
from physics_constant import *

## Saha equation to determine the ionization fraction of [SII]
## It MUST be divided by the electron density

print 'Computing the Saha equation for the ionization fraction...'

cs = 10.
T_gas = 1.e4*(cs/10.)**2.

species = 'NeII'
if (species == 'NeII'):
    g1 = 4.
    g2 = 2.
    Ipot = 21.56*eV
elif (species == 'SIIa'):
    g1 = 4.
    g2 = 2.
    Ipot = 10.36*eV
elif (species == 'SIIb'):
    g1 = 4.
    g2 = 2.
    Ipot = 10.36*eV
elif (species == 'SIIc'):
    g1 = 4.
    g2 = 6.
    Ipot = 10.36*eV

# Saha equation, but the number density of free electron is not included
N2overN1 = 2. *(2.*np.pi*m_e*k_b*T_gas)**1.5 / (h_planck)**3. * (g2/g1) * np.exp(-(Ipot / (k_b*T_gas)))
# X_II = N2overN1 / (1. + N2overN1)

# print X_II

# Let's use the ionized hydrogen fraction to derive the ionized oxygen fraction
# Since the ionisation potential is the same we can assume that the number density of free electron is also the same
Ipot_H = 13.6*eV
# This is the fraction of ionized hydrogen at T_gas=10^4 K
X_II_H = 0.7
n2n1_H = X_II_H / (1.-X_II_H)

n2n1 = n2n1_H * np.exp((Ipot_H-Ipot)/(k_b*T_gas))
X_II = n2n1 / (1. + n2n1)
X_I_O = 1. - X_II

print X_II
