import numpy as np

## ----- ASTRONOMICAL CONSTANTS ----- ##
au = 1.496e13  			                                                        #[cm]
G = 6.672e-8 			                                                        #[cm3 g-1 s-2]
Msun = 2.e33 		                                                            #[g]
Lsun = 3.826e33                                                                 #[erg/s]
year = 31536000.0		                                                        #[s]
Mstar = 1.*Msun 		                                                        #[g]
MJ = 1.898e30  		                                                            #[g]
eV = 1.60218e-12  		                                                        #[erg]
h_planck = 6.6261e-27
speed_light = 2.9979e10                                                         #[cm/s]

N_part = 1.191e57
cs = 1.0e6 			                                                            #[cm/s]
m_h = 1.6726e-24                                                                #[g]
m_e = 9.1094e-28                                                                #[g]
mu = 1.
CC = 0.14
Phi_star = 1.e41                                                                #[photons/s]
alphab = 2.60e-13                                                               #[cm3 s-1]
T_gas = 1.e4                                                                    #[K]
k_b=1.38e-16                                                                    #[erg K-1]

## ----- NeII CONSTANTS -----##
m_atom_ne = 20.                                                                 #[g]
Ab_ne = 1.e-4
A_ul_ne = 8.39e-3                                                               #[s-1]
lambda_ne = 12.81e-4
X_ion_ne = 0.75
T_ul_ne = 1122.8                                                                #[K]
n_cr_ne = 5.0e5                                                                 #[cm-3]

## ----- SIIa CONSTANTS -----##
m_atom_sa = 32.                                                                  #[g]
Ab_sa = 1.45e-5
A_ul_sa = 1.9e-1                                                                 #[s-1]
lambda_sa = 406.98e-7
X_ion_sa = 1.0
T_ul_sa = 35354.                                                                 #[K]
n_cr_sa = 2.6e6                                                                  #[cm-3]

## ----- SIIc CONSTANTS -----##
m_atom_sc = 32.                                                                  #[g]
Ab_sc = 1.45e-5
A_ul_sc = 2.0e-4                                                                 #[s-1]
lambda_sc = 671.83e-7
X_ion_sc = 1.0
T_ul_sc = 21416.                                                                 #[K]
n_cr_sc = 1.7e3                                                                  #[cm-3]

## ----- OI CONSTANTS -----##
m_atom_o = 16.                                                                  #[g]
Ab_o = 5.37e-4
A_ul_o = 5.6e-3                                                                 #[s-1]
lambda_o = 630.0e-7
X_ion_o = 1.0
T_ul_o = 22830.                                                                 #[K]
n_cr_o = 1.8e6                                                                  #[cm-3]

## ----- PHYSICS SCALING FACTORS ----- ##
Rg = G*Mstar / (cs**2.)	                                                        #[cm]
ng = CC * (3.*Phi_star / (4.*np.pi*alphab*((Rg)**3)))**0.5                      #[cm-3]
rhog = ng * m_h * mu 	                                                        #[g cm-3]
