import numpy as np

## ----- ASTRONOMICAL CONSTANTS ----- ##
au = 1.496e13  			                                                        #[cm]
G = 6.672e-8 			                                                        #[cm3 g-1 s-2]
Msun = 2.e33 		                                                            #[g]
Lsun = 3.826e33                                                                 #[erg/s]
year = 31536000.0		                                                        #[s]
Mstar = 1.*Msun 		                                                        #[g]
MJ = 1.898e30  		                                                            #[g]
eV=1.60218e-12  		                                                        #[erg]
h_planck=6.6261e-27

N_part = 1.191e57
cs = 1.0e6 			                                                            #[cm/s]
m_h = 1.6726e-24                                                                #[g]
m_e=9.1094e-28                                                                  #[g]
mu = 1.
CC = 0.14
Phi_star = 1.e41                                                                #[photons/s]
alphab = 2.60e-13                                                               #[cm3 s-1]
T = 1.e4                                                                        #[K]
k_b=1.38e-16                                                                    #[erg K-1]

## ----- NeII CONSTANTS -----##
m_ne = 20.*m_h                                                                  #[g]
Ab_ne = 1.e-4
A_ul = 8.39e-3                                                                  #[s-1]
A_hnu = 1.332e-15
X_II = 0.75
T_ul = 1122.8                                                                   #[K]
n_cr = 5.0e5                                                                    #[cm-3]

## ----- PHYSICS SCALING FACTORS ----- ##
Rg = G*Mstar / (cs**2.)	                                                        #[cm]
ng = CC * (3.*Phi_star / (4.*np.pi*alphab*((Rg)**3)))**0.5                      #[cm-3]
rhog = ng * m_h * mu 	                                                        #[g cm-3]
v_th = cs * ((m_h/m_ne)**0.5)                                                   #[cm/s]
