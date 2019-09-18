# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
#	Giulia Ballabio, gb258@leicester.ac.uk
#	Created on Jan, 2018
#
#	Reproduce the wind plot in Alexander et al. 2013.
#
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
from physics_constant import *

plt.style.use('default')
plt.rcParams['font.family'] = 'serif'
# plt.style.use('classic')
# plt.rcParams['font.family'] = 'sans'
# plt.rcParams['font.serif'] = 'Helvetica'

## ––––– create a polar grid ––––– ##
radius = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/grid_r.dat', 'r')]))
theta = np.arange(0., 1.3089, np.pi/600.) + (np.pi/12.)
## This theta is to plot data_hydro_midplane
# theta = np.arange(0., 1.5707, np.pi/600.)

dr = np.zeros(len(radius))
for i in range(len(radius)-1):
    dr[i] = radius[i+1] - radius[i]
dtheta = np.pi / len(theta)

## ––––– get the data from the files ––––– ##
# incl_deg = input("Insert the inclination angle used in the code (in degrees): ")
# b = input("Insert the value of b: ")
# r_in = input("Insert the inner radius: ")
# r_out = input("And the outer radius: ")
incl_deg = 90.0
b = 1.00
r_in = 0.1
r_out = 9.5
cs = 10.
str_cs = 10
mdot = 'nonorm'

## ––––– choose the species ––––– ##
species = 'NeII'
if (species == 'NeII'):
    m_atom = m_atom_ne
    Ab = Ab_ne
    A_ul = A_ul_ne
    lambda_ion = lambda_ne
    X_ion = X_ion_ne
    n_cr = n_cr_ne
    T_ul = T_ul_ne
if (species == 'SIIa'):
    m_atom = m_atom_sa
    Ab = Ab_sa
    A_ul = A_ul_sa
    lambda_ion = lambda_sa
    X_ion = X_ion_sa
    n_cr = n_cr_sa
    T_ul = T_ul_sa
if (species == 'SIIc'):
    m_atom = m_atom_sc
    Ab = Ab_sc
    A_ul = A_ul_sc
    lambda_ion = lambda_sc
    X_ion = X_ion_sc
    n_cr = n_cr_sc
    T_ul = T_ul_sc
else:
    m_atom = m_atom_o
    Ab = Ab_o
    A_ul = A_ul_o
    lambda_ion = lambda_o
    X_ion = X_ion_o
    n_cr = n_cr_o
    T_ul = T_ul_o

path_file = '../cs'+str(str_cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))

rho_mean = np.array(map(float, [lines.split()[0] for lines in open(str(path_file)+'/rho_grid.txt', 'r')]))
v_phi = np.array(map(float, [lines.split()[0] for lines in open(str(path_file)+'/v_grid.txt', 'r')]))

rho_hydro = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/rho_mean.dat', 'r')]))
v_r_hydro = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/v_r_mean.dat', 'r')]))
v_theta_hydro = -np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/v_th_mean.dat', 'r')]))
v_phi_hydro = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/v_phi_mean.dat', 'r')]))

## ––––– create a grid (r, theta) ––––– ##
grid_r, grid_theta = np.meshgrid(radius, theta, indexing='ij')
dr, dtheta = np.meshgrid(dr, dtheta, indexing='ij')

rho_2d = rho_mean.reshape(len(radius), len(theta))
v_phi_2d = v_phi.reshape(len(radius), len(theta))
rho_hydro_2d = rho_hydro.reshape(len(radius), len(theta))
v_phi_hydro_2d = v_phi_hydro.reshape(len(radius), len(theta))

rho_cr = n_cr * m_h * mu
rho_cr_2d = [[rho_cr for i in range(len(theta))] for j in range(len(radius))]

## ––––– plot the boundary condition ––––– ##
rho0 = []
rho0_hydro = []
for i in range(len(radius)):
    rho0.append(rho_2d[i][0]*rhog)
    rho0_hydro.append(rho_hydro_2d[i][0]*rhog)
plt.figure()
plt.loglog(radius*Rg/au, rho0, 'k')
plt.hlines(rho_cr, radius[0]*Rg/au, radius[len(radius)-1]*Rg/au, 'r', label=r'$\rho_{cr}$')
# plt.axis([1.e-2, 50., 1.e-3, 1.e3])
plt.title(r'Boundary condition', fontsize=15)
plt.xlabel(r'R / AU', fontsize=15)
plt.ylabel(r'$\rho_0$ / $\rho_{g}$', fontsize=15)
plt.legend(loc='best')
plt.savefig(str(path_file)+'/boundary_condition.png', format='png', dpi=300, bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/boundary_condition.png', format='png', dpi=300, bbox_inches='tight')
# plt.show()
plt.close()

## ––––– plot the density at a fixed theta as function of the radius ––––– ##
angle = 30.
angle_rad = angle * (np.pi/180.)
j = 1
while (theta[j] < angle_rad):
    j += 1
j_fixangle = j
rho_r = []
rho_hydro_r = []
for i in range(len(radius)):
    rho_r.append(rho_2d[i][j_fixangle])
    rho_hydro_r.append(rho_hydro_2d[i][j_fixangle])
plt.figure()
plt.loglog(radius, rho_r, 'r', label='model')
plt.loglog(radius, rho_hydro_r, 'b', label='hydro')
plt.axis([1.e-2, 50., 1.e-3, 1.e3])
plt.xlabel(r'R / R$_{g}$',fontsize=15)
plt.ylabel(r'$\rho$ / $\rho_{g}$', fontsize=15)
plt.legend(loc='best')
plt.savefig(str(path_file)+'/density_fixedtheta_'+str(angle)+'deg.png', format='png', bbox_inches='tight')
# plt.show()
plt.close()

## ––––– convert to a cartesian coordinate system ––––– ##
r = grid_r*np.cos(grid_theta)
z = grid_r*np.sin(grid_theta)

## ––––– formatting the colorbar ticks ––––– ##
def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${}$'.format(b)

plt.figure()
CS = plt.pcolormesh(r, z, rho_hydro_2d, cmap='hot', norm=LogNorm(), vmin=0.05, vmax=50.)
# plt.contour(r, z, cs, colors='cyan')
# plt.pcolormesh(r, -z, rho2d, cmap='viridis', norm=LogNorm(), vmin=0.05, vmax=20.)
# cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
cbar = plt.colorbar(CS)
plt.axis([0., 1.25, 0., 0.9])
# plt.axis([0.02, 0.04, 0, 0.02])
# plt.axis([0.02, 10., 0., 10.])
plt.xlabel(r'R / R$_{g}$',fontsize=15)
plt.ylabel(r'z / R$_{g}$',fontsize=15)
cbar.set_label(r'Log($\rho$ / $\rho_{g}$)')
# plt.savefig(str(path_file)+'/density_polar.png', format='png', bbox_inches='tight', dpi=300)
plt.savefig('../data_hydro/'+str(species)+'/density_polar_zoom.png', format='png', bbox_inches='tight', dpi=300)
plt.show()
plt.close()

plt.figure()
CS = plt.pcolormesh(r, z, v_phi_2d, cmap='winter', norm=LogNorm())#, vmin=0.25, vmax=10.)
cbar = plt.colorbar(CS)
# plt.axis([0., 1.25, 0, 0.9])
plt.axis([0.02, 10., 0., 10.])
plt.xlabel(r'R / R$_{g}$',fontsize=15)
plt.ylabel(r'z / R$_{g}$',fontsize=15)
cbar.set_label(r'Log($v_{K}$)')
plt.savefig(str(path_file)+'/keplerian_velocity.png', format='png', bbox_inches='tight')#, dpi=1000)
# plt.savefig('../data_hydro/'+str(species)+'/keplerian_velocity.png', format='png', bbox_inches='tight', dpi=300)
# plt.show()
plt.close()

## ––––– Creating a 2D map of the flux ––––– ##
v_th = cs / ((m_atom)**0.5)
Temp = T_gas * (cs/10.)**2.
nu = speed_light/lambda_ion
A_hnu = A_ul*h_planck*nu
constants = Ab*A_hnu*X_ion

dV = np.zeros((len(radius), len(theta)))
r_cr = np.zeros((len(theta)))
z_cr = np.zeros((len(theta)))
n_e = np.zeros((len(radius), len(theta)))
C_ul = np.zeros((len(radius), len(theta)))
P_u = np.zeros((len(radius), len(theta)))
cell_flux = np.zeros((len(radius), len(theta)))
for i in range(len(radius)):
    for j in range(len(theta)):
        ## N.B. n_cr is in cm^-3 so must be also n_e and dV
        dV[i][j] = (2.*np.pi) * radius[i]*Rg * radius[i]*Rg * np.sin(theta[j]) * dr[i]*Rg * dtheta[j]
        n_e[i][j] = rho_hydro_2d[i][j]*rhog * (ng/rhog)
        if (n_e[i][j] > 0.0):
            C_ul[i][j] = 1. + (n_cr/n_e[i][j])
            P_u[i][j] = 1. / ((2.*C_ul[i][j]*np.exp(T_ul/Temp)) + 1.)
            cell_flux[i][j] = constants * P_u[i][j] * n_e[i][j] * dV[i][j]
        else:
            cell_flux[j][j] = 0.0
# if (species == 'NeII'):
#     max_flux = np.amax(cell_flux)

plt.figure()
CS = plt.pcolormesh(r*Rg/au, z*Rg/au, cell_flux/np.amax(cell_flux), cmap='inferno', norm=LogNorm() , vmin=1.e-1, vmax=1.e0)
# plt.contour(r_cr*Rg/au, z_cr*Rg/au, n_cr_2d, color='k')
cbar = plt.colorbar(CS)
# plt.axis([0., 20., 0., 15.])
plt.xlabel(r'R / AU',fontsize=15)
plt.ylabel(r'z / AU',fontsize=15)
cbar.set_label(r'Log(L)') # / L_{\odot}
# plt.savefig(str(path_file)+'/line_flux.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('../data_hydro/'+str(species)+'/line_flux.png', format='png', bbox_inches='tight', dpi=300)
plt.savefig('../data_hydro/'+str(species)+'/line_flux.pdf', format='pdf', bbox_inches='tight', dpi=300)
plt.show()
plt.close()

ne_r = []
for i in range(len(radius)):
    ne_r.append(n_e[i][0])
# plt.figure()
plt.loglog(radius*Rg/au, ne_r, 'k')
plt.hlines(n_cr, radius[0]*Rg/au, radius[len(radius)-1]*Rg/au, 'r', label=r'$n_{cr}$')
# plt.axis([1.e-2, 50., 1.e-3, 1.e3])
plt.xlabel(r'R / AU',fontsize=15)
plt.ylabel(r'$n_e$', fontsize=15)
# plt.legend(loc='best')
plt.savefig(str(path_file)+'/numerical_density_midplane.png', format='png', dpi=300, bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/numerical_density_midplane.png', format='png', bbox_inches='tight')
plt.show()
plt.close()

flux_r = []
for i in range(len(radius)):
    flux_r.append(cell_flux[i][1]/np.amax(cell_flux))
plt.figure()
plt.loglog(radius*Rg/au, flux_r, 'k')
plt.plot(radius*Rg/au, (radius*Rg/au)**(-b), 'r', label='$R^{-b}$')
# plt.yscale('log')
# plt.axis([1.e-2, 50., 1.e-3, 1.e3])
plt.xlabel(r'R / AU',fontsize=15)
plt.ylabel(r'Flux(R)', fontsize=15)
plt.legend(loc='best')
plt.savefig(str(path_file)+'/flux_midplane.png', format='png', dpi=300, bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/flux_midplane.png', format='png', bbox_inches='tight')
plt.show()
plt.close()
