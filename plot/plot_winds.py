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
from matplotlib.colors import LogNorm
from physics_constant import *

#plt.style.use('seaborn-darkgrid')
#plt.rcParams["font.family"] = 'Montreal SF'
#plt.rcParams["font.family"] = 'DejaVu Sans'

## ––––– create a polar grid ––––– ##
radius = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/grid_r.dat', 'r')]))
theta = np.arange(0., 1.3089, np.pi/600.)+(np.pi/12.)

## ––––– get the data from the files ––––– ##
incl_deg = input("Insert the inclination angle used in the code (in degrees): ")
b = input("Insert the value of b: ")
r_in = input("Insert the inner radius: ")
rho_mean = np.array(map(float, [lines.split()[0] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/rho_grid.txt', 'r')]))
v_phi = np.array(map(float, [lines.split()[0] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/v_phi_grid.txt', 'r')]))

# rho_mean = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/rho_mean.dat', 'r')]))
# v_r = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/v_r_mean.dat', 'r')]))
# v_theta = -np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/v_th_mean.dat', 'r')]))
# v_phi = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/v_phi_mean.dat', 'r')]))

## ––––– create a grid (r, theta) ––––– ##
grid_r, grid_theta = np.meshgrid(radius, theta, indexing='ij')
rho_2d = rho_mean.reshape(len(radius), len(theta))

v_phi_2d = v_phi.reshape(len(radius), len(theta))

rho_cr = n_cr * m_h * mu / rhog
rho_cr_2d = [[rho_cr for i in range(len(theta))] for j in range(len(radius))]

## ––––– plot the boundary condition ––––– ##
rho0 = []
for i in range(len(radius)):
    rho0.append(rho_2d[i][0])
plt.figure()
plt.loglog(radius, rho0, 'r')
plt.axis([1.e-2, 50., 1.e-3, 1.e3])
#plt.title(r'Boundary condition', fontsize=15)
plt.xlabel(r'R / R$_{g}$',fontsize=15)
plt.ylabel(r'$\rho$ / $\rho_{g}$', fontsize=15)
plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/boundary_condition.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro/boundary_condition.png', format='png', bbox_inches='tight')
plt.show()
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
CS = plt.pcolormesh(r, z, rho_2d, cmap='hot', norm=LogNorm(), vmin=0.05, vmax=20.)
#plt.contour(rho_cr_2d, colors='r')
#plt.pcolormesh(r, -z, rho2d, cmap='viridis', norm=LogNorm(), vmin=0.05, vmax=20.)
#cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
cbar = plt.colorbar(CS)
#plt.axis([0., 1.25, 0, 0.9])
#plt.axis([0.02, 0.04, 0, 0.02])
plt.xlabel(r'R / R$_{g}$',fontsize=15)
plt.ylabel(r'z / R$_{g}$',fontsize=15)
cbar.set_label(r'Log($\rho$ / $\rho_{g}$)')
plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/density_polar.png', format='png', bbox_inches='tight')#, dpi=1000)
# plt.savefig('../data_hydro/density_polar.png', format='png', bbox_inches='tight')#, dpi=1000)
plt.show()
plt.close()

plt.figure()
CS = plt.pcolormesh(r, z, v_phi_2d, cmap='winter', norm=LogNorm(), vmin=0.05, vmax=20.)
cbar = plt.colorbar(CS)
#plt.axis([0., 1.25, 0, 0.9])
plt.xlabel(r'R / R$_{g}$',fontsize=15)
plt.ylabel(r'z / R$_{g}$',fontsize=15)
cbar.set_label(r'Log($v_{K}$)')
plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/keplerian_velocity.png', format='png', bbox_inches='tight')#, dpi=1000)
# plt.savefig('../data_hydro/keplerian_velocity.png', format='png', bbox_inches='tight')#, dpi=1000)
plt.show()
plt.close()

plt.figure()
plt.plot(z,v_phi_2d[0.5,:])
plt.show()
plt.close()

# x = np.array(map(float, [lines.split()[0] for lines in open('../../input_from_model/b1.0/streamline.txt', 'r')]))
# y = np.array(map(float, [lines.split()[1] for lines in open('../../input_from_model/b1.0/streamline.txt', 'r')]))
# #x = x / Rg
# #print Rg/au
# #y = y / (Rg / au)
# delta = 0.05
# x = x - 0.9
# plt.figure()
# CS = plt.pcolormesh(r, z, rho_2d, cmap='hot', norm=LogNorm(), vmin=0.05, vmax=20.)
# plt.contour(rho_cr_2d, colors='r')
# for i in range(0, int(5/delta)):
#     x_shift = x + delta*i
#     plt.plot(x_shift, y, color='white', linewidth=1)
# cbar = plt.colorbar(CS)
# plt.axis([0., 1.25, 0, 0.9])
# plt.xlabel(r'R / R$_{g}$',fontsize=15)
# plt.ylabel(r'z / R$_{g}$',fontsize=15)
# cbar.set_label(r'Log($\rho$ / $\rho_{g}$)')
# plt.savefig('./density_scalefree_b1.0.png', format='png', bbox_inches='tight')#, dpi=1000)
# #plt.show()
# plt.close()

## ––––– convert to physical units in cgs ––––– ##
#plt.figure()
#plt.pcolormesh(r*Rg, z*Rg, (rho_2d)*ng, cmap='viridis', norm=LogNorm())#, vmin=0.1, vmax=10.)
##cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#cbar = plt.colorbar()
##plt.axis([0., 1.25, 0, 0.9])
#plt.xlabel(r'R',fontsize=15)
#plt.ylabel(r'z',fontsize=15)
#cbar.set_label(r'Log($\rho$)')
#plt.savefig('./density_phys.png', format='png', bbox_inches='tight')#, dpi=500)
#plt.show()
#plt.close()
#
### ––––– make a 3D grid from 2D grid ––––– ##
### ––––– extend theta from 0 up to pi ––––– ##
#rev_theta = -theta[::-1]
#grid_r, grid_rev_theta = np.meshgrid(radius, rev_theta, indexing='ij')
#rev_z = grid_r*np.sin(grid_rev_theta)
#rev_rho2d = np.fliplr(rho_2d)
#rev_r = np.fliplr(r)
#
#rho = np.concatenate((rev_rho2d, rho_2d), axis=1)
#R = np.concatenate((rev_r, r), axis=1)
#Z = np.concatenate((rev_z, z), axis=1)
#
#plt.figure()
#plt.pcolormesh(R, Z, rho, cmap='magma', norm=LogNorm(), vmin=0.05, vmax=20.)
##cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
#cbar = plt.colorbar()
#plt.axis([0., 1.25, -0.9, 0.9])
#plt.xlabel(r'R / R$_{g}$',fontsize=15)
#plt.ylabel(r'z / R$_{g}$',fontsize=15)
#cbar.set_label(r'log$_{10}$($\rho$ / $\rho_{g}$)')
#plt.savefig('./density(0,pi).png', format='png', bbox_inches='tight')#, pdi=1000)
#plt.show()
#plt.close()

## ––––– rotate all the structure along the z-axis to make a 3d structure ––––– ##
## ––––– define the phi angle ––––– ##
# phi = np.arange(0., 2.*np.pi, np.pi/600.)
# grid_r_3d, grid_theta_3d, grid_phi_3d = np.meshgrid(radius, theta, phi, indexing='ij')

## ––––– velocity map ––––– ##
# plt.figure()
# plt.pcolormesh(r, z, vr_2d, cmap='hot', norm=LogNorm(), vmin=0.05, vmax=20.)
# #cbar = plt.colorbar(format=ticker.FuncFormatter(fmt))
# cbar = plt.colorbar()
# #plt.quiver(r, z, v_r, v_theta, units='width')
# plt.axis([0., 1.25, 0, 0.9])
# plt.xlabel(r'R / R$_{g}$',fontsize=15)
# plt.ylabel(r'z / R$_{g}$',fontsize=15)
# cbar.set_label(r'log$_{10}$($v_r$ / $c_s$)')
# plt.savefig('./velocity_polar.png', format='png', bbox_inches='tight')#, dpi=1000)
# plt.show()
# plt.close()
