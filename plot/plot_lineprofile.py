import numpy as np
import matplotlib.pyplot as plt

#plt.style.use('seaborn-darkgrid')
plt.rcParams["font.family"] = 'Verdana'

## ----- GET THE DATA FROM THE OUTPUT FILE FROM FORTRAN ----- ##
incl_deg = input("Insert the inclination angle used in the code (in degrees): ")
b = input("Insert the value of b: ")
v = np.array(map(float, [lines.split()[0] for lines in open('./b'+str(b)+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
line_flux = np.array(map(float, [lines.split()[1] for lines in open('./b'+str(b)+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))

#print np.amax(line_flux)
peak_edgeon = np.amax(line_flux)

#blueshift according to convention
if incl_deg == 0.0:
    v = -1. * v
    
    plt.figure()
    plt.plot(v, line_flux / np.amax(line_flux), 'r')
    plt.title(r'Line profile: i='+str(round(incl_deg, 2))+'$^\circ$', fontsize=15)
    plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
    plt.ylabel(r'Normalized L(v)', fontsize=15)
    plt.axis([-40., 40., 0., 1.2])
    plt.savefig('./b'+str(b)+'/line_profile_i'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
    plt.show()
else:
    plt.figure()
    plt.plot(v, line_flux / np.amax(line_flux), 'r')
    plt.title(r'Line profile: i='+str(round(incl_deg, 2))+'$^\circ$', fontsize=15)
    plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
    plt.ylabel(r'Normalized L(v)', fontsize=15)
    plt.axis([-40., 40., 0., 1.2])
    plt.savefig('./b'+str(b)+'/line_profile_i'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
    plt.show()

#r = np.array(map(float, [lines.split()[0] for lines in open('./b'+str(b)+'/bound_cond.txt', 'r')]))
#rho = np.array(map(float, [lines.split()[1] for lines in open('./b'+str(b)+'/bound_cond.txt', 'r')]))
#v_theta = np.array(map(float, [lines.split()[2] for lines in open('./b'+str(b)+'/bound_cond.txt', 'r')]))
#
#plt.figure()
#plt.loglog(r, rho, 'r')
#plt.title(r'Boundary condition', fontsize=15)
#plt.xlabel(r'r [au]', fontsize=15)
#plt.ylabel(r'$rho_{0}$', fontsize=15)
#plt.savefig('./b'+str(b)+'/boundary_condition.png', format='png', bbox_inches='tight')
#plt.show()
#
#plt.figure()
#plt.plot(r, -1.*v_theta, 'r')
#plt.axhline(y=3., color='k', linestyle='--')
#plt.xlabel(r'r [au]', fontsize=15)
#plt.ylabel(r'-$v_{theta}$ [km/s]', fontsize=15)
#plt.axis([0., 30., 0., 6.])
#plt.savefig('./b'+str(b)+'/vth_atbound_'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#plt.show()
#
#plt.figure()
#plt.plot(r, -1.*rho*v_theta, 'r')
#plt.axvline(x=1.3305, color='k', linestyle='--')
#plt.xlabel(r'r [au]', fontsize=15)
#plt.ylabel(r'$rho * v_{theta}$', fontsize=15)
#plt.savefig('./b'+str(b)+'/rho_vth_'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#plt.show()

#plt.figure()
#plt.plot(r, -1.*rho*v_theta, 'r')
#plt.axvline(x=1.3305, color='k', linestyle='--')
#plt.xlabel(r'r [au]', fontsize=15)
#plt.ylabel(r'$rho * v_{theta}$', fontsize=15)
#plt.axis([0., 7., -0.5e-11, 1.0e-11])
#plt.savefig('./b'+str(b)+'/rho_vth_zoom_'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#plt.show()
