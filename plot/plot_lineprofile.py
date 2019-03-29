import numpy as np
import matplotlib.pyplot as plt

## ----- GET THE DATA FROM THE OUTPUT FILE FROM FORTRAN ----- ##
# incl_deg = input("Insert the inclination angle used in the code (in degrees): ")
# b = input("Insert the value of b: ")
# r_in = input("Insert the inner radius: ")
# r_out = input("The outer radius: ")
# cs = input("And the sound speed: ")
b = 0.75 #[0.75, 1.00, 1.50, 2.00]
incl_deg = [0.0, 20.0, 45.0, 60.0, 90.0] #[0.0, 20.0, 45.0, 60.0, 90.0]
r_in = 0.1 #[0.03, 0.1, 1.0]
r_out = 9.5 #[5.0, 9.5]
cs = 10 #3 5
species = 'SII'

quantity = incl_deg
v = []
line_flux = []
plt.figure()
for i in range(len(quantity)):
    v = -np.array(map(float, [lines.split()[0] for lines in open('../cs'+str(cs)+'kms/'+str(species)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux = np.array(map(float, [lines.split()[1] for lines in open('../cs'+str(cs)+'kms/'+str(species)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))

    c = [float(i)/float(len(quantity)), 0.0, float(len(quantity)-i)/float(len(quantity))]
    # if incl_deg[i] == 90.0:
    #     v = -1. * v
    # else:
    #     v = v
    plt.plot(v, line_flux / np.amax(line_flux), color=c, label='i='+str(quantity[i]))
    # plt.vlines(0., 0., 1.3, 'k', linestyle='--')
    plt.title(r'Line profile', fontsize=15)
    plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
    plt.ylabel(r'Normalized L(v)', fontsize=15)
    # plt.axis([-40., 40., 0., 1.2])
    plt.legend(loc='best')
    # plt.savefig('./lineprofiles/'+str(species)+'/line_profile_i'+str(round(incl_deg, 2))+'_r'+str(r_in)+'_r'+str(r_out)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
    plt.savefig('./lineprofiles/'+str(species)+'/line_profile_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
    # plt.savefig('./lineprofiles/'+str(species)+'/line_profile_b'+str('{:.2f}'.format(round(b, 2)))+'_i'+str(round(incl_deg, 2))+'_r'+str(r_out)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
    # plt.savefig('./lineprofiles/'+str(species)+'/line_profile_b'+str('{:.2f}'.format(round(b, 2)))+'_i'+str(round(incl_deg, 2))+'_r'+str(r_in)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
plt.show()

# v = np.array(map(float, [lines.split()[0] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
# line_flux = np.array(map(float, [lines.split()[1] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))

#print np.amax(line_flux)
# peak_edgeon = np.amax(line_flux)

#blueshift according to convention
# if incl_deg == 0.0:
#     v = -1. * v
#
#     plt.figure()
#     plt.plot(v, line_flux / np.amax(line_flux), 'r')
#     plt.vlines(0., 0., 1.3, 'k', linestyle='--')
#     plt.title(r'Line profile: i='+str(round(incl_deg, 2))+'$^\circ$', fontsize=15)
#     plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
#     plt.ylabel(r'Normalized L(v)', fontsize=15)
#     plt.axis([-40., 40., 0., 1.2])
#     plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#     plt.show()
# else:
#     plt.figure()
#     plt.plot(v, line_flux / np.amax(line_flux), 'r')
#     plt.vlines(0., 0., 1.3, 'k', linestyle='--')
#     plt.title(r'Line profile: i='+str(round(incl_deg, 2))+'$^\circ$', fontsize=15)
#     plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
#     plt.ylabel(r'Normalized L(v)', fontsize=15)
#     plt.axis([-40., 40., 0., 1.2])
#     plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#     plt.show()
#
# r = np.array(map(float, [lines.split()[0] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/bound_cond.txt', 'r')]))
# rho = np.array(map(float, [lines.split()[1] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/bound_cond.txt', 'r')]))
# v_theta = np.array(map(float, [lines.split()[2] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/bound_cond.txt', 'r')]))
#
# # plt.figure()
# plt.loglog(r, rho, 'r')
# plt.title(r'Boundary condition', fontsize=15)
# plt.xlabel(r'r [au]', fontsize=15)
# plt.ylabel(r'$rho_{0}$', fontsize=15)
# plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/boundary_condition.png', format='png', bbox_inches='tight')
# plt.show()

#plt.figure()
#plt.plot(r, -1.*v_theta, 'r')
#plt.axhline(y=3., color='k', linestyle='--')
#plt.xlabel(r'r [au]', fontsize=15)
#plt.ylabel(r'-$v_{theta}$ [km/s]', fontsize=15)
#plt.axis([0., 30., 0., 6.])
#plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/vth_atbound_'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#plt.show()
#
#plt.figure()
#plt.plot(r, -1.*rho*v_theta, 'r')
#plt.axvline(x=1.3305, color='k', linestyle='--')
#plt.xlabel(r'r [au]', fontsize=15)
#plt.ylabel(r'$rho * v_{theta}$', fontsize=15)
#plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/rho_vth_'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#plt.show()

#plt.figure()
#plt.plot(r, -1.*rho*v_theta, 'r')
#plt.axvline(x=1.3305, color='k', linestyle='--')
#plt.xlabel(r'r [au]', fontsize=15)
#plt.ylabel(r'$rho * v_{theta}$', fontsize=15)
#plt.axis([0., 7., -0.5e-11, 1.0e-11])
#plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'/incl_'+str(round(incl_deg, 2))+'/rho_vth_zoom_'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight')
#plt.show()
