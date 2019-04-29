import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.style.use('classic')
plt.rcParams['font.family'] = 'sans' #Courier New'
plt.rcParams['font.serif'] = 'Helvetica'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12

## ----- GET THE DATA FROM THE OUTPUT FILE FROM FORTRAN ----- ##
# incl_deg = input("Insert the inclination angle used in the code (in degrees): ")
# b = input("Insert the value of b: ")
# r_in = input("Insert the inner radius: ")
# r_out = input("The outer radius: ")
# cs = input("And the sound speed: ")
b = 1.50 #[0.75, 1.00, 1.50, 2.00]
incl_deg = [0.0, 20.0, 45.0, 60.0, 90.0] #[0.0, 20.0, 45.0, 60.0, 90.0]
r_in = 0.1 #[0.03, 0.1, 1.0]
r_out = 9.5 #[5.0, 9.5]
cs = 10 #3 5
species = 'SIIa'

quantity = incl_deg
v = []
line_flux = []
path_file = []
plt.figure()
for i in range(len(quantity)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg[i], 2)))
    v = -np.array(map(float, [lines.split()[0] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux = np.array(map(float, [lines.split()[1] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))

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
# plt.show()

b = 1.00
incl_deg = [0.0, 20.0, 45.0, 60.0, 75.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
species = 'NeII'

v = []
line_flux = []
path_file = []
v_hydro = []
line_flux_hydro = []
path_file_hydro = []
for i in range(len(incl_deg)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg[i], 2)))
    v.append(map(float, [lines.split()[0] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux.append(map(float, [lines.split()[1] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    # v = -np.array(map(float, [lines.split()[0] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    # line_flux = np.array(map(float, [lines.split()[1] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    path_file_hydro.append('../data_hydro/incl_'+str(round(incl_deg[i], 2)))
    v_hydro.append(map(float, [lines.split()[0] for lines in open(str(path_file_hydro[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux_hydro.append(map(float, [lines.split()[1] for lines in open(str(path_file_hydro[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))

# fig = plt.figure()
# gs = gridspec.GridSpec(3, 2, figure=fig) #wspace=0.1, hspace=0.15,
# for j in range(0,3):
#     for k in range(0,2):
#         ax = plt.subplot(gs[j, k])
#         plt.plot(v[2*j+k], line_flux[2*j+k] / np.amax(line_flux[2*j+k]), color='red', label='$i='+str(incl_deg[2*j+k])+'$')
#         plt.plot(v_hydro[2*j+k], line_flux_hydro[2*j+k] / np.amax(line_flux_hydro[2*j+k]), color='blue')
#         ax.set_xlabel(r'$v [\frac{km}{s}]$', fontsize=10)
#         ax.set_ylabel(r'$Normalized \, L(v)$', fontsize=10)
#         fig.add_subplot(ax)
#         # plt.title(r'$i=$'+str(incl_deg[2*j+k]))
#         leg = ax.legend(loc='upper right', frameon=False, handlelength=0, handletextpad=0)
#         for item in leg.legendHandles:
#             item.set_visible(False)
# plt.savefig('./lineprofiles/'+str(species)+'/line_profile_comp.eps', format='eps', bbox_inches='tight', dpi=300)
# plt.show()

# fig = plt.figure()
fig, ax = plt.subplots(3, 2, sharex='col', sharey='row')
for row in range(3):
    for col in range(2):
        ax[row,col].plot(v[2*row+col], line_flux[2*row+col] / np.amax(line_flux[2*row+col]), color='red', label='$i='+str(incl_deg[2*row+col])+'$')
        ax[row,col].plot(v_hydro[2*row+col], line_flux_hydro[2*row+col] / np.amax(line_flux_hydro[2*row+col]), color='blue')
        leg = ax[row,col].legend(loc='upper right', frameon=False, handlelength=0, handletextpad=0)
        for item in leg.legendHandles:
            item.set_visible(False)
        ax[row,col].axis([-39.,39.,0.,1.1])
ax[2,0].set_xlabel(r'$v [\frac{km}{s}]$')
ax[2,1].set_xlabel(r'$v [\frac{km}{s}]$')
ax[0,0].set_ylabel(r'$Normalized \, L(v)$')
ax[1,0].set_ylabel(r'$Normalized \, L(v)$')
ax[2,0].set_ylabel(r'$Normalized \, L(v)$')
plt.subplots_adjust(hspace=0., wspace=0.)
plt.savefig('./lineprofiles/'+str(species)+'/line_profile_comp.eps', format='eps', bbox_inches='tight', dpi=300)
plt.show()

# plt.figure()
# plt.plot(v, line_flux / np.amax(line_flux), color='r', label='model')
# plt.plot(v_hydro, line_flux_hydro / np.amax(line_flux_hydro), color='b', label='hydro')
# # plt.title(r'Line profile', fontsize=15)
# plt.axis([-40.,40.,0.,1.1])
# plt.xlabel(r'$v [\frac{km}{s}]$', fontsize=15)
# plt.ylabel(r'$Normalized \, L(v)$', fontsize=15)
# # plt.legend(loc='best')
# # plt.savefig('./lineprofiles/'+str(species)+'/line_profile_comp_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'_cs'+str(cs)+'_i'+str(round(incl_deg, 2))+'.png', format='png', bbox_inches='tight', dpi=300)
# plt.savefig('./lineprofiles/'+str(species)+'/line_profile_comp_i'+str(round(incl_deg, 1))+'.eps', format='eps', bbox_inches='tight', dpi=300)
# plt.show()
