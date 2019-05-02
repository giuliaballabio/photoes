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
b = [0.75, 1.00, 1.50, 2.00]
incl_deg = 4.0 #[0.0, 20.0, 45.0, 60.0, 90.0] #[0.0, 20.0, 45.0, 60.0, 90.0]
r_in = 0.03 #[0.03, 0.1, 1.0]
r_out = 9.9 #[5.0, 9.5]
cs = 10 #3 5
species = 'TWHya_NeII'

quantity = b
v = []
line_flux = []
path_file = []
plt.figure()
for i in range(len(quantity)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/data_b'+str('{:.2f}'.format(round(b[i], 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2)))
    v = -np.array(map(float, [lines.split()[0] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
    line_flux = np.array(map(float, [lines.split()[1] for lines in open(str(path_file[i])+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))

    c = [float(i)/float(len(quantity)), 0.0, float(len(quantity)-i)/float(len(quantity))]
    # if incl_deg[i] == 90.0:
    #     v = -1. * v
    # else:
    #     v = v
    plt.plot(v, line_flux / np.amax(line_flux), color=c, label='b='+str(quantity[i]))
    # plt.vlines(0., 0., 1.3, 'k', linestyle='--')
    plt.title(r'Line profile', fontsize=15)
    plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
    plt.ylabel(r'Normalized L(v)', fontsize=15)
    plt.axis([-40., 40., 0., 1.2])
    plt.legend(loc='best')
    plt.savefig('./lineprofiles/'+str(species)+'/line_profile_i'+str(round(incl_deg, 2))+'_r'+str(r_in)+'_r'+str(r_out)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
    # plt.savefig('./lineprofiles/'+str(species)+'/line_profile_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
    # plt.savefig('./lineprofiles/'+str(species)+'/line_profile_b'+str('{:.2f}'.format(round(b, 2)))+'_i'+str(round(incl_deg, 2))+'_r'+str(r_out)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
    # plt.savefig('./lineprofiles/'+str(species)+'/line_profile_b'+str('{:.2f}'.format(round(b, 2)))+'_i'+str(round(incl_deg, 2))+'_r'+str(r_in)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
plt.show()

b = [0.75, 1.00, 1.50, 2.00]
incl_deg = [0.0, 20.0, 45.0, 60.0, 75.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
species = 'NeII'

path_file = []
path_file_hydro = []
for j in range(len(b)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/mdot10e-9/data_b'+str('{:.2f}'.format(round(b[j], 2)))+'_r'+str(r_in)+'_r'+str(r_out))

v1 = []
line_flux1 = []
v2 = []
line_flux2 = []
v3 = []
line_flux3 = []
v4 = []
line_flux4 = []
v_hydro = []
line_flux_hydro = []
for i in range(len(incl_deg)):
    v1.append(map(float, [lines.split()[0] for lines in open(str(path_file[0])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux1.append(map(float, [lines.split()[1] for lines in open(str(path_file[0])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    v2.append(map(float, [lines.split()[0] for lines in open(str(path_file[1])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux2.append(map(float, [lines.split()[1] for lines in open(str(path_file[1])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    v3.append(map(float, [lines.split()[0] for lines in open(str(path_file[2])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux3.append(map(float, [lines.split()[1] for lines in open(str(path_file[2])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    v4.append(map(float, [lines.split()[0] for lines in open(str(path_file[3])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux4.append(map(float, [lines.split()[1] for lines in open(str(path_file[3])+'/incl_'+str(round(incl_deg[i], 2))+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))

    path_file_hydro.append('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg[i], 2)))
    v_hydro.append(map(float, [lines.split()[0] for lines in open(str(path_file_hydro[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))
    line_flux_hydro.append(map(float, [lines.split()[1] for lines in open(str(path_file_hydro[i])+'/line_profile_i'+str(round(incl_deg[i], 2))+'.txt', 'r')]))

fig, ax = plt.subplots(3, 2, sharex='col', sharey='row')
for row in range(3):
    for col in range(2):
        ax[row,col].plot(v1[2*row+col], line_flux1[2*row+col] / np.amax(line_flux1[2*row+col]), color='#a50f15', label='$b='+str('{:.2f}'.format(round(b[0], 2)))+'$')
        ax[row,col].plot(v2[2*row+col], line_flux2[2*row+col] / np.amax(line_flux2[2*row+col]), color='#de2d26', label='$b='+str('{:.2f}'.format(round(b[1], 2)))+'$')
        ax[row,col].plot(v3[2*row+col], line_flux3[2*row+col] / np.amax(line_flux3[2*row+col]), color='#fb6a4a', label='$b='+str('{:.2f}'.format(round(b[2], 2)))+'$')
        # ax[row,col].plot(v4[2*row+col], line_flux4[2*row+col] / np.amax(line_flux4[2*row+col]), color='#fcae91', label='$b='+str('{:.2f}'.format(round(b[3], 2)))+'$')
        ax[row,col].plot(v_hydro[2*row+col], line_flux_hydro[2*row+col] / np.amax(line_flux_hydro[2*row+col]), color='blue', label='$hydro$') #, label='$i='+str(incl_deg[2*row+col])+'$')
        leg1 = ax[0,0].legend(loc='upper left')
        leg2 = ax[row,col].legend(['$i='+str(incl_deg[2*row+col])+'$'], loc='upper right', frameon=False, handlelength=0, handletextpad=0)
        for item in leg2.legendHandles:
            item.set_visible(False)
        ax[row,col].add_artist(leg2)
        ax[row,col].axis([-39.,39.,0.,1.1])
ax[2,0].set_xlabel(r'$v [\frac{km}{s}]$')
ax[2,1].set_xlabel(r'$v [\frac{km}{s}]$')
ax[0,0].set_ylabel(r'$Normalized \, L(v)$')
ax[1,0].set_ylabel(r'$Normalized \, L(v)$')
ax[2,0].set_ylabel(r'$Normalized \, L(v)$')
plt.subplots_adjust(hspace=0., wspace=0.)
plt.savefig('./lineprofiles/'+str(species)+'/line_profile_comp_mdot10e-9.eps', format='eps', bbox_inches='tight', dpi=300)
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
