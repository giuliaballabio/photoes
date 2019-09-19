import numpy as np
import matplotlib.pyplot as plt
import re

plt.style.use('classic')
plt.rcParams['font.family'] = 'serif'
# plt.rcParams['font.size'] = 10
# plt.rcParams['axes.labelsize'] = 10
# plt.rcParams['axes.titlesize'] = 10
# plt.rcParams['xtick.labelsize'] = 8
# plt.rcParams['ytick.labelsize'] = 8
# plt.rcParams['legend.fontsize'] = 10
# plt.rcParams['figure.titlesize'] = 12

b = [0.75, 1.00, 1.50] #, 2.00]
incl_deg = 90.0
r_in = 0.1
r_out = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
cs = 10
R = 3.e4
species = 'NeII'

## N.B. When you change to nonorm, remember to change also f.readlines()[24:]
## for mdot10e-8 etc is f.readlines()[33:]
mdot = 'mdot10e-9'

path_file = []
for j in range(len(b)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[j], 2)))+'_r'+str(r_in)+'_r'+str(r_out))


## ---------------------- PLOT FLUX AS A FUNCTION OF THE OUTER RADIUS FOR DIFFERENT b ----------------------

value = []
for i in range(len(r_out)):
    with open('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[0], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'/incl_'+str(round(incl_deg,2))+'/photoes_b'+str('{:.2f}'.format(round(b[0], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'_i'+str(round(incl_deg,2))+'.txt', 'r') as f:
        # lines = f.readlines()[24:]
        lines = f.readlines()[33:]
        value.append([x.split('\n')[0] for x in lines])
f.close()
string = []
flux1 = []
for i in range(len(r_out)):
    string.append(value[i][0])
    # flux1.append(string[i][18:24])
    flux1.append(float(string[i][18:40]))
flux1 = np.array(flux1)

value = []
for i in range(len(r_out)):
    with open('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[1], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'/incl_'+str(round(incl_deg,2))+'/photoes_b'+str('{:.2f}'.format(round(b[1], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'_i'+str(round(incl_deg,2))+'.txt', 'r') as f:
        # lines = f.readlines()[24:]
        lines = f.readlines()[33:]
        value.append([x.split('\n')[0] for x in lines])
f.close()
string = []
flux2 = []
for i in range(len(r_out)):
    string.append(value[i][0])
    # flux2.append(string[i][18:24])
    flux2.append(float(string[i][18:40]))
flux2 = np.array(flux2)

value = []
for i in range(len(r_out)):
    with open('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[2], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'/incl_'+str(round(incl_deg,2))+'/photoes_b'+str('{:.2f}'.format(round(b[2], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'_i'+str(round(incl_deg,2))+'.txt', 'r') as f:
        # lines = f.readlines()[24:]
        lines = f.readlines()[33:]
        value.append([x.split('\n')[0] for x in lines])
f.close()
string = []
flux3 = []
for i in range(len(r_out)):
    string.append(value[i][0])
    # flux3.append(string[i][18:24])
    flux3.append(float(string[i][18:40]))
flux3 = np.array(flux3)

plt.figure()
plt.plot(r_out, flux1, color='#addd8e', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#addd8e', label='$b=0.75$')
plt.plot(r_out, flux2, color='#31a354', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#31a354', label='$b=1.00$')
plt.plot(r_out, flux3, color='#006837', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#006837', label='$b=1.50$')
plt.xlabel(r'$R_{out} / R_g$', fontsize=15)
plt.ylabel(r'$L_{NeII} [L_{\odot}]$', fontsize=15)
plt.xticks(np.arange(min(r_out), max(r_out)+0.5, 0.5))
# plt.title('R$_{in}$ = '+str(r_in)+' au - i = '+str(round(incl_deg,2)))
# plt.axis([4.5, 10.0, 0.e-6, 6.e-6])
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/flux_cs'+str(cs)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/flux_cs'+str(cs)+'_'+str(mdot)+'.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.show()


## ---------------------- PLOT FLUX AS FUNCTION OF OUTER RADIUS FOR DIFFERENT DENSITY NORMALISATION FACTORS ----------------------

b = 1.00
incl_deg = 90.0
r_in = 0.1
r_out = [5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
cs = 10
R = 3.e4
species = 'NeII'
mdot = ['mdot10e-10', 'mdot10e-9', 'mdot10e-8']

path_file = []
for j in range(len(mdot)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot[j])+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out))

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg)):
    with open(str(path_file[0])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f1:
        lines = f1.readlines()[10:]
        v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f1.close()
v_peak2 = []
v_centr2 = []
fwhm2 = []
for i in range(len(incl_deg)):
    with open(str(path_file[1])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f2:
        lines = f2.readlines()[10:]
        v_peak2.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr2.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm2.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f2.close()
v_peak3 = []
v_centr3 = []
fwhm3 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f3:
        lines = f3.readlines()[10:]
        v_peak3.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr3.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm3.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f3.close()

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#fc9272', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fc9272', label='$\dot{M}(<25 \, au) = 10^{-10} M_{\odot}/yr$')
plt.plot(incl_deg, np.abs(v_centr2), color='#ef3b2c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#ef3b2c', label='$\dot{M}(<25 \, au) = 10^{-9} M_{\odot}/yr$')
plt.plot(incl_deg, np.abs(v_centr3), color='#99000d', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#99000d', label='$\dot{M}(<25 \, au) = 10^{-8} M_{\odot}/yr$')
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 6.5])
plt.legend(loc='best')
# plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vcentr_densitynorm_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_densitynorm_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()
