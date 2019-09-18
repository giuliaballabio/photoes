import numpy as np
import matplotlib.pyplot as plt
# import matplotlib as mpl
# print(mpl.font_manager.get_cachedir())

plt.style.use('classic')

## ---------------------- PLOT OBSERVABLES FOR DIFFERENT SPECIES ----------------------

b = 1.00
incl_deg = [0.0, 20.0, 45.0, 60.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
R = 3.e4
species = ['NeII', 'SIIa', 'SIIc', 'OI']

path_file = []
for j in range(len(species)):
    path_file.append('../data_hydro/'+str(species[j]))

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
v_peak4 = []
v_centr4 = []
fwhm4 = []
for i in range(len(incl_deg)):
    with open(str(path_file[3])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f4:
        lines = f4.readlines()[10:]
        v_peak4.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr4.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm4.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f4.close()

## UPLOAD THE DATA FROM RICHARD'S HYDRO SIMULATIONS
vpeak_hydro = np.array(map(float, [lines.split()[3] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))
incl_hydro = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))

## CONSIDER THE DATA FROM Sacco et al. 2012
## RX J1615.3-3255 IS TAKEN FROM de Boer et al. 2016
vpeak_data = [-10.5, -4.4, -7.5, -8.3, -10.5] # km/s
err_vpeak = [2.7, 2.1, 2.8, 2.7, 2.0]
incl_data = [75.0, 30.0, 47.0, 20.0, 35.0]
#incl_data = [75.0, 30.0, 5.0, 20.0, 35.0]
ID = ['T Cha', 'MP Mus', 'RX J1615.3-3255', 'SR 21', 'V4046 Sgr']

## WE PLOT ALSO THE DATA FROM Pascucci & Sterzik 2009
vpeak_data2 = [-6.2, -3.3, -4.7]
err_vpeak2 = [0.3, 0.7, 2.5]
incl_data2 = [4.0, 45.0, 75.0]
name = ['TW Hya', 'CS Cha', 'T Cha']

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#c6dbef', linestyle='dashed', marker='o', markeredgecolor='#c6dbef', label=str(species[0]))
plt.plot(incl_deg, np.abs(v_centr4), color='#6baed6', linestyle='dashed', marker='o', markeredgecolor='#6baed6', label=str(species[3]))
plt.plot(incl_deg, np.abs(v_centr2), color='#2171b5', linestyle='dashed', marker='o', markeredgecolor='#2171b5', label='SII 406.98nm')
plt.plot(incl_deg, np.abs(v_centr3), color='#08306b', linestyle='dashed', marker='o', markeredgecolor='#08306b', label='SII 671.83nm')
plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='hydro sim')
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='Sacco et al. (2012)')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='Pascucci & Sterzik (2009)')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg1), max(incl_deg1)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5., 95., -1., 14.])
plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/comparisons/vpeak_hydro_species.png', format='png', dpi=300, bbox_inches='tight')
plt.show()


plt.figure()
plt.plot(incl_deg, fwhm1, color='#b3cde3', linestyle='dashed', marker='o', markeredgecolor='#b3cde3', label=str(species[0]))
plt.plot(incl_deg, fwhm4, color='#8c96c6', linestyle='dashed', marker='o', markeredgecolor='#8c96c6', label=str(species[3]))
plt.plot(incl_deg, fwhm2, color='#8856a7', linestyle='dashed', marker='o', markeredgecolor='#8856a7', label='SII 406.98nm')
plt.plot(incl_deg, fwhm3, color='#810f7c', linestyle='dashed', marker='o', markeredgecolor='#810f7c', label='SII 671.83nm')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$FWHM$', fontsize=15)
plt.xticks(np.arange(min(incl_deg1), max(incl_deg1)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, 0.0, 30.0])
plt.legend(loc='best')
plt.savefig('./observables/comparisons/fwhm_hydro_species.png', format='png', dpi=300, bbox_inches='tight')
plt.show()


## ---------------- PLOT THE VELOCITY AT PEAK AND FWHM AS FUNCTIONS OF THE INCLINATION FOR A SINGLE SPECIE ----------------------

incl_deg = [0.0, 20.0, 45.0, 60.0, 90.0]
cs = 10
R = 3.e4
species = 'SIIa'

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg)):
    with open('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f1:
    # with open('../data_hydro_midplane/incl_'+str(round(incl_deg[i],2))+'/observables.txt', 'r') as f1:
        lines = f1.readlines()[10:]
        v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f1.close()


## UPLOAD THE DATA FROM RICHARD'S HYDRO SIMULATIONS
vpeak_hydro = np.array(map(float, [lines.split()[3] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))
incl_hydro = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))

## CONSIDER THE DATA FROM Sacco et al. 2012
## RX J1615.3-3255 IS TAKEN FROM de Boer et al. 2016
vpeak_data = [-10.5, -4.4, -7.5, -8.3, -10.5] # km/s
err_vpeak = [2.7, 2.1, 2.8, 2.7, 2.0]
incl_data = [75.0, 30.0, 47.0, 20.0, 35.0]
#incl_data = [75.0, 30.0, 5.0, 20.0, 35.0]
ID = ['T Cha', 'MP Mus', 'RX J1615.3-3255', 'SR 21', 'V4046 Sgr']

## WE PLOT ALSO THE DATA FROM Pascucci & Sterzik 2009
vpeak_data2 = [-6.2, -3.3, -4.7]
err_vpeak2 = [0.3, 0.7, 2.5]
incl_data2 = [4.0, 45.0, 75.0]
name = ['TW Hya', 'CS Cha', 'T Cha']

plt.figure()
# plt.plot(incl_deg, np.abs(v_peak1), color='#006837', linestyle='dashed', marker='o', markeredgecolor='#006837', label='hydro')
plt.plot(incl_deg, np.abs(v_centr1), color='#d7301f', linestyle='dashed', marker='o', markeredgecolor='#d7301f')#, label='hydro $v_{centr}$')
plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted')#, label='Alexander')
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', label='Sacco et al. (2012)')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', label='Pascucci & Sterzik (2009)')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.axis([-5.0, 95.0, -1.0, 14.0])
plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
# plt.savefig('./observables/'+str(species)+'/vpeak_hydro_midplane.png', format='png', bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/vcentr_hydro_R'+str(R)+'.eps', format='eps', bbox_inches='tight', dpi=300)
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#8856a7', linestyle='dashed', marker='o', markeredgecolor='#8856a7', label='hydro')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'FWHM', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.axis([-5.0, 95.0, 5.0, 16.0])
plt.legend(loc='best')
# plt.savefig('./observables/'+str(species)+'/fwhm_hydro_midplane.png', format='png', bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/fwhm_hydro_R'+str(R)+'.png', format='png', bbox_inches='tight', dpi=300)
plt.show()
