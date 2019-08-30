import numpy as np
import matplotlib.pyplot as plt
from physics_constant import *

plt.style.use('classic')
plt.rcParams['font.family'] = 'sans' #Courier New'
plt.rcParams['font.serif'] = 'Helvetica'
# plt.rcParams['font.monospace'] = 'Ubuntu Mono'
#plt.rcParams['font.size'] = 10
#plt.rcParams['axes.labelsize'] = 10
#plt.rcParams['axes.titlesize'] = 10
#plt.rcParams['xtick.labelsize'] = 8
#plt.rcParams['ytick.labelsize'] = 8
#plt.rcParams['legend.fontsize'] = 10
#plt.rcParams['figure.titlesize'] = 12

## ---------------------- PLOT OBSERVABLES FOR DIFFERENT SOUND SPEEDS ----------------------

b = 1.00
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = [3, 5, 10]
R = 3.e4
species = 'NeII'
mdot='mdot10e-8'

path_file = []
for j in range(len(cs)):
    path_file.append('../cs'+str(cs[j])+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out))

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

## UPLOAD THE DATA FROM RICHARD'S HYDRO SIMULATIONS
vpeak_hydro = np.array(map(float, [lines.split()[3] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))
fwhm_hydro = np.array(map(float, [lines.split()[4] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))
incl_hydro = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))

## CONSIDER THE DATA FROM Sacco et al. 2012
## RX J1615.3-3255 IS TAKEN FROM de Boer et al. 2016
vpeak_data = [-10.5, -4.4, -7.5, -8.3, -10.5] # km/s
err_vpeak = [2.7, 2.1, 2.8, 2.7, 2.0]
incl_data = [75.0, 30.0, 47.0, 20.0, 35.0]
#incl_data = [75.0, 30.0, 5.0, 20.0, 35.0]
ID = ['T Cha', 'MP Mus', 'RX J1615.3-3255', 'SR 21', 'V4046 Sgr']

## WE PLOT ALSO THE DATA FROM Pascucci & Sterzik 2009
vpeak_data2 = [-6.2, -3.3, -4.]
err_vpeak2 = [0.3, 0.7, 2.]
incl_data2 = [4.0, 45.0, 75.0]
name = ['TW Hya', 'CS Cha', 'T Cha']

## MORE DATA FROM Baldovic-Saavedra et al. 2012
vpeak_data3 = [-2.13, -3.6, -2.9]
err_vpeak3 = [6.5, 1.5, 0.7]
incl_data3 = [60.0, 87.0, 35.0]
name3 = ['V892 Tau', 'CoKu Tau 1', 'FS Tau A']

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, -1.0, 6.0])
#plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize='small')
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/vcentr_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peak} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, -1.0, 6.0])
#plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize='small')
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/vpeak_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peak} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1.0, 91.0, -1.0, 16.0])
plt.legend(loc='upper center', bbox_to_anchor=(1., 1.05), fontsize='small')
#plt.legend(loc='upper left')
plt.savefig('./observables/'+str(species)+'/vpeak_soundspeed_data_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, fwhm2, color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, fwhm3, color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
# plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$FWHM$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, 5.0, 20.0])
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

## ---------------------- PLOT THE FWHM FOR DIFFERENT RESOLUTIONS ---------------------- ##

R = [3.e4, 10.e4]

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[0])+'.txt', 'r') as f1:
        lines = f1.readlines()[10:]
        v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f1.close()
v_peak2 = []
v_centr2 = []
fwhm2 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[1])+'.txt', 'r') as f2:
        lines = f2.readlines()[10:]
        v_peak2.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr2.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm2.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f2.close()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#8c96c6', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8c96c6', label='$R = '+str(R[0])+'$')
plt.plot(incl_deg, fwhm2, color='#8856a7', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8856a7', label='$R = '+str(R[1])+'$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$FWHM$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, 0.0, 30.0])
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_resolution_b'+str(b)+'_cs'+str(cs[2])+'_mdot10e-8.png', format='png', dpi=300, bbox_inches='tight')
plt.show()


## ---------------------- PLOT OBSERVABLES FOR DIFFERENT SPECIES ----------------------

b = 1.00
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
R = 3.e4
species = ['NeII', 'SIIa', 'SIIc', 'OI']
mdot='mdot10e-8'

path_file = []
for j in range(len(species)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species[j])+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out))

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

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#b3cde3', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#b3cde3', label='$['+str(species[0])+']$')
plt.plot(incl_deg, np.abs(v_centr4), color='#8c96c6', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8c96c6', label='$['+str(species[3])+']$')
plt.plot(incl_deg, np.abs(v_centr2), color='#8856a7', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8856a7', label='$[SII]\,406.98\,nm$')
plt.plot(incl_deg, np.abs(v_centr3), color='#810f7c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#810f7c', label='$[SII]\,671.83\,nm$')
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='hydro sim')
# plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='Sacco et al. (2012)')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
# plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='Pascucci & Sterzik (2009)')
# for i in range(len(name)):
#     plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5., 95., -1., 7.])
plt.legend(loc='best')
# plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/comparisons/vcentr_species_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#b3cde3', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#b3cde3', label='$['+str(species[0])+']$')
plt.plot(incl_deg, np.abs(v_peak4), color='#8c96c6', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8c96c6', label='$['+str(species[3])+']$')
plt.plot(incl_deg, np.abs(v_peak2), color='#8856a7', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8856a7', label='$[SII]\,406.98\,nm$')
plt.plot(incl_deg, np.abs(v_peak3), color='#810f7c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#810f7c', label='$[SII]\,671.83\,nm$')
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
# plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='Sacco et al. (2012)')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
# plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='Pascucci & Sterzik (2009)')
# for i in range(len(name)):
#     plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peak} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5., 95., -1., 7.0])
plt.legend(loc='best')
# plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/comparisons/vpeak_species_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#b3cde3', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#b3cde3', label='$['+str(species[0])+']$')
plt.plot(incl_deg, fwhm4, color='#8c96c6', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8c96c6', label='$['+str(species[3])+']$')
plt.plot(incl_deg, fwhm2, color='#8856a7', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8856a7', label='$[SII]\,406.98\,nm$')
plt.plot(incl_deg, fwhm3, color='#810f7c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#810f7c', label='$[SII]\,671.83\,nm$')
# plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$FWHM$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, 5.0, 20.0])
plt.legend(loc='best')
plt.savefig('./observables/comparisons/fwhm_species_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

## ---------------------- PLOT THE LAUNCHING RADIUS DERIVED FROM THE FWHM ----------------------
rlaunch1 = []
rlaunch2 = []
rlaunch3 = []
rlaunch4 = []
rlaunch1 = np.array(rlaunch1)
rlaunch2 = np.array(rlaunch2)
rlaunch3 = np.array(rlaunch3)
rlaunch4 = np.array(rlaunch4)
fwhm1 = np.array(fwhm1)
fwhm2 = np.array(fwhm2)
fwhm3 = np.array(fwhm3)
fwhm4 = np.array(fwhm4)

const = 4.*G*Mstar/au
## The units of rlaunch are cm, but the velocity is in km/s
rlaunch1 = const/(fwhm1*fwhm1*10.**5.*10.**5.)
rlaunch2 = const/(fwhm2*fwhm2*10.**5.*10.**5.)
rlaunch3 = const/(fwhm3*fwhm3*10.**5.*10.**5.)
rlaunch4 = const/(fwhm4*fwhm4*10.**5.*10.**5.)

plt.figure()
plt.plot(incl_deg, rlaunch1, color='#b3cde3', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#b3cde3', label='$['+str(species[0])+']$')
plt.plot(incl_deg, rlaunch2, color='#8c96c6', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8c96c6', label='$['+str(species[3])+']$')
plt.plot(incl_deg, rlaunch3, color='#8856a7', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8856a7', label='$[SII]\,406.98\,nm$')
plt.plot(incl_deg, rlaunch4, color='#810f7c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#810f7c', label='$[SII]\,671.83\,nm$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$R_0 \, [au]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
# plt.axis([-5.0, 95.0, 5.0, 20.0])
plt.legend(loc='best')
plt.savefig('./observables/comparisons/rlaunch_species_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

## ---------------------- PLOT THE VELOCITY AT PEAK OF THE NON CONVOLVED LINES FOR DIFFERENT SPECIES ----------------------

b = 1.00
incl_deg = [0.0, 20.0, 45.0, 60.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
R = 3.e4
species = ['NeII', 'SIIa', 'SIIc', 'OI']

path_file = []
for j in range(len(species)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species[j])+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out))

v_peak_line1 = []
for i in range(len(incl_deg)):
    with open(str(path_file[0])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f1:
        line = f1.readlines()[5:6]
        v_peak_line1.append(map(float, [x.split('\t\t\t')[0] for x in line]))
f1.close()
v_peak_line2 = []
for i in range(len(incl_deg)):
    with open(str(path_file[1])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f2:
        line = f2.readlines()[5:6]
        v_peak_line2.append(map(float, [x.split('\t\t\t')[0] for x in line]))
f2.close()
v_peak_line3 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f3:
        line = f3.readlines()[5:6]
        v_peak_line3.append(map(float, [x.split('\t\t\t')[0] for x in line]))
f3.close()
v_peak_line4 = []
for i in range(len(incl_deg)):
    with open(str(path_file[3])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f4:
        line = f4.readlines()[5:6]
        v_peak_line4.append(map(float, [x.split('\t\t\t')[0] for x in line]))
f4.close()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak_line1), color='#b3cde3', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#b3cde3', label='$['+str(species[0])+']$')
plt.plot(incl_deg, np.abs(v_peak_line4), color='#8c96c6', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8c96c6', label='$['+str(species[3])+']$')
plt.plot(incl_deg, np.abs(v_peak_line2), color='#8856a7', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#8856a7', label='$SII 406.98nm$')
plt.plot(incl_deg, np.abs(v_peak_line3), color='#810f7c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#810f7c', label='$SII 671.83nm$')
plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco et al. \, (2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='$Pascucci & Sterzik \, (2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peakline} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5., 95., -1., 14.])
plt.legend(loc='lower right', bbox_to_anchor=(1.26, 0.05), fontsize = 'small')
plt.savefig('./observables/comparisons/vpeakline_species_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

## ---------------------- PLOT OBSERVABLES FOR DIFFERENT DENSITY NORMALISATION FACTORS ----------------------

b = 1.00
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
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
# plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='Sacco et al. (2012)')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
# plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='Pascucci & Sterzik (2009)')
# for i in range(len(name)):
#     plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5., 95., -1., 7.])
plt.legend(loc='best')
# plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vcentr_densitynorm_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#fc9272', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fc9272', label='$\dot{M}(<25 \, au) = 10^{-10} M_{\odot}/yr$')
plt.plot(incl_deg, np.abs(v_peak2), color='#ef3b2c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#ef3b2c', label='$\dot{M}(<25 \, au) = 10^{-9} M_{\odot}/yr$')
plt.plot(incl_deg, np.abs(v_peak3), color='#99000d', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#99000d', label='$\dot{M}(<25 \, au) = 10^{-8} M_{\odot}/yr$')
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
# plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='Sacco et al. (2012)')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
# plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='Pascucci & Sterzik (2009)')
# for i in range(len(name)):
#     plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peak} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5., 95., -1., 7.0])
plt.legend(loc='best')
# plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vpeak_densitynorm_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#fc9272', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fc9272', label='$\dot{M}(<25 \, au) = 10^{-10} M_{\odot}/yr$')
plt.plot(incl_deg, fwhm2, color='#ef3b2c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#ef3b2c', label='$\dot{M}(<25 \, au) = 10^{-9} M_{\odot}/yr$')
plt.plot(incl_deg, fwhm3, color='#99000d', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#99000d', label='$\dot{M}(<25 \, au) = 10^{-8} M_{\odot}/yr$')
# plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$FWHM$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, 5.0, 20.0])
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_densitynorm_b'+str(b)+'_cs'+str(cs)+'_R'+str(R)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()


## ---------------------- PLOT OBSERVABLES FOR DIFFERENT b ----------------------

b = [0.75, 1.00, 1.50, 2.00]
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
R = 3.e4
species = 'NeII'
mdot='mdot10e-8'

path_file = []
for j in range(len(b)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[j], 2)))+'_r'+str(r_in)+'_r'+str(r_out))

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

#vpeak_modelhydro = []
#vcentr_modelhydro = []
#fwhm_modelhydro = []
#for i in range(len(incl_deg)):
#    with open('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f5:
#        lines = f5.readlines()[10:]
#        vpeak_modelhydro.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
#        vcentr_modelhydro.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
#        fwhm_modelhydro.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
#f5.close()

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$') #a50f15
plt.plot(incl_deg, np.abs(v_centr2), color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$') #de2d26
plt.plot(incl_deg, np.abs(v_centr3), color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$') #fb6a4a
plt.plot(incl_deg, np.abs(v_centr4), color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$') #fcae91
# plt.plot(incl_deg, np.abs(vcentr_modelhydro), color='#31a354', linestyle='dashed', marker='o', markeredgecolor='#31a354', label='$hydro$')
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
# plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='Sacco et al. (2012)')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
# plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='Pascucci & Sterzik (2009)')
# for i in range(len(name)):
#     plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-2., 92., -0.5, 6.5])
plt.legend(loc='best') #'upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$') #a50f15
plt.plot(incl_deg, np.abs(v_peak2), color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$') #de2d26
plt.plot(incl_deg, np.abs(v_peak3), color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$') #fb6a4a
plt.plot(incl_deg, np.abs(v_peak4), color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$') #fcae91
# plt.plot(incl_deg, np.abs(vpeak_modelhydro), color='#31a354', linestyle='dashed', marker='o', markeredgecolor='#31a354', label='$hydro$')
plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
# plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='Sacco et al. (2012)')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
# plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', capsize=3, label='Pascucci & Sterzik (2009)')
# for i in range(len(name)):
#     plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peak} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-2., 92., -0.5, 6.5])
plt.legend(loc='best') #'upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vpeak_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$')
plt.plot(incl_deg, fwhm2, color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$')
plt.plot(incl_deg, fwhm3, color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$')
plt.plot(incl_deg, fwhm4, color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$')
#plt.plot(incl_deg, fwhm_modelhydro, color='#31a354', linestyle='dashed', marker='o', markeredgecolor='#31a354', label='$hydro$')
#plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$FWHM$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, 0.0, 40.0])
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.show()

## ---------------------- PLOT THE FWHM FOR DIFFERENT Rout ---------------------- ##

# b = [0.75, 1.00, 1.50, 2.00]
# incl_deg = 90.0
# r_in = 0.1
# # r_out = [0.9, 1.0, 1.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5]
# r_out = [0.9, 1.0, 1.5, 5.0, 9.5]
# cs = 10
# R = 3.e4
# species = 'NeII'
# mdot='mdot10e-8'
#
# path_file = []
# for j in range(len(b)):
#     path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[j], 2)))+'_r'+str(r_in)+'_r'+str(r_out))
#
# v_peak1 = []
# v_centr1 = []
# fwhm1 = []
# for i in range(len(r_out)):
#     with open(str(path_file[0])+'/incl_'+str(round(r_out,2))+'/observables_R'+str(R)+'.txt', 'r') as f1:
#         lines = f1.readlines()[10:]
#         v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
#         v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
#         fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
# f1.close()
# v_peak2 = []
# v_centr2 = []
# fwhm2 = []
# for i in range(len(r_out)):
#     with open(str(path_file[1])+'/incl_'+str(round(r_out[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f2:
#         lines = f2.readlines()[10:]
#         v_peak2.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
#         v_centr2.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
#         fwhm2.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
# f2.close()
# v_peak3 = []
# v_centr3 = []
# fwhm3 = []
# for i in range(len(r_out)):
#     with open(str(path_file[2])+'/incl_'+str(round(r_out[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f3:
#         lines = f3.readlines()[10:]
#         v_peak3.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
#         v_centr3.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
#         fwhm3.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
# f3.close()
# v_peak4 = []
# v_centr4 = []
# fwhm4 = []
# for i in range(len(r_out)):
#     with open(str(path_file[2])+'/incl_'+str(round(r_out[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f4:
#         lines = f4.readlines()[10:]
#         v_peak4.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
#         v_centr4.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
#         fwhm4.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
# f4.close()
#
# plt.figure()
# plt.plot(r_out, fwhm1, color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$')
# plt.plot(r_out, fwhm2, color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$')
# plt.plot(r_out, fwhm3, color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$')
# plt.plot(r_out, fwhm4, color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$')
# plt.xlabel(r'$R_{out} [R_g]$', fontsize=15)
# plt.ylabel(r'$FWHM$', fontsize=15)
# plt.xticks(np.arange(min(r_out), max(r_out)+0.5, 0.5))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
# plt.axis([4.5, 10.0, 0.0, 30.0])
# plt.legend(loc='best')
# plt.savefig('./observables/'+str(species)+'/fwhm_rout_b'+str(b)+'_cs'+str(cs)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
# plt.show()
