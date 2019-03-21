import numpy as np
import matplotlib.pyplot as plt
import re


incl_deg = [0.0, 10.0, 20.0, 45.0, 60.0, 75.0, 90.0]
cs = 10

## ---------------- PLOT THE VELOCITY AT PEAK AND FWHM AS FUNCTIONS OF THE INCLINATION ----------------------

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg)):
    with open('../data_hydro/incl_'+str(round(incl_deg[i],2))+'/observables.txt', 'r') as f1:
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
plt.plot(incl_deg, np.abs(v_peak1), color='#006837', linestyle='dashed', marker='o', markeredgecolor='#006837', label='hydro')
plt.plot(incl_deg, np.abs(v_centr1), color='#d7301f', linestyle='dashed', marker='o', markeredgecolor='#d7301f', label='hydro $v_{centr}$')
plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='Alexander')
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', label='Sacco et al. (2012)')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', label='Pascucci & Sterzik (2009)')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peak} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.axis([-5.0, 95.0, -1.0, 14.0])
plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/vpeak_hydro.png', format='png', bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#8856a7', linestyle='dashed', marker='o', markeredgecolor='#8856a7', label='hydro')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'FWHM', fontsize=15)
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.axis([-5.0, 95.0, 5.0, 30.0])
plt.legend(loc='best')
plt.savefig('./observables/fwhm_hydro.png', format='png', bbox_inches='tight')
plt.show()
