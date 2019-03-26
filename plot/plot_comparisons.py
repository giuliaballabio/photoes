import numpy as np
import matplotlib.pyplot as plt

b = 1.00
incl_deg1 = [0.0, 1.0, 5.0, 10.0, 20.0, 35.0, 45.0, 50.0, 60.0, 75.0, 90.0]
incl_deg2 = [0.0, 10.0, 20.0, 35.0, 45.0, 60.0, 70.0, 80.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = [3, 5, 10]
R = 3.e4

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg2)):
    with open('../cs'+str(cs[0])+'kms/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg2[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f1:
        lines = f1.readlines()[10:]

        v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f1.close()
v_peak2 = []
v_centr2 = []
fwhm2 = []
for i in range(len(incl_deg2)):
    with open('../cs'+str(cs[1])+'kms/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg2[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f2:
        lines = f2.readlines()[10:]

        v_peak2.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr2.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm2.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f2.close()
v_peak3 = []
v_centr3 = []
fwhm3 = []
for i in range(len(incl_deg1)):
    with open('../cs'+str(cs[2])+'kms/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg1[i],2))+'/observables_R'+str(R)+'.txt', 'r') as f3:
        lines = f3.readlines()[10:]

        v_peak3.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr3.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm3.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f3.close()

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
plt.plot(incl_deg2, np.abs(v_centr1), color='#c6dbef', linestyle='dashed', marker='o', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg2, np.abs(v_centr2), color='#2171b5', linestyle='dashed', marker='o', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg1, np.abs(v_centr3), color='#08306b', linestyle='dashed', marker='o', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='hydro sim')
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', label='Sacco et al. (2012)')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', label='Pascucci & Sterzik (2009)')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{centroid} \, [km/s]$', fontsize=15)
plt.xticks(np.arange(min(incl_deg1), max(incl_deg1)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, -1.0, 14.0])
plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/vpeak_b'+str(b)+'_r'+str(r_in)+'_r'+str(r_out)+'_R'+str(R)+'.png', format='png', bbox_inches='tight')
plt.show()


cs = 10
R = [3.e4, 10.e4]

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg1)):
    with open('../cs'+str(cs)+'kms/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg1[i],2))+'/observables_R'+str(R[0])+'.txt', 'r') as f1:
        lines = f1.readlines()[10:]

        v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f1.close()
v_peak2 = []
v_centr2 = []
fwhm2 = []
for i in range(len(incl_deg1)):
    with open('../cs'+str(cs)+'kms/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg1[i],2))+'/observables_R'+str(R[1])+'.txt', 'r') as f2:
        lines = f2.readlines()[10:]

        v_peak2.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr2.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm2.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f2.close()

plt.figure()
plt.plot(incl_deg1, fwhm1, color='#8c96c6', linestyle='dashed', marker='o', markeredgecolor='#8c96c6', label='R = '+str(R[0]))
plt.plot(incl_deg1, fwhm2, color='#8856a7', linestyle='dashed', marker='o', markeredgecolor='#8856a7', label='R = '+str(R[1]))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'FWHM', fontsize=15)
plt.xticks(np.arange(min(incl_deg1), max(incl_deg1)+10., 10.0))
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-5.0, 95.0, 0.0, 30.0])
plt.legend(loc='best')
plt.savefig('./observables/fwhm_b'+str(b)+'_r'+str(r_in)+'_r'+str(r_out)+'_cs'+str(cs)+'.png', format='png', bbox_inches='tight')
plt.show()
