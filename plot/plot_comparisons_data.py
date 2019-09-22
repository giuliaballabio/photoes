import numpy as np
import matplotlib.pyplot as plt
import re
from physics_constant import *

plt.style.use('classic')
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='medium')
plt.rc('ytick', labelsize='medium')
plt.rc('axes', titlesize='xx-large')
plt.rc('axes', labelsize='xx-large')
plt.rc('legend', fontsize='large')

## ---------------- COMPARE WITH DATA OF [NeII] LINE FOR DIFFERENT b ----------------------

b = [0.75, 1.00, 1.50, 2.00]
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
R = 3.e4
species = 'NeII'
mdot = 'mdot10e-9'

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

# UPLOAD THE DATA FROM RICHARD'S HYDRO SIMULATIONS
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
vpeak_data2 = [-6.2, -3.3, -4.7]
err_vpeak2 = [0.3, 0.7, 2.5]
incl_data2 = [4.0, 45.0, 75.0]
name = ['TW Hya', 'CS Cha', 'T Cha']

## MORE DATA FROM Baldovic-Saavedra et al. 2012
vpeak_data3 = [-2.13, -3.6, -2.9]
err_vpeak3 = [6.5, 1.5, 0.7]
incl_data3 = [60.0, 87.0, 35.0]
name3 = ['V892 Tau', 'CoKu Tau 1', 'FS Tau A']

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$') #a50f15
plt.plot(incl_deg, np.abs(v_centr2), color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$') #de2d26
plt.plot(incl_deg, np.abs(v_centr3), color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$') #fb6a4a
plt.plot(incl_deg, np.abs(v_centr4), color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$') #fcae91
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 14.])
plt.tight_layout()
plt.legend(loc='best') #'upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$') #a50f15
plt.plot(incl_deg, np.abs(v_peak2), color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$') #de2d26
plt.plot(incl_deg, np.abs(v_peak3), color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$') #fb6a4a
plt.plot(incl_deg, np.abs(v_peak4), color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$') #fcae91
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{peak} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 14.0])
plt.tight_layout()
plt.legend(loc='best') #'upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vpeak_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vpeak_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

## ---------------- COMPARE WITH DATA OF [NeII] LINE FOR DIFFERENT cs ----------------------

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

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{peak} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[NeII] \, 12.81 \mu m$', loc='right')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 17.])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0.1, 1.2), loc='upper left')
# plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/vpeak_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vpeak_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[NeII] \, 12.81 \mu m$')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 14.])
plt.tight_layout()
# plt.legend(loc='upper center', bbox_to_anchor=(1., 1.05), fontsize='small')
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/vcentr_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

## ---------------- COMPARE WITH DATA OF [OI] LINE FOR DIFFERENT b ----------------------

b = [0.75, 1.00, 1.50, 2.00]
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = 10
R = 3.e4
species = 'OI'
mdot = 'mdot10e-9'

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

## CONSIDER THE DATA FROM Banzatti et al. 2019
## Narrow Component
vcentr_data = [0.3, -0.4, -4., 0.9, -2.6, -1.3, -11., -5.3, -0.3, 0.4, -2.5, 5., -4., -12., 0.5, -0.9, -4.6]
err_vcentr = [1.1, 0.5, 1., 0.9, 0.5, 5., 0.7, 0.6, 5., 5., 3.6, 5.1, 1., 5., 0.6, 1.1, 5.]
fwhm_data = [26., 14.5, 25., 18., 24., 13., 18., 25.4, 12., 25., 29., 67., 47., 15., 26., 18., 18.] # km/s
err_fwhm = [2., 0.4, 1.1, 2., 0.6, 1.2, 0.5, 0.9, 1., 0.7, 1.5, 4., 2., 1., 1., 1., 2.]
## Broad Component
fwhm_data1 = [96., 54., 120., 111., 80., 67., 42., 177., 53., 86., 162., 223., 288., 44., 84., 79., 36.] # km/s
err_fwhm1 = [4.6, 2., 9., 6., 2., 2., 2., 7., 5., 3., 10., 17., 28., 6., 6., 8., 6.]

incl_data = [71., 20., 18., 39., 65., 60., 32., 26., 9., 55., 38., 66., 37., 35., 54., 54., 50.]
ID = ['AATau', 'AS205N', 'AS353A', 'BPTau', 'CWTau', 'DFTau', 'DGTau', 'DKTau', 'DRTau', 'FMTau', 'FZTau', 'ITTau', 'RNO90', 'RULup', 'RXJ1842', 'V853Oph', 'VVCrAS']

## CONSIDER THE DATA FROM Rigliaco et al. 2013
fwhm_data2 = [14.1, 42.6, 47.1, 8.2, 42.8, 44.1, 40.5, 34.8, 55.7, 45.1, 57.6, 28.3, 51.1, 28.9, 38.8] # km/s
err_fwhm2 = [2.7, 13.6, 8.4, 0.6, 20.0, 1.5, 6.6, 2., 6.3, 1.1, 6.9, 7.0, 5.6, 8.5, 5.9]
incl_data2 = [37., 30., 35., 4., 60., 80., 52., 63., 78., 35., 70., 37., 23., 45., 42.]
ID2 = ['DR Tau', 'MP Mus', 'V4046 Sgr', 'TW Hya', 'AS 353A', 'CW Tau', 'CY Tau', 'DF Tau', 'DN Tau', 'DQ Tau', 'FM Tau', 'GG Tau', 'GK Tau', 'GM Aur', 'UY Aur']

plt.figure()
plt.plot(incl_deg, fwhm1, color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$')
plt.plot(incl_deg, fwhm2, color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$')
plt.plot(incl_deg, fwhm3, color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$')
# plt.plot(incl_deg, fwhm4, color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$')
#plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(fwhm_data[i])+0.3), color='k')
plt.errorbar(incl_data, np.abs(fwhm_data1), yerr=err_fwhm1, color='#969696', markeredgecolor='None', linestyle='None', marker='*', capsize=3, label='$BC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(fwhm_data1[i])+0.3), color='#969696')
# plt.errorbar(incl_data2, np.abs(fwhm_data2), yerr=err_fwhm2, color='k', linestyle='None', marker='d', capsize=3, label='$Rigliaco\,et\,al.\,(2013)$')
# for i in range(len(ID2)):
#     plt.annotate(ID2[i], (incl_data2[i]+0.3, np.abs(fwhm_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$FWHM$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1.0, 91.0, 0.0, 300.0])
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/fwhm_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$')
plt.plot(incl_deg, fwhm2, color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$')
plt.plot(incl_deg, fwhm3, color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$')
# plt.plot(incl_deg, fwhm4, color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$')
#plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(fwhm_data[i])+0.5), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$FWHM$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., 5., 32.])
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_dataNC.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/fwhm_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_dataNC.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$') #a50f15
plt.plot(incl_deg, np.abs(v_centr2), color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$') #de2d26
plt.plot(incl_deg, np.abs(v_centr3), color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$') #fb6a4a
# plt.plot(incl_deg, np.abs(v_centr4), color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$') #fcae91
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.errorbar(incl_data, np.abs(vcentr_data), yerr=err_vcentr, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(vcentr_data[i])+0.5), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 14.])
plt.tight_layout()
plt.legend(loc='best') #'upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/'+str(species)+'/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

## ---------------- COMPARE WITH DATA OF [OI] LINE FOR DIFFERENT cs ----------------------

b = 1.00
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = [3, 5, 10]
R = 3.e4
species = 'OI'
mdot='mdot10e-9'

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

## CONSIDER THE DATA FROM Banzatti et al. 2019
## Narrow Component
vcentr_data = [0.3, -0.4, -4., 0.9, -2.6, -1.3, -11., -5.3, -0.3, 0.4, -2.5, 5., -4., -12., 0.5, -0.9, -4.6]
err_vcentr = [1.1, 0.5, 1., 0.9, 0.5, 5., 0.7, 0.6, 5., 5., 3.6, 5.1, 1., 5., 0.6, 1.1, 5.]
fwhm_data = [26., 14.5, 25., 18., 24., 13., 18., 25.4, 12., 25., 29., 67., 47., 15., 26., 18., 18.] # km/s
err_fwhm = [2., 0.4, 1.1, 2., 0.6, 1.2, 0.5, 0.9, 1., 0.7, 1.5, 4., 2., 1., 1., 1., 2.]
## Broad Component
fwhm_data1 = [96., 54., 120., 111., 80., 67., 42., 177., 53., 86., 162., 223., 288., 44., 84., 79., 36.] # km/s
err_fwhm1 = [4.6, 2., 9., 6., 2., 2., 2., 7., 5., 3., 10., 17., 28., 6., 6., 8., 6.]

incl_data = [71., 20., 18., 39., 65., 60., 32., 26., 9., 55., 38., 66., 37., 35., 54., 54., 50.]
ID = ['AATau', 'AS205N', 'AS353A', 'BPTau', 'CWTau', 'DFTau', 'DGTau', 'DKTau', 'DRTau', 'FMTau', 'FZTau', 'ITTau', 'RNO90', 'RULup', 'RXJ1842', 'V853Oph', 'VVCrAS']

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data, np.abs(vcentr_data), yerr=err_vcentr, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(vcentr_data[i])+0.5), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 14.])
plt.tight_layout()
# plt.legend(loc='upper center', bbox_to_anchor=(1., 1.05), fontsize='small')
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/vcentr_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(fwhm1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(fwhm2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(fwhm3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(fwhm_data[i])+0.5), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., 5., 32.])
plt.tight_layout()
# plt.legend(loc='upper center', bbox_to_anchor=(1., 1.05), fontsize='small')
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/fwhm_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()
