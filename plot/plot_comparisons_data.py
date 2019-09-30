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

## WE NEED TO PLOT THE FWHM, NOT HALF
fwhm1 = np.array(fwhm1)*2.
fwhm2 = np.array(fwhm2)*2.
fwhm3 = np.array(fwhm3)*2.
fwhm4 = np.array(fwhm4)*2.

## UPLOAD THE DATA FROM RICHARD'S HYDRO SIMULATIONS
vpeak_hydro = np.array(map(float, [lines.split()[3] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))
fwhm_hydro = np.array(map(float, [lines.split()[4] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))
incl_hydro = np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/NeII_inclination.dat', 'r')]))

## CONSIDER THE DATA FROM Sacco et al. 2012
## RX J1615.3-3255 IS TAKEN FROM de Boer et al. 2016
vpeak_data = [-10.5, -4.4, -7.5, -8.3, -10.5] # km/s
err_vpeak = [2.7, 2.1, 2.8, 2.7, 2.0]
fwhm_data = [44.9, 15.9, 20.5, 15.1, 22.5]
err_fwhm = [3.2, 1.4, 2.7, 1.2, 0.5]
incl_data = [75.0, 30.0, 47.0, 20.0, 35.0]
#incl_data = [75.0, 30.0, 5.0, 20.0, 35.0]
ID = ['T Cha', 'MP Mus', 'RX J1615.3-3255', 'SR 21', 'V4046 Sgr']

## WE PLOT ALSO THE DATA FROM Pascucci & Sterzik 2009
vpeak_data2 = [-6.2, -3.3, -4.7]
err_vpeak2 = [0.3, 0.7, 2.5]
fwhm_data2 = [14.6, 27., 42.]
err_fwhm2 = [0.7, 2., 4.]
incl_data2 = [4.0, 45.0, 75.0]
name = ['TW Hya', 'CS Cha', 'T Cha']

## MORE DATA FROM Baldovic-Saavedra et al. 2012
vpeak_data3 = [-2.13, -3.6, -2.9]
err_vpeak3 = [6.5, 1.5, 0.7]
fwhm_data3 = [26., 55.2, 26.8]
err_fwhm3 = [3.4, 3.3, 1.7]
incl_data3 = [60.0, 87.0, 35.0]
name3 = ['V892 Tau', 'CoKu Tau 1', 'FS Tau A']

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$') #a50f15
plt.plot(incl_deg, np.abs(v_centr2), color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$') #de2d26
plt.plot(incl_deg, np.abs(v_centr3), color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$') #fb6a4a
# plt.plot(incl_deg, np.abs(v_centr4), color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$') #fcae91
# plt.plot(incl_hydro, np.abs(vpeak_hydro), color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]-12.0, np.abs(vpeak_data3[i])-0.1))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 17.5])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0., 1.), loc='upper left', ncol=2)
plt.savefig('./observables/'+str(species)+'/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$') #a50f15
plt.plot(incl_deg, np.abs(v_peak2), color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$') #de2d26
plt.plot(incl_deg, np.abs(v_peak3), color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$') #fb6a4a
# plt.plot(incl_deg, np.abs(v_peak4), color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$') #fcae91
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]-12.0, np.abs(vpeak_data3[i])-0.1))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{peak} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 17.5])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0., 1.), loc='upper left', ncol=2)
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

## WE NEED TO PLOT THE FWHM, NOT HALF
fwhm1 = np.array(fwhm1)*2.
fwhm2 = np.array(fwhm2)*2.
fwhm3 = np.array(fwhm3)*2.

plt.figure()
plt.plot(incl_deg, np.abs(v_peak1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_peak3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(vpeak_data3), yerr=err_vpeak3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]-12.0, np.abs(vpeak_data3[i])-0.1))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{peak} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[NeII] \, 12.81 \mu m$')#, loc='right')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 17.5])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0., 1.), loc='upper left', ncol=2)
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
    plt.annotate(name3[i], (incl_data3[i]-12.0, np.abs(vpeak_data3[i])-0.1))
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[NeII] \, 12.81 \mu m$')
# plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
plt.axis([-1., 91., -0.5, 17.5])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0., 1.), loc='upper left', ncol=2)
plt.savefig('./observables/'+str(species)+'/vcentr_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, fwhm2, color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, fwhm3, color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data2, np.abs(fwhm_data2), yerr=err_fwhm2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(fwhm_data2[i])+0.3))
plt.errorbar(incl_data3, np.abs(fwhm_data3), yerr=err_fwhm3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
for i in range(len(name3)):
    plt.annotate(name3[i], (incl_data3[i]-12.0, np.abs(fwhm_data3[i])-0.1))
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(fwhm_data[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$FWHM \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[NeII] \, 12.81 \mu m$')
plt.axis([-1., 91., 9., 59.])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0., 1.), loc='upper left', ncol=2)
plt.savefig('./observables/'+str(species)+'/fwhm_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/fwhm_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

# plt.plot(np.abs(v_centr1), fwhm1, color='#c6dbef', linestyle='None', linewidth=2.5, marker='o', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
# plt.plot(np.abs(v_centr2), fwhm2, color='#2171b5', linestyle='None', linewidth=2.5, marker='o', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
# plt.plot(np.abs(v_centr3), fwhm3, color='#08306b', linestyle='None', linewidth=2.5, marker='o', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
# plt.errorbar(np.abs(vpeak_data2), np.abs(fwhm_data2), xerr=err_vpeak2, yerr=err_fwhm2, color='k', linestyle='None', marker='*', capsize=3, label='$Pascucci\,&\,Sterzik\,(2009)$')
# for i in range(len(name)):
#     plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(fwhm_data2[i])+0.3))
# plt.errorbar(np.abs(vpeak_data3), np.abs(fwhm_data3), xerr=err_vpeak3, yerr=err_fwhm3, color='k', linestyle='None', marker='d', capsize=3, label='$Baldovin-Saavedra\,(2012)$')
# for i in range(len(name3)):
#     plt.annotate(name3[i], (incl_data3[i]-12.0, np.abs(fwhm_data3[i])-0.1))
# plt.errorbar(np.abs(vpeak_data), np.abs(fwhm_data), xerr=err_vpeak, yerr=err_fwhm, color='k', linestyle='None', marker='o', capsize=3, label='$Sacco\,et\,al.\,(2012)$')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(fwhm_data[i])+0.3))
# plt.xlabel(r'$- v_{centroid} \, [km/s]$')
# plt.ylabel(r'$FWHM \, [km/s]$')
# # plt.axis([-1., 91., 5., 17.])
# plt.xscale('log')
# plt.yscale('log')
# plt.tight_layout()
# # plt.legend(loc='best')
# plt.savefig('./observables/'+str(species)+'/obs_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
# plt.savefig('./observables/'+str(species)+'/eps/obs_soundspeed_b'+str(b)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
# plt.show()

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

## WE NEED TO PLOT THE FWHM, NOT HALF
fwhm1 = np.array(fwhm1)*2.
fwhm2 = np.array(fwhm2)*2.
fwhm3 = np.array(fwhm3)*2.
fwhm4 = np.array(fwhm4)*2.

## CONSIDER THE DATA FROM Banzatti et al. 2019
## Narrow Component
vcentr_data = [0.3, -0.4, 0.9, -1.3, -5.3, -0.3, 0.4, -2.5, 5., -4., -12., -4.6]
err_vcentr = [1.1, 0.5, 0.9, 5., 0.6, 5., 5., 3.6, 5.1, 1., 5., 5.]
fwhm_data = [26., 14.5, 18., 13., 25.4, 12., 25., 29., 67., 47., 15., 18.] # km/s
err_fwhm = [2., 0.4, 2., 1.2, 0.9, 1., 0.7, 1.5, 4., 2., 1., 2.]
## Broad Component
fwhm_data1 = [96., 54., 111., 67., 177., 53., 86., 162., 223., 288., 44., 36.] # km/s
err_fwhm1 = [4.6, 2., 6., 2., 7., 5., 3., 10., 17., 28., 6., 6.]

incl_data = [71., 20., 39., 60., 26., 9., 55., 38., 66., 37., 35., 50.]
ID = ['AATau', 'AS205N', 'BPTau', 'DFTau', 'DKTau', 'DRTau', 'FMTau', 'FZTau', 'ITTau', 'RNO90', 'RULup', 'VVCrAS']

## Divide the data in order to plot the ID labels clearly
vcentr_sx = [-4., -11., 0.5]
err_vcentr_sx = [1.,0.7, 0.6]
fwhm_sx = [25., 18., 26.]
err_fwhm_sx =[1.1, 0.5, 1.]
fwhm1_sx = [120., 42., 84.]
err_fwhm1_sx = [9., 2., 6.]
incl_sx = [18., 32., 54.]
ID_sx = ['AS353A', '  DGTau', 'RXJ1842']

vcentr_down = [-2.6, -0.9]
err_vcentr_down = [0.5, 1.1]
fwhm_down = [24., 18.]
err_fwhm_down =[0.6, 1.]
fwhm1_down =[80., 79.]
err_fwhm1_down = [2., 8.]
incl_down = [65., 54.]
ID_down = ['CWtau', 'V853Oph']

## CONSIDER THE DATA FROM Rigliaco et al. 2013
# fwhm_data2 = [14.1, 42.6, 47.1, 8.2, 42.8, 44.1, 40.5, 34.8, 55.7, 45.1, 57.6, 28.3, 51.1, 28.9, 38.8] # km/s
# err_fwhm2 = [2.7, 13.6, 8.4, 0.6, 20.0, 1.5, 6.6, 2., 6.3, 1.1, 6.9, 7.0, 5.6, 8.5, 5.9]
# incl_data2 = [37., 30., 35., 4., 60., 80., 52., 63., 78., 35., 70., 37., 23., 45., 42.]
# ID2 = ['DR Tau', 'MP Mus', 'V4046 Sgr', 'TW Hya', 'AS 353A', 'CW Tau', 'CY Tau', 'DF Tau', 'DN Tau', 'DQ Tau', 'FM Tau', 'GG Tau', 'GK Tau', 'GM Aur', 'UY Aur']

# plt.figure()
# plt.plot(incl_deg, fwhm1, color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$')
# plt.plot(incl_deg, fwhm2, color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$')
# plt.plot(incl_deg, fwhm3, color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$')
# # plt.plot(incl_deg, fwhm4, color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$')
# #plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
# plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(fwhm_data[i])+0.3), color='k')
# plt.errorbar(incl_sx, np.abs(fwhm_sx), yerr=err_fwhm_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
# for i in range(len(ID_sx)):
#     plt.annotate(ID_sx[i], (incl_sx[i]+0.3, np.abs(fwhm_sx[i])+0.3), color='k')
# plt.errorbar(incl_down, np.abs(fwhm_down), yerr=err_fwhm_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
# for i in range(len(ID_down)):
#     plt.annotate(ID_down[i], (incl_down[i]+0.3, np.abs(fwhm_down[i])+0.3), color='k')
# plt.errorbar(incl_data, np.abs(fwhm_data1), yerr=err_fwhm1, color='#969696', markeredgecolor='None', linestyle='None', marker='*', capsize=3, label='$BC$')
# for i in range(len(ID)):
#     plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(fwhm_data1[i])+0.3), color='#969696')
# plt.errorbar(incl_sx, np.abs(fwhm1_sx), yerr=err_fwhm1_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
# for i in range(len(ID_sx)):
#     plt.annotate(ID_sx[i], (incl_sx[i]+0.3, np.abs(fwhm1_sx[i])+0.3), color='k')
# plt.errorbar(incl_down, np.abs(fwhm1_down), yerr=err_fwhm1_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
# for i in range(len(ID_down)):
#     plt.annotate(ID_down[i], (incl_down[i]+0.3, np.abs(fwhm1_down[i])+0.3), color='k')
# # plt.errorbar(incl_data2, np.abs(fwhm_data2), yerr=err_fwhm2, color='k', linestyle='None', marker='d', capsize=3, label='$Rigliaco\,et\,al.\,(2013)$')
# # for i in range(len(ID2)):
# #     plt.annotate(ID2[i], (incl_data2[i]+0.3, np.abs(fwhm_data2[i])+0.3))
# plt.xlabel(r'$i \, [^{\circ}]$')
# plt.ylabel(r'$FWHM$')
# plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
# plt.title('$[OI] \, 6300 \AA$')
# # plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg')
# plt.axis([-1.0, 91.0, 0.0, 300.0])
# plt.tight_layout()
# plt.legend(loc='best')
# plt.savefig('./observables/'+str(species)+'/fwhm_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
# plt.savefig('./observables/'+str(species)+'/eps/fwhm_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
# plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#fecc5c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fecc5c', label='$b='+str(b[0])+'$')
plt.plot(incl_deg, fwhm2, color='#fd8d3c', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#fd8d3c', label='$b='+str(b[1])+'$')
plt.plot(incl_deg, fwhm3, color='#f03b20', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#f03b20', label='$b='+str(b[2])+'$')
# plt.plot(incl_deg, fwhm4, color='#bd0026', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#bd0026', label='$b='+str(b[3])+'$')
#plt.plot(incl_hydro, fwhm_hydro, color='k', linestyle='dotted', label='$Alexander \, (2008)$')
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(fwhm_data[i])+0.5), color='k')
plt.errorbar(incl_sx, np.abs(fwhm_sx), yerr=err_fwhm_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_sx)):
    plt.annotate(ID_sx[i], (incl_sx[i]-10.0, np.abs(fwhm_sx[i])+0.5), color='k')
plt.errorbar(incl_down, np.abs(fwhm_down), yerr=err_fwhm_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_down)):
    plt.annotate(ID_down[i], (incl_down[i]+0.5, np.abs(fwhm_down[i])-1.2), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$FWHM \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
plt.axis([-1., 91., 5., 51.])
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
    plt.annotate(ID[i], (incl_data[i]+0.8, np.abs(vcentr_data[i])-0.1), color='k')
plt.errorbar(incl_sx, np.abs(vcentr_sx), yerr=err_vcentr_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_sx)):
    plt.annotate(ID_sx[i], (incl_sx[i]-10., np.abs(vcentr_sx[i])-0.3), color='k')
plt.errorbar(incl_down, np.abs(vcentr_down), yerr=err_vcentr_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_down)):
    plt.annotate(ID_down[i], (incl_down[i]+0.5, np.abs(vcentr_down[i])-0.1), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
plt.axis([-1., 91., -0.5, 14.])
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_b_cs'+str(cs)+'_R'+str(R)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

## ---------------------- PLOT THE [OI] FWHM FOR DIFFERENT RESOLUTIONS ---------------------- ##

R = [30.e3, 21.e3, 15.e3, 10.e3]
species = 'OI'
mdot = 'mdot10e-9'

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
v_peak3 = []
v_centr3 = []
fwhm3 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[2])+'.txt', 'r') as f3:
        lines = f3.readlines()[10:]
        v_peak3.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr3.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm3.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f3.close()
v_peak4 = []
v_centr4 = []
fwhm4 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[3])+'.txt', 'r') as f4:
        lines = f4.readlines()[10:]
        v_peak4.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr4.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm4.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f4.close()

## WE NEED TO PLOT THE FWHM, NOT HALF
fwhm1 = np.array(fwhm1)*2.
fwhm2 = np.array(fwhm2)*2.
fwhm3 = np.array(fwhm3)*2.
fwhm4 = np.array(fwhm4)*2.

plt.figure()
plt.plot(incl_deg, fwhm1, color='#edf8b1', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#edf8b1', label='$R = '+str(R[0])+'$')
plt.plot(incl_deg, fwhm2, color='#a1dab4', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#a1dab4', label='$R = '+str(R[1])+'$')
plt.plot(incl_deg, fwhm3, color='#41b6c4', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#41b6c4', label='$R = '+str(R[2])+'$')
plt.plot(incl_deg, fwhm4, color='#225ea8', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#225ea8', label='$R = '+str(R[3])+'$')
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(fwhm_data[i])+0.5), color='k')
plt.errorbar(incl_sx, np.abs(fwhm_sx), yerr=err_fwhm_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_sx)):
    plt.annotate(ID_sx[i], (incl_sx[i]-10.0, np.abs(fwhm_sx[i])+0.5), color='k')
plt.errorbar(incl_down, np.abs(fwhm_down), yerr=err_fwhm_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_down)):
    plt.annotate(ID_down[i], (incl_down[i]+0.5, np.abs(fwhm_down[i])-1.2), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$FWHM \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
plt.axis([-1., 91., 5., 51.])
plt.tight_layout()
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/fwhm_resolution_b'+str(b[2])+'_cs'+str(cs)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/fwhm_resolution_b'+str(b[2])+'_cs'+str(cs)+'_'+str(mdot)+'.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

## ---------------- COMPARE WITH DATA OF [OI] LINE FOR DIFFERENT cs ---------------------- ##

b = 1.00
incl_deg = [0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
r_in = 0.1
r_out = 9.5
cs = [3, 5, 10]
species = 'OI'
mdot='mdot10e-9'

path_file = []
for j in range(len(cs)):
    path_file.append('../cs'+str(cs[j])+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out))

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg)):
    with open(str(path_file[0])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[0])+'.txt', 'r') as f1:
        lines = f1.readlines()[10:]
        v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f1.close()
v_peak2 = []
v_centr2 = []
fwhm2 = []
for i in range(len(incl_deg)):
    with open(str(path_file[1])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[0])+'.txt', 'r') as f2:
        lines = f2.readlines()[10:]
        v_peak2.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr2.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm2.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f2.close()
v_peak3 = []
v_centr3 = []
fwhm3 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[0])+'.txt', 'r') as f3:
        lines = f3.readlines()[10:]
        v_peak3.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr3.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm3.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f3.close()

## WE NEED TO PLOT THE FWHM, NOT HALF
fwhm1 = np.array(fwhm1)*2.
fwhm2 = np.array(fwhm2)*2.
fwhm3 = np.array(fwhm3)*2.

plt.figure()
plt.plot(incl_deg, np.abs(v_centr1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(v_centr3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data, np.abs(vcentr_data), yerr=err_vcentr, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.8, np.abs(vcentr_data[i])-0.1), color='k')
plt.errorbar(incl_sx, np.abs(vcentr_sx), yerr=err_vcentr_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_sx)):
    plt.annotate(ID_sx[i], (incl_sx[i]-10., np.abs(vcentr_sx[i])-0.3), color='k')
plt.errorbar(incl_down, np.abs(vcentr_down), yerr=err_vcentr_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_down)):
    plt.annotate(ID_down[i], (incl_down[i]+0.5, np.abs(vcentr_down[i])-0.1), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$- v_{centroid} \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
plt.axis([-1., 91., -0.5, 14.])
plt.tight_layout()
plt.legend(loc='upper right')
plt.savefig('./observables/'+str(species)+'/vcentr_soundspeed_b'+str(b)+'_R'+str(R[0])+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/vcentr_soundspeed_b'+str(b)+'_R'+str(R[0])+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, np.abs(fwhm1), color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, np.abs(fwhm2), color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, np.abs(fwhm3), color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(fwhm_data[i])+0.5), color='k')
plt.errorbar(incl_sx, np.abs(fwhm_sx), yerr=err_fwhm_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_sx)):
    plt.annotate(ID_sx[i], (incl_sx[i]-10.0, np.abs(fwhm_sx[i])+0.5), color='k')
plt.errorbar(incl_down, np.abs(fwhm_down), yerr=err_fwhm_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_down)):
    plt.annotate(ID_down[i], (incl_down[i]+0.5, np.abs(fwhm_down[i])-1.2), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$FWHM \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
plt.axis([-1., 91., 9., 32.])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0., 1.1), loc='upper left')
plt.savefig('./observables/'+str(species)+'/fwhm_soundspeed_b'+str(b)+'_R'+str(R[0])+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/fwhm_soundspeed_b'+str(b)+'_R'+str(R[0])+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()

## ---------------- COMPARE WITH DATA OF [OI] LINE FOR DIFFERENT cs and R ---------------------- ##

v_peak4 = []
v_centr4 = []
fwhm4 = []
for i in range(len(incl_deg)):
    with open(str(path_file[0])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[1])+'.txt', 'r') as f4:
        lines = f4.readlines()[10:]
        v_peak4.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr4.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm4.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f4.close()
v_peak5 = []
v_centr5 = []
fwhm5 = []
for i in range(len(incl_deg)):
    with open(str(path_file[1])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[1])+'.txt', 'r') as f5:
        lines = f5.readlines()[10:]
        v_peak5.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr5.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm5.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f5.close()
v_peak6 = []
v_centr6 = []
fwhm6 = []
for i in range(len(incl_deg)):
    with open(str(path_file[2])+'/incl_'+str(round(incl_deg[i],2))+'/observables_R'+str(R[1])+'.txt', 'r') as f6:
        lines = f6.readlines()[10:]
        v_peak6.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr6.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm6.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f6.close()

## WE NEED TO PLOT THE FWHM, NOT HALF
fwhm4 = np.array(fwhm4)*2.
fwhm5 = np.array(fwhm5)*2.
fwhm6 = np.array(fwhm6)*2.

plt.figure()
plt.plot(incl_deg, fwhm1, color='#c6dbef', linestyle='--', linewidth=2.5, marker='None', markeredgecolor='#c6dbef')
plt.plot(incl_deg, fwhm2, color='#2171b5', linestyle='--', linewidth=2.5, marker='None', markeredgecolor='#2171b5')
plt.plot(incl_deg, fwhm3, color='#08306b', linestyle='--', linewidth=2.5, marker='None', markeredgecolor='#08306b')
plt.plot(incl_deg, fwhm4, color='#c6dbef', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#c6dbef', label='$c_{s} = 3 \, km/s$')
plt.plot(incl_deg, fwhm5, color='#2171b5', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#2171b5', label='$c_{s} = 5 \, km/s$')
plt.plot(incl_deg, fwhm6, color='#08306b', linestyle='-', linewidth=2.5, marker='None', markeredgecolor='#08306b', label='$c_{s} = 10 \, km/s$')
plt.plot([], [], color='k', linestyle='-', linewidth=2.5, label='$R = '+str(R[1])+'$')
plt.plot([], [], color='k', linestyle='--', linewidth=2.5, label='$R = '+str(R[0])+'$')
plt.errorbar(incl_data, np.abs(fwhm_data), yerr=err_fwhm, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3, label='$NC$')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.5, np.abs(fwhm_data[i])+0.5), color='k')
plt.errorbar(incl_sx, np.abs(fwhm_sx), yerr=err_fwhm_sx, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_sx)):
    plt.annotate(ID_sx[i], (incl_sx[i]-10.0, np.abs(fwhm_sx[i])+0.5), color='k')
plt.errorbar(incl_down, np.abs(fwhm_down), yerr=err_fwhm_down, color='k', markeredgecolor='None', linestyle='None', marker='o', capsize=3)
for i in range(len(ID_down)):
    plt.annotate(ID_down[i], (incl_down[i]+0.5, np.abs(fwhm_down[i])-1.2), color='k')
plt.xlabel(r'$i \, [^{\circ}]$')
plt.ylabel(r'$FWHM \, [km/s]$')
plt.xticks(np.arange(min(incl_deg), max(incl_deg)+10., 10.0))
plt.title('$[OI] \, 6300 \AA$')
plt.axis([-1., 91., 9., 39.])
plt.tight_layout()
plt.legend(bbox_to_anchor=(0., 1.), loc='upper left', ncol=2)
plt.savefig('./observables/'+str(species)+'/fwhm_soundspeed_res_b'+str(b)+'_'+str(mdot)+'_data.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/eps/fwhm_soundspeed_res_b'+str(b)+'_'+str(mdot)+'_data.eps', format='eps', dpi=300, bbox_inches='tight')
plt.show()
