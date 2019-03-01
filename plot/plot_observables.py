import numpy as np
import matplotlib.pyplot as plt
import re


b = [0.75, 1.00, 1.50] #, 2.00]
incl_deg = [0.0, 45.0, 90.0]
r_in = [0.03, 0.1, 1.0]
r_out = [5.0, 9.5]

## ---------------- PLOT THE VELOCITY AT PEAK AND FWHM AS FUNCTIONS OF THE INCLINATION ----------------------

v_peak1 = []
v_centr1 = []
fwhm1 = []
for i in range(len(incl_deg)):
    with open('../data_b'+str('{:.2f}'.format(round(b[0], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[1])+'/incl_'+str(round(incl_deg[i],2))+'/observables.txt', 'r') as f1:
        lines = f1.readlines()[10:]

        v_peak1.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr1.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm1.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f1.close()
v_peak2 = []
v_centr2 = []
fwhm2 = []
for i in range(len(incl_deg)):
    with open('../data_b'+str('{:.2f}'.format(round(b[1], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[1])+'/incl_'+str(round(incl_deg[i],2))+'/observables.txt', 'r') as f2:
        lines = f2.readlines()[10:]

        v_peak2.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr2.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm2.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f2.close()
v_peak3 = []
v_centr3 = []
fwhm3 = []
for i in range(len(incl_deg)):
    with open('../data_b'+str('{:.2f}'.format(round(b[2], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[1])+'/incl_'+str(round(incl_deg[i],2))+'/observables.txt', 'r') as f3:
        lines = f3.readlines()[10:]

        v_peak3.append(map(float, [x.split('\t\t\t')[0] for x in lines]))
        v_centr3.append(map(float, [x.split('\t\t\t')[1] for x in lines]))
        fwhm3.append(map(float, [x.split('\t\t\t')[2] for x in lines]))
f3.close()

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
plt.plot(incl_deg, np.abs(v_peak1), color='#addd8e', linestyle='dashed', marker='o', markeredgecolor='#addd8e', label='b=0.75')
plt.plot(incl_deg, np.abs(v_peak2), color='#31a354', linestyle='dashed', marker='o', markeredgecolor='#31a354', label='b=1.00')
plt.plot(incl_deg, np.abs(v_peak3), color='#006837', linestyle='dashed', marker='o', markeredgecolor='#006837', label='b=1.50')
plt.plot(incl_deg, np.abs(v_centr1), color='#fdcc8a', linestyle='dashed', marker='o', markeredgecolor='#fdcc8a', label='b=0.75 $v_{centr}$')
plt.plot(incl_deg, np.abs(v_centr2), color='#fc8d59', linestyle='dashed', marker='o', markeredgecolor='#fc8d59', label='b=1.00 $v_{centr}$')
plt.plot(incl_deg, np.abs(v_centr3), color='#d7301f', linestyle='dashed', marker='o', markeredgecolor='#d7301f', label='b=1.50 $v_{centr}$')
# plt.plot(incl_data, np.abs(vpeak_data), color='k', linestyle='None', marker='o', markerfacecolor=None)
plt.errorbar(incl_data, np.abs(vpeak_data), yerr=err_vpeak, color='k', linestyle='None', marker='o', label='Sacco et al. (2012)')
for i in range(len(ID)):
    plt.annotate(ID[i], (incl_data[i]+0.3, np.abs(vpeak_data[i])+0.3))
plt.errorbar(incl_data2, np.abs(vpeak_data2), yerr=err_vpeak2, color='k', linestyle='None', marker='o', markerfacecolor='None', label='Pascucci & Sterzik (2009)')
for i in range(len(name)):
    plt.annotate(name[i], (incl_data2[i]+0.3, np.abs(vpeak_data2[i])+0.3))
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'$- v_{peak} \, [km/s]$', fontsize=15)
plt.title('R$_{in}$ = '+str(r_in[1])+' au - R$_{out}$ = '+str(r_out[1]))
plt.axis([-5.0, 95.0, -1.0, 14.0])
plt.legend(loc='upper right', bbox_to_anchor=(1.26, 1.05), fontsize = 'small')
plt.savefig('./observables/vpeak_r'+str(r_in[1])+'_r'+str(r_out[1])+'.png', format='png', bbox_inches='tight')
plt.show()

plt.figure()
plt.plot(incl_deg, fwhm1, color='#b3cde3', linestyle='dashed', marker='o', markeredgecolor='#b3cde3', label='b=0.75')
plt.plot(incl_deg, fwhm2, color='#8c96c6', linestyle='dashed', marker='o', markeredgecolor='#8c96c6', label='b=1.00')
plt.plot(incl_deg, fwhm3, color='#8856a7', linestyle='dashed', marker='o', markeredgecolor='#8856a7', label='b=1.50')
plt.xlabel(r'$i \, [^{\circ}]$', fontsize=15)
plt.ylabel(r'FWHM', fontsize=15)
plt.title('R$_{in}$ = '+str(r_in[1])+' au - R$_{out}$ = '+str(r_out[1]))
plt.axis([-5.0, 95.0, 5.0, 30.0])
plt.legend(loc='best')
plt.savefig('./observables/fwhm_r'+str(r_in[1])+'_r'+str(r_out[1])+'.png', format='png', bbox_inches='tight')
# plt.show()

## ---------------------- PLOT THE FLUX AS A FUNCTION OF THE OUTER RADIUS ----------------------

value = []
for i in range(len(r_out)):
    with open('../data_b'+str('{:.2f}'.format(round(b[0], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[i])+'/incl_'+str(round(incl_deg[2],2))+'/photoes_b'+str('{:.2f}'.format(round(b[0], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[i])+'_i'+str(incl_deg[2])+'.txt', 'r') as f:
        lines = f.readlines()[24:]
        value.append([x.split('\n')[0] for x in lines])
f.close()
str1 = value[0][0]
str2 = value[1][0]
## This array contains all the numbers in the line 'Total Flux etc'
## array[0] = integer part
## array[1] = decimal part
## array[2] = exponent
array1 = map(float, re.findall('\d+', str1))
array2 = map(float, re.findall('\d+', str2))
## We need to move the decimal point to get the right number
d1 = np.floor(np.log10(np.abs(array1[1])))+1
array1[1] = array1[1]*(10.**(-d1))
d2 = np.floor(np.log10(np.abs(array2[1])))+1
array2[1] = array2[1]*(10.**(-d2))
## Now we build the value of the flux
flux1 = [(array1[0]+array1[1])*(10.**(-array1[2])), (array2[0]+array2[1])*(10.**(-array2[2]))]

value = []
for i in range(len(r_out)):
    with open('../data_b'+str('{:.2f}'.format(round(b[1], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[i])+'/incl_'+str(round(incl_deg[2],2))+'/photoes_b'+str('{:.2f}'.format(round(b[1], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[i])+'_i'+str(incl_deg[2])+'.txt', 'r') as f:
        lines = f.readlines()[24:]
        value.append([x.split('\n')[0] for x in lines])
f.close()
str1 = value[0][0]
str2 = value[1][0]
## This array contains all the numbers in the line 'Total Flux etc'
## array[0] = integer part
## array[1] = decimal part
## array[2] = exponent
array1 = map(float, re.findall('\d+', str1))
array2 = map(float, re.findall('\d+', str2))
## We need to move the decimal point to get the right number
d1 = np.floor(np.log10(np.abs(array1[1])))+1
array1[1] = array1[1]*(10.**(-d1))
d2 = np.floor(np.log10(np.abs(array2[1])))+1
array2[1] = array2[1]*(10.**(-d2))
## Now we build the value of the flux
flux2 = [(array1[0]+array1[1])*(10.**(-array1[2])), (array2[0]+array2[1])*(10.**(-array2[2]))]

value = []
for i in range(len(r_out)):
    with open('../data_b'+str('{:.2f}'.format(round(b[2], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[i])+'/incl_'+str(round(incl_deg[2],2))+'/photoes_b'+str('{:.2f}'.format(round(b[2], 2)))+'_r'+str(r_in[1])+'_r'+str(r_out[i])+'_i'+str(incl_deg[2])+'.txt', 'r') as f:
        lines = f.readlines()[24:]
        value.append([x.split('\n')[0] for x in lines])
f.close()
str1 = value[0][0]
str2 = value[1][0]
## This array contains all the numbers in the line 'Total Flux etc'
## array[0] = integer part
## array[1] = decimal part
## array[2] = exponent
array1 = map(float, re.findall('\d+', str1))
array2 = map(float, re.findall('\d+', str2))
## We need to move the decimal point to get the right number
d1 = np.floor(np.log10(np.abs(array1[1])))+1
array1[1] = array1[1]*(10.**(-d1))
d2 = np.floor(np.log10(np.abs(array2[1])))+1
array2[1] = array2[1]*(10.**(-d2))
## Now we build the value of the flux
flux3 = [(array1[0]+array1[1])*(10.**(-array1[2])), (array2[0]+array2[1])*(10.**(-array2[2]))]

plt.figure()
plt.plot(r_out, flux1, color='#41b6c4', linestyle='dashed', marker='o', markeredgecolor='#41b6c4', label='b=0.75')
plt.plot(r_out, flux2, color='#225ea8', linestyle='dashed', marker='o', markeredgecolor='#225ea8', label='b=1.00')
plt.plot(r_out, flux3, color='#081d58', linestyle='dashed', marker='o', markeredgecolor='#081d58', label='b=1.50')
plt.xlabel(r'$R_{out}$', fontsize=15)
plt.ylabel(r'$L_{NeII}$', fontsize=15)
plt.title('R$_{in}$ = '+str(r_in[1])+' au - i = '+str(round(incl_deg[2],2)))
plt.axis([4.0, 10.0, 0.e-6, 6.e-6])
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.legend(loc='best')
plt.savefig('./observables/fluxNeII_r'+str(r_in[1])+'_i'+str(round(incl_deg[2],2))+'.png', format='png', bbox_inches='tight')
plt.show()
