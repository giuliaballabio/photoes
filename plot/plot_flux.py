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
mdot = 'nonorm'

path_file = []
for j in range(len(b)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[j], 2)))+'_r'+str(r_in)+'_r'+str(r_out))


## ---------------------- PLOT THE FLUX AS A FUNCTION OF THE OUTER RADIUS ----------------------

value = []
for i in range(len(r_out)):
    with open('../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b[0], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'/incl_'+str(round(incl_deg,2))+'/photoes_b'+str('{:.2f}'.format(round(b[0], 2)))+'_r'+str(r_in)+'_r'+str(r_out[i])+'_i'+str(round(incl_deg,2))+'.txt', 'r') as f:
        lines = f.readlines()[24:]
        # lines = f.readlines()[33:]
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
        lines = f.readlines()[24:]
        # lines = f.readlines()[33:]
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
        lines = f.readlines()[24:]
        # lines = f.readlines()[33:]
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
plt.plot(r_out, flux1, color='#addd8e', linestyle='None', marker='o', markeredgecolor='#addd8e', label='$b=0.75$')
plt.plot(r_out, flux2, color='#31a354', linestyle='None', marker='o', markeredgecolor='#31a354', label='$b=1.00$')
plt.plot(r_out, flux3, color='#006837', linestyle='None', marker='o', markeredgecolor='#006837', label='$b=1.50$')
plt.xlabel(r'$R_{out}$', fontsize=15)
plt.ylabel(r'$L_{NeII}$', fontsize=15)
plt.xticks(np.arange(min(r_out), max(r_out)+0.5, 0.5))
# plt.title('R$_{in}$ = '+str(r_in)+' au - i = '+str(round(incl_deg,2)))
plt.axis([4.5, 10.0, 0.e-6, 6.e-6])
plt.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
plt.legend(loc='best')
plt.savefig('./observables/'+str(species)+'/flux_cs'+str(cs)+'_'+str(mdot)+'.png', format='png', dpi=300, bbox_inches='tight')
plt.savefig('./observables/'+str(species)+'/flux_cs'+str(cs)+'_'+str(mdot)+'.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.show()
