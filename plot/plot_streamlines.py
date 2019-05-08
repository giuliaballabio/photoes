import numpy as np
import matplotlib.pyplot as plt

plt.style.use('classic')

b = [0.75, 1.00, 1.50, 2.00]
incl_deg = 0.0
r_in = 0.1
r_out = 9.5
cs = 10
species = 'NeII'
path_file = []
for j in range(len(b)):
    path_file.append('../cs'+str(cs)+'kms/'+str(species)+'/data_b'+str('{:.2f}'.format(round(b[j], 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg,2)))

x1 = np.array(map(float, [lines.split()[0] for lines in open(str(path_file[0])+'/streamline_cartcoord.txt', 'r')]))
y1 = np.array(map(float, [lines.split()[1] for lines in open(str(path_file[0])+'/streamline_cartcoord.txt', 'r')]))
x2 = np.array(map(float, [lines.split()[0] for lines in open(str(path_file[1])+'/streamline_cartcoord.txt', 'r')]))
y2 = np.array(map(float, [lines.split()[1] for lines in open(str(path_file[1])+'/streamline_cartcoord.txt', 'r')]))
x3 = np.array(map(float, [lines.split()[0] for lines in open(str(path_file[2])+'/streamline_cartcoord.txt', 'r')]))
y3 = np.array(map(float, [lines.split()[1] for lines in open(str(path_file[2])+'/streamline_cartcoord.txt', 'r')]))
x4 = np.array(map(float, [lines.split()[0] for lines in open(str(path_file[3])+'/streamline_cartcoord.txt', 'r')]))
y4 = np.array(map(float, [lines.split()[1] for lines in open(str(path_file[3])+'/streamline_cartcoord.txt', 'r')]))
xtrue1 = np.array(map(float, [lines.split()[0] for lines in open('../../Cathie_solution/b0.75.dat', 'r')]))
ytrue1 = np.array(map(float, [lines.split()[1] for lines in open('../../Cathie_solution/b0.75.dat', 'r')]))
xtrue2 = np.array(map(float, [lines.split()[0] for lines in open('../../Cathie_solution/b1.dat', 'r')]))
ytrue2 = np.array(map(float, [lines.split()[1] for lines in open('../../Cathie_solution/b1.dat', 'r')]))
xtrue3 = np.array(map(float, [lines.split()[0] for lines in open('../../Cathie_solution/b1.5.dat', 'r')]))
ytrue3 = np.array(map(float, [lines.split()[1] for lines in open('../../Cathie_solution/b1.5.dat', 'r')]))

# plt.figure()
# #plt.plot(xtrue1, ytrue1, color='darkgreen', linewidth=0.5)#, label='b=0.75')
# plt.plot(xtrue2, ytrue2, color='mediumseagreen', linewidth=0.5)#, label='b=1.0')
# #plt.plot(xtrue3, ytrue3, color='limegreen', linewidth=0.5)#, label='b=1.5')
# #plt.plot(x1, y1, color='orangered', linewidth=1, label='b=0.75')
# plt.plot(x2, y2, color='darkorange', linewidth=1, label='b=1.0')
# #plt.plot(x3, y3, color='orange', linewidth=1, label='b=1.5')
# #plt.plot(x4, y4, color='gold', linewidth=1, label='b=2.0')
# #plt.title(r'Streamline topology - b = '+str(b), fontsize=15)
# plt.xlabel(r'$R / R_g$', fontsize=15)
# plt.ylabel(r'$z / R_g$', fontsize=15)
# #plt.axis([0, 5, 0, 5])
# plt.legend()
# plt.savefig('./fortran_version/plots/selfsimilar_solutions_b'+str(b)+'.png', format='png', bbox_inches='tight')
# #plt.show()

plt.figure()
# plt.plot(xtrue1, ytrue1, color='#addd8e', linewidth=1.3, linestyle='--', label='$b=0.75, u_b=0.85$')
# plt.plot(xtrue2+0.5, ytrue2, color='#31a354', linewidth=1.3, linestyle='--', label='$b=1.0, u_b=0.77$')
# plt.plot(xtrue3+1., ytrue3, color='#006837', linewidth=1.3, linestyle='--', label='$b=1.5, u_b=0.56$')
plt.plot(x1, y1, color='#fecc5c', linewidth=1.5, label='$b=0.75, \, u_b=0.85$')
plt.plot(x2+0.5, y2, color='#fd8d3c', linewidth=1.5, label='$b=1.00, \, u_b=0.77$')
plt.plot(x3+1., y3, color='#f03b20', linewidth=1.5, label='$b=1.50, \, u_b=0.56$')
plt.plot(x4+1.5, y4, color='#bd0026', linewidth=1.5, label='$b=2.00, \, u_b=0.29$')
plt.plot(1.06, 0.33, 'kx', linewidth=2.)
plt.plot(1.09+0.5, 0.35, 'kx', linewidth=2.)
plt.plot(1.17+1., 0.30, 'kx', linewidth=2.)
plt.plot(1.23+1.5, 0.16, 'kx', linewidth=2.)
# plt.title(r'Streamline topology',fontsize=15)
plt.xlabel(r'$R / R_g$',fontsize=15)
plt.ylabel(r'$z / R_g$',fontsize=15)
plt.axis([0.,5.,0.,2.0])
plt.legend(loc='upper left', bbox_to_anchor=(0.0, 1.1))
plt.savefig('./streamlines/selfsimilar_solutions_2.png', format='png', dpi=300, bbox_inches='tight')
plt.show()
