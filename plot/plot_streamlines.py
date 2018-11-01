import numpy as np
import matplotlib.pyplot as plt

#plt.style.use('seaborn-darkgrid')
#plt.rcParams["font.family"] = 'Montreal SF'

# CATHIE'S DATA
b = 1.0
x1 = np.array(map(float, [lines.split()[0] for lines in open('../input_from_model/b0.75/streamline.txt', 'r')]))
y1 = np.array(map(float, [lines.split()[1] for lines in open('../input_from_model/b0.75/streamline.txt', 'r')]))
x2 = np.array(map(float, [lines.split()[0] for lines in open('../input_from_model/b1.0/streamline.txt', 'r')]))
y2 = np.array(map(float, [lines.split()[1] for lines in open('../input_from_model/b1.0/streamline.txt', 'r')]))
x3 = np.array(map(float, [lines.split()[0] for lines in open('../input_from_model/b1.5/streamline.txt', 'r')]))
y3 = np.array(map(float, [lines.split()[1] for lines in open('../input_from_model/b1.5/streamline.txt', 'r')]))
x4 = np.array(map(float, [lines.split()[0] for lines in open('../input_from_model/b2.0/streamline.txt', 'r')]))
y4 = np.array(map(float, [lines.split()[1] for lines in open('../input_from_model/b2.0/streamline.txt', 'r')]))
X = np.array(map(float, [lines.split()[0] for lines in open('./Cathie_solutions/output.txt', 'r')]))
Y = np.array(map(float, [lines.split()[1] for lines in open('./Cathie_solutions/output.txt', 'r')]))
xtrue1 = np.array(map(float, [lines.split()[0] for lines in open('./Cathie_solutions/b0.75.dat', 'r')]))
ytrue1 = np.array(map(float, [lines.split()[1] for lines in open('./Cathie_solutions/b0.75.dat', 'r')]))
xtrue2 = np.array(map(float, [lines.split()[0] for lines in open('./Cathie_solutions/b1.dat', 'r')]))
ytrue2 = np.array(map(float, [lines.split()[1] for lines in open('./Cathie_solutions/b1.dat', 'r')]))
xtrue3 = np.array(map(float, [lines.split()[0] for lines in open('./Cathie_solutions/b1.5.dat', 'r')]))
ytrue3 = np.array(map(float, [lines.split()[1] for lines in open('./Cathie_solutions/b1.5.dat', 'r')]))

plt.figure()
plt.plot(xtrue1, ytrue1, color='darkgreen', linewidth=0.5)#, label='b=0.75')
plt.plot(xtrue2, ytrue2, color='mediumseagreen', linewidth=0.5)#, label='b=1.0')
plt.plot(xtrue3, ytrue3, color='limegreen', linewidth=0.5)#, label='b=1.5')
#plt.plot(x1, y1, color='orangered', linewidth=1, label='b=0.75')
plt.plot(x2, y2, color='darkorange', linewidth=1, label='b=1.0')
#plt.plot(x3, y3, color='orange', linewidth=1, label='b=1.5')
#plt.plot(x4, y4, color='gold', linewidth=1, label='b=2.0')
#plt.plot(X, Y, color='navy', linewidth=1, label='Cathie')
#plt.title(r'Streamline topology - b = '+str(b), fontsize=15)
plt.xlabel(r'$R / R_g$', fontsize=15)
plt.ylabel(r'$z / R_g$', fontsize=15)
#plt.axis([0, 5, 0, 5])
plt.legend()
plt.savefig('./fortran_version/plots/selfsimilar_solutions_b'+str(b)+'.png', format='png', bbox_inches='tight')
#plt.show()

plt.figure()
plt.plot(xtrue1, ytrue1, color='darkgreen', linewidth=2, linestyle='--', label='$b=0.75, u_b=0.85$')
plt.plot(xtrue2+0.5, ytrue2, color='mediumseagreen', linewidth=2, linestyle='--', label='$b=1.0, u_b=0.77$')
plt.plot(xtrue3+1., ytrue3, color='limegreen', linewidth=2, linestyle='--', label='$b=1.5, u_b=0.56$')
plt.step(x1, y1, color='orangered', linewidth=2)#, label='b=0.75')
plt.step(x2+0.5, y2, color='darkorange', linewidth=2)#, label='b=1.0')
plt.step(x3+1., y3, color='orange', linewidth=2)#, label='b=1.5')
plt.plot(1.06, 0.33, 'kx', linewidth=2.)
plt.plot(1.09+0.5, 0.35, 'kx', linewidth=2.)
plt.plot(1.17+1., 0.30, 'kx', linewidth=2.)
#plt.title(r'Streamline topology',fontsize=15,**csfont)
plt.xlabel(r'$R / R_g$',fontsize=15)
plt.ylabel(r'$z / R_g$',fontsize=15)
plt.axis([0,5,0,5])
plt.legend()
plt.savefig('./fortran_version/plots/selfsimilar_solutions_real.png', format='png', bbox_inches='tight')
plt.show()

x = x2
y = y2

delta = 0.1
x = x - 0.9
plt.figure()
for i in range(1, int(5/delta)):
    x_shift = x + delta*i
    plt.plot(x_shift, y, color='orange', linewidth=1)
plt.xlabel(r'$R / R_g$', fontsize=15)
plt.ylabel(r'$z / R_g$', fontsize=15)
plt.axis([0, 5, 0, 2.5])
#plt.legend()
plt.savefig('./fortran_version/plots/scalefree_b'+str(b)+'.png', format='png', bbox_inches='tight')
plt.show()


#r = np.sqrt(x4*x4 + y4*y4)
#f = open('grid_r.txt', 'w+')
#f.write("\n".join(str(elem) for elem in x5))
#
#f.close()






