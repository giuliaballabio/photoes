import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.style.use('classic')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 10
plt.rcParams['xtick.labelsize'] = 8
plt.rcParams['ytick.labelsize'] = 8
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.titlesize'] = 12

speed_light = 299792.458                     #km/s
cs = 10
species = 'OI'
mdot = 'mdot10e-9'

## GET THE DATA FROM THE OUTPUT FILE FROM FORTRAN ##
# incl_deg = 90.0
# b = input("Insert the value of b: ")
# r_in = input("Insert the inner radius: ")
# r_out = input("And the outer radius: ")
incl_deg = 90.0
b_input = 2.00
r_inner = 1.0
r_outer = 9.5

b = b_input
r_in = r_inner
r_out = r_outer

path_file = '../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))

## DATA FROM MODEL
v = -1.*np.array(map(float, [lines.split()[0] for lines in open(str(path_file)+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
line_flux = np.array(map(float, [lines.split()[1] for lines in open(str(path_file)+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
## DATA FROM HYDRO SIMULATIONS
# v = -1.*np.array(map(float, [lines.split()[0] for lines in open('../data_hydro_midplane/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
# line_flux = np.array(map(float, [lines.split()[1] for lines in open('../data_hydro_midplane/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
# v = -1.*np.array(map(float, [lines.split()[0] for lines in open('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
# line_flux = np.array(map(float, [lines.split()[1] for lines in open('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))

## CALCULATE THE VELOCITY AT THE PEAK OF THE LINE
peak_flux = np.amax(line_flux)
i = 0
while (line_flux[i] < peak_flux):
    i += 1
i_max = i
v_peak = v[i_max]
# print 'v_peak = '+str(v_peak)+' km/s'

## CONVOLUTION WITH A GAUSSIAN FUNCTION
def gaussian(x,norm,mean,sigma):
    # norm = 1./np.sqrt(2. * np.pi * sigma**2.)
    return norm*np.exp(-(x-mean)**2/(2.*sigma**2.))

## PROPERTIES OF THE TELESCOPE BEAM ##
# Telescope spectral resolution
# VISIR R=30000
# MIKE R=19000 25000
R = 15000.
delta_v = speed_light / R
sigma_telescope = delta_v / 2.
norm = 1./np.sqrt(2. * np.pi * sigma_telescope**2.)
gauss_telescope = gaussian(v,norm,0.,sigma_telescope)

## CONVOLUTION ##
convolution = np.convolve(line_flux, gauss_telescope, mode="same")
convolution = convolution / np.amax(convolution)

plt.figure()
plt.plot(v, convolution, color='b', label='convolution')
plt.plot(v, line_flux / np.amax(line_flux), color='r', label='line')
plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
plt.ylabel(r'Normalized L(v)', fontsize=15)
plt.legend(loc='best')
plt.savefig(str(path_file)+'/convolution_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro_midplane/incl_'+str(round(incl_deg, 2))+'/convolution_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/convolution_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.show()

## FIT THE CONVOLUTION WITH A GAUSSIAN
## FIND THE VELOCITY AT THE PEAK OF THE CONVOLUTION
peak_conv = np.amax(convolution)
i = 0
while (convolution[i] < peak_conv):
    i += 1
i_max_conv = i
v_peak_conv = v[i_max_conv]

## FIND THE FWHM OF THE CONVOLUTION
# st_dev = np.std(v)
# sigma_conv = 2.*np.sqrt(2.*np.log(2.)) * st_dev
def FWHM(X,Y):
    half_max = max(Y) / 2.
    #find when function crosses line half_max (when sign of diff flips)
    #take the 'derivative' of signum(half_max - Y[])
    d = np.sign(half_max - np.array(Y[0:-1])) - np.sign(half_max - np.array(Y[1:]))
    #plot(X[0:len(d)],d) #if you are interested
    #find the left and right most indexes
    left_idx = np.where(d > 0)[0]
    right_idx = np.where(d < 0)[-1]
    return np.abs(X[right_idx] - X[left_idx]) #return the difference (full width)

sigma_conv = FWHM(v,convolution)/2.

## FIND THE PARAMETERS OF THE GAUSSIAN FIT
popt,pcov = curve_fit(gaussian,v,convolution,p0=[1.,v_peak_conv,sigma_conv[0]])

plt.figure()
plt.plot(v, line_flux / np.amax(line_flux), color='r', label='Model')
plt.plot(v, convolution,'b',label='R = '+str(R))
plt.plot(v, gaussian(v,*popt),'k--',label='Gaussian fit')
plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
plt.ylabel(r'Normalized L(v)', fontsize=15)
plt.axis([-40., 40., 0., 1.2])
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' Rg - R$_{out}$ = '+str(r_out)+' Rg - i = '+str(incl_deg))
plt.legend(loc='best')
plt.savefig(str(path_file)+'/gaussian_fit_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro_midplane/incl_'+str(round(incl_deg, 2))+'/gaussian_fit_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/gaussian_fit_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.show()

## CALCULATE THE CUMULATIVE FUNCTION
## PYTHON MODULE TO CALCULATE THE CUMULATIVE FUNCTION
L_cum = np.cumsum(line_flux)
## MY IMPLEMENTATION OF THE CUMULATIVE FUNCTION
# def cumulative(list):
#     total = 0
#     for x in line_flux:
#         total += x
#         yield total
# L_cum = list(cumulative(line_flux))
## NORMALISE THE CUMULATIVE FUNCTION
L_cum = L_cum / L_cum[len(L_cum)-1]

## CALCULATE THE CENTROID VELOCITY
L_centr = 0.5
i = 0
while (L_cum[i] < L_centr):
    i = i+1
i_centr = i
v_centr = v[i_centr]

plt.figure()
plt.plot(-v, L_cum, color='r')
plt.vlines(v_centr, 0., 2.5, 'k', linestyle='--')
plt.hlines(L_centr, -50., 50., 'k', linestyle='--')
plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
plt.ylabel(r'Cumulative flux', fontsize=15)
plt.axis([-40., 40., 0., 1.2])
plt.savefig(str(path_file)+'/cumulative.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro_midplane/incl_'+str(round(incl_deg, 2))+'/cumulative.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/cumulative.png', format='png', bbox_inches='tight')
# plt.show()

## WRITE THE OBSERVABLES INTO A FILE
f = open(str(path_file)+'/observables_R'+str(R)+'.txt', 'w+')
# f = open('../data_hydro_midplane/incl_'+str(round(incl_deg, 2))+'/observables_R'+str(R)+'.txt', 'w+')
# f = open('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/observables_R'+str(R)+'.txt', 'w+')
f.write('b = '+str(b)+' - R_in = '+str(r_in)+' Rg - R_out = '+str(r_out)+' Rg - i = '+str(incl_deg)+'\n')
f.write('\n')
f.write('------------------------------------------------------------------------------- \n')
f.write('PROPERTIES OF THE CONVOLUTION \n')
f.write('Velocity at peak NOT convolved [km/s] \t Velocity at peak [km/s] \t FWHM \n')
f.write(str(v_peak)+'\t\t\t'+str(v_peak_conv)+'\t\t\t'+str(sigma_conv[0])+'\n')
f.write('\n')
f.write('------------------------------------------------------------------------------- \n')
f.write('PROPERTIES OF THE GAUSSIAN FIT \n')
f.write('Velocity at peak [km/s] \t Centroid velocity [km/s] \t FWHM \n')
f.write(str(popt[1])+'\t\t\t'+str(v_centr)+'\t\t\t'+str(popt[2]))
f.close()
