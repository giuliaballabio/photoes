import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

speed_light=299792.458                     #km/s

## GET THE DATA FROM THE OUTPUT FILE FROM FORTRAN ##
# incl_deg = 20.0
# b = input("Insert the value of b: ")
# r_in = input("Insert the inner radius: ")
# r_out = input("And the outer radius: ")
incl_deg = 20.0
b_input = 0.75
r_inner = 0.1
r_outer = 5.0

b = b_input
r_in = r_inner
r_out = r_outer

v = np.array(map(float, [lines.split()[0] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
line_flux = np.array(map(float, [lines.split()[1] for lines in open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))

if incl_deg == 90.0:
    v = v
else:
    v = -1. * v

## PROPERTIES OF THE GAUSSIAN CONVOLUTION ##
R = 3.e4 # Telescope resolution
delta_v = speed_light / R
sigma_telescope = delta_v / 2.

## CALCULATE THE VELOCITY AT THE PEAK OF THE LINE
peak_flux = np.amax(line_flux)
i=0
while (line_flux[i] < peak_flux):
    i += 1
i_max = i
v_peak = v[i_max]
# print 'v_peak = '+str(v_peak)+' km/s'

## CONVOLUTION WITH A GAUSSIAN FUNCTION
# convoluted = []
# total = []
# for i in range(len(v)):
#     for j in range(len(v)):
#         convoluted = ((1./np.sqrt(2. * np.pi * sigma**2.))*np.exp(-((v[i]-v[j])/sigma)**2./2.) * line_flux[j])
#     total.append(np.sum(convoluted))
#
# convoluted = np.array(convoluted)
# total = np.array(total)
# print convoluted.shape
# print total.shape
def gaussian(x,norm,mean,sigma):
    # norm = 1./np.sqrt(2. * np.pi * sigma**2.)
    return norm*np.exp(-(x-mean)**2/(2.*sigma**2.))

norm = 1./np.sqrt(2. * np.pi * sigma_telescope**2.)
gauss_telescope = gaussian(v,norm,0.,sigma_telescope)
convolution = np.convolve(line_flux, gauss_telescope, mode="same")
convolution = convolution / np.amax(convolution)

plt.figure()
plt.plot(v, convolution, color='b', label='convolution')
plt.plot(v, line_flux / np.amax(line_flux), color='r', label='line')
plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
plt.ylabel(r'Normalized L(v)', fontsize=15)
plt.legend(loc='best')
plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/convolution.png', format='png', bbox_inches='tight')
# plt.show()

## FIT THE CONVOLUTION WITH A GAUSSIAN
## FIND THE MEAN OF THE CONVOLUTION
peak_conv = np.amax(convolution)
i=0
while (line_flux[i] < peak_flux):
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
# print sigma_conv[0]

## FIND THE PARAMETERS OF THE GAUSSIAN FIT
popt,pcov = curve_fit(gaussian,v,convolution,p0=[1.,v_peak_conv,sigma_conv])

plt.figure()
plt.plot(v, line_flux / np.amax(line_flux), color='r', label='Model')
plt.plot(v, convolution,'b',label='R = 30,000')
plt.plot(v, gaussian(v,*popt),'k--',label='Gaussian fit')
plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
plt.ylabel(r'Normalized L(v)', fontsize=15)
plt.axis([-40., 40., 0., 1.2])
plt.title('b = '+str(b)+' - R$_{in}$ = '+str(r_in)+' au - R$_{out}$ = '+str(r_out)+'au - i = '+str(incl_deg))
plt.legend(loc='best')
plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/gaussian_fit.png', format='png', bbox_inches='tight')
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
# print 'v_centr = '+str(v_centr)

# if incl_deg == 90.0:
#     v = v
# else:
#     v = -1. * v

plt.figure()
plt.plot(v, L_cum, color='r')
plt.vlines(v_centr, 0., 2.5, 'k', linestyle='--')
plt.hlines(L_centr, -50., 50., 'k', linestyle='--')
plt.xlabel(r'v [$\frac{km}{s}$]', fontsize=15)
plt.ylabel(r'Cumulative flux', fontsize=15)
plt.axis([-40., 40., 0., 1.2])
plt.savefig('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/cumulative.png', format='png', bbox_inches='tight')
# plt.show()

## WRITE THE OBSERVABLES INTO A FILE
f = open('../data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))+'/observables.txt', 'w+')
f.write('b = '+str(b)+' - R_in = '+str(r_in)+' au - R_out = '+str(r_out)+'au - i = '+str(incl_deg)+'\n')
f.write('\n')
f.write('------------------------------------------------------------------------------- \n')
f.write('PROPERTIES OF THE CONVOLUTION \n')
f.write('Velocity at peak [km/s] \t Centroid velocity [km/s] \t FWHM \n')
f.write(str(v_peak_conv)+'\t\t\t'+str(v_centr)+'\t\t\t'+str(sigma_conv[0])+'\n')
f.write('\n')
f.write('------------------------------------------------------------------------------- \n')
f.write('PROPERTIES OF THE GAUSSIAN FIT \n')
f.write('Velocity at peak [km/s] \t Centroid velocity [km/s] \t FWHM \n')
f.write(str(popt[1])+'\t\t\t'+str(v_centr)+'\t\t\t'+str(popt[2]))
f.close()
