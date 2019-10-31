import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.cm as cm
from scipy.optimize import curve_fit
from scipy import signal
import scipy.integrate as integrate
from physics_constant import *

plt.style.use('classic')
plt.rc('font', family='serif')
plt.rc('xtick', labelsize='medium')
plt.rc('ytick', labelsize='medium')
plt.rc('axes', titlesize='xx-large')
plt.rc('axes', labelsize='xx-large')
plt.rc('legend', fontsize='large')

cs = '10'
species = 'NeII'
mdot = 'mdot10e-9'

## GET THE DATA FROM THE OUTPUT FILE FROM FORTRAN ##
# incl_deg = 90.0
# b = input("Insert the value of b: ")
# r_in = input("Insert the inner radius: ")
# r_out = input("And the outer radius: ")
incl_deg = 0.0
b_input = 0.75
r_inner = 0.1
r_outer = 9.5

b = b_input
r_in = r_inner
r_out = r_outer

path_file = '../cs'+str(cs)+'kms/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))
# path_file = '../cs'+str(cs)+'/'+str(species)+'/'+str(mdot)+'/data_b'+str('{:.2f}'.format(round(b, 2)))+'_r'+str(r_in)+'_r'+str(r_out)+'/incl_'+str(round(incl_deg, 2))

## DATA FROM MODEL
v = -1.*np.array(map(float, [lines.split()[0] for lines in open(str(path_file)+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
line_flux = np.array(map(float, [lines.split()[1] for lines in open(str(path_file)+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
## DATA FROM HYDRO SIMULATIONS
# v = -1.*np.array(map(float, [lines.split()[0] for lines in open('../data_hydro_midplane/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
# line_flux = np.array(map(float, [lines.split()[1] for lines in open('../data_hydro_midplane/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/line_profile_i'+str(round(incl_deg, 2))+'.txt', 'r')]))
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
    return norm*np.exp(-(x-mean)**2./(2.*sigma**2.))

## PROPERTIES OF THE TELESCOPE BEAM ##
# Telescope spectral resolution
# VISIR R=30000
# MIKE R=19000 25000
R = 30000.
speed_light = 299792.458                     #km/s
delta_v = speed_light/R
sigma_telescope = delta_v/2.
norm = (2.*np.pi*sigma_telescope*sigma_telescope)**(-0.5)
gauss_telescope = gaussian(v,norm,0.,sigma_telescope)

## NORMALIZATION TO THE INTEGRAL OF THE LINE
lineflux_norm = np.abs(integrate.simps(line_flux, v))
line_flux = line_flux / lineflux_norm

# CONVOLUTION
convolution = signal.convolve(line_flux, gauss_telescope, mode='same') / sum(gauss_telescope)

## TEST THE CONVOLUTION PYTHON LIBRARY
## THE INTEGRAL OF THE LINE MUST BE THE SAME
# print integrate.simps(convolution,v)
# print integrate.simps(line_flux,v)

plt.figure()
plt.plot(v, convolution, color='b', label='convolution')
plt.plot(v, line_flux, color='r', label='line')
plt.xlabel(r'v [$\frac{km}{s}$]')
plt.ylabel(r'Normalized L(v)')
plt.tight_layout()
plt.axis([-40., 40., 0., np.max(line_flux)+(np.max(line_flux)*0.25)])
plt.legend(loc='best')
plt.savefig(str(path_file)+'/convolution_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro_midplane/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/convolution_R'+str(R)+'.png', format='png', bbox_inches='tight')
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
## p0=[1,mean,sigma]
## popt is the optimal values for the parameters so that the sum of the squared residuals is minimized
## pcov is the estimated covariance of popt; the diagonals provide the variance of the parameter estimate
popt, pcov = curve_fit(gaussian,v,convolution,p0=[1.,v_peak_conv,sigma_conv[0]])

plt.figure()
plt.plot(v, line_flux, color='r', label='Model')
plt.plot(v, convolution,'b',label='R = '+str(R))
plt.plot(v, gaussian(v,*popt),'k--',label='Gaussian fit')
plt.xlabel(r'v [$\frac{km}{s}$]')
plt.ylabel(r'Normalized L(v)')
plt.tight_layout()
plt.axis([-40., 40., 0., np.max(line_flux)+(np.max(line_flux)*0.25)])
plt.legend(loc='best')
plt.savefig(str(path_file)+'/gaussian_fit_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro_midplane/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/gaussian_fit_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/gaussian_fit_R'+str(R)+'.png', format='png', bbox_inches='tight')
# plt.show()

## CALCULATE THE ERROR ON THE PARAMETERS FROM THE FUNCTION
## To compute one standard deviation of the parameters use:
perr = np.sqrt(np.diag(pcov))

## CALCULATE THE CHI-SQUARE TEST FOR THE FIT
def gaussfunc(x,mean,sigma):
    return 1./(sigma*np.sqrt(2.*np.pi))*np.exp(-(x-mean)**2./(2.*sigma**2.))

mean = np.arange(-10., 10., 0.1)
sigma = np.arange(1., 20., 0.1)
mean_grid, sigma_grid = np.meshgrid(mean, sigma, indexing='ij')

def chi2reduced(x,y):
    chisq = []
    for j in range(len(x)):
        for k in range(len(y)):
            gauss = []
            add = []
            for i in range(len(v)):
                gauss.append(gaussfunc(v[i],x[j],y[k]))
                add.append((gauss[i]-convolution[i])**2./(0.07*np.max(convolution))**2.)
            # Reduced chi-squared
            chisq.append(np.sum(add)/(len(convolution)-2))
    chisq = np.array(chisq)
    chisq = np.reshape(chisq,(len(x),len(y)))
    return chisq

chisq = chi2reduced(mean,sigma)
# chisq = ((mean_grid[:]-np.abs(popt[1]))**2./np.abs(popt[1])) + ((sigma_grid[:]-np.abs(popt[2]))**2./np.abs(popt[2]))
chisq = np.where(chisq <= 15.,chisq,15.)
# chisq = np.where(chisq <= 9.,chisq,9.)
# mask_white = np.ma.array(chisq, mask=chisq<=15.)
mask_1chisq = np.ma.array(chisq, mask=chisq>=1.5)

plt.figure()
# plt.plot(mean_grid, sigma_grid, marker='.', color='k', linestyle='none')
plt.contourf(mean_grid, 2.*sigma_grid, chisq, cmap='viridis')
plt.plot([min(mean),max(mean)], [2.*popt[2],2.*popt[2]], 'k')
plt.plot([popt[1],popt[1]], [2.*min(sigma),2.*max(sigma)], 'k')
# plt.pcolormesh(mean_grid, 2.*sigma_grid, chisq, cmap='viridis', vmin=0., vmax=7.)
cbar = plt.colorbar()
# plt.contourf(mean_grid, sigma_grid, mask_white, colors='white')
CS = plt.contour(mean_grid, 2.*sigma_grid, mask_1chisq, levels=[1.0], colors='k')
plt.xlabel(r'<$v_{peak}$>')
plt.ylabel(r'<$FWHM$>')
plt.tight_layout()
# h,_ = CS.legend_elements()
# plt.legend([h[0]], ['$\chi^2 = 1.0$'])
# plt.axis([0.,3.,5.,20.])
plt.savefig(str(path_file)+'/chisquare_jk.png', format='png', bbox_inches='tight')
# plt.show()

## FIND THE ERROR BARS ON v_peak AND FWHM
k = 0
while(sigma[k] <= popt[2]):
    k += 1
k_bestsigma = k
idx_mean = np.argwhere(np.diff(np.sign(chisq[:,k_bestsigma] - 1.))).flatten()
plt.figure()
plt.plot(mean, chisq[:,k_bestsigma], 'r')
plt.plot([min(mean),max(mean)], [1.,1.], 'k')
plt.plot([popt[1],popt[1]], [min(chisq[:,k_bestsigma]),max(chisq[:,k_bestsigma])], 'k')
plt.plot(mean[idx_mean], chisq[idx_mean,k_bestsigma], 'bo')
plt.xlabel(r'$v_{peak}$')
plt.ylabel(r'$Reduced\,\chi^2$')
plt.savefig(str(path_file)+'/chisquare_vpeak.png', format='png', bbox_inches='tight')
# plt.show()

err_mean_inf = np.abs(popt[1]-mean[idx_mean[0]])
err_mean_sup = np.abs(popt[1]-mean[idx_mean[1]])

j = 0
while(mean[j] <= popt[1]):
    j += 1
j_bestmean = j
idx_sigma = np.argwhere(np.diff(np.sign(chisq[j_bestmean,:] - 1.))).flatten()
plt.figure()
plt.plot(2.*sigma, chisq[j_bestmean,:], 'r')
plt.plot([2.*min(sigma),2.*max(sigma)], [1.,1.], 'k')
plt.plot([2.*popt[2],2.*popt[2]], [min(chisq[j_bestmean,:]),max(chisq[j_bestmean,:])], 'k')
plt.plot(2.*sigma[idx_sigma], chisq[j_bestmean,idx_sigma], 'bo')
plt.xlabel(r'$FWHM$')
plt.ylabel(r'$Reduced\,\chi^2$')
plt.savefig(str(path_file)+'/chisquare_width.png', format='png', bbox_inches='tight')
# plt.show()

err_sigma_inf = 2.*np.abs(popt[2]-sigma[idx_sigma[0]])
err_sigma_sup = 2.*np.abs(popt[2]-sigma[idx_sigma[1]])

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
plt.xlabel(r'v [$\frac{km}{s}$]')
plt.ylabel(r'Cumulative flux')
plt.tight_layout()
plt.axis([-40., 40., 0., 1.2])
plt.savefig(str(path_file)+'/cumulative.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro_midplane/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/cumulative.png', format='png', bbox_inches='tight')
# plt.savefig('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/cumulative.png', format='png', bbox_inches='tight')
# plt.show()

## WRITE THE OBSERVABLES INTO A FILE
f = open(str(path_file)+'/observables_R'+str(R)+'.txt', 'w+')
# f = open('../data_hydro_midplane/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/observables_R'+str(R)+'.txt', 'w+')
# f = open('../data_hydro/'+str(species)+'/incl_'+str(round(incl_deg, 2))+'/observables_R'+str(R)+'.txt', 'w+')
f.write('b = '+str(b)+' - R_in = '+str(r_in)+' Rg - R_out = '+str(r_out)+' Rg - i = '+str(incl_deg)+'\n')
f.write('\n')
f.write('------------------------------------------------------------------------------- \n')
f.write('PROPERTIES OF THE CONVOLUTION \n')
f.write('Velocity at peak NOT convolved [km/s] \t Velocity at peak [km/s] \t Half FWHM \n')
f.write(str(v_peak)+'\t\t\t'+str(v_peak_conv)+'\t\t\t'+str(sigma_conv[0])+'\n')
f.write('\n')
f.write('------------------------------------------------------------------------------- \n')
f.write('PROPERTIES OF THE GAUSSIAN FIT \n')
f.write('Velocity at peak [km/s] \t Centroid velocity [km/s] \t Half FWHM \t Error on the mean \t Error on the FWHM \n')
f.write(str(popt[1])+'\t\t\t'+str(v_centr)+'\t\t\t'+str(popt[2])+'\t\t\t'+str(perr[1])+'\t\t\t'+str(perr[2])+'\n')
f.write('\n')
f.write('------------------------------------------------------------------------------- \n')
f.write('PROPERTIES OF THE CHI SQUARED \n')
f.write('Err_mean_inf \t Err_mean_sup \t Err_sigma_inf \t Err_sigma_sup \n')
f.write(str(err_mean_inf)+'\t\t\t'+str(err_mean_sup)+'\t\t\t'+str(err_sigma_inf)+'\t\t\t'+str(err_sigma_sup))
f.close()
