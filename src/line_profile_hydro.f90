! ---------------------------------------------------------------------------------------
!
!          Giulia Ballabio,gb258@leicester.ac.uk
!          Created on March 2018
!
!           Compute forbidden lines for photoevaporative winds.
!                Input file -> hydro simulations.
!
!        N.B. run the code with these flags:
!        gfortran -Wunused-variable -Wextra -ffpe-trap=invalid,zero,overflow
!           -finit-real=snan -pedantic -fbounds-check -g -fopenmp -o output program.f90
!
! ---------------------------------------------------------------------------------------


program lineprofile

use omp_lib

implicit none

integer                                          :: i,j,k,l
integer,parameter                                :: n_r=1113,n_theta0=300,n_theta=2*300,n_phi=4*300,n_v=800
double precision,dimension(1:n_r)                :: r,r_in,r_out,dr
!double precision,dimension(1:n_r-1)             :: dr
double precision,dimension(1:n_theta)            :: theta,sinth,costh
double precision,dimension(1:n_theta)            :: dA,dmass
double precision,dimension(1:n_phi)              :: phi,sinphi,cosphi
double precision                                 :: ratio_r,dtheta,dphi
double precision                                 :: incl_deg,incl_rad,sinincl,cosincl,tot_flux,Mdot
double precision,dimension(1:n_r,1:n_theta0)     :: rho2d,v_r2d,v_theta2d,v_phi2d
double precision,dimension(1:n_r,1:n_theta)      :: rho,n_e,v_r,v_theta,v_phi
double precision,dimension(1:n_r,1:n_theta)      :: dV,C,cell_flux,v_los
double precision,dimension(1:n_v)                :: v,line_flux
double precision                                 :: Rg,ng,rhog,vth,vel_convert,nu,A_hnu,constants
double precision                                 :: v_los_r,v_los_th,v_los_phi
logical,save                                     :: init=.false.
character(len=5)                                 :: str_i

!! PHYSICAL CONSTANTS !!
double precision,parameter                       :: au=1.496d13,year=31536000.0,G=6.672d-8
double precision,parameter                       :: km=6.6846d-9,s=3.171d-8
double precision,parameter                       :: Msun=1.989d33,Lsun=3.826d33,Mstar=1.*Msun,MJ=1.898d30
double precision,parameter                       :: pi=3.14159,cs=1.0d6,m_h=1.6726d-24,mu=1.
double precision,parameter                       :: h_planck=6.6261d-27,speed_light=2.9979d10
double precision,parameter                       :: CC=0.14,Phi_star=0.75d41,alphab=2.60d-13,T=1.d4

!! NEII CONSTANTS !!
double precision,parameter                       :: m_atom=20.,Ab_ne=1.d-4,A_ul=8.39d-3,lambda_ne=12.81d-4
double precision,parameter                       :: X_II=0.75,n_cr=5.0d5,T_ul=1122.8

!! PHYSICS SCALING FACTORS !!
Rg=(G*Mstar)/(cs**2)                                                      ! [cm]
ng=CC*sqrt((3.0*Phi_star)/(4.0*pi*alphab*(Rg*Rg*Rg)))                   ! [cm^-3]
rhog=mu*m_h*ng                                                          ! [g/cm^3]
vth=cs/sqrt(m_atom)                                                         ! [cm/s]
!! CONVERSION: [km/s] in [au/years]
vel_convert=km/s

print *,'-----------------------------------------------------------'
print *,' '
print *,'    Compute forbidden lines for photoevaporative winds'
print *,' '
print *,'-----------------------------------------------------------'
print *,'Scaling factors:'
print *,'Rg=',Rg,'cm -',Rg/au,'au'
print *,'ng=',ng,'cm^-3'
print *,'rhog=',rhog,'g/cm^3'
print *,'ng/rhog=',ng/rhog

!! READ GRID FILE AND CREATE A GRID AT THE BOUNDARY OF THE CELL !!
print *,'Creating the 2D grid...'
open(unit=100,file='../grid_r.dat')
do i=1,n_r
		read(100,*) r(i)
enddo
close(100)
ratio_r=sqrt(r(2)/r(1))
do i=1,n_r
		r_in(i)=r(i)/ratio_r
		r_out(i)=r(i)*ratio_r
		dr(i)=r_out(i)-r_in(i)
enddo
open(unit=110,file='../grid_th.dat')
do j=1,n_theta0
		read(110,*) theta(j)
enddo
close(110)
! dtheta=pi/dble(n_theta)
! do j=1,n_theta
! 		theta(j)=(j+0.5)*dtheta
! enddo
dphi=2.*pi/dble(n_phi)
do k=1,n_phi
		phi(k)=(k+0.5)*dphi
enddo

!! UNCOMMENT THESE LINES TO DEFINE AN ANGULAR GRID AT THE CENTRE OF THE CELL !!
!! REMEMBER TO CHANGE THE DEFINITION OF dr !!
!do i=1,n_r-1
!    dr(i)=r(i+1)-r(i)
!enddo
!dtheta=pi/n_theta
!do j=1,n_theta
!    theta(j)=j*dtheta
!enddo
!dphi=2.*pi/n_phi
!do k=1,n_phi
!    phi(k)=k*dphi
!enddo

!! USEFUL VARIABLES TO MAKE THE COMPUTATION FASTER !!
do j=1,n_theta
		sinth(j)=sin(theta(j))
		costh(j)=cos(theta(j))
enddo
do k=1,n_phi
		sinphi(k)=sin(phi(k))
		cosphi(k)=cos(phi(k))
enddo

!! GET THE DATA FROM THE HYDRO-SIMULATIONS AND CREATE 2D ARRAYS !!
print *,'Reading data from files...'
open(unit=200,file='../rho_mean.dat')
open(unit=210,file='../v_r_mean.dat')
open(unit=220,file='../v_th_mean.dat')
open(unit=230,file='../v_phi_mean.dat')
do i=1,n_r
		do j=1,n_theta0
				read(200,*) rho2d(i,j)
				read(210,*) v_r2d(i,j)
				read(220,*) v_theta2d(i,j)
				read(230,*) v_phi2d(i,j)
		enddo
enddo
close(200)
close(210)
close(220)
close(230)

!! REVERSE ALONG THETA AXIS !!
print *,'Building the 3D field...'
do i=1,n_r
		do j=1,n_theta0
				!! DISC ZONE FOR z>0 !!
				rho(i,j)=rho2d(i,n_theta0+1-j)
				v_r(i,j)=v_r2d(i,n_theta0+1-j)
				v_theta(i,j)=v_theta2d(i,n_theta0+1-j)
				v_phi(i,j)=v_phi2d(i,n_theta0+1-j)
				!! DISC ZONE FOR z<0 !!
				! rho(i,n_theta/2+50+j)=rho2d(i,j)
				! v_r(i,n_theta/2+50+j)=v_r2d(i,j)
				! v_theta(i,n_theta/2+50+j)=v_theta2d(i,j)
				! v_phi(i,n_theta/2+50+j)=v_phi2d(i,j)
				rho(i,n_theta0+j)=rho2d(i,j)
				v_r(i,n_theta0+j)=v_r2d(i,j)
				v_theta(i,n_theta0+j)=v_theta2d(i,j)
				v_phi(i,n_theta0+j)=v_phi2d(i,j)
		enddo
		!! DISC ZONE ON THE MIDPLANE !!
		! rho(i,n_theta0+1:n_theta/2+50)=0.d0 !rho2d(i,n_theta0)
		! v_r(i,n_theta0+1:n_theta/2+50)=0.d0 !v_r2d(i,n_theta0)
		! v_theta(i,n_theta0+1:n_theta/2+50)=0.d0 !v_theta2d(i,n_theta0)
		! v_phi(i,n_theta0+1:n_theta/2+50)=0.d0 !v_phi2d(i,n_theta0)
enddo

!! CONVERT TO PHYSICAL UNITS !!
print *,'Converting to physical units...'
!! CONVERSION: code units -> cm -> au !!
r(:)=r(:)*Rg/au
dr(:)=dr(:)*Rg/au
!! CONVERSION: code units -> g/cm**3 -> Msun/au**3 !!
rho(:,:)=rho(:,:)*rhog/(Msun/(au**3))
!! CONVERSION: code units -> cm/s -> km/s !!
v_r(:,:)=v_r(:,:)*(cs/(2.0*pi))*1.e-5 !*(year/au)/vel_convert
v_theta(:,:)=v_theta(:,:)*(cs/(2.0*pi))*1.e-5 !*(year/au)/vel_convert
v_phi(:,:)=v_phi(:,:)*(cs/(2.0*pi))*1.e-5 !*(year/au)/vel_convert

!! CONVERSION: cm/s -> km/s !!
vth=vth*1.e-5

!! WRITE THE DATA INTO A FILE TO PLOT THE BOUNDARY CONDITION !!
if(.not.init) then
		open(unit=300,file='bound_cond.txt')
else
		open(unit=300,file='bound_cond.txt',status='old',position='append')
endif
do i=1,n_r
		write(300,'(4(es18.10,1X))') r(i),rho(i,299),v_theta(i,299),v_phi(i,299)
enddo
close(300)

!! COMPUTE THE MASS FLUX !!
print *,'Calculating the mass flux...'
do j=1,n_theta
		dA(j)=2.0*r(1113)*r(1113)*sinth(j)*dtheta
		dmass(j)=rho(1113,j)*v_r(1113,j)*dA(j)
enddo
Mdot=sum(dmass)
print *,'-----------------------------------------------------------'
print *,'   Total mass flux=',Mdot,'Msun/yr'
print *,'-----------------------------------------------------------'

!! COMPUTE THE FLUX FOR A SINGLE CELL !!
!! N.B. THE CONSTANTS ARE ALL IN CGS: CONVERT QUANTITIES IN CGS !!
print *,'Calculating the flux for a single cell...'
nu=speed_light/lambda_ne
A_hnu=A_ul*h_planck*nu
constants=Ab_ne*A_hnu*X_II
do i=1,n_r
		do j=1,n_theta
				!! CONVERSION: volume [au**3] -> [cm**3] !!
				dV(i,j)=r(i)*r(i)*sinth(j)*dr(i)*dtheta*dphi*(au**3)
				!! CONVERSION: mass density [Msun/au**3] -> numerical density [particles/cm**3] !!
				n_e(i,j)=rho(i,j)*(ng/rhog)*(Msun/(au**3))
				if(n_e(i,j) > 0.0) then
						C(i,j)=1.d0+(n_cr/n_e(i,j))
						cell_flux(i,j)=constants/((2.d0*C(i,j)*exp(-1.0*T_ul/T))+1.d0)*n_e(i,j)*dV(i,j)
				else
						cell_flux(i,j)=0.d0
				endif
		enddo
enddo
!! WRITE THE DATA INTO A FILE TO PLOT A MAP OF THE FLUX !!
if(.not.init) then
    open(unit=400,file='./cellflux.txt')
else
    open(unit=400,file='./cellflux.txt',status='old',position='append')
endif
do i=1,n_r
    do j=1,n_theta
        write(400,'(1(es20.10e5,1X))') cell_flux(i,j)
    enddo
enddo
close(400)

tot_flux=sum(cell_flux)*n_phi
print *,'-----------------------------------------------------------'
print *,'   Total flux=',tot_flux,'erg/s' !=',tot_flux/Lsun,'Lsun'
print *,'-----------------------------------------------------------'

!! CREATE VELOCITY ARRAY !!
do l=1,n_v
v(l)=l*0.1-40.
enddo

!! DEFINE THE INCLINATION ANGLE !!
!print *,'Please enter the inclination angle of the disc in degrees: '
!read(*,*) incl_deg
incl_deg=90.0
str_i='90.0'
incl_rad=incl_deg*(pi/180.)

!! COMPUTE THE LINE PROFILE !!
print *,'Computing the line profile...'
line_flux(:)=0.d0

!! USEFUL VARIABLES TO MAKE THE COMPUTATION FASTER !!
sinincl=sin(incl_rad)
cosincl=cos(incl_rad)

!$OMP PARALLEL &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(i,j,l,v_los_r,v_los_th,v_los_phi,v_los) &
!$OMP REDUCTION(+: line_flux)
!$OMP DO SCHEDULE(runtime)
do k=1,n_phi
		do j=1,n_theta
				do i=1,n_r
						v_los_r=(costh(j)*cosphi(k)*sinincl+sinth(j)*cosincl)*v_r(i,j)
						v_los_th=(-sinth(j)*cosphi(k)*sinincl+costh(j)*cosincl)*v_theta(i,j)
						v_los_phi=(-sinphi(k)*sinincl)*v_phi(i,j)
						v_los(i,j)=v_los_r+v_los_th+v_los_phi
						do l=1,n_v
								line_flux(l)=line_flux(l)+(cell_flux(i,j)/(sqrt(2.d0*pi)*vth)) &
									*exp((-1.d0*(v(l)-v_los(i,j))*(v(l)-v_los(i,j)))/(2.d0*(vth*vth)))
						enddo
				enddo
		enddo
enddo
!$OMP END DO
!$OMP END PARALLEL

!! WRITE THE DATA INTO A FILE TO PLOT THE LINE PROFILE !!
if(.not.init) then
		open(unit=500,file='line_profile_i'//trim(str_i)//'.txt')
else
		open(unit=500,file='line_profile_i'//trim(str_i)//'.txt',status='old',position='append')
endif
do l=1,n_v
		write(500,'(2(es18.10,1X))') v(l),line_flux(l)
enddo
close(500)

print *,'-----------------------------------------------------------'

end program
