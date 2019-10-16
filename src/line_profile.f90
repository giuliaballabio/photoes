
! ----------------------------------------------------------------------------------------------
!
!          Giulia Ballabio,gb258@leicester.ac.uk
!          Created on March 2018
!
!           Compute forbidden lines for photoevaporative winds.
!                Input file -> semianalythical solutions.
!
!        N.B. run the code with these flags:
!        > gfortran -Wunused-variable -Wextra -ffpe-trap=invalid,zero,overflow -pedantic
!           -finit-real=snan -fdefault-real-8 -fdefault-double-8 -fbounds-check -g
!           -fopenmp -o output line_profile.f90
!        or alternatively with:
!        > ifort -g -check all -fpe0 -warn -traceback -debug extended -qopenmp -o output
!           line_profile.f90
!
! ----------------------------------------------------------------------------------------------


program lineprofile

use omp_lib

implicit none

integer                              :: i,j,k,l,npoints,l_in,l_out,l_25
integer,parameter                    :: n_r=1324,n_theta0=250,n_theta=2*300,n_phi=4*300,n_v=1600,n=5d7
real,dimension(1:n_r)                :: r,r_in,r_out,dr,centre_r
!real,dimension(1:n_r-1)             :: dr
real,dimension(1:n)                  :: r_stream,theta_stream,x_stream,y_stream
real,dimension(1:n)                  :: rho_stream,v_r_stream,v_theta_stream
real,dimension(1:n_theta)            :: theta,sinth,costh,centre_theta
! real,dimension(1:n_theta)            :: dA,dmass,dmass_ng
real,dimension(1:n_r)                :: dA,dmass,dmass_ng
real,dimension(1:n_phi)              :: phi,sinphi,cosphi
real                                 :: ratio_r,dtheta,dphi,r_inner,r_outer,b,b_input,ub,theta_max
real                                 :: incl_deg,incl_rad,sinincl,cosincl,tot_flux,Mdot,Mdot_25
!real                                 :: t_in,t_fin
real,dimension(1:n_r,1:n_theta0)     :: rho2d,v_r2d,v_theta2d,v_phi2d,Rb
real,dimension(1:n_r,1:n_theta)      :: rho,rho_ng,n_e,v_r,v_theta,v_phi
real,dimension(1:n_r,1:n_theta)      :: dV,C,cell_flux,v_los
real,dimension(1:n_v)                :: v,line_flux
real                                 :: Rg,ng,rhog,vth,vel_convert,nu,A_hnu,constants,Temp
real                                 :: v_los_r,v_los_th,v_los_phi
logical,save                         :: init=.false.
character(len=6)                     :: str_i,species_flag
real                                 :: m_atom,Ab,A_ul,lambda,X_ion,n_cr,T_ul

!! PHYSICAL CONSTANTS IN CGS !!
real,parameter                       :: au=1.496d13,year=31536000.0,G=6.672d-8
real,parameter                       :: km=6.6846d-9,s=3.171d-8,eV=1.60218d-12
real,parameter                       :: Msun=1.989d33,Lsun=3.826d33,Mstar=1.*Msun,MJ=1.898d30
real,parameter                       :: pi=3.14159,m_h=1.6726d-24,mu=1.,m_e=9.1094d-28
real,parameter                       :: cs=10.0d5
real,parameter                       :: h_planck=6.6261d-27,speed_light=2.9979d10
real,parameter                       :: CC=0.14,Phi_star=1.0d41,alphab=2.60d-13,T0=1.d4,k_b=1.38d-16

! call cpu_time(t_in)

species_flag='NeII'
call which_species(species_flag,m_atom,Ab,A_ul,T_ul,n_cr,X_ion,lambda)

!! PHYSICS SCALING FACTORS !!
Rg=(G*Mstar)/(cs**2)                                                      ! [cm]
ng=CC*sqrt((3.0*Phi_star)/(4.0*pi*alphab*(Rg*Rg*Rg)))                     ! [cm^-3]
rhog=mu*m_h*ng                                                            ! [g/cm^3]
vth=cs/sqrt(m_atom)                                                       ! [cm/s]
!! CONVERSION: [km/s] in [au/years]
vel_convert=km/s

write(*,*) '-----------------------------------------------------------'
write(*,*) ' '
write(*,*) '    Compute forbidden lines for photoevaporative winds'
write(*,*) ' '
write(*,*) '-----------------------------------------------------------'
write(*,*) 'Scaling factors:'
write(*,*) 'Rg =',Rg,'cm =',Rg/au,'au'
write(*,*) 'ng =',ng,'cm^-3'
write(*,*) 'rhog =',rhog,'g/cm^3'
write(*,*) 'ng/rhog =',ng/rhog

!! READ GRID FILE AND CREATE A GRID AT THE BOUNDARY OF THE CELL !!
!! The values of radii are in units of Rg
write(*,*) 'Creating the 2D grid from the hydro simulations...'
open(unit=112,file='../../../../../src/grid_r_30.dat')
do i=1,n_r
    read(112,*) r(i)
enddo
close(112)
ratio_r=sqrt(r(2)/r(1))
do i=1,n_r
    r_in(i)=r(i)/ratio_r
    r_out(i)=r(i)*ratio_r
    dr(i)=r_out(i)-r_in(i)
    !write(*,*) i,r_in(i),r(i),r_out(i),dr(i),ratio_r
enddo
dtheta=pi/dble(n_theta)
do j=1,n_theta
    theta(j)=(j-1)*dtheta
enddo
dphi=2.*pi/dble(n_phi)
do k=1,n_phi
    phi(k)=(k-1)*dphi
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

!! CALCULATE THE NUMBER OF POINTS OF THE STREAMLINE !!
write(*,*) 'Counting the number of points of the streamline...'
open(unit=123,file='./streamline_polarcoord.txt')
npoints=0
do
    read(123,*,end=100)
    npoints=npoints+1
enddo
100 close(123)

!! GET r AND theta FOR THE FIRST STREAMLINE !!
open(unit=134,file='./streamline_polarcoord.txt')
do i=1,npoints
    read(134,*) r_stream(i),theta_stream(i)
enddo
close(134)
!! GET x AND y FOR THE FIRST STREAMLINE !!
open(unit=145,file='./streamline_cartcoord.txt')
do i=1,npoints
    read(145,*) x_stream(i),y_stream(i)
enddo
close(145)

!! DEFINING THE WIND LAUNCHING REGION !!
write(*,*) 'Setting the wind launching region...'
r_inner=0.1
r_outer=30.0
!! FIND THE INDEX THAT CORRESPONDS TO THE INNER AND OUTER RADII !!
l=1
do while (r(l).le.r_inner)
    l=l+1
enddo
l_in=l
l=1
do while (r(l).le.r_outer)
    l=l+1
enddo
l_out=l

!! CALCULATING THE CENTRE OF EACH CELL !!
!! These quantites are in units of Rg !!
write(*,*) 'Calculating the centre of the cell...'
do i=1,n_r
    centre_r(i)=r(i)+(dr(i)/2.)
enddo
do j=1,n_theta-1
    centre_theta(j)=theta(j)+(dtheta/2.)
enddo

!! GET THE DATA FROM THE SCLAE-FREE STREAMLINE !!
write(*,*) 'Reading data from files...'
open(unit=156,file='./rhov_fields.txt',status='old')
do i=1,npoints
    read(156,*) rho_stream(i),v_r_stream(i),v_theta_stream(i)
enddo
close(156)

!! DEFINE A THETA MAX FOR THE FORBIDDEN REGION !!
theta_max=atan(y_stream(npoints)/x_stream(npoints))
! j=1
! do while (centre_theta(j).le.theta_max)
!     j=j+1
! enddo
! j_max=j

!! NORMALISATION FACTOR FOR THE DENSITY DETERMINED AT THE FLOW BASE !!
!! rho(R=Rg) = rhog !!
!! N.B. THE CONVERSION IN PHYSICAL UNITS IS DONE LATER !!
write(*,*) 'Normalizing the streamlines...'
b_input=1.50
ub=0.56
b=b_input

rho2d(:,:)=0.0 !1.5e-17
v_r2d(:,:)=0.0
v_theta2d(:,:)=0.0
v_phi2d(:,:)=0.0 !0001
Rb(:,:)=0.0
do i=l_in,l_out
    do j=1,n_theta0
        k=1
        if (centre_theta(j)<theta_max) then
            do while (theta_stream(k)<centre_theta(j))
                k=k+1
            enddo
            Rb(i,j)=centre_r(i)/r_stream(k)
            rho2d(i,j)=(rho_stream(k))*((Rb(i,j))**(-b)) !*((Rb(i,j)/(Rg/au))**(-b))
            v_r2d(i,j)=ub*v_r_stream(k)
            v_theta2d(i,j)=ub*v_theta_stream(k)
            !! V_PHI: KEPLERIAN ANGULAR MOMENTUM CONSERVATION
            !v_phi2d(i,j)=((r_stream(k)/(Rg/au))*Rb(i,j))**(-0.5) !*(Mstar/Msun)**0.5
            !! V_PHI: KEPLERIAN VELOCITY
            v_phi2d(i,j)=((x_stream(k)/(Rg/au))*Rb(i,j))**(-0.5) !*(Mstar/Msun)**0.5
        elseif (centre_theta(j)>theta_max) then
            rho2d(i,j)=0.d0 !1.5e-15
            v_r2d(i,j)=0.d0 !5.e-1
            v_theta2d(i,j)=0.d0 !0.5
            v_phi2d(i,j)=0.d0 !0.05
        endif
    enddo
enddo

!! SET THE DENSITY AT THE BASE EQUAL TO ZERO FOR r > 10Rg !!
do i=l_in,l_out
    if (centre_r(i)>10.) then
        rho2d(i,250)=0.d0
    endif
enddo

if(.not.init) then
    open(unit=178,file='./rho_grid.txt')
else
    open(unit=178,file='./rho_grid.txt',status='old',position='append')
endif
if(.not.init) then
    open(unit=189,file='./v_grid.txt')
else
    open(unit=189,file='./v_grid.txt',status='old',position='append')
endif
do i=1,n_r
    do j=1,n_theta0
        write(178,'(1(es18.10,1X))') rho2d(i,j)
        write(189,'(3(es18.10,1X))') v_phi2d(i,j),v_r2d(i,j),v_theta2d(i,j)
    enddo
enddo
close(178)
close(189)

!! REVERSE ALONG THETA AXIS !!
!! THE STREAMLINES START FROM THE MIDPLANE !!
write(*,*) 'Building the 3D field...'
!! WIND STARTING FROM THE MIDPLANE
do i=1,n_r
    do j=1,n_theta0
        !! DISC ZONE FOR z>0 UP TO 250 !!
        rho(i,50+j)=rho2d(i,n_theta0+1-j)
        v_r(i,50+j)=v_r2d(i,n_theta0+1-j)
        v_theta(i,50+j)=v_theta2d(i,n_theta0+1-j)
        v_phi(i,50+j)=v_phi2d(i,n_theta0+1-j)
        !! DISC ZONE FOR z<0 DOWN TO 250 !!
        rho(i,n_theta/2+j)=rho2d(i,j) !(i,n_theta0+1-j)
        v_r(i,n_theta/2+j)=v_r2d(i,j) !(i,n_theta0+1-j)
        v_theta(i,n_theta/2+j)=v_theta2d(i,j) !(i,n_theta0+1-j)
        v_phi(i,n_theta/2+j)=v_phi2d(i,j) !(i,n_theta0+1-j)
    enddo
    !! REGIONS NEAR THE Z AXIS !!
    rho(i,1:50)=0.d0
    v_r(i,1:50)=0.d0
    v_theta(i,1:50)=0.d0
    v_phi(i,1:50)=0.d0
    rho(i,n_theta/2+n_theta0:n_theta)=0.d0
    v_r(i,n_theta/2+n_theta0:n_theta)=0.d0
    v_theta(i,n_theta/2+n_theta0:n_theta)=0.d0
    v_phi(i,n_theta/2+n_theta0:n_theta)=0.d0
enddo
!! WIND STARTING AT A FIXED THETA, AS THE HYDRO SIMULATIONS
! do i=1,n_r
!     do j=1,n_theta0
!         !! DISC ZONE FOR z>0 !!
!         rho(i,j)=rho2d(i,n_theta0+1-j)
!         v_r(i,j)=v_r2d(i,n_theta0+1-j)
!         v_theta(i,j)=v_theta2d(i,n_theta0+1-j)
!         v_phi(i,j)=v_phi2d(i,n_theta0+1-j)
!         !! DISC ZONE FOR z<0 !!
!         rho(i,n_theta/2+50+j)=rho2d(i,j)
!         v_r(i,n_theta/2+50+j)=v_r2d(i,j)
!         v_theta(i,n_theta/2+50+j)=v_theta2d(i,j)
!         v_phi(i,n_theta/2+50+j)=v_phi2d(i,j)
!     enddo
!     !! DISC ZONE ON THE MIDPLANE !!
!     rho(i,n_theta0+1:n_theta/2+50)=0.d0
!     v_r(i,n_theta0+1:n_theta/2+50)=0.d0
!     v_theta(i,n_theta0+1:n_theta/2+50)=0.d0
!     v_phi(i,n_theta0+1:n_theta/2+50)=0.d0
! enddo

if(.not.init) then
    open(unit=198,file='./rho.txt')
else
    open(unit=198,file='./rho.txt',status='old',position='append')
endif
if(.not.init) then
    open(unit=200,file='./velocities.txt')
else
    open(unit=200,file='./velocities.txt',status='old',position='append')
endif
do i=1,n_r
    do j=1,n_theta
        write(198,'(1(es18.10,1X))') rho(i,j)
        write(200,'(3(es18.10,1X))') v_phi(i,j), v_r(i,j), v_theta(i,j)
    enddo
enddo
close(198)
close(200)

!! CONVERT TO PHYSICAL UNITS !!
write(*,*) 'Converting to physical units...'
!! CONVERSION: code units -> cm -> au !!
r(:)=r(:)*Rg/au
dr(:)=dr(:)*Rg/au
!! CONVERSION: code units -> g/cm**3 -> Msun/au**3 !!
rho(:,:)=rho(:,:)*rhog/(Msun/(au**3.))
!! CONVERSION: code units -> cm/s -> km/s !!
v_r(:,:)=v_r(:,:)*cs*1.e-5 !(cs/(2.0*pi))*1.e-5
v_theta(:,:)=v_theta(:,:)*cs*1.e-5 !(cs/(2.0*pi))*1.e-5
v_phi(:,:)=v_phi(:,:)*cs*1.e-5 !(cs/(2.0*pi))*1.e-5
!! CONVERSION: cm/s -> km/s !!
vth=vth*1.e-5

!! WRITE THE DATA INTO A FILE TO PLOT THE BOUNDARY CONDITION !!
if(.not.init) then
    open(unit=210,file='./bound_cond.txt')
else
    open(unit=210,file='./bound_cond.txt',status='old',position='append')
endif
do i=1,n_r
    write(210,'(4(es18.10,1X))') r(i),rho(i,250),v_theta(i,250),v_phi(i,250)
    !write(11,'(4(es18.10,1X))') r(i),rho(i,249),v_theta(i,249),v_phi(i,250)
enddo
close(210)

!! COMPUTE THE MASS FLUX FROM THE ANALYTHICAL MODEL !!
! write(*,*) 'Calculating the mass flux at the outer radius...'
! do j=1,n_theta
!     dA(j)=r_out(l_out)*r_out(l_out)*sinth(j)*dtheta
!     dmass(j)=rho(l_out,j)*v_r(l_out,j)*dA(j)
! enddo
! Mdot=sum(dmass) !*n_phi
! write(*,*) '-----------------------------------------------------------'
! write(*,*) '   Total mass flux =',Mdot,'Msun/yr'
! write(*,*) '-----------------------------------------------------------'

write(*,*) 'Calculating the mass flux at the base of the wind...'
do i=1,n_r
    dA(i)=r(i)*dr(i)
    dmass(i)=rho(i,250)*v_theta(i,250)*dA(i)
enddo
Mdot=sum(dmass)
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Total mass flux =',Mdot,'Msun/yr'
write(*,*) '-----------------------------------------------------------'
l=1
do while (r(l).le.25.)
    l=l+1
enddo
l_25=l
Mdot_25=0.0
dmass=0.0
do i=1,l_25
    dA(i)=r(i)*dr(i)
    dmass(i)=rho(i,250)*v_theta(i,250)*dA(i)
enddo
Mdot_25=sum(dmass)
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Mass flux (r<25 au) =',Mdot_25,'Msun/yr'
write(*,*) '-----------------------------------------------------------'

!! NORMALIZE THE DENSITY SUCH AS Mdot(r<25au) = 10^-10 Msun/yr
rho_ng(:,:)=rho(:,:)/ng
rho(:,:)=rho(:,:)*(1.d-10/Mdot_25)
Mdot_25=0.0
dmass=0.0
do i=1,l_25
    dA(i)=r(i)*dr(i)
    dmass(i)=rho(i,250)*v_theta(i,250)*dA(i)
enddo
Mdot_25=sum(dmass)
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Normalised mass flux (r<25 au) =',Mdot_25,'Msun/yr'
write(*,*) '-----------------------------------------------------------'
!! Define the density without the correction for ng
Mdot=0.0
dmass=0.0
do i=1,n_r
    dA(i)=r(i)*dr(i)
    dmass(i)=rho(i,250)*v_theta(i,250)*dA(i)
    dmass_ng(i)=rho_ng(i,250)*v_theta(i,250)*dA(i)
enddo
Mdot=sum(dmass)
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Normalised total mass flux =',Mdot,'Msun/yr'
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Normalised ng =',Mdot/sum(dmass_ng)
write(*,*) '-----------------------------------------------------------'

!! COMPUTE THE FLUX FOR A SINGLE CELL !!
!! N.B. THE CONSTANTS ARE ALL IN CGS: CONVERT QUANTITIES IN CGS !!
write(*,*) 'Calculating the flux for a single cell...'
Temp=T0*(cs/10.0d5)**2.
nu=speed_light/lambda
A_hnu=A_ul*h_planck*nu
constants=Ab*A_hnu*X_ion
do i=1,n_r
    do j=1,n_theta
        !! CONVERSION: volume [au**3] -> [cm**3] !!
        dV(i,j)=r(i)*r(i)*sinth(j)*dr(i)*dtheta*dphi*(au**3)
        !! CONVERSION: mass density [Msun/au**3] -> numerical density [particles/cm**3] !!
        n_e(i,j)=rho(i,j)*(ng/rhog)*(Msun/(au**3))
        if(n_e(i,j) > 0.0) then
            C(i,j)=1.d0+(n_cr/n_e(i,j))
            cell_flux(i,j)=constants/((2.d0*C(i,j)*exp(T_ul/Temp))+1.d0)*n_e(i,j)*dV(i,j)
        else
            cell_flux(i,j)=0.d0
        endif
    enddo
enddo
!! WRITE THE DATA INTO A FILE TO PLOT A MAP OF THE FLUX !!
if(.not.init) then
    open(unit=219,file='./cellflux.txt')
else
    open(unit=219,file='./cellflux.txt',status='old',position='append')
endif
do i=1,n_r
    do j=1,n_theta
        write(219,'(1(es20.10e5,1X))') cell_flux(i,j)
    enddo
enddo
close(219)
tot_flux=sum(cell_flux)*n_phi
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Total flux =',tot_flux/Lsun,'Lsun'
write(*,*) '-----------------------------------------------------------'

!! CREATE VELOCITY ARRAY !!
do l=1,n_v
    v(l)=l*0.1-80.
enddo

!! DEFINE THE INCLINATION ANGLE !!
!write(*,*) 'Please enter the inclination angle of the disc in degrees: '
!read(*,*) incl_deg
!write(*,*) 'Write the value of i in the format for the name of the file: '
!read(*,*) str_i
incl_deg=0.0
str_i='0.0'
incl_rad=incl_deg*(pi/180.)

!! USEFUL VARIABLES TO MAKE THE COMPUTATION FASTER !!
sinincl=sin(incl_rad)
cosincl=cos(incl_rad)

!! COMPUTE THE LINE PROFILE !!
write(*,*) 'Computing the line profile...'
line_flux(:)=0.d0
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
    open(unit=223,file='./line_profile_i'//trim(str_i)//'.txt')
else
    open(unit=223,file='./line_profile_i'//trim(str_i)//'.txt',status='old',position='append')
endif
do l=1,n_v
    write(223,'(2(es20.10e5,1X))') v(l),line_flux(l)
enddo
close(223)

! call cpu_time(t_fin)
!
! write(*,*) t_fin-t_in,'seconds'
! write(*,*) (t_fin-t_in)/60.,'minutes'
! write(*,*) (t_fin-t_in)/3600.,'hours'
write(*,*) '-----------------------------------------------------------'

contains
!! CHOOSE THE SPECIES !!
subroutine which_species(species_flag,m_atom_x,Ab_x,A_ul_x,T_ul_x,n_cr_x,X_ion_x,lambda_x)
    implicit none
    character(len=6),intent(in)          :: species_flag
    real,intent(out)                     :: m_atom_x,Ab_x,A_ul_x,lambda_x,X_ion_x,n_cr_x,T_ul_x

    ! !! [NEII] CONSTANTS !!
    ! real,parameter                       :: m_atom_ne=20.,Ab_ne=1.d-4,A_ul_ne=8.39d-3,lambda_ne=12.81d-4
    ! real,parameter                       :: X_ion_ne=0.75,n_cr_ne=5.0d5,T_ul_ne=1122.8
    !
    ! !! [SII] CONSTANTS !!
    ! real,parameter                       :: m_atom_s=32.,Ab_s=1.45d-5,A_ul_s=1.9d-1,lambda_s=406.98d-7
    ! real,parameter                       :: X_ion_s=1.0,n_cr_s=2.6d6,T_ul_s=35354.
    ! real,parameter                       :: Ipot_s=10.36 !value in eV
    !
    ! !! [OI] CONSTANTS !!
    ! real,parameter                       :: m_atom_o=16.,Ab_o=5.37d-4,A_ul_o=5.6d-3,lambda_o=630.0d-7
    ! real,parameter                       :: X_ion_o=1.0,n_cr_o=1.8d6,T_ul_o=22830.

    if (species_flag=='NeII') then
        m_atom_x=20.
        Ab_x=1.d-4
        A_ul_x=8.39d-3
        lambda_x=12.81d-4
        X_ion_x=0.75
        n_cr_x=5.0d5
        T_ul_x=1122.8
    else if (species_flag=='SIIa') then
        m_atom_x=32.
        Ab_x=1.45d-5
        A_ul_x=1.9d-1
        lambda_x=406.98d-7
        X_ion_x=1.0
        n_cr_x=2.6d6
        T_ul_x=35354.
    else if (species_flag=='SIIb') then
        m_atom_x=32.
        Ab_x=1.45d-5
        A_ul_x=7.7d-2
        lambda_x=407.75d-7
        X_ion_x=1.0
        n_cr_x=1.9e6
        T_ul_x=35430.
    else if (species_flag=='SIIc') then
        m_atom_x=32.
        Ab_x=1.45d-5
        A_ul_x=2.0d-4
        lambda_x=671.83d-7
        X_ion_x=1.0
        n_cr_x=1.7d3
        T_ul_x=21416.
    else
        m_atom_x=16.
        Ab_x=5.37d-4
        A_ul_x=5.6d-3
        lambda_x=630.0d-7
        X_ion_x=0.3
        n_cr_x=1.8d6
        T_ul_x=22830.
    endif

end subroutine which_species

end program
