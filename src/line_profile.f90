
! ---------------------------------------------------------------------------------------
!
!          Giulia Ballabio,gb258@leicester.ac.uk
!          Created on March 2018
!
!           Compute forbidden lines for photoevaporative winds.
!                Input file -> semianalythical solutions.
!
!        N.B. run the code with these flags:
!        > gfortran -Wunused-variable -Wextra -ffpe-trap=invalid,zero,overflow
!           -finit-real=snan -pedantic -fbounds-check -g -fopenmp -o output program.f90
!        or alternatively with:
!        > ifort -g -check all -fpe0 -warn -traceback -debug extended -qopenmp -o output
!           program.f90
!
! ---------------------------------------------------------------------------------------


program lineprofile

use omp_lib

implicit none

integer                                          :: i,j,k,l,index_i,index_j,npoints,l_in,l_out
integer,parameter                                :: n_r=1113,n_theta0=250,n_theta=2*300,n_phi=4*300,n_v=800,n=1d7
double precision,dimension(1:n_r)                :: r,r_in,r_out,dr
!double precision,dimension(1:n_r-1)             :: dr
double precision,dimension(1:n)                  :: r_stream,theta_stream,x_stream,y_stream
double precision,dimension(1:n)                  :: rho_stream,v_r_stream,v_theta_stream,v_phi_stream
double precision,dimension(1:n_theta)            :: theta,sinth,costh
double precision,dimension(1:n_theta)            :: dA,dmass
double precision,dimension(1:n_phi)              :: phi,sinphi,cosphi
double precision                                 :: ratio_r,dtheta,dphi,r_inner,r_outer,r1,x1,b,b_input
double precision                                 :: incl_deg,incl_rad,sinincl,cosincl,tot_flux,Mdot
double precision                                 :: t_in,t_fin
integer,dimension(1:n_r,1:n_theta0)              :: ncount
double precision,dimension(1:n_r,1:n_theta0)     :: rho2d,v_r2d,v_theta2d,v_phi2d
double precision,dimension(1:n_r,1:n_theta0)     :: sum_rho,sum_vr,sum_vth,sum_vphi
!double precision,dimension(1:n_r,1:n_theta0)     :: rho_mean2d,v_r_mean2d,v_th_mean2d,v_phi_mean2d
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

call cpu_time(t_in)

!! PHYSICS SCALING FACTORS !!
Rg=(G*Mstar)/(cs**2)                                                      ! [cm]
ng=CC*sqrt((3.0*Phi_star)/(4.0*pi*alphab*(Rg*Rg*Rg)))                   ! [cm^-3]
rhog=mu*m_h*ng                                                          ! [g/cm^3]
vth=cs/sqrt(m_atom)                                                         ! [cm/s]
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
write(*,*) 'Creating the 2D grid from the hydro simulations...'
open(unit=112,file='../../data_hydro/grid_r.dat')
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

!! GET THE DATA FROM THE HYDRO-SIMULATIONS AND CREATE 2D ARRAYS !!
!write(*,*) 'Reading data from files...'
!open(unit=2,file='../data_hydro/rho_mean.dat')
!open(unit=3,file='../data_hydro/v_r_mean.dat')
!open(unit=4,file='../data_hydro/v_th_mean.dat')
!open(unit=5,file='../data_hydro/v_phi_mean.dat')
!do i=1,n_r
!    do j=1,n_theta0
!        read(2,*) rho2d(i,j)
!        read(3,*) v_r2d(i,j)
!        read(4,*) v_theta2d(i,j)
!        read(5,*) v_phi2d(i,j)
!    enddo
!enddo
!close(2)
!close(3)
!close(4)
!close(5)

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

!! SHIFT THE STREAMLINE AT THE INNER RADIUS OF THE LAUNCHING REGION !!
write(*,*) 'Setting the wind launching region...'
r_inner=0.03
r_outer=5.
! First value to shift the streamline at zero
x1=x_stream(1)
r1=r_stream(1)
do i=1,npoints
    r_stream(i)=r_stream(i)-r1+r_inner
    x_stream(i)=x_stream(i)-x1+r_inner
enddo

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

!! GET THE DATA FROM THE FIRST STREAMLINE !!
write(*,*) 'Reading data from files...'
open(unit=156,file='./rhov_fields.txt',status='old')
do i=1,npoints
    read(156,*) rho_stream(i),v_r_stream(i),v_theta_stream(i) !!,v_phi_stream(i)
enddo
close(156)

!! CALCULATE THE KEPLERIAN VELOCITY v_phi FOR THE FIRST STREAMLINE !!
do i=1,npoints
    v_phi_stream(i)=x_stream(i)**(-0.5)
enddo

!! MAP THE STREAMLINE INTO THE GRID !!
!! AT EACH STEP THE STREAMLINE IS SHIFTED AT THE CENTRE OF THE NEXT CELL !!
!! THEN WE RUN ALONG THE STREAMLINE AND CHECK AT WHICH CELL EACH POINT BELONGS TO !!
write(*,*) 'Binning the streamline into the grid...'
ncount(:,:)=0
rho2d(:,:)=0.
v_r2d(:,:)=0.
v_theta2d(:,:)=0.
v_phi2d(:,:)=0.
sum_rho(:,:)=0.
sum_vr(:,:)=0.
sum_vth(:,:)=0.
sum_vphi(:,:)=0.
!!$OMP PARALLEL &
!!$OMP DEFAULT(SHARED) &
!!$OMP PRIVATE(i,j,k,l,index_i,index_j,r_stream) &
!!$OMP REDUCTION(+: sum_rho,sum_vr,sum_vth,sum_vphi,ncount)
!!$OMP DO SCHEDULE(runtime)
do l=l_in,l_out
    do k=1,npoints
        r_stream(k)=r_stream(k)-r(l)+r(l+1)
        x_stream(k)=x_stream(k)-r(l)+r(l+1)
        write(*,*) r_stream(k)
        write(*,*) x_stream(k)
        v_phi_stream(k)=x_stream(k)**(-0.5)
        do i=1,n_r-1
            if (r(i).le.r_stream(k).and.r_stream(k).lt.r(i+1))then
                index_i=i
                exit
            endif
        enddo
        do j=1,n_theta0-1
            if (theta(j).le.theta_stream(k).and.theta_stream(k).lt.theta(j+1))then
                index_j=j
                exit
            endif
        enddo
        if (index_i.le.0.or.index_j.le.0) then
            write(*,*) 'Warning! i-index or j-index out of boundary'
            write(*,*) x_stream(k),r(1),r(n_r)
            write(*,*) theta_stream(k),theta(1),theta(n_theta0)
        endif
        sum_rho(index_i,index_j)=sum_rho(index_i,index_j)+rho_stream(k)
        sum_vr(index_i,index_j)=sum_vr(index_i,index_j)+v_r_stream(k)
        sum_vth(index_i,index_j)=sum_vth(index_i,index_j)+v_theta_stream(k)
        sum_vphi(index_i,index_j)=sum_vphi(index_i,index_j)+v_phi_stream(k)
        ncount(index_i,index_j)=ncount(index_i,index_j)+1
    enddo
enddo
!!$OMP END DO
!!$OMP END PARALLEL

write(*,*)
if(.not.init) then
    open(unit=167,file='./coordinates.txt')
else
    open(unit=167,file='./coordinates.txt',status='old',position='append')
endif
do k=1,npoints
    write(167,'(4(es18.10,1X))') x_stream(k),r_stream(k)
enddo
close(167)

!! CALCULATE THE MEAN OF THE DENSITY AND VELOCITY IN EACH CELL !!
write(*,*) 'Calculating the average density and velocity in each cell...'
do i=1,n_r
    do j=1,n_theta0
        if(ncount(i,j)>0)then
            rho2d(i,j)=sum_rho(i,j)/ncount(i,j)
            v_r2d(i,j)=sum_vr(i,j)/ncount(i,j)
            v_theta2d(i,j)=sum_vth(i,j)/ncount(i,j)
            v_phi2d(i,j)=sum_vphi(i,j)/ncount(i,j)
        endif
    enddo
enddo

!! NORMALISATION FACTOR FOR THE DENSITY !!
!! rho(R=Rg) = rhog !!
b_input=0.75
b=b_input
do i=1,n_r
    rho(i,:)=rho(i,:)*((r(i))**(-b))
enddo

if(.not.init) then
    open(unit=178,file='./rho_grid.txt')
else
    open(unit=178,file='./rho_grid.txt',status='old',position='append')
endif
if(.not.init) then
    open(unit=189,file='./v_phi_grid.txt')
else
    open(unit=189,file='./v_phi_grid.txt',status='old',position='append')
endif
do i=1,n_r
    do j=1,n_theta0
        write(178,'(4(es18.10,1X))') rho2d(i,j)
        write(189,'(4(es18.10,1X))') v_phi2d(i,j) !v_r2d(i,j),v_theta2d(i,j),
    enddo
enddo
close(178)
close(189)

!! REVERSE ALONG THETA AXIS !!
write(*,*) 'Building the 3D field...'
do i=1,n_r
    do j=1,n_theta0
        !! DISC ZONE FOR z>0 !!
        rho(i,j)=rho2d(i,n_theta0+1-j)
        v_r(i,j)=v_r2d(i,n_theta0+1-j)
        v_theta(i,j)=-1.*v_theta2d(i,n_theta0+1-j)
        v_phi(i,j)=v_phi2d(i,n_theta0+1-j)
        !! DISC ZONE FOR z<0 !!
        rho(i,n_theta/2+50+j)=rho2d(i,j)
        v_r(i,n_theta/2+50+j)=v_r2d(i,j)
        v_theta(i,n_theta/2+50+j)=v_theta2d(i,j)
        v_phi(i,n_theta/2+50+j)=v_phi2d(i,j)
    enddo
    !! DISC ZONE ON THE MIDPLANE !!
    rho(i,n_theta0+1:n_theta/2+50)=0.d0 !rho2d(i,n_theta0)
    v_r(i,n_theta0+1:n_theta/2+50)=0.d0 !v_r2d(i,n_theta0)
    v_theta(i,n_theta0+1:n_theta/2+50)=0.d0 !v_theta2d(i,n_theta0)
    v_phi(i,n_theta0+1:n_theta/2+50)=0.d0 !v_phi2d(i,n_theta0)
enddo

!! CONVERT TO PHYSICAL UNITS !!
write(*,*) 'Converting to physical units...'
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
    open(unit=210,file='./bound_cond.txt')
else
    open(unit=210,file='./bound_cond.txt',status='old',position='append')
endif
do i=1,n_r
    write(210,'(4(es18.10,1X))') r(i),rho(i,250),v_theta(i,250),v_phi(i,250)
    !write(11,'(4(es18.10,1X))') r(i),rho(i,249),v_theta(i,249),v_phi(i,250)
enddo
close(210)

!! COMPUTE THE MASS FLUX FROM THE HYDRO DATA !!
!write(*,*) 'Calculating the mass flux...'
!do j=1,n_theta
!    dA(j)=2.0*r_out(1113)*r_out(1113)*sinth(j)*dtheta
!    dmass(j)=rho(1113,j)*v_r(1113,j)*dA(j)
!enddo
!! COMPUTE THE MASS FLUX FROM THE ANALYTHICAL MODEL !!
write(*,*) 'Calculating the mass flux...'
do j=1,n_theta
    dA(j)=2.0*r_out(l_out)*r_out(l_out)*sinth(j)*dtheta
    dmass(j)=rho(l_out,j)*v_r(l_out,j)*dA(j)
enddo
Mdot=sum(dmass)
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Total mass flux =',Mdot,'Msun/yr'
write(*,*) '-----------------------------------------------------------'

!! COMPUTE THE FLUX FOR A SINGLE CELL !!
!! N.B. THE CONSTANTS ARE ALL IN CGS: CONVERT QUANTITIES IN CGS !!
write(*,*) 'Calculating the flux for a single cell...'
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
tot_flux=sum(cell_flux)*n_phi
write(*,*) '-----------------------------------------------------------'
write(*,*) '   Total flux =',tot_flux,'Lsun'
write(*,*) '-----------------------------------------------------------'

!! CREATE VELOCITY ARRAY !!
do l=1,n_v
    v(l)=l*0.1-40.
enddo

!! DEFINE THE INCLINATION ANGLE !!
!write(*,*) 'Please enter the inclination angle of the disc in degrees: '
!read(*,*) incl_deg
!write(*,*) 'Write the value of i in the format for the name of the file: '
!read(*,*) str_i
incl_deg=90.0
str_i='90.0'
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
    write(223,'(2(es18.10,1X))') v(l),line_flux(l)
enddo
close(223)

call cpu_time(t_fin)

write(*,*) t_fin-t_in,'seconds'
write(*,*) (t_fin-t_in)/60.,'minutes'
write(*,*) (t_fin-t_in)/3600.,'hours'
write(*,*) '-----------------------------------------------------------'

end program
