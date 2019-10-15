
! ---------------------------------------------------------------------------------------
!
!          Giulia Ballabio, gb258@leicester.ac.uk
!          Created on June 2018
!
!        Compute the semi-analythical solution for photoevaporative winds
!            from Clarke and Alexander, 2016
!            This version has been created following Cathie's code ssbasegb.f
!            I also implemented a RK4 method for integration stability
!
!         N.B. run the code with these flags:
!        > gfortran -Wunused-variable -Wextra -ffpe-trap=invalid,zero,overflow
!           -finit-real=snan -fno-range-check -pedantic -fbounds-check -g -o output program.f90
!        or alternatively with:
!        > ifort -g -check all -fpe0 -warn -traceback -debug extended -qopenmp -o output
!           program.f90
!
! ---------------------------------------------------------------------------------------

program selfsimilar

implicit none
integer                                          :: i
integer,parameter                                :: n=10e12
double precision,parameter                       :: cs=10.0d5
double precision,parameter                       :: pi=3.141592  !!,G=6.672d-8,Msun=1.989d33,Mstar=1.*Msun
double precision                                 :: b,ub,reff0,gm,b_input
double precision                                 :: x0,y0,u0,dx,dy,rho,r,phi,theta,flux !,uph
double precision                                 :: x1,y1,u1,tanth1,a1,phi1,theta1
double precision                                 :: x2,y2,u2,ur2,uth2,tanth2,a2,phi2,theta2
double precision                                 :: dudy,ddx1,ddx2,ddxdy
double precision                                 :: rk4,dphi,dphidy,ft

!write (*,*) 'Insert a value for b: '
!read (*,*) b
!write (*,*) 'Now, insert the corrisponding value for ub: '
!read (*,*) ub

!b_array=(0.5,0.75,1.,1.5,2.)
!ub_array=(0.92,0.85,0.77,0.56,0.29)

b_input=0.75
ub=0.85

!! INITIAL CONDITIONS
x0=1.
y0=0.
u0=ub
b=b_input
theta=0.5*pi
phi=datan(y0/x0)

!! I don't know what this parameter does,but it's very important!
!! It should be the gravity, therefore set to zero according to the model
gm=0.0 !1.d-6
reff0=u0*u0*x0/b
dx=1.d-6
dy=1.d-6
x1=x0+dx
y1=(2.*reff0*dx-dx*dx)**0.5
! dy=y1-y0
u1=u0
theta1=datan((reff0-dx)/y1)
phi1=datan(y1/x1)
tanth1=dtan(theta1)
a1=(x1*x1*tanth1-x1*y1)/(1.+tanth1*tanth1)**0.5
do i=1,n
    ft=(1.+tanth1*tanth1)**0.5
    a1=(x1*x1*tanth1-x1*y1)/ft
    ddx1=a1*u0*((u1**(-2.)-1.)*dudy(x1,y1,u1,tanth1,gm,b,u0)+(dphidy(x1,y1,tanth1,gm)/u1))*EXP(0.5*(u1**2.-u0**2.) &
        +dphi(x1,y1,gm))*ft**3./(tanth1**2.*x1*(x1+y1*tanth1))
    ddx2=ft**2.*(x1*tanth1-y1)/(tanth1**3.*x1*(x1+y1*tanth1))
    ddxdy=ddx1+ddx2
    rho=dexp(-0.5*(u1**2.-u0**2.)-dphi(x1,y1,gm))
    flux=rho*u1*a1
    a2=(x1**2.*tanth1-x1*y1)/(1.+tanth1**2.)**0.5
    x2=x1+dy/tanth1+0.5*ddxdy*dy**2.
    y2=y1+dy
    !u2=u1+dudy(x1,y1,u1,tanth1,gm,b,u0)*dy
    u2=u1+rk4(x1,y1,u1,dy,tanth1,gm,b,u0)
    ur2=u1*dcos(theta1-phi1)
    uth2=u1*dsin(theta1-phi1)
    tanth2=(1./tanth1+ddxdy*dy)**(-1.)
    if (tanth2<y2/x2) stop
    phi2=datan(y2/x2)
    theta2=datan(tanth2)
    r=(x1*x1+y1*y1)**0.5
    !!uph=(1/x1)**0.5
    open(unit=1,file='./streamline_cartcoord.txt')
    write(1,'(3(es18.10,1X))') x1,y1
    open(unit=2,file='./streamline_polarcoord.txt')
    write(2,'(2(es18.10,1X))') r,phi1
    !open(unit=3,file='./velocity.txt')
    !write(3,'(1(es18.10,1X))') u1
    open(unit=4,file='./rhov_fields.txt')
    write(4,'(4(es18.10,1X))') rho, ur2, uth2 !!, uph
    x1=x2
    y1=y2
    u1=u2
    tanth1=tanth2
    phi1=phi2
    theta1=theta2
    a1=a2
enddo

close(1)
close(2)
!close(3)
close(4)
end program

!!--------------------------------------------------------------------------
!!   FUNCTIONS
!!--------------------------------------------------------------------------

!! FUNCTION FOR THE EULER METHOD
double precision function dudy(x,y,u,tanth,gm,b,u0)

double precision                                 :: x,y,u,tanth
double precision                                 :: ft,a1,rhs1,rhs2,rhs3,rh,rhs,lhs1,lhs2,lhs
double precision                                 :: b,u0,gm
double precision                                 :: dphi,dphidy,gefft

ft=(1.+tanth*tanth)**0.5
a1=(x*x*tanth-x*y)/ft
rhs1=b*ft/(x*tanth-y)
rhs2=-u**2.*(x*tanth-y)/(ft*x*(x+y*tanth))
rhs3=u*u0*tanth*EXP(0.5*(u**2.-u0**2.)+dphi(x,y,gm))*dphidy(x,y,tanth,gm)/(x*(x+y*tanth))
rh=-tanth*((x+y*tanth)*dphidy(x,y,tanth,gm))/(ft*(x*tanth-y))
!write (*,*) rh,rhs3
rhs=rhs1+rhs2+rhs3+rh+gefft(x,y,tanth,gm)
lhs1=a1*u0*(1.-u**2.)*tanth*EXP(0.5*(u**2.-u0**2.)+dphi(x,y,gm))/(x*(x+y*tanth))
lhs2=u*tanth*(x+y*tanth)/(ft*(x*tanth-y))
lhs=lhs1+lhs2
dudy=rhs/lhs
return
end function dudy

!! FUNCTION FOR THE RUNGE-KUTTA METHOD
double precision function rk4(x,y,u,dy,tanth,gm,b,u0)

double precision                                 :: x,y,u,tanth,dy,dudy
double precision                                 :: gm,b,u0
double precision                                 :: k1,k2,k3,k4

k1=dudy(x,y,u,tanth,gm,b,u0)
!half_step=y1+0.5*dy
k2=dudy(x,y+0.5*dy,u+0.5*dy*k1,tanth,gm,b,u0)
k3=dudy(x,y+0.5*dy,u+0.5*dy*k2,tanth,gm,b,u0)
k4=dudy(x,y+dy,u+dy*k3,tanth,gm,b,u0)
rk4=(k1+(2.*k2)+(2.*k3)+k4)*dy/6.
return
end function rk4

!! These functions are zero in absence of gravity (gm=0)
double precision function dphi(x,y,gm)

double precision                                 :: x,y,r
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
dphi=gm*(-r**(-1.)+0.5*(x)**(-2.)+0.5)
return
end function dphi

double precision function dphidy(x,y,tanth,gm)

double precision                                 :: x,y,tanth,r
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
dphidy=gm*(((x+y*tanth)/(tanth*r**3.))-(1./(tanth*x**3.)))
return
end function dphidy

double precision function gefft(x,y,tanth,gm)

double precision                                 :: x,y,tanth,r,ft
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
ft=(1.+tanth**2.)**0.5
gefft=gm*((tanth/(ft*x**3.))-((x*tanth-y)/(ft*r**3.)))
return
end function gefft
