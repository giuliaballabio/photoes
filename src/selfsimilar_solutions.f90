
! ---------------------------------------------------------------------------------------
!
!          Giulia Ballabio, gb258@leicester.ac.uk
!          Created on June 2018
!
!        Compute the semi-analythical solution for photoevaporative winds
!            from Clarke and Alexander, 2016
!            This version has been created following Cathie's code ssbasegb.f
!
!         N.B. run the code with these flags:
!        gfortran -Wunused-variable -Wextra -ffpe-trap=invalid,zero,overflow
!           -finit-real=snan -pedantic -fbounds-check -g -o output program.f90
!
! ---------------------------------------------------------------------------------------

program selfsimilar

implicit none

integer                                          :: i
integer,parameter                                :: n=5000000
double precision,parameter                       :: cs=1.d6
double precision                                 :: b,reff0,gm
double precision                                 :: x0,y0,u0,dx,dy
double precision                                 :: x1,y1,u1,tanth1,a1,x2,y2,u2,tanth2,a2
double precision                                 :: ft,rhs1,rhs2,rhs3,rh,rhs,lhs1,lhs2,lhs
double precision                                 :: dudy,ddx1,ddx2,ddxdy
double precision                                 :: dphi,dphidy,gefft
external                                         :: derivs

!ub=(0.92,0.85,0.77,0.56,0.29)
!b=(0.5,0.75,1.,1.5,2.)

!! INITIAL CONDITIONS
x0=1.
y0=1.
u0=0.85
b=0.75
! I don't know what this parameter does,but it's very important!
! It should be the gravity set to zero
gm=0.0 !1.d-6

reff0=u0*u0*x0/b

dx=1.d-8
dy=1.d-6

x1=x0+dx
y1=(2.*reff0*dx-dx*dx)**0.5
u1=u0

tanth1=(reff0-dx)/y1
!a1=(x1*x1*tanth1-x1*y1)/(1.+tanth1*tanth1)**0.5

do i=1,n
ft=(1.+tanth1*tanth1)**0.5
a1=(x1*x1*tanth1-x1*y1)/(1.+tanth1*tanth1)**0.5
rhs1=b*ft/(x1*tanth1-y1)
rhs2=-u1**2.*(x1*tanth1-y1)/(ft*x1*(x1+y1*tanth1))
rhs3=u1*u0*tanth1*EXP(0.5*(u1**2.-u0**2.)+dphi(x1,y1,gm))*dphidy(x1,y1,tanth1,gm) &
/(x1*(x1+y1*tanth1))
rh=-tanth1*((x1+y1*tanth1)*dphidy(x1,y1,tanth1,gm))/(ft*(x1*tanth1-y1))
rhs=rhs1+rhs2+rhs3+rh+gefft(x1,y1,tanth1,gm)
lhs1=a1*u0*(1.-u1**2.)*tanth1*EXP(0.5*(u1**2.-u0**2.) &
+dphi(x1,y1,gm))/(x1*(x1+y1*tanth1))
lhs2=u1*tanth1*(x1+y1*tanth1)/(ft*(x1*tanth1-y1))
lhs=lhs1+lhs2
dudy=rhs/lhs
ddx1=a1*u0*((u1**(-2.)-1.)*dudy+(dphidy(x1,y1,tanth1,gm)/u1))*EXP(0.5*(u1**2.-u0**2.) &
+dphi(x1,y1,gm))*ft**3./(tanth1**2.*x1*(x1+y1*tanth1))
ddx2=ft**2.*(x1*tanth1-y1)/(tanth1**3.*x1*(x1+y1*tanth1))
ddxdy=ddx1+ddx2
a2=(x1**2.*tanth1-x1*y1)/(1.+tanth1**2.)**0.5
x2=x1+dy/tanth1+0.5*ddxdy*dy**2.
y2=y1+dy
u2=u1+dudy*dy
tanth2=(1./tanth1+ddxdy*dy)**(-1.)
open(unit=8,file='b0.75.txt')
write(8,'(2(es18.10,1X))') x2,y2
x1=x2
y1=y2
u1=u2
tanth1=tanth2
a1=a2
enddo

end program

!--------------------------------------------------------------------------
! Functions
!--------------------------------------------------------------------------

double precision function dphi(x,y,gm)

integer,parameter                                :: n=12000000
double precision                                 :: x,y,r
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
dphi=gm*(-r**(-1.)+0.5*(x)**(-2.)+0.5)
return
end function dphi

double precision function dphidy(x,y,tanth,gm)

integer,parameter                                :: n=12000000
double precision                                 :: x,y,tanth,r
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
dphidy=gm*(((x+y*tanth)/(tanth*r**3.))-(1./(tanth*x**3.)))
return
end function dphidy

double precision function gefft(x,y,tanth,gm)

integer,parameter                                :: n=12000000
double precision                                 :: x,y,tanth,r,ft
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
ft=(1.+tanth**2.)**0.5
gefft=gm*((tanth/(ft*x**3.))-((x*tanth-y)/(ft*r**3.)))
return
end function gefft


