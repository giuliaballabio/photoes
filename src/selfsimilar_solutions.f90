
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
integer,parameter                                :: n=1200000
double precision,parameter                       :: cs=1.d6
double precision                                 :: b,reff0,gm
double precision                                 :: x0,y0,u0,dx,dy
double precision                                 :: x(n),y(n),u(n),tanth(n),a(n),dudy(n) !,x2,y2,u2,tanth2,a2
double precision                                 :: ft !,rhs1,rhs2,rhs3,rh,rhs,lhs1,lhs2,lhs
double precision                                 :: ddx1,ddx2,ddxdy
real                                             :: dphi,dphidy,gefft
external                                         :: deriv

!ub=(0.92,0.85,0.77,0.56,0.29)
!b=(0.5,0.75,1.,1.5,2.)

!! INITIAL CONDITIONS
x0=1.
y0=0.
u0=0.85
b=0.75
! I don't know what this parameter does,but it's very important!
! It should be the gravity set to zero
gm=0.0 !1.d-6

reff0=u0*u0*x0/b

dx=1.d-8
dy=1.d-6

x(1)=x0+dx
y(1)=(2.*reff0*dx-dx*dx)**0.5
u(1)=u0

tanth(1)=(reff0-dx)/y(1)
a(1)=(x(1)*x(1)*tanth(1)-x(1)*y(1))/(1.+tanth(1)*tanth(1))**0.5

do i=1,n
    ft=(1.+tanth(i)*tanth(i))**0.5
    call deriv(y(i),u(i),dudy(i))
    print *,dudy(i)
    ddx1=a(i)*u0*((u(i)**(-2.)-1.)*dudy(i)+(dphidy(x(i),y(i),tanth(i),gm) &
        /u(i)))*EXP(0.5*(u(i)**2.-u0**2.)+dphi(x(i),y(i),gm))*ft**3./(tanth(i)**2.*x(i)*(x(i)+y(i)*tanth(i)))
    ddx2=ft**2.*(x(i)*tanth(i)-y(i))/(tanth(i)**3.*x(i)*(x(i)+y(i)*tanth(i)))
    ddxdy=ddx1+ddx2
    a(i+1)=(x(i)**2.*tanth(i)-x(i)*y(i))/(1.+tanth(i)**2.)**0.5
    x(i+1)=x(i)+dy/tanth(i)+0.5*ddxdy*dy**2.
    y(i+1)=y(i)+dy
    !call derivs(y(i),u(i),dudy(i))
    !call rungekutta(u(i),dudy(i),n,y(i),dy,u(i+1),derivs)
    !call rungekutta(y(i),u(i),dy,dudy(x(i),y(i),u(i),a(i),tanth(i)))
    call rk4(y(i),u(i),dy,deriv,u(i+1))
    !u(i+1)=u(i)+dudy(i)*dy
    tanth(i+1)=(1./tanth(i)+ddxdy*dy)**(-1.)
    open(unit=8,file='b0.75.txt')
    write(8,'(2(es18.10,1X))') x(i+1),y(i+1)
    x(i)=x(i+1)
    y(i)=y(i+1)
    u(i)=u(i+1)
    tanth(i)=tanth(i+1)
    a(i)=a(i+1)
enddo

end program

!--------------------------------------------------------------------------
! Functions
!--------------------------------------------------------------------------

real function dphi(x,y,gm)

implicit none
double precision                                 :: x,y,r
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
dphi=gm*(-r**(-1.)+0.5*(x)**(-2.)+0.5)
return
end function dphi

real function dphidy(x,y,tanth,gm)

implicit none
double precision                                 :: x,y,tanth,r
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
dphidy=gm*(((x+y*tanth)/(tanth*r**3.))-(1./(tanth*x**3.)))
return
end function dphidy

real function gefft(x,y,tanth,gm)

implicit none
double precision                                 :: x,y,tanth,r,ft
double precision                                 :: gm

r=(x**2.+y**2.)**0.5
ft=(1.+tanth**2.)**0.5
gefft=gm*((tanth/(ft*x**3.))-((x*tanth-y)/(ft*r**3.)))
return
end function gefft

!--------------------------------------------------------------------------
! Subroutines
!--------------------------------------------------------------------------

subroutine deriv(y,u,uprime)

double precision                                 :: x,y,u,a,tanth
double precision                                 :: ft,rhs1,rhs2,rhs3,rh,rhs,lhs1,lhs2,lhs
double precision                                 :: b,u0,gm,uprime

b=0.75
gm=0.0
dx=1.d-8
reff0=u0*u0/b
tanth=(reff0-dx)/y0
a=(x*x*tanth-x*y)/(1.+tanth*tanth)**0.5

ft=(1.+tanth*tanth)**0.5
rhs1=b*ft/(x*tanth-y)
rhs2=-u**2.*(x*tanth-y)/(ft*x*(x+y*tanth))
rhs3=u*u0*tanth*EXP(0.5*(u**2.-u0**2.)+dphi(x,y,gm))*dphidy(x,y,tanth,gm)/(x*(x+y*tanth))
rh=-tanth*((x+y*tanth)*dphidy(x,y,tanth,gm))/(ft*(x*tanth-y))
rhs=rhs1+rhs2+rhs3+rh+gefft(x,y,tanth,gm)
lhs1=a*u0*(1.-u**2.)*tanth*EXP(0.5*(u**2.-u0**2.) &
+dphi(x,y,gm))/(x*(x+y*tanth))
lhs2=u*tanth*(x+y*tanth)/(ft*(x*tanth-y))
lhs=lhs1+lhs2
uprime=rhs/lhs
return
end

subroutine rk4(t0,u0,dt,f,u)

!*****************************************************************************
!
!! RK4 takes one Runge-Kutta step for a scalar ODE.
!
!  Discussion:
!
!    It is assumed that an initial value problem, of the form
!
!      du/dt=f(t,u)
!      u(t0)=u0
!
!    is being solved.
!
!    If the user can supply current values of t, u, a stepsize dt, and a
!    function to evaluate the derivative, this function can compute the
!    fourth-order Runge Kutta estimate to the solution at time t+dt.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real (kind=8) T0, the current time.
!
!    Input, real (kind=8) U0, the solution estimate at the current time.
!
!    Input, real (kind=8) DT, the time step.
!
!    Input, external F, a subroutine of the form
!      subroutine f(t,u,uprime)
!    which evaluates the derivative uprime given the time T and
!    solution vector U.
!
!    Output, real (kind=8) U, the fourth-order Runge-Kutta solution
!    estimate at time T0+DT.
!
implicit none

real (kind=8) dt
external f
real (kind=4) f0
real (kind=4) f1
real (kind=4) f2
real (kind=4) f3
real (kind=8) t0
real (kind=8) t1
real (kind=8) t2
real (kind=8) t3
real (kind=8) u
real (kind=8) u0
real (kind=8) u1
real (kind=8) u2
real (kind=8) u3
!
!  Get four sample values of the derivative.
!
call f(t0,u0,f0)

t1=t0+dt/2.0D+00
u1=u0+dt*f0/2.0D+00
call f(t1,u1,f1)

t2=t0+dt/2.0D+00
u2=u0+dt*f1/2.0D+00
call f(t2,u2,f2)

t3=t0+dt
u3=u0+dt*f2
call f(t3,u3,f3)
!
!  Combine them to estimate the solution U at time T1.
!
u=u0+dt*(f0+2.0D+00*f1+2.0D+00*f2+f3)/6.0D+00

return
end

!subroutine derivs(x,func,dfuncdx)
!
!double precision                                 :: x,dfuncdx,func,dx
!
!dfuncdx=(func(x+dx)-func(x-dx))/(2.*dx)
!return
!end subroutine
!
!subroutine rungekutta(y,dydx,n,x,h,yout,derivs)
!
!! Given the values for the variables y(1:n) and their derivatives dydx(1:n) known at x,use the RK4 method to advance the
!! solution over an interval h and return the incremented variables as yout(1:n),which need not be a distinct array from y.
!! The user supplies the subroutine derivs(x,y,dydx),which returns derivatives dydx at x.
!
!integer                                          :: n,j
!integer,parameter                                :: nmax=1
!double precision                                 :: h,x,dydx(n),y(n),yout(n)
!external                                         :: derivs
!double precision                                 :: h6,hh,xh,dym(nmax),dyt(nmax),yt(nmax)
!
!hh=h*0.5
!h6=h/6.
!xh=x+hh
!
!! First step
!do j=1,n
!    yt(j)=y(j)+hh*dydx(j)
!enddo
!! Second step
!call derivs(xh,yt,dyt)
!do j=1,n
!    yt(j)=y(j)+hh*dyt(j)
!enddo
!! Third step
!call derivs(xh,yt,dym)
!do j=1,n
!    yt(j)=y(j)+h*dym(j)
!    dym(j)=dyt(j)+dym(j)
!enddo
!! Fourth step
!call derivs(x+h,yt,dyt)
!do j=1,n
!    yout(j)=y(j)+h6*(dydx(j)+dyt(j)+2.*dym(j))
!enddo
!return
!end subroutine
!
!subroutine rungekutta(x,y,h,f)
!
!implicit none
!integer                                          :: i
!integer,parameter                                :: n=1200000
!double precision                                 :: x,y,h,m1,m2,m3,m4
!real                                             :: f
!
!do i=1,n
!m1=h*f(x,y)
!m2=h*f(x+0.5*h,y+0.5*m1)
!m3=h*f(x+0.5*h,y+0.5*m2)
!m4=h*f(x+h,y+m3)
!x=x+h
!y=y+(m1+2.0*m2+2.0*m3+m4)*(1/6.0)
!enddo
!return
!end subroutine
