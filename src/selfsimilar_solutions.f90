program selfsimilarity

implicit none

integer                          :: i,n
real, dimension(1:50001)         :: x,y,u,dx,du,dA,ddx
real                             :: du_dy,dA_dy,ddx_dy
!! BIN WIDTH
real, parameter                  :: dy=0.0001,ymax=5.
!! DIMENSIONAL VARIABLES
real, parameter                  :: cs=1.e5,ub=0.85*cs,b=0.75
integer                          :: infile
logical,save                     :: init=.false.

n=int(ymax/dy)

!! INITIAL CONDITIONS AT P
x(1)=1.
y(1)=0.0
u(1)=ub/cs
dx(1)=0.

do i=1,n
    du(i)=du_dy(x(i),y(i),u(i),dx(i))
    y(i+1)=y(i)+dy
    !! EULER METHOD TO CALCULATE u
    u(i+1)=u(i)+dy*du(i)
    !! CALCULTE dA AND ddx FROM EQS. 14 & 15
    dA(i)=dA_dy(x(i),y(i),u(i),dx(i))
    ddx(i)=ddx_dy(x(i),y(i),u(i),dx(i),dA(i))
    !! CALCULATE x AND dx AT THE SUBSEQUENT POINT
    dx(i+1)=dx(i)+ddx(i)*dy
    x(i+1)=x(i)+dx(i)*dy+0.5*ddx(i)*dy*dy
enddo

if (.not.init) then
    open(infile,file='output.txt')
    !write(infile,"('# ',2('[',i2.2,1x,a12,']',1x))") 1,'x',2,'y'
else
    open(infile,file='output.txt',status='old',position='append')
endif
do i=1,n
    write(infile,'(2(es18.10,1X))') x(i),y(i)
enddo
close(infile)

end program

!--------------------------------------------------------------------------
! Useful Functions
!--------------------------------------------------------------------------

real function f(x,y,u,dx)

real                             :: x,y,u,dx
real                             :: c1,c2,e,d
real, parameter                  :: cs=1.e5,ub=0.85*cs

c1=-(ub**2*u/cs)**2
e=0.5*(ub/cs)**2*(u**2-1.)
c2=u*ub**2
d=dx*(x*dx+y)

if (d.eq.0) then
    f=(ub*ub/cs)**2
else
    f=(c1/d)*EXP(e)*(1.-(cs/(ub*u))**2)+c2*(x*dx+y)/((1.+dx**2)**0.5*(x-y*dx))
endif

end function f

real function g(x,y,u,dx)

real                             :: x,y,u,dx
real                             :: c1,c2,d
real, parameter                  :: cs=1.e5,ub=0.85*cs,b=0.75

c1=b*cs**2
c2=(ub*u)**2
d=(1+x**2)**0.5*x*(x*dx+y)

if (d.eq.0) then
    g=b*(cs**2)
else
    g=c1*(1.+dx**2)**0.5/(x-y*dx)-c2*dx*(x-y*dx)/d
endif

end function g

real function du_dy(x,y,u,dx)

real                             :: x,y,u,dx

du_dy=g(x,y,u,dx)/f(x,y,u,dx)

end function du_dy

real function dA_dy(x,y,u,dx)

real                             :: x,y,u,dx
real                             :: a,e
real, parameter                  :: cs=1.e5,ub=0.85*cs,b=0.75

a=(ub/cs)**2-(1/u)**2
e=0.5*(ub/cs)**2*(u**2-1.)

dA_dy=a*du_dy(x,y,u,dx)*EXP(e)

end function dA_dy

real function ddx_dy(x,y,u,dx,dA)

real                             :: x,y,u,dx,dA
real                             :: d

d=x*y*(1.+(dx**2))+dx !x*(y+x*dx)

if (d.eq.0) then
    ddx_dy=1/x     
    !487.089
    !0.859871e4
else
    ddx_dy=((1.+dx)**2*dx*(x-y*dx)-(1.+dx**2)**1.5*dA)/d
endif

end function ddx_dy
