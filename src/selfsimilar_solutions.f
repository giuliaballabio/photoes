          implicit real*8(a-h,o-z)
c
          b_input=1.00
          gm=0.0
          n=0
          pi=4.*datan(1.d0)
          x0=1.
          y0=0.
          ub=0.77
          u0=ub
          phi=datan(y0/x0)
          th=0.5*pi
          reff0=u0**2.*x0/b
          dx=1.d-7
          u1=u0
          x1=x0+dx
          y1=(2.*reff0*dx-dx*dx)**0.5
          th1=datan((reff0-dx)/y1)
          dy=1.d-5
          tanth1=dtan(th1)
          al=x1**2.*tanth1-x1*y1
          al=al/(1.+tanth1**2.)**0.5
          a0=al
1         continue
c         return here each timestep
          n=n+1
          ft1=(1.+tanth1**2.)**0.5
          rhs=b*ft1/(x1*tanth1-y1)
          rhs1=rhs
          rhs=rhs-u1**2.*(x1*tanth1-y1)/(ft1*x1*(x1+y1*tanth1))
          rhs2=u1**2.*(x1*tanth1-y1)/(ft1*x1*(x1+y1*tanth1))
          rh=(x1+y1*tanth1)*dphidy(gm,tanth1,x1,y1)/(ft1*(x1*tanth1-y1))
          rh=-tanth1*rh
          rhs=rhs+rh
          rhs3=u1*tanth1*u0/(x1*(x1+y1*tanth1))
          rhs3=rhs3*dexp(dphi(gm,x1,y1)+0.5*(u1**2.-u0**2.))
          rhs3=rhs3*dphidy(gm,tanth1,x1,y1)
          rhs=rhs+rhs3+gefft(gm,tanth1,x1,y1)
           alhs=a0*u0*(1.-u1**2.)*tanth1/(x1*(x1+y1*tanth1))
           alhs=alhs*dexp(0.5*(u1**2.-u0**2.)+dphi(gm,x1,y1))
           alhs1=alhs
           alhs=alhs+u1*tanth1*(x1+y1*tanth1)/(ft1*(x1*tanth1-y1))
           alhs2=u1*tanth1*(x1+y1*tanth1)/(ft1*(x1*tanth1-y1))
           ratio=alhs1/alhs2
           dudy=rhs/alhs
           xdd=(u1**(-2.)-1.)*dudy
           xdd=xdd+(dphidy(gm,tanth1,x1,y1))/u1
           xdd=xdd*dexp(0.5*(u1**2.-u0**2.)+dphi(gm,x1,y1))
           xdd=xdd*a0*u0
c           write(6,*) dphidy(gm,tanth1,x1,y1),dphi(gm,x1,y1)
           xdd0=xdd
           xdd=xdd*ft1**3./(tanth1**2.*x1*(x1+y1*tanth1))
           xdd=xdd+ft1**2.*(x1*tanth1-y1)/(tanth1**3.*x1*(x1+y1*tanth1))
           x2=x1+dy/tanth1+0.5*xdd*dy*dy
           y2=y1+dy
           tanth2=(1./tanth1+xdd*dy)**(-1.)
           u2=u1+dudy*dy
           a=x1**2.*tanth1-x1*y1
           a=a/ft1
           dady=(a-al)/dy
           dadych=a0*u0*(1.-u1**(-2.))*dudy
           dadych=dadych+a0*(u0/u1)*dphidy(gm,tanth1,x1,y1)
           dadych=dadych*dexp(0.5*(u1**2.-u0**2.)+dphi(gm,x1,y1))
c           write(48,*) dady,dadych
           rho=dexp(-0.5*(u1**2.-u0**2.)-dphi(gm,x1,y1))
           flux=rho*u1*a
           be=0.5*u1*u1+dphi(gm,x1,y1)-dlog(u1*a)
           if(dfloat(n)/100..eq.dfloat(n/100)) then
c           write(492,100) x2,y2,u2,flux,be,alhs
           open(unit=1,file='./streamline_cartcoord.txt')
           write(1,'(3(es18.10,1X))') x2,y2
           ratiow=(x1+y1*tanth1)/(x1*tanth1-y1)
           ratiow2=ratiow*tanth1/ft1
           comp=(tanth1*u1*u0*dexp(0.5*u1*u1))/(x1*(x1+y1*tanth1))
           usq=u2*u2
c           write(493,100) x2,y2,usq,ratiow,ratiow2,comp
            if(dabs(alhs).lt.1.d-2.or.alhs.lt.0.) then
c           write(494,*) x2,y2,u2,flux,rhs,alhs
            stop
            end if
           end if
           if(tanth2.lt.y2/x2) stop
           phi2=datan(y2/x2)
           theta2=datan(tanth2)
           duds=(u2-u1)*(dsin(theta2)/dy)
           r=(x2**2.+y2**2.)**0.5
            rmom = b/(r*dsin(theta2-phi2))
           rmom=rmom-u2*duds/dtan(theta2-phi2)
           ur2=u1*dcos(theta2-phi2)
           uth2=u1*dsin(theta2-phi2)
           reff=(1.+tanth2**2.)**1.5/(tanth2**3.*xdd)
           almom=u2**2./reff
           u1=u2
           x1=x2
           y1=y2
           al=a
           tanth1=tanth2
           open(unit=2,file='./streamline_polarcoord.txt')
           write(2,'(2(es18.10,1X))') r,phi2
           open(unit=4,file='./rhov_fields.txt')
           write(4,'(4(es18.10,1X))') rho, ur2, uth2 !!, uph
           if(dfloat(n)/1..eq.dfloat(n/1)) then
c           write(482,*)almom,rmom
           end if
           if(n.lt.12000000) goto 1
           stop
100        format(6(d11.3))
           end
           function gefft(gm,tanth,x,y)
           implicit real*8(a-h,o-z)
           r=(x*x+y*y)**0.5
           ft=(1.+tanth*tanth)**0.5
           gefft=gm*tanth/(ft*x**3.)
           gefft=gefft-(gm/(ft*r**3.))*(x*tanth-y)
           return
           end
           function dphidy(gm,tanth,x,y)
           implicit real*8(a-h,o-z)
           dphidy=(x+y*tanth)*(x**2.+y**2.)**(-1.5)
           dphidy=dphidy/tanth-(tanth*x**3.)**(-1.)
           dphidy=dphidy*gm
           return
           end
           function dphi(gm,x,y)
           implicit real*8(a-h,o-z)
           r=(x**2.+y**2.)**0.5
           dphi=-r**(-1.)+0.5*(x**(-2.))+0.5
           dphi=dphi*gm
           return
           end
