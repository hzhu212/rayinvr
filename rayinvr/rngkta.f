c                 
c     version 1.3  Aug 1992
c
c     Runge-Kutta routine for RAYINVR
c
c     ----------------------------------------------------------------
c                 
      subroutine rngkta(x,z,y,f,h,hmin,e,func,g,s,t) 
c                 
c     routine to solve a 2x2 system of first order o.d.e.'s 
c     using the runge-kutta method with error control
c                 
      real y(2),f(2),t(2),s(2),g(2) 
      logical*1 be,bh,br,bx 
      logical bxx
      common /rkc$/ bxx
c                 
      bh=.true.   
      br=.true.   
      bx=.true.   
      bxx=.true.
      if(z.gt.x.and.h.lt.0.) h=-h
      if(z.lt.x.and.h.gt.0.) h=-h
10    xs=x        
c                 
      g(1)=y(1) 
      g(2)=y(2) 
c                 
20    hs=h        
      q=x+h-z     
      be=.true.   
      if((h.gt.0..and.q.ge.0.).or.(h.lt.0..and.q.le.0.)) then 
        h=z-x       
        br=.false.  
      end if
      h3=h/3.     
c                 
      call func(x,y,f) 
c              
      do 210 i=1,2 
         q=h3*f(i) 
         t(i)=q 
         r=q   
         y(i)=g(i)+r 
210   continue 
      x=x+h3   
c
      call func(x,y,f) 
c                 
      do 220 i=1,2 
         q=h3*f(i) 
         r=.5*(q+t(i)) 
         y(i)=g(i)+r 
220   continue 
c                 
      call func(x,y,f) 
c                 
      do 230 i=1,2 
         q=h3*f(i) 
         r=3.*q 
         s(i)=r 
         r=.375*(r+t(i)) 
         y(i)=g(i)+r 
230   continue 
      x=x+.5*h3 
c
      call func(x,y,f) 
c                 
      do 240 i=1,2 
         q=h3*f(i) 
         r=t(i)+4.*q 
         t(i)=r 
         r=1.5*(r-s(i)) 
         y(i)=g(i)+r 
240   continue 
      x=x+.5*h 
c
      call func(x,y,f) 
c                 
      do 250 i=1,2 
         q=h3*f(i) 
         r=.5*(q+t(i)) 
         q=abs(r+r-1.5*(q+s(i))) 
         y(i)=g(i)+r 
         r=abs(y(i)) 
         if(r.lt..001) then
           r=e  
         else
           r=e*r 
         end if
         if(q.ge.r.and.bx) then
           br=.true. 
           bh=.false. 
           h=.5*h 
           if(abs(h).lt.hmin) then
             sigh=1. 
             if(h.lt.0.) sigh=-1. 
             h=sigh*hmin 
             bx=.false. 
             bxx=.false.
           end if
           y(1)=g(1) 
           y(2)=g(2) 
           x=xs  
           go to 20 
         end if
         if(q.ge..03125*r) be=.false. 
250   continue 
c                 
      if(be.and.bh.and.br) then 
        h=h+h       
        bx=.true.   
      end if
      bh=.true.   
      if(br) go to 10 
      h=hs        
c
      return      
      end         
c
c     ----------------------------------------------------------------
c
      subroutine rkdumb(y,x1,x2,derivs)
c                 
c     routine to solve a 2x2 system of first order o.d.e.'s 
c     using a 4th-order runge-kutta method without error control
c                 
      dimension y(2),dydx(2),yt(2),dyt(2),dym(2)
      x=x1
      h=x2-x1
      if(x+h.eq.x) return
c
      call derivs(x,y,dydx)
c
      hh=h*0.5
      h6=h/6.
      xh=x+hh
      do 11 i=1,2
        yt(i)=y(i)+hh*dydx(i)
11    continue
c
      call derivs(xh,yt,dyt)
c
      do 12 i=1,2
        yt(i)=y(i)+hh*dyt(i)
12    continue
c
      call derivs(xh,yt,dym)
c
      do 13 i=1,2
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
c
      call derivs(x+h,yt,dyt)
c
      do 14 i=1,2
        y(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
14    continue
c
      return
      end
c        
c     ----------------------------------------------------------------
c        
      subroutine odexfi(x,y,f)
c      
c     pair of first order o.d.e.'s solved by runge kutta method
c     with x as independent variable
c      
      include 'rayinvr.par'
      real y(2),f(2)
      include 'rayinvr.com'
c 
      sa=sign(1.,y(2))
      n1=int(sa*y(2)*factan)+1
      f(1)=mcotan(n1)*y(2)+sa*bcotan(n1)
      term1=c(layer,iblk,3)+c(layer,iblk,4)*x
      term2=x*(c(layer,iblk,1)+c(layer,iblk,2)*x)+
     +       term1*y(1)+c(layer,iblk,5)
      vxv=(x*(c(layer,iblk,8)+c(layer,iblk,9)*x)+
     +    c(layer,iblk,10)*y(1)+c(layer,iblk,11))/
     +    (term2*(c(layer,iblk,6)*x+c(layer,iblk,7)))
      vzv=term1/term2
      f(2)=vzv-vxv*f(1)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine odezfi(x,y,f)
c
c     pair of first order o.d.e.'s solved by runge kutta method
c     with z as independent variable
c
      include 'rayinvr.par'
      real y(2),f(2)
      include 'rayinvr.com'
c 
      sa=sign(1.,y(2))
      n1=int(sa*y(2)*factan)+1
      f(1)=mtan(n1)*y(2)+sa*btan(n1)
      term1=c(layer,iblk,3)+c(layer,iblk,4)*y(1)
      term2=y(1)*(c(layer,iblk,1)+c(layer,iblk,2)*y(1))+
     +       term1*x+c(layer,iblk,5)
      vxv=(y(1)*(c(layer,iblk,8)+c(layer,iblk,9)*y(1))+
     +    c(layer,iblk,10)*x+c(layer,iblk,11))/
     +    (term2*(c(layer,iblk,6)*y(1)+c(layer,iblk,7)))
      vzv=term1/term2
      f(2)=vzv*f(1)-vxv
c
      return
      end
