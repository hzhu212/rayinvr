c
c     version 1.3  Aug 1992
c     
c     Ray tracing routines for TRAMP
c                 
c     ----------------------------------------------------------------
c                 
      subroutine trace(npt,ifam,ir,iturn,iamp,ihamp,iflag,idl,idr)
c                 
c     trace a single ray through the model
c                 
      include 'tramp.par'
      external odex,odez,odexfi,odezfi
      real y(2),f(2),w1(2),w2(2),w3(2),left
      include 'tramp.com'
c
      ntray=ntray+1
      dstepf=1.   
      n2=0        
      n3=0        
      iflag=0     
      if(ivg(layer,iblk).gt.0) then
        hn=dstep(xr(npt),zr(npt))/hdenom
      else        
        hn=smax/hdenom
      end if      
      hminn=hn*hmin
1000  if(idump.eq.1) write(12,5) ifam,ir,npt,xr(npt),zr(npt),
     +  ar(npt,1)*pi18,ar(npt,2)*pi18,vr(npt,1),vr(npt,2),
     +  layer,iblk,id,iwave
5     format(i2,i3,i4,2f8.3,2f8.2,2f7.2,4i3)
      if(iwave.eq.-1.and.vr(npt,2).le..001) then
        vr(npt,2)=0.0 
        iflag=1   
        if(ir.ne.0) write(11,25)
25      format('***  ray stopped - s-wave cannot propagate  ***')
        go to 900 
      end if      
      if(npt.gt.1) then
        if(((ar(npt,2)-fid*pi2)*(ar(npt-1,2)-fid*pi2)).le.0.) then
          if(iturn.ne.0) go to 900
          if(iwave.eq.-1) iep=iep+1
        end if    
      end if      
      if(ihamp.eq.1.and.npt.ne.1) then
        s1=s(layer,iblk,1)
        s2=s(layer,iblk,2)
        z1=s1*xr(npt)+b(layer,iblk,1)
        z2=s2*xr(npt)+b(layer,iblk,2) 
        sh(npt-1)=(s1*(z2-zr(npt))+s2*(zr(npt)-z1))/(z2-z1)
      end if      
      isrkc=1     
      if((fid*ar(npt,2)).ge.pi4.and.(fid*ar(npt,2)).le.pi34) then
c                 
c       solve o.d.e.'s w.r.t. x
c                 
        x=xr(npt) 
        y(1)=zr(npt)
        y(2)=ar(npt,2)
c                 
        if(ivg(layer,iblk).eq.0) then
c                 
          call strait(x,y,npt,0)
c                 
        else      
c
          z=x+fid*dstep(xr(npt),zr(npt))/dstepf
c
          call check(0,z,zr(npt))
c                 
          if(ifast.eq.0) then
            call rngkta(x,z,y,f,hn,hminn,tol,odex,w1,w2,w3)
          else
            call rkdumb(y,x,z,odexfi)
            x=z
          end if
c                 
        end if    
        if(npt.eq.ppray) go to 999
        npt=npt+1 
        xr(npt)=x 
        zr(npt)=y(1)
      else        
c                 
c       solve o.d.e.'s w.r.t. z
c                 
        x=zr(npt) 
        y(1)=xr(npt)
        y(2)=ar(npt,2)
c                 
        if(ivg(layer,iblk).eq.0) then
c                 
          call strait(x,y,npt,1)
c                 
        else      
c
          if((fid*ar(npt,2)).le.pi2) then
            z=x+dstep(xr(npt),zr(npt))/dstepf
          else    
            z=x-dstep(xr(npt),zr(npt))/dstepf
          end if  
c
          call check(1,xr(npt),z)
c                 
          if(ifast.eq.0) then
            call rngkta(x,z,y,f,hn,hminn,tol,odez,w1,w2,w3)
          else
            call rkdumb(y,x,z,odezfi)
            x=z
          end if
c                 
        end if    
        if(npt.eq.ppray) go to 999
        npt=npt+1 
        xr(npt)=y(1)
        zr(npt)=x 
      end if      
      ar(npt,1)=y(2)
      ar(npt,2)=y(2)
      vp(npt,1)=vel(xr(npt),zr(npt))
      vs(npt,1)=vp(npt,1)*vsvp(layer,iblk)
      if(iwave.eq.1) then
        vr(npt,1)=vp(npt,1)
      else        
        vr(npt,1)=vs(npt,1) 
      end if      
      vr(npt,2)=vr(npt,1)
      vp(npt,2)=vp(npt,1)
      vs(npt,2)=vs(npt,1)
      qr(npt)=q(layer,iblk,(3-iwave)/2)
      top=s(layer,iblk,1)*xr(npt)+b(layer,iblk,1)
      bottom=s(layer,iblk,2)*xr(npt)+b(layer,iblk,2)
      left=xbnd(layer,iblk,1)
      right=xbnd(layer,iblk,2)
      if(zr(npt).gt.top.and.zr(npt).le.bottom.and.xr(npt).gt.left.
     +   and.xr(npt).le.right) then
        dstepf=1.
        go to 1000
      end if
c                 
      call adjpt(npt,top,bottom,left,right,ifam,ir,iturn,iamp,
     +           iflag,idl,idr)
c                 
      qr(npt)=q(layer,iblk,(3-iwave)/2)
      if(iflag.eq.0) go to 1000
      go to 900   
999   vr(500,2)=0.0
      iflag=1     
      if(ir.ne.0) write(11,15)
15    format('***  ray stopped - consists of too many points  ***')
900   ntpts=ntpts+npt
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine xzpt(xpt,zpt,layers,iblks,iflag)
c                 
c     determine layer and block of (x,y) point in model
c      
      include 'tramp.par'
      real left 
      include 'tramp.com'           
c
      iflag=0     
      do 10 i=1,nlayer
         do 20 j=1,nblk(i)
            top=s(i,j,1)*xpt+b(i,j,1)
            bottom=s(i,j,2)*xpt+b(i,j,2)
            left=xbnd(i,j,1)
            right=xbnd(i,j,2)
            if(zpt.ge.top.and.zpt.le.bottom.and.xpt.ge.left.
     +         and.xpt.le.right) then
              layers=i
              iblks=j
              return
            end if
20       continue 
10    continue    
      iflag=1     
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine block(x,layer1,iblk1)
c                 
c     determine block of point x in layer
c      
      include 'tramp.par'
      include 'tramp.com'           
c
      do 10 i=1,nblk(layer1)
         if(x.ge.xbnd(layer1,i,1).and.x.le.xbnd(layer1,i,2)) then
           iblk1=i
           return 
         end if   
10    continue    
c
      write(6,5) layer1,x 
5     format(/'***  block undetermined  ***'/
     +'layer=',i3,'  x=',f10.3/)
c
      stop        
      end         
c                
c     ----------------------------------------------------------------                
c                 
      subroutine check(iflag,x,z)
c                 
c     check that the new ray point will still be inside the current
c     trapezoid; if not, adjust the step length so it falls 1 meter
c     outside the trapezoid
c        
      include 'tramp.par'
      include 'tramp.com'
c   
      if(iflag.eq.0) then
        if(x.lt.xbnd(layer,iblk,1)) x=xbnd(layer,iblk,1)-.001
        if(x.gt.xbnd(layer,iblk,2)) x=xbnd(layer,iblk,2)+.001
      else
        z1=s(layer,iblk,1)*xbnd(layer,iblk,1)+b(layer,iblk,1)
        z2=s(layer,iblk,1)*xbnd(layer,iblk,2)+b(layer,iblk,1)
        zt=s(layer,iblk,1)*x+b(layer,iblk,1)
        z3=s(layer,iblk,2)*xbnd(layer,iblk,1)+b(layer,iblk,2)
        z4=s(layer,iblk,2)*xbnd(layer,iblk,2)+b(layer,iblk,2)
        zb=s(layer,iblk,2)*x+b(layer,iblk,2)
        if(id.lt.0) then
          zt=amin1(z1,zt)
          zb=amax1(z3,zb)
        else
          zt=amin1(z2,zt)
          zb=amax1(z4,zb)
        end if
        if(z.lt.zt) z=zt-.001
        if(z.gt.zb) z=zb+.001
      end if
c
      return
      end
c                 
c     ----------------------------------------------------------------
c                 
      function dstep(x,z)
c                 
c     calculate ray step length - step length is a function of the
c     absolute value of the velocity divided by the sum of the
c     absolute values of the velocity gradients in the x and z 
c     directions  
c      
      include 'tramp.par'
      include 'tramp.com' 
c   
      vx=(c(layer,iblk,8)*x+c(layer,iblk,9)*x**2+c(layer,iblk,10)*z+
     +    c(layer,iblk,11))/(c(layer,iblk,6)*x+c(layer,iblk,7))**2
      vz=(c(layer,iblk,3)+c(layer,iblk,4)*x)/
     +   (c(layer,iblk,6)*x+c(layer,iblk,7))
      dstep=step*vel(x,z)/(abs(vx)+abs(vz))
      if(dstep.lt.smin) dstep=smin
      if(dstep.gt.smax) dstep=smax 
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      function vel(x,z)
c                 
c     calculate p-wave velocity at point (x,z) in model
c      
      include 'tramp.par'
      include 'tramp.com' 
c            

      vel=(c(layer,iblk,1)*x+c(layer,iblk,2)*x**2+c(layer,iblk,3)*z
     +   +c(layer,iblk,4)*x*z+c(layer,iblk,5))/(c(layer,iblk,6)*x
     +   +c(layer,iblk,7))
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine adhoc(x1,z1,a1,x2,z2,a2,mb,bb,it,x,z,a)
c                 
c     calculate intersection point (x,z) along boundary given by:
c                 
c        it=0 --  slope mb, intercept bb
c        it=1 --  x=bb
c                 
c     between points (x1,z1) and (x2,z2) using a straight line path
c     between (x1,z1) and (x2,z2). also calculate angle of 
c     incidence (a) at boundary
c                 
      real mb,mr
      if((x1-x2).eq.0..and.(z1-z2).eq.0.) then
        x=x1      
        z=z1      
        a=(a1+a2)/2.
      else        
        if((10.*abs(x2-x1)).gt.abs(z2-z1)) then
          mr=(z2-z1)/(x2-x1)
          br=z1-mr*x1
          if(it.eq.0) then
            if((mb-mr).eq.0.) then
              x=(x1+x2)/2.
            else  
              x=(br-bb)/(mb-mr)
            end if
          else    
            x=bb  
          end if  
          z=mr*x+br
        else      
          mr=(x2-x1)/(z2-z1)
          br=x1-mr*z1
          if(it.eq.0) then
            if((1.-mr*mb).eq.0.) then
              z=(z1+z2)/2.
            else  
              z=(br*mb+bb)/(1.-mr*mb)
            end if
          else    
            if(abs(mr).gt..000001) then
              z=(bb-br)/mr
            else  
              z=(z1+z2)/2. 
            end if 
          end if  
          x=mr*z+br
        end if    
        d1=sqrt((x1-x)**2+(z1-z)**2)
        d2=sqrt((x2-x)**2+(z2-z)**2)
        if((d1+d2).eq.0.) then
          a=(a1+a2)/2.
        else      
          a=(d1*a2+d2*a1)/(d1+d2)
        end if    
      end if      
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine corner(x,z,x1,z1,a1,i1,x2,z2,a2,i2,xn,zn,an,icase)
c                 
c     determine which point (x1,z1) or (x2,z2) is nearest to (x,z)
c     and assign new point (xn,zn) accordingly as well as an and icase
c                 
      d1=sqrt((x-x1)**2+(z-z1)**2)
      d2=sqrt((x-x2)**2+(z-z2)**2)
      if(d1.lt.d2) then
        icase=i1  
        xn=x1     
        zn=z1     
        an=a1     
      else        
        icase=i2  
        xn=x2     
        zn=z2     
        an=a2     
      end if      
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine strait(x,y,npt,iflag)
c                 
c     if the upper and lower velocities in a layer are nearly the same
c     (or equal) the ray path is essentially a straight line segment.
c     therefore, no need to use runge kutta routine. this routine
c     calculates new ray coordinates and tangent for a zero gradient
c     layer.      
c                 
c                  iflag=0 --  x is present x coordinate
c                  iflag=1 --  x is present z coordinate
c      
      include 'tramp.par'
      real y(2)
      include 'tramp.com'           
c
      if(iflag.eq.0) then
        x=xr(npt)+fid*smax/dstepf
        y(1)=zr(npt)+(x-xr(npt))/tan(ar(npt,2))
        y(2)=ar(npt,2)
      else        
        if((fid*ar(npt,2)).le.pi2) then
          x=x+smax/dstepf
        else      
          x=x-smax/dstepf
        end if    
        y(1)=xr(npt)+tan(ar(npt,2))*(x-zr(npt))
        y(2)=ar(npt,2)
      end if      
      return      
      end         
c        
c     ----------------------------------------------------------------
c        
      subroutine odex(x,y,f)
c      
c     pair of first order o.d.e.'s solved by runge kutta method
c     with x as independent variable
c      
      include 'tramp.par'
      real y(2),f(2)
      logical ok
      include 'tramp.com'
      common /rkc$/ ok
c 
      f(1)=1./tan(y(2))
      term1=c(layer,iblk,3)+c(layer,iblk,4)*x
      term2=c(layer,iblk,1)*x+c(layer,iblk,2)*x*x+
     +       term1*y(1)+c(layer,iblk,5)
      vxv=(c(layer,iblk,8)*x+c(layer,iblk,9)*x*x+
     +    c(layer,iblk,10)*y(1)+c(layer,iblk,11))/
     +    (term2*(c(layer,iblk,6)*x+c(layer,iblk,7)))
      vzv=term1/term2
      f(2)=vzv-vxv*f(1)
c
      if(.not.ok) then
        if(isrkc.eq.1.and.idump.eq.1) then
          write(12,5) vel(x,y(1))*(vxv**2+vzv**2)**.5
5         format('***  possible inaccuracies in rngkta w.r.t. x  ***'/
     +           '|del(v)| = ',f12.5)
        end if
        irkc=1
        isrkc=0
      end if
      return
      end
c        
c     ----------------------------------------------------------------
c        
      subroutine odez(x,y,f)
c      
c     pair of first order o.d.e.'s solved by runge kutta method
c     with z as independent variable
c      
      include 'tramp.par'
      real y(2),f(2)
      logical ok
      include 'tramp.com'
      common /rkc$/ ok
c 
      f(1)=tan(y(2))
      term1=c(layer,iblk,3)+c(layer,iblk,4)*y(1)
      term2=c(layer,iblk,1)*y(1)+c(layer,iblk,2)*y(1)*y(1)+
     +       term1*x+c(layer,iblk,5)
      vxv=(c(layer,iblk,8)*y(1)+c(layer,iblk,9)*y(1)*y(1)+
     +    c(layer,iblk,10)*x+c(layer,iblk,11))/
     +    (term2*(c(layer,iblk,6)*y(1)+c(layer,iblk,7)))
      vzv=term1/term2
      f(2)=vzv*f(1)-vxv
c
      if(.not.ok) then
        if(isrkc.eq.1.and.idump.eq.1) then
          write(12,5) vel(y(1),x)*(vxv**2+vzv**2)**.5
5         format('***  possible inaccuracies in rngkta w.r.t. z  ***'/
     +           '|del(v)| = ',f12.5)
        end if
        irkc=1
        isrkc=0
      end if
      return
      end
