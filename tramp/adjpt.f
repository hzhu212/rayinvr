c
c     version 1.3  Aug 1992
c
c     Adjpt routine for TRAMP
c                 
c     ----------------------------------------------------------------
c                 
      subroutine adjpt(n,top,bot,lef,rig,ifam,ir,iturn,iamp,iflag,
     +                 idl,idr)
c                 
c     ray has intersected a model boundary so must determine correct
c     (x,z) coordinates and angle at boundary 
c                 
      include 'tramp.par'
      real lef
      integer icasec(4)
      include 'tramp.com'
c      
      data icasec/3,4,1,2/
c           
      iconvt=0    
      if(nbnd.eq.piray) then
        xr(n)=xr(n-1) 
        zr(n)=zr(n-1)
        ar(n,1)=ar(n-1,2)
        vr(n,2)=0.0
        iflag=1   
        if(ir.ne.0) write(11,5)
5     format('***  ray stopped - intersected too many boundaries ***')
        return    
      end if      
c                 
c     make appropriate adjustment if point outside current model 
c     block is above top boundary and below bottom boundary
c                 
100   if(zr(n).le.top.and.zr(n).gt.bot) then
        xmmm=(xr(n)+xr(n-1))/2.
        zmmm=(zr(n)+zr(n-1))/2.
        am=(ar(n,1)+ar(n-1,2))/2.
500     top=s(layer,iblk,1)*xmmm+b(layer,iblk,1)
        bot=s(layer,iblk,2)*xmmm+b(layer,iblk,2)
        if(zmmm.le.top.or.zmmm.gt.bot.or.xmmm.le.lef.or.xmmm.gt.rig) 
     +  then
          xr(n)=xmmm 
          zr(n)=zmmm 
          ar(n,1)=am
          go to 100
        else      
          xmmm=(xr(n)+xmmm)/2.
          zmmm=(zr(n)+zmmm)/2.
          am=(ar(n,1)+am)/2.
          go to 500
        end if    
      end if      
c                 
      if(zr(n).le.top.and.xr(n).gt.lef.and.xr(n).le.rig) then
c                 
c      ray intersects upper boundary
c                 
       icase=1    
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            s(layer,iblk,1),b(layer,iblk,1),0,xn,zn,an)
       go to 1001
      end if      
      if(zr(n).le.top.and.xr(n).gt.rig) then
c                 
c      ray intersects upper right corner
c                 
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            s(layer,iblk,1),b(layer,iblk,1),0,xt,zt,at)
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            0.,rig,1,xri,zri,ari)
       call corner(xr(n-1),zr(n-1),xt,zt,at,1,xri,zri,ari,2,
     +             xn,zn,an,icase)
       if(icase.eq.2.and.iblk.lt.nblk(layer)) then
         if(ivg(layer,iblk+1).eq.-1) icase=1
       end if
       go to 1001
      end if      
      if(zr(n).gt.bot.and.xr(n).gt.lef.and.xr(n).le.rig) then
c                 
c      ray intersects lower boundary
c                 
       icase=3    
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            s(layer,iblk,2),b(layer,iblk,2),0,xn,zn,an)
       go to 1001
      end if      
      if(zr(n).gt.bot.and.xr(n).gt.rig) then
c                 
c      ray intersects lower right corner
c                 
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            s(layer,iblk,2),b(layer,iblk,2),0,xb,zb,ab)
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            0.,rig,1,xri,zri,ari)
       call corner(xr(n-1),zr(n-1),xb,zb,ab,3,xri,zri,ari,2,
     +             xn,zn,an,icase)
       if(icase.eq.2.and.iblk.lt.nblk(layer)) then
         if(ivg(layer,iblk+1).eq.-1) icase=3
       end if
       go to 1001
      end if      
      if(xr(n).gt.rig.and.zr(n).gt.top.and.zr(n).le.bot) then
c                 
c      ray intersects right boundary 
c                 
       icase=2    
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            0.,rig,1,xn,zn,an)
       go to 1001
      end if      
      if(zr(n).le.top.and.xr(n).le.lef) then
c                 
c      ray intersects upper left corner
c                 
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            s(layer,iblk,1),b(layer,iblk,1),0,xt,zt,at)
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            0.,lef,1,xl,zl,al)
       call corner(xr(n-1),zr(n-1),xt,zt,at,1,xl,zl,al,4,xn,zn,an,icase)
       if(icase.eq.4.and.iblk.gt.1) then
         if(ivg(layer,iblk-1).eq.-1) icase=1
       end if
       go to 1001
      end if      
      if(xr(n).le.lef.and.zr(n).gt.top.and.zr(n).le.bot) then
c                 
c      ray intersects left boundary
c                 
       icase=4    
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            0.,lef,1,xn,zn,an)
       go to 1001
      end if      
      if(zr(n).gt.bot.and.xr(n).le.lef) then
c                 
c      ray intersects lower left corner
c                 
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            s(layer,iblk,2),b(layer,iblk,2),0,xb,zb,ab)
       call adhoc(xr(n-1),zr(n-1),ar(n-1,2),xr(n),zr(n),ar(n,1),
     +            0.,lef,1,xl,zl,al)
       call corner(xr(n-1),zr(n-1),xb,zb,ab,3,xl,zl,al,4,xn,zn,an,icase)
       if(icase.eq.4.and.iblk.gt.1) then
         if(ivg(layer,iblk-1).eq.-1) icase=3
       end if
       go to 1001
      end if      
c                 
c     check to see if ray has crossed the same boundary in two
c     successive points
c                 
1001  continue
      if(nbnd.gt.0) then
        if(nptamp(nbnd).eq.n-1) then
          if(icasec(icase).eq.icasel) then
            if(n2.gt.nstepr) then
              n2=0 
              icasel=icase
            else  
              n2=n2+1
              n=n-1
              dstepf=dstepf*2.
              return
            end if
          else    
            n2=0  
            icasel=icase
          end if  
        else      
          n2=0    
          icasel=icase
        end if    
      else        
        if(n.eq.2.and.layer.eq.1.and.icase.eq.1) then
          if(n2.gt.nstepr) then 
            n2=0 
            icasel=icase 
          else    
            n2=n2+1
            n=1   
            dstepf=dstepf*2.
            return
          end if  
        else      
          n2=0    
          icasel=icase
        end if    
      end if      
      nbnd=nbnd+1 
      nptamp(nbnd)=n
c                 
      ar(n,1)=an   
      if((fid*ar(n,1)).lt.0.) then
        id=-id    
        fid=float(id)
      end if      
      if((fid*ar(n,1)).gt.pi) then
        ar(n,1)=fid*pit2+ar(n,1)
        id=-id    
        fid=float(id)
      end if      
      xr(n)=xn     
      zr(n)=zn     
      vp(n,1)=vel(xn,zn)
      vs(n,1)=vp(n,1)*vsvp(layer,iblk)
      if(iwave.eq.1) then
        vr(n,1)=vp(n,1) 
      else        
        vr(n,1)=vs(n,1)
      end if      
      vpamp(nbnd,1)=vp(n,1)
      vsamp(nbnd,1)=vs(n,1) 
      irtamp(nbnd,1)=iwave
c                 
      go to (1000,2000,3000,4000), icase
c                 
c     ray intersects upper boundary
c                 
1000  lstart=layer    
      istart=iblk
      if(ibsmth.eq.0) then
        alpha=atan(s(layer,iblk,1))
      else        
        npl=ifix((xr(n)-xmin-.001)/xsinc)+1
        npr=npl+1 
        xpl=xmin+float(npl-1)*xsinc
        alpha=(cosmth(layer,npr)-cosmth(layer,npl))*(xr(n)-xpl)
     +        /xsinc+cosmth(layer,npl)
      end if      
      aamp(nbnd,1)=pi-fid*(ar(n,1)+alpha)
      if(abs(aamp(nbnd,1)).ge.pi2.and.n3.le.nstepr.
     +and.dstepf.lt.1.e6) go to 999
      n3=0        
      dstepf=1.   
c                 
c     check to see if ray is at top of model
c                 
      if(layer.gt.1) go to 1300
      if(refll(ircbnd).ne.-1) then
        vr(n,2)=1.0
        iflag=1   
        return    
      else        
        vpamp(nbnd,2)=.001
        vsamp(nbnd,2)=.001
        vo=.001   
c       vpamp(nbnd,2)=3.3
c       vsamp(nbnd,2)=1.8 
c       vo=3.3   
      end if      
1300  nccbnd=nccbnd+1
      if(icbnd(iccbnd).eq.nccbnd) then
        iconvt=1  
        iwave=-iwave
        iccbnd=iccbnd+1
      end if      
      irtamp(nbnd,2)=iwave
      vi=vel(xr(n),zr(n))
      if(iwave.eq.-1) vi=vi*vsvp(layer,iblk)
      if(layer.gt.1) then
1320    i1=layer-1    
        call block(xn,i1,ib)
        vpamp(nbnd,2)=(vm(i1,ib,4)-vm(i1,ib,3))*(xn-xbnd(i1,ib,1))/
     +                (xbnd(i1,ib,2)-xbnd(i1,ib,1))+vm(i1,ib,3)
        if(ivg(i1,ib).ne.-1) go to 1310
        layer=layer-1     
        if(refll(ircbnd).eq.-layer) then
          if(ir.ne.0) write(11,55) 
55        format( 
     +    '***  ray stopped - reflected from pinchout  ***')
          vr(n,2)=0.0 
          iflag=1 
          return  
        end if    
        if(layer.gt.1) go to 1320
        if(refll(ircbnd).ne.-1) then
          vr(n,2)=1.0
          iflag=1 
          return  
        else      
          vpamp(nbnd,2)=.001
          vsamp(nbnd,2)=.001
          vo=.001 
          go to 1330
        end if    
1310    vsamp(nbnd,2)=vpamp(nbnd,2)*vsvp(i1,ib)
        if(iwave.eq.1) then
           vo=vpamp(nbnd,2)
        else      
           vo=vsamp(nbnd,2)
        end if    
      end if      
1330  if(refll(ircbnd).eq.-lstart) then
        layer=lstart  
        ircbnd=ircbnd+1
        if(iconvt.eq.0) then
          aamp(nbnd,2)=aamp(nbnd,1)
          go to 1110
        end if    
        at1=(vi/vr(n,1))*sin(aamp(nbnd,1))
        if(abs(at1).lt.1.) go to 1100
        if(ir.ne.0) write(11,35)
35      format(   
     + '***  ray stopped - converted ray cannot reflect/refract  ***')
        vr(n,2)=0.0
        iflag=1   
        return    
      end if      
      at2=(vo/vr(n,1))*sin(aamp(nbnd,1))
      if(abs(at2).lt.1.) go to 1200
c                 
c     ray reflects off upper boundary
c                 
      if(istop.gt.0.or.ir.eq.0) then
        vr(n,2)=0.0
        iflag=1   
        if(ir.ne.0) write(11,15)
15      format('***  ray stopped - illegal reflection  ***')
        return    
      end if      
      if(iconvt.eq.0) then
        aamp(nbnd,2)=aamp(nbnd,1)
        go to 1110
      end if      
      at1=(vi/vr(n,1))*sin(aamp(nbnd,1))
      if(abs(at1).ge.1.) then
        if(ir.ne.0) write(11,35)
        vr(n,2)=0.0
        iflag=1   
        return    
      end if      
1100  aamp(nbnd,2)=asin(at1)
1110  irtamp(nbnd,3)=1
      ar(n,2)=fid*aamp(nbnd,2)-alpha
      vr(n,2)=vi   
      vp(n,2)=vpamp(nbnd,1)
      vs(n,2)=vsamp(nbnd,1)
      if((fid*ar(n,2)).lt.0.) then
        id=-id    
        fid=float(id)
      end if      
      if(vm(layer,iblk,1).eq.0.) then
        layer=lstart
        iblk=istart
      end if
      return      
c                 
c       ray refracts through upper boundary
c                 
1200  layer=layer-1       
      iblk=ib        
      irtamp(nbnd,3)=2
      aamp(nbnd,2)=asin(at2)
      ar(n,2)=fid*(pi-aamp(nbnd,2))-alpha
      vr(n,2)=vo   
      vp(n,2)=vpamp(nbnd,2)
      vs(n,2)=vsamp(nbnd,2)
      if(iwave.eq.-1.and.vr(n,2).le..001) then
        vr(n,2)=0.0 
        iflag=1   
        if(ir.ne.0) write(11,45)
45      format('***  ray stopped - s-wave cannot propagate  ***')
        return    
      end if      
      if((fid*ar(n,2)).gt.pi) then
        ar(n,2)=fid*pit2+ar(n,2)
        id=-id    
        fid=float(id)
      end if      
      return      
c                 
c     ray intersects lower boundary
c                 
3000  lstart=layer    
      istart=iblk
      if(layer.eq.nlayer) then
        vr(n,2)=0.0
        iflag=1   
        return    
      end if      
      if(ibsmth.eq.0) then
        alpha=atan(s(layer,iblk,2))
      else        
        npl=ifix((xr(n)-xmin-.001)/xsinc)+1
        npr=npl+1 
        xpl=xmin+float(npl-1)*xsinc
        alpha=(cosmth(layer+1,npr)-cosmth(layer+1,npl))*(xr(n)-xpl)
     +        /xsinc+cosmth(layer+1,npl)
      end if      
      aamp(nbnd,1)=fid*(ar(n,1)+alpha)
      if(abs(aamp(nbnd,1)).ge.pi2.and.n3.le.nstepr.
     +and.dstepf.lt.1.e6) go to 999
      n3=0        
      dstepf=1.   
      if(difbnd.eq.layer) then
        ar(n,2)=fid*(pi-aamp(nbnd,1))-alpha
        vr(n,2)=vr(n,1)
        iflag=3   
        return    
      end if      
      nccbnd=nccbnd+1
      if(icbnd(iccbnd).eq.nccbnd) then
        iconvt=1  
        iwave=-iwave 
        iccbnd=iccbnd+1
      end if      
      irtamp(nbnd,2)=iwave
      vi=vel(xr(n),zr(n))
      if(iwave.eq.-1) vi=vi*vsvp(layer,iblk)
3020  i1=layer+1      
      call block(xn,i1,ib)
      vpamp(nbnd,2)=(vm(i1,ib,2)-vm(i1,ib,1))*(xn-xbnd(i1,ib,1))/
     +              (xbnd(i1,ib,2)-xbnd(i1,ib,1))+vm(i1,ib,1)
      if(ivg(i1,ib).ne.-1) go to 3010
      layer=layer+1       
      if(refll(ircbnd).eq.layer) then
        if(ir.ne.0) then
          write(11,55) 
          vr(n,2)=0.0 
        end if    
        if(idray(1).lt.layer) idray(1)=layer
        idray(2)=2
        iflag=1   
        return    
      end if      
      if(layer.eq.nlayer) then
        vr(n,2)=0.0
        iflag=1   
        return    
      end if      
      if(difbnd.eq.layer) then
        ar(n,2)=fid*(pi-aamp(nbnd,1))-alpha
        vr(n,2)=vr(n,1)
        iflag=3   
        return    
      end if      
      go to 3020  
3010  vsamp(nbnd,2)=vpamp(nbnd,2)*vsvp(i1,ib)
      if(iwave.eq.1) then
        vo=vpamp(nbnd,2)
      else        
        vo=vsamp(nbnd,2) 
      end if      
      if(refll(ircbnd).eq.lstart) then
        layer=lstart  
        ircbnd=ircbnd+1
        if(iconvt.eq.0) then
          aamp(nbnd,2)=aamp(nbnd,1)
          go to 3110 
        end if    
        at1=(vi/vr(n,1))*sin(aamp(nbnd,1))
        if(abs(at1).lt.1.) go to 3100
        if(ir.ne.0) write(11,35)
        vr(n,2)=0.0
        iflag=1   
        return    
      end if      
      if(iheadf(layer).eq.1.and.vr(n,1).lt.vo) then
        crita=asin(vr(n,1)/vo)
        diff=abs(aamp(nbnd,1)-crita)
        if(diff.le.crit) go to 3300
      end if      
      at2=(vo/vr(n,1))*sin(aamp(nbnd,1))
      if(abs(at2).lt.1.) go to 3200
c                 
c     ray reflects off lower boundary
c                 
      if(istop.gt.0.or.ir.eq.0) then
        vr(n,2)=0.0
        iflag=1   
        if(ir.ne.0) write(11,15)
        idray(2)=2
        return    
      end if      
      if(iconvt.eq.0) then
        aamp(nbnd,2)=aamp(nbnd,1)
        go to 3110
      end if      
      at1=(vi/vr(n,1))*sin(aamp(nbnd,1))
      if(abs(at1).ge.1.) then
        if(ir.ne.0) write(11,35)
        vr(n,2)=0.0
        iflag=1   
        return    
      end if      
3100  aamp(nbnd,2)=asin(at1)
3110  irtamp(nbnd,3)=1
      ar(n,2)=fid*(pi-aamp(nbnd,2))-alpha
      vr(n,2)=vi   
      vp(n,2)=vpamp(nbnd,1)
      vs(n,2)=vsamp(nbnd,1)
      if((fid*ar(n,2)).lt.0.) then
        id=-id    
        fid=float(id)
      end if      
      if((fid*ar(n,2)).gt.pi) then
        ar(n,2)=fid*pit2+ar(n,2)
        id=-id    
        fid=float(id)
      end if      
      if(ir.ne.0.or.idr.ne.1) idray(2)=2  
      if(iturn.eq.1.and.idr.eq.2) iflag=1
      if(vm(layer,iblk,1).eq.0.) then
        layer=lstart
        iblk=istart
      end if
      return      
c                 
c     ray refracts through lower boundary
c                 
3200  if(istop.gt.0.and.idray(2).eq.3) then
        vr(n,2)=0.
        iflag=1
        if(ir.ne.0) write(11,115)
115     format('***  ray stopped - illegal headwave refraction  ***')
        return
      end if
      if(ir.ne.0.and.istop.eq.2.and.(layer+1).gt.idl) then
        vr(n,2)=0.
        iflag=1
        if(ir.ne.0) write(11,125)
125     format('***  ray stopped - illegal refraction (istop=2)  ***')
        return
      end if
      layer=layer+1       
      iblk=ib        
      irtamp(nbnd,3)=2
      aamp(nbnd,2)=asin(at2)
      ar(n,2)=fid*aamp(nbnd,2)-alpha
      vr(n,2)=vo   
      vp(n,2)=vpamp(nbnd,2)
      vs(n,2)=vsamp(nbnd,2)
      if(iwave.eq.-1.and.vr(n,2).le..001) then
        vr(n,2)=0.0 
        iflag=1   
        write(11,45)
        return    
      end if      
      if((fid*ar(n,2)).lt.0.) then
        id=-id    
        fid=float(id)
      end if      
      if(idray(1).lt.layer) idray(1)=layer
      return      
c                 
c     head waves to be generated along lower boundary
c                 
3300  ar(n,2)=fid*(pi-crita)-alpha
      vr(n,2)=vo   
      vp(n,2)=vpamp(nbnd,2)
      vs(n,2)=vsamp(nbnd,2)
      iflag=2     
      idray(1)=layer  
      return      
c                 
c     ray intersects right boundary 
c                 
2000  jt=nblk(layer)  
      ja=1        
      go to 5000  
c                 
c       ray intersects left boundary
c                 
4000  jt=1        
      ja=-1       
c                 
5000  if(iblk.eq.jt) then
        vr(n,2)=0.0
        iflag=1   
        return    
      end if      
      irtamp(nbnd,2)=iwave
      iblk=iblk+ja      
      vp(n,2)=vel(xn,zn)
      vs(n,2)=vp(n,2)*vsvp(layer,iblk)
      if(iwave.eq.1) then
        vr(n,2)=vp(n,2)
      else        
        vr(n,2)=vs(n,2)
      end if      
      vo=vr(n,2)   
      vpamp(nbnd,2)=vp(n,2)
      vsamp(nbnd,2)=vs(n,2) 
      aamp(nbnd,1)=pi2-fid*ar(n,1)
      if(abs(aamp(nbnd,1)).ge.pi2.and.n3.le.nstepr.
     +and.dstepf.lt.1.e6) go to 999
      n3=0        
      dstepf=1.   
      if((fid*ar(n,1)).gt.pi2) aamp(nbnd,1)=-aamp(nbnd,1)
      at=(vo/vr(n,1))*sin(aamp(nbnd,1))
      if(abs(at).ge.1.) then
c                 
c       ray reflects off vertical boundary
c                 
        if(istop.gt.0.or.ir.eq.0) then
          vr(n,2)=0.0
          iflag=1 
          if(ir.ne.0) write(11,15)
          return  
        end if    
        irtamp(nbnd,3)=1
        aamp(nbnd,2)=aamp(nbnd,1)
        iblk=iblk-ja    
        ar(n,2)=-ar(n,1)
        vr(n,2)=vr(n,1)
        vp(n,2)=vp(n,1)
        vs(n,2)=vs(n,1)
        id=-id    
        fid=float(id)
      else        
c                 
c       ray refracts through vertical boundary
c                 
        irtamp(nbnd,3)=2
        aamp(nbnd,2)=asin(at)
        if((fid*ar(n,1)).le.pi2) then
          ar(n,2)=fid*(pi2-aamp(nbnd,2))
        else      
          ar(n,2)=fid*(pi2+aamp(nbnd,2))
        end if    
        if(iwave.eq.-1.and.vr(n,2).le..001) then
          vr(n,2)=0.0 
          iflag=1 
          write(11,45)
          return  
        end if    
      end if      
      return      
c                 
c     ray has intersected a boundary with an angle of incidence
c     beyond 90 degrees
c                 
999   n3=n3+1    
      n=n-1      
      nbnd=nbnd-1
      dstepf=dstepf*2.
c
      return      
      end         
