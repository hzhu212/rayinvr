c
c     version 1.3  Aug 1992
c
c     Head wave and diffraction routines for TRAMP
c                 
c     ----------------------------------------------------------------
c                 
      subroutine hdwave(is,npt,ifam,nr,xshot,a1,nskip,itt,l,iszero,
     +    idot,iray,iamp,caust,spamp,nrskip,irayp,xsmax,isw,ampsmt,
     +    istep,anglew)
c                 
c     generate head waves at a spacing of hws along bottom of layer
c     and calculate corresponding amplitudes if iamp>0
c                 
      include 'tramp.par'
      real real,ls1,vave(ppray),amps(pnrayf),
     +     ls2,ls21,len,lenmin,lh,vpave(ppray),vsave(ppray),ab(ppray/2)
      complex q1,q2,z,zc,z1
      integer w,itt(1)
      include 'tramp.com'
c      
      nhdw=1      
      w=iwave     
      on2=-omega/2.
      tb=0.       
      nbnd=0      
      iturn=0     
      idh=id      
      fidh=fid    
      l=layer     
      ib=iblk     
      ilr=(id+3)/2
      ihamp=0     
      irhbnd=ircbnd
      ichbnd=iccbnd
      nchbnd=nccbnd+1
      vpk=vp(npt,1)
      vsk=vs(npt,1)
      vk=vr(npt,1)
c                 
c     make preliminary head wave amplitude calculations 
c                 
      if(iamp.gt.0) then 
        if(idump.eq.1) write(14,5) ifam,nr
5       format(/'head wave:   ifam=',i3,'  ir=',i5//
     +'irt iw1 iw2  vp1     vs1     vp2     vs2    angle     zc',
     +5x,'phase') 
        ihamp=1   
        vp1=vp(1,2)
        vs1=vs(1,2)
        if(w.eq.1) then 
          v1=vp1  
        else      
          v1=vs1  
        end if    
        ls1=sqrt(vk/v1)
        vpave(1)=(vp(1,2)+vp(2,1))/2.
        vsave(1)=(vs(1,2)+vs(2,1))/2.
        if(w.eq.1) then
          vave(1)=vpave(1)
        else      
          vave(1)=vsave(1) 
        end if    
        if(npt.gt.2) then
          do 60 i=2,npt-1
             vpave(i)=(vp(i,2)+vp(i+1,1))/2.
             vsave(i)=(vs(i,2)+vs(i+1,1))/2.
             if(w.eq.1) then
               vave(i)=vpave(i)
             else 
               vave(i)=vsave(i) 
             end if
             ab(i)=atan(sh(i-1))
60        continue
          do 10 i=2,npt-1
             az1=ar(i,1)+ab(i)
             az2=ar(i,2)+ab(i)
             ls1=ls1*cos(az1)/cos(az2)
10        continue
        end if    
c                 
        ls21=0.   
        do 20 i=1,npt-1
           ls21=ls21+tr(i)*vave(i)
20      continue  
c                 
        q1=cmplx(0.,1.)
        if(npt.gt.2) then
          do 30 i=2,npt-1
             az=abs(ar(i,1)+ab(i))
             if(az.gt.pi2) az=pi-az
             if(az.gt.1.5704) az=1.5704
             call zprtz(2,w,w,vpave(i-1),vpave(i),
     +                  vsave(i-1),vsave(i),az,zc)
             if(idump.eq.1) write(14,25) 2,w,w,vpave(i-1),
     +         vsave(i-1),vpave(i),vsave(i),az*pi18,cabs(zc),
     +         atan2(aimag(zc),real(zc))*pi18
25           format(i2,2i4,6f8.3,f8.2)
             q1=q1*cabs(zc)
30        continue
        end if    
c                 
        call block(xr(npt),l+1,ib1)
        vpbnd=(vm(l+1,ib1,2)-vm(l+1,ib1,1))*(xr(npt)-xbnd(l+1,ib1,1))/
     +        (xbnd(l+1,ib1,2)-xbnd(l+1,ib1,1))+vm(l+1,ib1,1)
        vsbnd=vpbnd*vsvp(l+1,ib1)
        if(w.eq.1) then
          vbnd=vpbnd
        else      
          vbnd=vsbnd
        end if    
        if((vk/vbnd).ge.1.) then
          tc=pi2
        else
          tc=asin(vk/vbnd)
        end if
        call zprtz(2,w,w,vpk,vpbnd,vsk,vsbnd,tc,zc)
        if(idump.eq.1) write(14,25) 2,w,w,vpk,vsk,vpbnd,vsbnd,tc*pi18,
     +                 cabs(zc),atan2(aimag(zc),real(zc))*pi18
        zl=.5*cabs(zc)*tan(tc)/omega
        z1=zl*q1/ls1
        len=0.    
        vhw=vbnd  
        if(idump.eq.1) write(14,1)
1       format(' ') 
      end if      
c                 
      call ttime(isw,xshot,npt,nr,a1,ifam,itt,iszero,iamp,2,idump,
     +           omega)
      qh=qa       
      if(iray.ne.0) call pltray(npt,nskip,idot,irayps,0,0.)
c                 
      t1=time     
      l1=layer+1  
      xr(1)=xr(npt)
      zr(1)=zr(npt)
      call block(xr(1),l1,ib1)
      vpbnd=(vm(l1,ib1,2)-vm(l1,ib1,1))*(xr(1)-xbnd(l1,ib1,1))/
     +      (xbnd(l1,ib1,2)-xbnd(l1,ib1,1))+vm(l1,ib1,1)
      vsbnd=vpbnd*vsvp(l1,ib1)
      if(w.eq.1) then
        vbnd=vpbnd
      else        
        vbnd=vsbnd
      end if      
      vr(1,1)=0.0 
1020  vp(1,2)=(vm(layer,iblk,4)-vm(layer,iblk,3))*(xr(1)-
     +        xbnd(layer,iblk,1))/(xbnd(layer,iblk,2)-
     +        xbnd(layer,iblk,1))+vm(layer,iblk,3) 
      if(vp(1,2).gt.0.) go to 1010
      layer=layer-1
      call block(xr(1),layer,iblk)
      go to 1020  
1010  vs(1,2)=vp(1,2)*vsvp(layer,iblk)
      if(w.eq.1) then
        vr(1,2)=vp(1,2)
      else        
        vr(1,2)=vs(1,2)
      end if      
      qr(1)=q(layer,iblk,(3-iwave)/2)
      iflag2=-2   
      abnd=atan(s(layer,iblk,2))
1000  if(ibsmth.eq.0) then
        alpha=atan(s(layer,iblk,2))
      else
        npl=ifix((xr(1)-xmin-.001)/xsinc)+1
        npr=npl+1
        xpl=xmin+float(npl-1)*xsinc
        alpha=(cosmth(layer+1,npr)-cosmth(layer+1,npl))*(xr(1)-xpl)
     +        /xsinc+cosmth(layer+1,npl)
      end if
      alpha2=atan(s(layer,iblk,2))
      time=t1+tb  
      nbnd=0      
      ircbnd=irhbnd
      iccbnd=ichbnd
      nccbnd=nchbnd
      iwave=w
      if(icbnd(iccbnd).eq.nccbnd) then
        iwave=-iwave
        iccbnd=iccbnd+1
        if(iwave.eq.1) then
          vr(1,2)=vp(1,2)
        else
          vr(1,2)=vs(1,2)
        end if
      end if
      if((vr(1,2)/vbnd).ge.1.) then
        a=pi2
      else
        a=pi-asin(vr(1,2)/vbnd)
      end if
      if(abs(abnd-alpha).le.angbnd) then
        idray(1)=l
        idray(2)=3
        ar(1,2)=fid*a-alpha
        npt=1     
        call trace(npt,ifam,nr,iturn,iamp,ihamp,iflag)
c                 
        qa=qh     
        call ttime(isw,xshot,npt,nr,a1,ifam,itt,iszero,iamp,
     +             iflag2,idump,omega)
        if((iray.eq.1.or.(iray.eq.2.and.vr(npt,2).gt.0.)).
     +    and.mod(nhdw-1,nrskip).eq.0) 
     +    call pltray(npt,0,idot,irayps,istep,anglew)
c                 
c       complete head wave amplitude calculations
c                 
        if(iamp.gt.0.and.vr(npt,2).ne.0.) then
          if(len.lt.001) then
            xhmin=abs(xshot-xr(npt))
            lenmin=hedcut*xhmin
          end if  
          if(len.lt.lenmin) then
            lh=lenmin
          else    
            lh=len
          end if  
c                 
          do 70 i=1,npt-1
             vpave(i)=(vp(i,2)+vp(i+1,1))/2.
             vsave(i)=(vs(i,2)+vs(i+1,1))/2. 
             if(w.eq.1) then
               vave(i)=vpave(i)
             else 
               vave(i)=vsave(i)
             end if 
70        continue
          ls2=0.  
          do 40 i=1,npt-1
             ls2=ls2+tr(i)*vave(i)
40        continue
          ls2=sqrt((lh*vhw+ls21)/v1)
c                 
          q2=cmplx(1.,0.)
          if(npt.gt.2) then
            do 50 i=2,npt-1
               az=abs(pi-fid1*(ar(i,1)+atan(sh(i-1))))
               if(az.gt.pi2) az=pi-az
               if(az.gt.1.5704) az=1.5704
               call zprtz(2,w,w,vpave(i-1),vpave(i),vsave(i-1),
     +                    vsave(i),az,zc)
               if(idump.eq.1) write(14,25) 2,w,w,vpave(i-1),
     +           vsave(i-1),vpave(i),vsave(i),az*pi18,cabs(zc),
     +           atan2(aimag(zc),real(zc))*pi18
               q2=q2*cabs(zc)
50          continue
          end if  
c                 
          call zprtz(-2,w,w,vpbnd,vp(1,2),vsbnd,vs(1,2),1.570796,zc)
          if(idump.eq.1) write(14,25) 2,w,w,vpbnd,vsbnd,vp(1,2),
     +      vs(1,2),90.,cabs(zc),atan2(aimag(zc),real(zc))*pi18
          q2=q2*zc
c                 
          if(surcon.eq.1) then
            if(aamp(nbnd,1).gt.1.5704) aamp(nbnd,1)=1.5704
            call zprtz(3,w,icomp,vp(npt,1),.001,vs(npt,1),.001,
     +                 aamp(nbnd,1),zc)
            if(idump.eq.1) write(14,25) 3,w,w,vp(npt,1),vs(npt,1),
     +        0.,0.,aamp(nbnd,1)*pi18,cabs(zc),atan2(aimag(zc),
     +        real(zc))*pi18
            q2=q2*zc
          end if  
c                 
          z=z1*(vhw/lh)**1.5*q2/(sqrt(vk)*ls2)
          if(surcon.eq.0) then
            if(icomp.eq.1) then
              if(w.eq.1) then
                z=z*abs(cos(ar(npt,1)))
              else 
                z=z*abs(sin(ar(npt,1))) 
              end if
            else  
              if(w.eq.1) then
                z=z*abs(sin(ar(npt,1)))
              else 
                z=z*abs(cos(ar(npt,1))) 
              end if
            end if
          end if  
          z=z*qa  
          if(icbnd(1).eq.0) then
            z=z*spamp
          end if  
          amp(ntt-1)=cabs(z)
          zimag=aimag(z)
          zreal=real(z)
          phase(ntt-1)=(caust+1)*pi2
          if(idump.eq.1) write(14,35) lh,vhw,qa
35        format(/'hdw bnd length=',f10.3,'  vhw=',f10.3,
     +            '  q fctr=',f10.4/)
        end if    
      end if      
      id=idh      
      fid=fidh    
      hwsl=hws    
500   xrp=xr(1)+fid*hwsl*cos(alpha2)
      if(xrp.gt.xmin.and.xrp.lt.xmax) then
        vh1=vbnd
        la=l1     
        if((fid*xrp).gt.(fid*xbnd(la,ib1,ilr))) then
          hl=abs(xbnd(la,ib1,ilr)-xr(1))/cos(alpha2)
          hwsl=hwsl-hl
          vpbnd=vm(la,ib1,ilr)
          vsbnd=vpbnd*vsvp(la,ib1)
          if(w.eq.1) then
            vbnd=vpbnd 
          else    
            vbnd=vsbnd
          end if  
          vh2=vbnd
          if(abs(vh2-vh1).lt..001) then
            tb=tb+hl/vh1
          else
            tb=tb+hl*alog(vh2/vh1)/(vh2-vh1)
          end if 
          if(iamp.gt.0) then
            vhw=(len*vhw+hl*(vh2+vh1)/2.)/(len+hl)
            len=len+hl
          end if  
          xr(1)=xbnd(la,ib1,ilr)
          ib1=ib1+id 
530       if(vm(la,ib1,1).eq.0.) then
            la=la+1
            if(la.gt.nlayer) go to 990
            call block(xr(1)+.001*fid,la,ib1)
            go to 530
          end if  
          alpha=atan(s(la,ib1,1))
          vbnd=vm(la,ib1,-(id-3)/2)
          if(iq.eq.1) qh=qh*exp(on2*hl/
     +                   (((vh2+vh1)/2.)*q(la,ib1,(3-w)/2)))
          go to 500
        else      
          xr(1)=xrp
          zr(1)=s(la,ib1,1)*xr(1)+b(la,ib1,1)
          if(xsmax.gt.0..and.abs(xshot-xr(1)).gt.xsmax) go to 990
          vr(1,1)=0.
          ls=l    
520       call block(xr(1),ls,ib)
          if(vm(ls,ib,2).le.0.) then 
            ls=ls-1
            if(ls.le.0) go to 990
            go to 520
          end if  
510       vp(1,2)=(vm(ls,ib,4)-vm(ls,ib,3))*(xr(1)-xbnd(ls,ib,1))/
     +            (xbnd(ls,ib,2)-xbnd(ls,ib,1))+vm(ls,ib,3)
          vs(1,2)=vp(1,2)*vsvp(ls,ib)
          if(w.eq.1) then
            vr(1,2)=vp(1,2)
          else    
            vr(1,2)=vs(1,2) 
          end if  
          la=l1   
540       call block(xr(1),la,ib1)
          if(vm(la,ib1,1).eq.0.) then
            la=la+1
            if(la.gt.nlayer) go to 990
            go to 540
          end if  
          vpbnd=(vm(la,ib1,2)-vm(la,ib1,1))*(xr(1)-xbnd(la,ib1,1))/
     +          (xbnd(la,ib1,2)-xbnd(la,ib1,1))+vm(la,ib1,1)
          vsbnd=vpbnd*vsvp(la,ib1)
          if(w.eq.1) then
            vbnd=vpbnd
          else    
            vbnd=vsbnd 
          end if  
          if(vbnd.le..001) then
            write(11,45)
45          format('***  ray stopped - s-wave cannot propagate  ***')
            go to 990
          end if  
          vh2=vbnd
          if(abs(vh2-vh1).lt..001) then
            tb=tb+hwsl/vh1
          else
            tb=tb+hwsl*alog(vh2/vh1)/(vh2-vh1)
          end if 
          if(iamp.gt.0) then
            vhw=(len*vhw+hwsl*(vh2+vh1)/2.)/(len+hwsl)
            len=len+hwsl
          end if  
          vratio=vr(1,2)/vbnd
          if(vratio.ge.1.) then
            a=pi2 
          else    
            a=pi-asin(vratio)
          end if  
          if(iq.eq.1) qh=qh*exp(on2*hwsl/
     +                   (((vh2+vh1)/2.)*q(la,ib1,(3-w)/2)))
          layer=ls
          iblk=ib 
          iflag2=-1
          nhdw=nhdw+1
          if(nhdw.gt.pnrayf) go to 990
          go to 1000
        end if    
      end if      
c
990   if(iamp.gt.0.and.itt(ifam).gt.0) then
        if(ampsmt.gt.0.) then
          nsmth=nint(ampsmt*(itt(ifam)-2))
          if(nsmth.gt.0) then
            do 9110 i=1,itt(ifam)
               if(amp(ntt-itt(ifam)-1+i).gt.0.) then
                 amps(i)=amp(ntt-itt(ifam)-1+i)
               else
                 amps(i)=0.
               end if
9110        continue
            do 9120 i=1,nsmth
               call smooth(amps,itt(ifam))
9120        continue
            do 9130 i=1,itt(ifam)
               amp(ntt-itt(ifam)-1+i)=amps(i)
9130        continue 
          end if
        end if
c
        do 9140 i=1,itt(ifam)
           ipos=ntt-itt(ifam)-1+i
           write(13,15) isw,nr,range(ipos),amp(ipos),phase(ipos)*pi18,
     +                  rayid(ipos)
15         format(2i5,f10.3,e12.3,f10.2,f6.1)
9140    continue
      end if
c
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine difrct(is,npt,ifam,nr,xshot,a1,nskip,itt,iszero,
     +           idot,iray,nrskip,irayps,xsmax,isw,istep,anglew)
c                 
c     generate damped head waves at a spacing of hws along bottom of
c     layer to simulate diffracted arrivals associated with the top
c     of a low-velocity zone - amplitudes cannot be calculated
c                 
      include 'tramp.par'
      integer w,itt(1)
      include 'tramp.com'
c      
      nhdw=1      
      w=iwave     
      difbnd=0    
      tb=0.       
      nbnd=0      
      iturn=0     
      idh=id      
      fidh=fid    
      l=layer     
      ib=iblk     
      ilr=(id+3)/2
      ihamp=0     
      irhbnd=ircbnd
      ichbnd=iccbnd
c                 
      call ttime(isw,xshot,npt,nr,a1,ifam,itt,iszero,0,2,idump,
     +           omega)
      if(iray.ne.0) call pltray(npt,nskip,idot,irayps,0,0.)
c                 
      t1=time     
      l1=layer+1  
      xr(1)=xr(npt)
      zr(1)=zr(npt)
      call block(xr(1),l1,ib1)
      vpbnd=(vm(l1,ib1,2)-vm(l1,ib1,1))*(xr(1)-xbnd(l1,ib1,1))/
     +      (xbnd(l1,ib1,2)-xbnd(l1,ib1,1))+vm(l1,ib1,1)
      vsbnd=vpbnd*vsvp(l1,ib1)
      if(w.eq.1) then
        vbnd=vpbnd
      else        
        vbnd=vsbnd
      end if      
      abnd=atan(s(layer,iblk,2))
      a=pi-abs(ar(npt,1)+abnd)
      vr(1,1)=0.0 
1020  vp(1,2)=(vm(layer,iblk,4)-vm(layer,iblk,3))*(xr(1)-
     +        xbnd(layer,iblk,1))/(xbnd(layer,iblk,2)-
     +        xbnd(layer,iblk,1))+vm(layer,iblk,3)
      if(vp(1,2).gt.0.) go to 1010
      layer=layer-1
      call block(xr(1),layer,iblk)
      go to 1020  
1010  vs(1,2)=vp(1,2)*vsvp(layer,iblk)
      if(w.eq.1) then
        vr(1,2)=vp(1,2)
      else        
        vr(1,2)=vs(1,2)
      end if      
      iflag2=-2   
1000  alpha=atan(s(layer,iblk,2))
      time=t1+tb  
      nbnd=0      
      ircbnd=irhbnd
      iccbnd=ichbnd
      if(abs(abnd-alpha).le.angbnd) then
        idray(1)=l
        idray(2)=3
        ar(1,2)=fid*a-alpha
        npt=1     
        call trace(npt,ifam,nr,iturn,iamp,ihamp,iflag)
c                 
        call ttime(isw,xshot,npt,nr,a1,ifam,itt,iszero,0,iflag2,
     +             idump,omega)
        if((iray.eq.1.or.(iray.eq.2.and.vr(npt,2).gt.0.)).
     +    and.mod(nhdw-1,nrskip).eq.0) 
     +    call pltray(npt,0,idot,irayps,istep,anglew)
c                 
        amp(ntt-1)=-1.
      end if      
      id=idh      
      fid=fidh    
      hwsl=hws    
500   xrp=xr(1)+fid*hwsl*cos(alpha2)
      if(xrp.gt.xmin.and.xrp.lt.xmax) then
        vh1=vbnd
        la=l1     
        if((fid*xrp).gt.(fid*xbnd(la,ib1,ilr))) then
          hl=abs(xbnd(la,ib1,ilr)-xr(1))/cos(alpha2)
          hwsl=hwsl-hl
          vpbnd=vm(la,ib1,ilr)
          vsbnd=vpbnd*vsvp(la,ib1)
          if(w.eq.1) then 
            vbnd=vpbnd
          else    
            vbnd=vsbnd
          end if  
          vh2=vbnd
          if(abs(vh2-vh1).lt..001) then
            tb=tb+hl/vh1
          else
            tb=tb+hl*alog(vh2/vh1)/(vh2-vh1)
          end if 
          xr(1)=xbnd(la,ib1,ilr)
          ib1=ib1+id 
530       if(vm(la,ib1,1).eq.0.) then
            la=la+1
            if(la.gt.nlayer) return
            call block(xr(1),la,ib1)
            go to 530
          end if  
          alpha=atan(s(la,ib1,1))
          vbnd=vm(la,ib1,-(id-3)/2)
          go to 500
        else      
          xr(1)=xrp
          zr(1)=s(la,ib1,1)*xr(1)+b(la,ib1,1)
          if(xsmax.gt.0..and.abs(xshot-xr(1)).gt.xsmax) return
          vr(1,1)=0.
          ls=l    
520       call block(xr(1),ls,ib)
          if(vm(ls,ib,2).le.0.) then
            ls=ls-1 
            if(ls.le.0) return
            go to 520
          end if  
510       vp(1,2)=(vm(ls,ib,4)-vm(ls,ib,3))*(xr(1)-xbnd(ls,ib,1))/
     +            (xbnd(ls,ib,2)-xbnd(ls,ib,1))+vm(ls,ib,3)
          vs(1,2)=vp(1,2)*vsvp(ls,ib) 
          if(w.eq.1) then
            vr(1,2)=vp(1,2) 
          else    
            vr(1,2)=vs(1,2) 
          end if  
          la=l1   
540       call block(xr(1),la,ib1)
          if(vm(la,ib1,1).le.0.) then 
            la=la+1 
            if(la.gt.nlayer) return
            go to 540
          end if  
          vpbnd=(vm(la,ib1,2)-vm(la,ib1,1))*(xr(1)-xbnd(la,ib1,1))/
     +          (xbnd(la,ib1,2)-xbnd(la,ib1,1))+vm(la,ib1,1)
          vsbnd=vpbnd*vsvp(la,ib1)
          if(w.eq.1) then 
            vbnd=vpbnd
          else    
            vbnd=vsbnd
          end if  
          if(vbnd.le..001) then
            write(11,45)
45          format('***  ray stopped - s-wave cannot propagate  ***')
            return
          end if  
          vh2=vbnd
          if(abs(vh2-vh1).lt..001) then
            tb=tb+hwsl/vh1
          else
            tb=tb+hwsl*alog(vh2/vh1)/(vh2-vh1)
          end if 
          layer=ls
          iblk=ib 
          iflag2=-1
          nhdw=nhdw+1
          if(nhdw.gt.pnrayf) return
          go to 1000
        end if    
      end if      
      return      
      end         
