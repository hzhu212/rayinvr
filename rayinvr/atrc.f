c                 
c     version 1.3  Aug 1992
c
c     Take-off angle search mode routines for RAYINVR
c
c     ----------------------------------------------------------------
c                 
      subroutine autotr(ang,layer1,iblk1,xshot,zshot,ifam,iturn,
     +           npt,iflag2,irays,nskip,idot,idr,irayps,istep)
c                 
c     trace a single ray through model to its turning point
c     (iturn=1) or to completion (iturn=0)
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c
      ircbnd=1    
      iccbnd=1    
      iwave=1     
      if(icbnd(1).eq.0) then
        iwave=-iwave
        iccbnd=2  
      end if      
      nptbnd=0
      nbnd=0
      npskp=npskip
      nccbnd=0     
      id=nint(fid1)
      fid=fid1    
      ir=0        
      invr=0
      angle=fid*(90.-ang)/pi18
      layer=layer1 
      iblk=iblk1  
      npt=1       
      xr(1)=xshot 
      zr(1)=zshot 
      ar(1,1)=0.0 
      ar(1,2)=angle
      vr(1,1)=0.  
      vp(1,1)=0.0 
      vs(1,1)=0.0 
      vp(1,2)=vel(xshot,zshot)
      vs(1,2)=vp(1,2)*vsvp(layer1,iblk1)
      if(iwave.eq.1) then
        vr(1,2)=vp(1,2) 
      else        
        vr(1,2)=vs(1,2)
      end if      
      idray(1)=layer1
      idray(2)=1  
      ifcbnd=0
c
      call trace(npt,ifam,ir,iturn,invr,0.,iflag2,0,idr,0,0,0,0)
c
      if(irays.eq.1) call pltray(npt,nskip,idot,irayps,istep,ang)
c
      if(idump.eq.1) write(12,5) ifam,ir,npt,xr(npt),zr(npt),
     +  ar(npt,1)*pi18,ar(npt,2)*pi18,vr(npt,1),vr(npt,2),
     +  layer,iblk,id,iwave
5     format(i2,i3,i4,2f8.3,2f8.2,2f7.2,4i3)
c
      return      
      end         
                 
c     ----------------------------------------------------------------
c                 
      subroutine auto(xshot,zshot,if,ifam,l,idr,amin,amax,aamin,aamax,
     +       layer1,iblk1,aainc,aaimin,nsmax,iflag,it,amina,amaxa,
     +       ia0,stol,irays,nskip,idot,irayps,xsmax,istep,nsmin)
c                 
c     determine min. and max. take-off angles for a specific ray code
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c      
      iflag=0     
      if(idr.eq.0) then
c                 
c       ray take-off angles input by user
c                 
        if(abs(amina).gt.180..or.abs(amaxa).gt.180.) then
          write(11,5)
5         format(/'***  error in specification of amin or amax  ***'/)
          iflag=1 
          return  
        end if    
        amin=amina
        amax=amaxa
        ia0=ia0+1 
        return    
      end if      
c                 
c     check to see if layer of ray code is above shot or below model
c                 
      if(l.lt.layer1.or.l.gt.nlayer.or.(l.ge.nlayer.and.idr.gt.1))
     +then        
        iflag=1   
        return    
      end if      
c
      if(idr.eq.1) then
c                 
c       search for refracted rays
c                 
        call calvel(xshot,zshot,layer1,iblk1,l,vshot,vtop,vbotom)
c                 
        if(l.gt.layer1.and.abs(vtop-vbotom).le..000001) then
          call block(xshot+fid*.0005,l,ib)
          if(id.eq.1) then
            i1=ib
            i2=nblk(l)
          else    
            i1=1  
            i2=ib
          end if  
          do 103 i=i1,i2
            if((vm(l,i,1).le.vm(l,i,3).or.vm(l,i,2).le.vm(l,i,4)).
     +      and.vm(l,i,1).ne.0.) go to 105
103       continue
          iflag=1 
          return  
105       vratio=vshot/vtop
          if(vratio.gt..99999) then
            ang=aamin
          else    
            ang=90.-asin(vratio)*pi18
          end if  
          ainc=1. 
          a1=0.   
          amin=-999999.
          amax=-999999.
          if(isrch.eq.0.and.tang(l,1).lt.100.) then
            amin=tang(l,1)
          else    
            do 101 n=1,nsmax
               call autotr(ang,layer1,iblk1,xshot,zshot,ifam,it,npt,
     +              iflag2,irays,nskip,idot,idr,irayps,istep)
               if(idray(1).eq.l.and.idray(2).eq.1.and.vr(npt,2).ne.
     +           0.) then
                 amin=ang
                 tang(l,1)=amin 
                 if(n.ge.nsmin) go to 104
                 if(stol.gt.0..and.n.gt.1) then
                   dp=sqrt((xr(npt)-xmmm)**2+(zr(npt)-zmmm)**2)
                   if(dp.lt.stol) go to 104
                 end if
                 if(amax.lt.-999998.) then
                   amax=amin
                 else
                   if(ang.gt.amax) amax=ang
                 end if
                 if(a1.eq.0.) then
                   ang=ang-ainc
                 else
                   ang=(a1+ang)/2.
                 end if
               else
                 if(idray(1).ge.l) then
                   if(a1.eq.0.) then
                     if(amin.lt.-999998.) then
                       ang=ang-ainc
                     else
                       ang=amin-ainc
                     end if
                   else
                     if(amin.lt.-999998.) then
                       ang=(a1+ang)/2.
                     else
                       ang=(a1+amin)/2.
                     end if
                   end if
                 end if
                 if(idray(1).lt.l) then
                   a1=ang
                   if(amin.lt.-999998.) then
                     ang=ang+ainc
                   else 
                     ang=(amin+ang)/2.
                   end if
                 end if
               end if
               xmmm=xr(npt)
               zmmm=zr(npt)
101         continue
          end if  
104       if(amin.lt.-999998.) then
            iflag=1
          else    
            if(isrch.eq.0.and.tang(l,2).lt.100.) then
              amax=tang(l,2) 
            else  
              a2=0.
              ang=amax+ainc
              do 102 n=1,nsmax
                 call autotr(ang,layer1,iblk1,xshot,zshot,ifam,it,npt,
     +                iflag2,irays,nskip,idot,idr,irayps,istep)
                 if(idray(1).eq.l.and.idray(2).eq.1.and.vr(npt,2).ne.
     +             0.) then 
                   amax=ang 
                   tang(l,2)=amax 
                   if(n.ge.nsmin) return
                   if(stol.gt.0..and.n.gt.1) then
                     dp=sqrt((xr(npt)-xmmm)**2+(zr(npt)-zmmm)**2)
                     if(dp.lt.stol) return
                   end if
                   if(a2.eq.0.) then
                     ang=ang+ainc*float(n)
                   else
                     ang=(a2+ang)/2.
                   end if 
                 else
                   a2=ang
                   ang=(amax+ang)/2.
                 end if
                 xmmm=xr(npt) 
                 zmmm=zr(npt)
102           continue 
            end if
          end if  
        else      
          vratio=vshot/vtop 
          if(vratio.gt..99999) then
            amin=aamin 
          else    
            amin=90.-asin(vratio)*pi18
          end if  
          vratio=vshot/vbotom
          if(vratio.gt..99999) then
            amax=aamin
          else    
            amax=90.-asin(vshot/vbotom)*pi18
          end if  
          ainc=(amax-amin)*aainc
          if(ainc.lt.aaimin) ainc=aaimin
          if(l.ne.layer1) then
            if(isrch.eq.0.and.tang(l,1).lt.100.) then
              amin=tang(l,1)
            else
              ang=amin
              a1=0.
              a2=0.
              amin=-999999.
              do 110 n=1,nsmax
                 call autotr(ang,layer1,iblk1,xshot,zshot,ifam,
     +           it,npt,iflag2,irays,nskip,idot,idr,irayps,istep)
                 if(idray(1).ge.l) then
                   a2=ang
                   if(idray(1).eq.l.and.idray(2).eq.1.and.vr(npt,2).
     +               ne.0.) then
                     amin=ang
                     tang(l,1)=amin
                     if(n.ge.nsmin) go to 111
                     if(stol.gt.0..and.n.gt.1) then
                       dp=sqrt((xr(npt)-xmmm)**2+(zr(npt)-zmmm)**2)
                       if(dp.lt.stol) go to 111
                     end if
                   end if
                   if(a1.eq.0.) then
                     ang=ang-ainc
                   else 
                     ang=(a1+a2)/2.
                   end if 
                 else
                   a1=ang
                   if(a2.eq.0.) then
                     ang=ang+ainc
                   else
                     ang=(a1+a2)/2.
                   end if
                 end if
                 xmmm=xr(npt)
                 zmmm=zr(npt)
110           continue
              if(amin.lt.-999998.) iflag=1
            end if
          end if
111       if(isrch.eq.0.and.tang(l,2).lt.100.) then
            amax=tang(l,2)
          else  
            ang=amax
            a1=0.
            a2=0.
            amax=-999999.
            do 130 n=1,nsmax
               call autotr(ang,layer1,iblk1,xshot,zshot,ifam,it,npt,
     +              iflag2,irays,nskip,idot,idr,irayps,istep)
               if((idray(1).lt.l).or.(idray(1).eq.l.and.idray(2).
     +         eq.1.and.vr(npt,2).ne.0.)) then
                 a1=ang
                 if(idray(1).eq.l.and.idray(2).eq.1.and.vr(npt,2).
     +             ne.0.) then
                   amax=ang
                   tang(l,2)=amax 
                   if(n.ge.nsmin) go to 131
                   if(stol.gt.0..and.n.gt.1) then
                     dp=sqrt((xr(npt)-xmmm)**2+(zr(npt)-zmmm)**2)
                     if(dp.lt.stol) go to 131
                   end if
                 end if
                 if(a2.eq.0.) then
                   ang=ang+ainc
                 else
                   ang=(a1+a2)/2.
                 end if
               else
                 a2=ang
                 if(a1.eq.0.) then
                   ang=ang-ainc
                 else 
                   ang=(a1+a2)/2.
                 end if
               end if
               xmmm=xr(npt) 
               zmmm=zr(npt) 
130         continue
          end if 
131       if(amax.lt.amin) then
            aminh=amin
            amin=amax
            amax=aminh
            iflag=0
            return
          end if
          if(amax.lt.-999998.) then
            if(iflag.ne.1) then
              amax=amin 
              iflag=0
            end if
          else  
            if(iflag.ne.0) then 
              amin=amax
              iflag=0
            end if
          end if
        end if    
      end if      
c
      if(idr.eq.2) then
c                 
c       search for reflected ray
c                 
        if(isrch.eq.0.and.tang(l,3).lt.100.) then
          amin=tang(l,3)
          amax=aamax
        else      
          call calvel(xshot,zshot,layer1,iblk1,l,vshot,vtop,vbotom)
          vratio=vshot/vbotom
          if(vratio.gt..99999) then
            amin=aamin 
          else    
            amin=90.-asin(vratio)*pi18
          end if  
          amax=aamax
          a1=0. 
          a2=0. 
          ainc=(amax-amin)*aainc
          if(ainc.lt.aaimin) ainc=aaimin
          ang=amin
          amin=-999999.
          do 210 n=1,nsmax
             call autotr(ang,layer1,iblk1,xshot,zshot,ifam,it,npt,
     +            iflag2,irays,nskip,idot,idr,irayps,istep)
               if((idray(1).lt.l).or.(idray(1).eq.l.and.idray(2).eq.1).
     +           or.(idray(1).eq.l.and.idray(2).eq.2.and.vr(npt,2).
     +           eq.0.).or.(xsmax.gt.0..and.it.eq.1.and.idray(1).eq.l.
     +           and.idray(2).eq.2.and.2.*abs(xshot-xr(npt)).gt.xsmax).
     +           or.(xsmax.gt.0..and.it.eq.0.and.idray(1).eq.l.and.
     +           idray(2).eq.2.and.abs(xshot-xr(npt)).gt.xsmax)) then
               a1=ang
               if(a2.eq.0.) then
                 ang=ang+ainc
               else
                 ang=(a1+a2)/2.
               end if
             else
               a2=ang
               if(idray(1).eq.l.and.idray(2).eq.2) then
                 amin=ang
                 tang(l,3)=amin
                 if(n.ge.nsmin) return
                 if(stol.gt.0..and.n.gt.1) then
                   dp=sqrt((xr(npt)-xmmm)**2+(zr(npt)-zmmm)**2)
                   if(dp.lt.stol) return
                 end if
               end if
               if(a1.eq.0.) then
                 ang=ang-ainc
               else
                 ang=(a1+a2)/2.
               end if
             end if
             xmmm=xr(npt)
             zmmm=zr(npt)
210       continue
          if(amin.lt.-999998.) iflag=1
        end if    
      end if      
c
      if(idr.eq.3) then
c
c       search for head wave
c
        if(isrch.eq.0.and.tang(l,4).lt.100.) then
          amin=tang(l,4)
          amax=tang(l,4)
        else      
          call calvel(xshot,zshot,layer1,iblk1,l+1,vshot,vtop,vbotom)
          vratio=vshot/vtop
          if(vratio.gt..99999) then 
            ang=aamin
          else    
            ang=90.-asin(vratio)*pi18
          end if  
          a1=0.   
          a2=0.   
          call calvel(xshot,zshot,layer1,iblk1,l,vshot,vtop,vbotom)
          if(vtop.lt.vbotom) then
            if(vtop.le.vshot) then 
               ainc=(90.-asin(vratio)*pi18-aamin)*aainc 
            else  
               ainc=(asin(vshot/vtop)-asin(vshot/vbotom))*pi18
            end if 
          else    
            ainc=1.
          end if  
          if(ainc.lt.aaimin) ainc=aaimin
          angm=ang
          do 1100 n=1,nsmax
             call autotr(ang,layer1,iblk1,xshot,zshot,ifam,it,npt,
     +            iflag2,irays,nskip,idot,idr,irayps,istep)
             if(idray(1).eq.l.and.iflag2.eq.2) then
               amax=ang
               amin=ang
               tang(l,4)=ang 
               return
             end if
             if(idray(1).gt.l) then
               a2=ang
               if(a1.eq.0.) then
                 ang=ang-ainc
               else
                 ang=(a1+a2)/2.
               end if
             else 
               a1=ang
               if(a2.eq.0.) then
                 ang=ang+ainc
               else
                 ang=(a1+a2)/2.
               end if
             end if
1100      continue
          ang=angm
          a1=0.
          a2=0.
          do 1200 n=1,nsmax
             call autotr(ang,layer1,iblk1,xshot,zshot,ifam,it,npt,
     +            iflag2,irays,nskip,idot,idr,irayps,istep)
             if(idray(1).eq.l.and.iflag2.eq.2) then
               amax=ang
               amin=ang
               tang(l,4)=ang
               return
             end if
             if(idray(1).le.l) then
               a2=ang
               if(a1.eq.0.) then
                 ang=ang-ainc
               else
                 ang=(a1+a2)/2.
               end if
             else
               a1=ang
               if(a2.eq.0.) then
                 ang=ang+ainc
               else
                 ang=(a1+a2)/2.
               end if
             end if
1200      continue
          if(idiff.eq.2.and.idray(1).ge.l) then
            amax=a2
            amin=a2
            tang(l,4)=ang 
            idifff=1
            return
          end if
          iflag=1 
        end if    
      end if      
c
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine calvel(xshot,zshot,ls,ibs,l,vshot,vtop,vbotom)
c                 
c     calculate the velocity at the shot point and the 
c     maximum velocity between the shot and the top and 
c     bottom of the layer for which rays are being searched 
c     for directly beneath the shot point
c 
      include 'rayinvr.par'
      include 'rayinvr.com'                
c                 
      layer=ls    
      iblk=ibs    
      iwave=1     
      iccbnd=1    
      vshot=vel(xshot,zshot)
      if(icbnd(1).eq.0) then
        vshot=vshot*vsvp(layer,iblk)
        iwave=-iwave 
        iccbnd=2  
      end if      
      vlmax=vshot  
      xfctr=(xshot-xbnd(layer,iblk,1))/
     +      (xbnd(layer,iblk,2)-xbnd(layer,iblk,1))
      vb1=(vm(layer,iblk,4)-vm(layer,iblk,3))*xfctr+vm(layer,iblk,3)
      if(iwave.eq.-1) vb1=vb1*vsvp(layer,iblk)
      if(vb1.gt.vlmax) vlmax=vb1
      nloop=l-layer
      if(nloop.eq.0) then
        vtop=vshot
        vbotom=vlmax 
        return    
      end if      
      do 10 k=1,nloop
         i=layer+k 
         do 20 j=1,nblk(i)
            if(xbnd(i,j,1).le.xshot.and.xbnd(i,j,2).ge.xshot) then
              vti=(vm(i,j,2)-vm(i,j,1))*xfctr+vm(i,j,1)
              vbi=(vm(i,j,4)-vm(i,j,3))*xfctr+vm(i,j,3)
              if(icbnd(iccbnd).eq.k) then
                iwave=-iwave
                iccbnd=iccbnd+1
              end if
              if(iwave.eq.-1) then
                vti=vti*vsvp(i,j)
                vbi=vbi*vsvp(i,j)
              end if
              if(i.lt.l) then
                if(vti.gt.vlmax) vlmax=vti
                if(vbi.gt.vlmax) vlmax=vbi
              else
                if(vti.gt.vlmax) then
                  vtop=vti
                  vlmax=vti
                else
                  vtop=vlmax
                end if
                if(vbi.gt.vlmax) then
                  vbotom=vbi
                else
                  vbotom=vlmax
                end if
              end if
            end if
20       continue 
10    continue    
c                 
      return      
      end         
