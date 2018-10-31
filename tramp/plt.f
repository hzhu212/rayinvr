c
c     version 1.3  Aug 1992
c                 
c     Plotting routines for TRAMP
c
c     ----------------------------------------------------------------
c                 
      subroutine pltmod(ncont,ibnd,imod,iaxlab,ivel,velht,idash,
     +                  idata,iroute)
c                 
c     plot the velocity model 
c                 
      include 'tramp.par'
      real xpc(pnsmth),xcp(pnsmth),zcp(pnsmth)
      integer ic(plcntr)
      include 'tramp.com'
c        
      data ic/plcntr*0/      
      if(iplots.eq.0) then 
        if(isep.eq.0.and.(itx.gt.0.or.idata.ne.0)) sep=sep+orig+tmm
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      end if
      call erase
      if(iaxlab.eq.1) then
        call axtick(xmin,xmax,xtmin,xtmax,ntickx,ndecix)
        call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
        call axis(orig,sep+zmm,xmin,xmax,xmm,xscale,0.,-1,
     +       xtmin,xtmax,ntickx,ndecix,'DISTANCE (km)',13,albht)
        call axis(orig,sep,zmin,zmax,zmm,zscale,90.,1,
     +       ztmin,ztmax,ntickz,ndeciz,'DEPTH (km)',10,albht) 
      end if      
      call box(orig,sep,orig+xmm,sep+zmm)     
c                 
      if(imod.eq.1) then
        if(icntr.eq.0) then 
          do 10 i=1,ncont
             nptsc=nzed(i)
             if(nptsc.eq.1) then
               nptsc=2
               xcp(1)=orig
               zcp(1)=(zm(i,1)-zmax)/zscale+sep
               xcp(2)=xmm+orig
               zcp(2)=zcp(1)
             else
               do 20 j=1,nptsc
                  xcp(j)=(xm(i,j)-xmin)/xscale+orig
                  zcp(j)=(zm(i,j)-zmax)/zscale+sep
20             continue
             end if
             if(idash.ne.1) then
               call line(xcp,zcp,nptsc)
             else
               dash=xmm/150.
               call dashln(xcp,zcp,nptsc,dash)
             end if
10        continue
          if(ibnd.eq.1) then
c                 
            do 30 i=1,nlayer
               if(nblk(i).gt.1) then
                 do 40 j=1,nblk(i)-1
                    xcp(1)=(xbnd(i,j,2)-xmin)/xscale+orig
                    xcp(2)=xcp(1)
                    zp1=s(i,j,1)*xbnd(i,j,2)+b(i,j,1)
                    zp2=s(i,j,2)*xbnd(i,j,2)+b(i,j,2)
                    if(abs(zp2-zp1).gt..001) then
                      zcp(1)=(zp1-zmax)/zscale+sep
                      zcp(2)=(zp2-zmax)/zscale+sep
                      if(idash.ne.1) then
                        call line(xcp,zcp,2)
                      else
                        call dashln(xcp,zcp,2,dash)
                      end if
                    end if
40               continue
               end if 
30          continue
          end if  
          if(ibsmth.eq.2) then
            do 600 i=1,npbnd 
               x=xmin+float(i-1)*xsinc
               xcp(i)=(x-xmin)/xscale+orig
600         continue 
            do 610 i=1,nlayer+1
               do 620 j=1,npbnd
                  zcp(j)=(cosmth(i,j)-zmax)/zscale+sep
620            continue
               call line(xcp,zcp,npbnd)
610         continue
          end if  
        else      
          ncx=int((xmax-xmin)/xcinc+.5)+1
          if(ncx.gt.pnsmth) ncx=pnsmth
          xcinc=(xmax-xmin-.02)/float(ncx-1)
          do 70 i=1,ncx 
             xcp(i)=xmin+.01+float(i-1)*xcinc
             xpc(i)=(xcp(i)-xmin)/xscale+orig
70        continue
          nl=0    
          if(xclab(1).lt.-9998.) xclab(1)=(xmin+xmax)/2.
          do 150 i=1,plcntr
             if(xclab(i).gt.-9998.) then
               nl=nl+1 
               do 160 j=1,ncx-1
                  if(xclab(i).ge.xcp(j).and.xclab(i).le.xcp(j+1)) then
                    ic(i)=j
                    go to 150
                  end if
160            continue
             end if
150       continue
          ht=albht*velht
          do 80 k=1,pcntr
             if(vcntr(k).lt.-998.) go to 80
             do 90 l=1,ncx
                vc=abs(vcntr(k)) 
                do 100 i=1,nlayer
                   do 110 j=1,nblk(i)
                      if(xbnd(i,j,1).le.xcp(l).and.xbnd(i,j,2).ge.
     +                  xcp(l)) then
                        if(vm(i,j,1).eq.0.) go to 100
                        vu=(vm(i,j,2)-vm(i,j,1))*(xcp(l)-xbnd(i,j,1))/
     +                     (xbnd(i,j,2)-xbnd(i,j,1))+vm(i,j,1)
                        if(icntr.lt.0) vu=vu*vsvp(i,j)
                        vl=(vm(i,j,4)-vm(i,j,3))*(xcp(l)-xbnd(i,j,1))/
     +                     (xbnd(i,j,2)-xbnd(i,j,1))+vm(i,j,3)
                        if(icntr.lt.0) vl=vl*vsvp(i,j)
                        if(vc.lt.vu) then
                          zcp(l)=s(i,j,1)*xcp(l)+b(i,j,1)
                          go to 90
                        end if
                        if((vc.ge.vu.and.vc.le.vl).or.
     +                     (vc.le.vu.and.vc.ge.vl)) then
                          z1=s(i,j,1)*xcp(l)+b(i,j,1) 
                          z2=s(i,j,2)*xcp(l)+b(i,j,2)
                          if(vu.ne.vl) then
                            zcp(l)=z1+(vc-vu)/(vl-vu)*(z2-z1)
                          else
                            zcp(l)=z1
                          end if
                          go to 90
                        end if
                      end if
110                continue
100             continue
                zcp(l)=zmax
90           continue
             if(ncsmth.gt.0) then
               do 120 i=1,ncsmth 
                  call smooth(zcp,ncx)
120            continue 
             end if
             if(idump.eq.1) write(23,555) (xcp(i),i=1,ncx)
             if(idump.eq.1) write(23,555) (zcp(i),i=1,ncx)
555          format(10(50f7.2))
             do 130 i=1,ncx
                zcp(i)=(zcp(i)-zmax)/zscale+sep
130          continue
             if(idash.ne.1) then
               call line(xpc,zcp,ncx)
             else 
               dash=xmm/150.
               call dashln(xpc,zcp,ncx,dash)
             end if
             if(vcntr(k).gt.0..and.nl.gt.0) then
               do 140 i=1,nl 
                  if(ic(i).gt.0) then
                    call number(xpc(ic(i)),
     +              zcp(ic(i))+.1*ht,ht,abs(vcntr(k)),0.,2)
                  end if
140            continue
             end if
80        continue
        end if    
      end if      
      if(ivel.ne.0.and.icntr.eq.0) then
        ht=albht*velht
        do 50 i=1,nlayer
           do 60 j=1,nblk(i)
              htt=ht
              if(ivg(i,j).ne.-1) then
165             xp=xbnd(i,j,1)
                z1=s(i,j,1)*xp+b(i,j,1)
                zp1=(z1-zmax)/zscale-1.5*htt+sep
                xp1=(xp-xmin)/xscale+0.5*htt+orig
                xp=xbnd(i,j,2)
                z2=s(i,j,1)*xp+b(i,j,1)
                zp2=(z2-zmax)/zscale-1.5*htt+sep
                xp2=(xp-xmin)/xscale-4.0*htt+orig
                xp=xbnd(i,j,1)
                z3=s(i,j,2)*xp+b(i,j,2)
                zp3=(z3-zmax)/zscale+0.5*htt+sep
                xp3=xp1
                xp=xbnd(i,j,2)
                z4=s(i,j,2)*xp+b(i,j,2)
                zp4=(z4-zmax)/zscale+0.5*htt+sep
                xp4=xp2
c
                if(abs(z1-z3).gt..005.and.abs(z2-z4).gt..005) then
                  xhmin=4.*htt
                else
                  xhmin=8.*htt
                end if
c
                if(xp2-xp1.lt.xhmin.or.(abs(z1-z3).gt..005.and.
     +            zp1-zp3.le.1.5*htt).or.(abs(z2-z4).gt..005.and.
     +            zp2-zp4.le.1.5*htt)) then
                  htt=htt/1.25
                  go to 165
                end if
c
                v1=vm(i,j,1)
                if(ivel.lt.0) v1=v1*vsvp(i,j)
                if(abs(z1-z3).gt..005) then
                  call number(xp1,zp1,htt,v1,0.,2)
                  v3=vm(i,j,3)
                  if(ivel.lt.0) v3=v3*vsvp(i,j)
                  call number(xp3,zp3,htt,v3,0.,2)
                else
                  call number(xp1+4.*htt,zp1,htt,v1,0.,2)
                end if
                v2=vm(i,j,2)
                if(ivel.lt.0) v2=v2*vsvp(i,j)
                if(abs(z2-z4).gt..005) then
                  call number(xp2,zp2,htt,v2,0.,2)
                  v4=vm(i,j,4)
                  if(ivel.lt.0) v4=v4*vsvp(i,j)
                  call number(xp4,zp4,htt,v4,0.,2)
                else
                  call number(xp2-4.*htt,zp2,htt,v2,0.,2)
                end if
              end if
60         continue 
50      continue  
      end if      
c
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine plttx(ifam,npts,iszero,idata,iaxlab,xshot,
     +                 idr,nshot,itxout,ibrka,ttunc,itrev,iroute)
c                 
c     plot the travel time curves for all rays reaching the surface -
c     for each group of rays a separate curve is drawn for each
c     ray code    
c                 
      include 'tramp.par'
      real xh(pnrayf),th(pnrayf),rh(pnrayf),xs(pnrayf),f(pnrayf),
     +     x(pnrayf),t(pnrayf),xshot(1),xsh(pnrayf),fh(pnrayf)
      character tlab*12 
      integer npts(1),idr(1),ibrka(1)
      include 'tramp.com'
c                 
      if(vred.ne.0.) then
        rvred=1./vred
      else
        rvred=0.
      end if
      n=0
      ib=0                
      if(iplots.eq.0) then 
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      end if
      if(isep.gt.0) call erase
      if(itrev.ne.1) then
        tadj=tmin
      else
        tadj=tmax
      end if
c
      if(iaxlab.eq.1) then
        call axtick(xmint,xmaxt,xtmint,xtmaxt,ntckxt,ndecxt)
        call axis(orig,orig,xmint,xmaxt,xmmt,xscalt,0.,1,
     +       xtmint,xtmaxt,ntckxt,ndecxt,'DISTANCE (km)',13,albht)
        if(vred.eq.0.) then
          nchart=8
          tlab='TIME (s)'
        else     
          i1=int(vred)
          id1=int((vred-float(i1))*10+.05)
          id2=int((vred-float(i1)-float(id1)/10.)*100.+.5)
          if(id1.eq.0.and.id2.eq.0) then
            nchart=9
            tlab='T-D/  (s)'
            tlab(5:5)=char(i1+48)
          else   
            if(id2.eq.0) then
              nchart=11
              tlab='T-D/    (s)'
              tlab(5:7)=char(i1+48)//'.'//char(id1+48)
            else 
              nchart=12
              tlab='T-D/     (s)'
              tlab(5:8)=char(i1+48)//'.'//char(id1+48)//char(id2+48)
            end if
          end if 
        end if
        call axtick(tmin,tmax,ttmin,ttmax,ntickt,ndecit)
        call axis(orig,orig,tmin,tmax,tmm,tscale,90.,1,
     +       ttmin,ttmax,ntickt,ndecit,tlab,nchart,albht)
      end if      
      call box(orig,orig,orig+xmmt,orig+tmm)           
c                 
      if(idata.eq.1) call pltdat(iszero,xshot,idr,nshot,tadj)
c                 
      xshoth=-99999.
      fidh=0.
      do 10 i=1,ifam
         nh1=n    
         if(npts(i).gt.0) then
           do 20 j=1,npts(i)
              xh(j)=range(j+nh1)
              th(j)=tt(j+nh1)
              rh(j)=rayid(j+nh1) 
              xsh(j)=xshtar(j+nh1)
              fh(j)=fidarr(j+nh1)
              n=n+1
20         continue
           nh2=0  
1000       if(nh2.lt.npts(i)) then
             k=0  
100          k=k+1
             if((k+nh2).le.npts(i)) then
               x(k)=(xh(k+nh2)-xmint)/xscalt+orig
               t(k)=(th(k+nh2)-tadj)/tscale+orig
               xs(k)=xsh(k+nh2)
               f(k)=fh(k+nh2)
               if(k.eq.1) then
                 rhc=rh(k+nh2)
                 go to 100
               else
                 if(rh(k+nh2).eq.rhc) go to 100
               end if
             end if 
             npt=k-1
             nh2=nh2+npt
             if(ibrka(i).eq.0.or.itx.gt.1.or.npt.lt.3) then
               if(itxout.gt.0) then
                 if(xs(1).ne.xshoth.or.f(1).ne.fidh) then
                   ib=0
                   write(17,5) xs(1),f(1),0.,0
                   xshoth=xs(1)
                   fidh=f(1)
                 end if
                 ib=ib+1
                 do 40 j=1,npt
                    xw=(x(j)-orig)*xscalt+xmint
                    tw=(t(j)-orig)*tscale+tadj+abs(xshoth-xw)*rvred
                    write(17,5) xw,tw,ttunc,ib
5                   format(3f10.3,i10)
40               continue
               end if
               if(npt.gt.1) then
                 if(itx.lt.2) then
                   call line(x,t,npt)
                 else
                   if(symht.gt.0.) then
                     do 41 j=1,npt
                        call ssymbl(x(j),t(j),symht,1)
41                   continue
                   end if
                 end if
               else
                 if(symht.gt.0.) call ssymbl(x(1),t(1),symht,1)
               end if
             else 
               ih=0
               if(itxout.gt.0) then
                 if(xs(1).ne.xshoth.or.f(1).ne.fidh) then
                   ib=0
                   write(17,5) xs(1),f(1),0.,0
                   xshoth=xs(1)
                   fidh=f(1)
                 end if
                 ib=ib+1 
                 xw=(x(1)-orig)*xscalt+xmint
                 tw=(t(1)-orig)*tscale+tadj+abs(xshoth-xw)*rvred
                 write(17,5) xw,tw,ttunc,ib
                 xw=(x(2)-orig)*xscalt+xmint
                 tw=(t(2)-orig)*tscale+tadj+abs(xshoth-xw)*rvred
                 write(17,5) xw,tw,ttunc,ib
               end if
               call plot(x(1),t(1),3) 
               call plot(x(2),t(2),2)
               do 30 j=3,npt
                  if(((x(j)-x(j-1))*(x(j-1)-x(j-2))).ge.0.or.
     +              ih.eq.(j-1)) then
                    xw=(x(j)-orig)*xscalt+xmint
                    tw=(t(j)-orig)*tscale+tadj+abs(xshoth-xw)*rvred
                    if(itxout.gt.0) write(17,5) xw,tw,ttunc,ib
                    call plot(x(j),t(j),2)
                  else
                    if(itxout.gt.0) then
                      ib=ib+1 
                      xw=(x(j)-orig)*xscalt+xmint
                      tw=(t(j)-orig)*tscale+tadj+abs(xshoth-xw)*rvred
                      write(17,5) xw,tw,ttunc,ib
                    end if
                    call plot(x(j),t(j),3)
                    ih=j 
                  end if
30             continue
               if(ih.eq.npt.and.symht.gt.0.) 
     +           call ssymbl(x(npt),t(npt),symht,1)
             end if 
             go to 1000
           end if 
         end if   
10    continue    
      if(itxout.gt.0) write(17,5) 0.,0.,0.,-1
c
      call empty
c
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine pltdat(iszero,xshot,idr,nshot,tadj)
c                 
c     plot observed travel times 
c                 
      include 'tramp.par'
      real xshot(1)
      integer idr(1)
      include 'tramp.com'
      nsfc=1
      isf=ilshot(nsfc)
c                
100   xp=xpf(isf)
      tp=tpf(isf)
      up=upf(isf)
      ip=ipf(isf)
c
      if(ip.lt.0) go to 999
      if(ip.eq.0) then
        xsp=xp   
        idp=sign(1.,tp)
        do 10 i=1,nshot
           xdiff=abs(xshot(i)-xsp)
           if(xdiff.lt..001.and.idr(i).eq.idp) then
             iplt=1
             isf=isf+1
             nsfc=nsfc+1
             icshot=i
             go to 100
           end if
10      continue
        iplt=0
        nsfc=nsfc+1
        isf=ilshot(nsfc)
        go to 100
      else
30      if(iplt.eq.1) then
          if(iszero.eq.0) then
            xplot=(xp-xmint)/xscalt+orig
          else
            xplot=((xp-xsp)*float(idp)-xmint)/xscalt+orig
          end if
          tplot=(tp-tadj)/tscale+orig
          if(itcol.ne.0) then
            if(itcol.ne.2) then
              ipcol=colour(mod(ip-1,10)+1)
            else
              ipcol=colour(mod(icshot-1,10)+1)
            end if
            call pcolor(ipcol)
          end if
          call plot(xplot,tplot-up/tscale,3)
          call plot(xplot,tplot+up/tscale,2)
        end if
      end if
      isf=isf+1
      go to 100
c          
999   continue
      if(itcol.ne.0) call pcolor(ifcol)
      return
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine pltray(npt,nskip,idot,irayps,istep,anglew)
c                 
c     plot one ray
c                 
      include 'tramp.par'
      real x(ppray),z(ppray),xa(ppray),za(ppray)
      character reply*1
      include 'tramp.com'
c                 
      npts=npt-nskip
      if(npts.lt.2) return    
c
      do 10 i=1,npts
         x(i)=(xr(i+nskip)-xmin)/xscale+orig
         z(i)=(zr(i+nskip)-zmax)/zscale+sep
10    continue    
c
      if(istep.eq.1) write(6,5) anglew
5     format('take-off angle: ',f10.5$)
c
      if(irayps.ne.1) then
        call line(x,z,npts)
      else
        dash=xmm/250.
        na=0
        if(vr(1+nskip,2).eq.vp(1+nskip,2)) then
          ips=1
        else
          ips=-1
        end if
        do 30 i=1,npts-1
           if(vr(i+nskip,2).eq.vp(i+nskip,2)) then
             if(ips.ne.1) then
               na=na+1
               xa(na)=x(i)
               za(na)=z(i)
               call dashln(xa,za,na,dash)
               na=0
               ips=1
             end if
             na=na+1
             xa(na)=x(i)
             za(na)=z(i)
           else
             if(ips.ne.-1) then
               na=na+1
               xa(na)=x(i)
               za(na)=z(i)
               call line(xa,za,na)
               na=0
               ips=-1
             end if
             na=na+1
             xa(na)=x(i)
             za(na)=z(i)
           end if
30      continue
        na=na+1
        xa(na)=x(npts)
        za(na)=z(npts)
        if(ips.eq.1) then
          call line(xa,za,na)
        else
          call dashln(xa,za,na,dash)
        end if
      end if
c
      if(idot.eq.1) then
        do 20 i=1,npts
           call ssymbl(x(i),z(i),symht,4)
20      continue  
      end if      
c
      if(istep.eq.1) then
        call empty
        read(5,15) reply
15      format(a1)
        if(reply(1:1).eq.'s') then
          call plotnd(1)
          stop
        end if
      end if
c
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine pltamp(ifam,npts,iaxlab,ibrka,iamout,ampunc,iroute,
     +                  amp1,amp2)
c                 
c     plot the amplitude-distance curves for all rays reaching the
c     surface - for each group of rays a separate curve is drawn 
c     for each ray code
c                 
      include 'tramp.par'
      real xh(pnrayf),ah(pnrayf),rh(pnrayf),xs(pnrayf),f(pnrayf),
     +     x(pnrayf),a(pnrayf),xsh(pnrayf),fh(pnrayf)
      integer npts(1),ibrka(1)
      include 'tramp.com'
c        
      ampmin=amp1
      ampmax=amp2
      ampunc=ampunc/100.
      n=0         
      ib=0
      if(iplots.eq.0) then
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      else
        call aldone
        call erase
      end if
c
      if(iaxlab.eq.1) then
        call axtick(xmint,xmaxt,xtmint,xtmaxt,ntckxt,ndecxt)
        call axtick(ampmin,ampmax,atmin,atmax,nticka,ndecia)
        call axis(orig,orig,xmint,xmaxt,xmmt,xscalt,0.,1,
     +       xtmint,xtmaxt,ntckxt,ndecxt,'DISTANCE (km)',13,albht)
        call axis(orig,orig,ampmin,ampmax,amm,ascale,90.,1,
     +       atmin,atmax,nticka,ndecia,'Log(Amplitude)',14,albht)
      end if      
      call box(orig,orig,orig+xmmt,orig+amm)           
c                 
      xshoth=-99999.
      fidh=0.
      do 10 i=1,ifam
         nh1=n    
         if(npts(i).gt.0) then
           nptsa=0
           do 20 j=1,npts(i)
              if(amp(j+nh1).gt.0.) then
                nptsa=nptsa+1
                xh(nptsa)=range(j+nh1)
                ah(nptsa)=amp(j+nh1)
                rh(nptsa)=rayid(j+nh1) 
                xsh(nptsa)=xshtar(j+nh1)
                fh(nptsa)=fidarr(j+nh1)
              end if
              n=n+1
20         continue
           if(nptsa.gt.0) then
             nh2=0
1000         if(nh2.lt.nptsa) then
               k=0
100            k=k+1
               if((k+nh2).le.nptsa) then
                 x(k)=(xh(k+nh2)-xmint)/xscalt+orig
                 a(k)=(alog10(ah(k+nh2))-ampmin)/ascale+orig
                 xs(k)=xsh(k+nh2)
                 f(k)=fh(k+nh2)
                 if(k.eq.1) then
                   rhc=rh(k+nh2)
                   go to 100
                 else
                   if(rh(k+nh2).eq.rhc) go to 100
                 end if
               end if 
               npt=k-1
               nh2=nh2+npt
               if(ibrka(i).eq.0.or.itx.gt.1.or.npt.lt.3) then
                 if(iamout.eq.1) then
                   if(xs(1).ne.xshoth.or.f(1).ne.fidh) then
                     ib=0
                     write(24,5) xs(1),f(1),0.,0
                     xshoth=xs(1)
                     fidh=f(1)
                   end if
                   ib=ib+1
                   do 40 j=1,npt
                      ampf=10.**((a(j)-orig)*ascale+ampmin)
                      write(24,5) (x(j)-orig)*xscalt+xmint,ampf,
     +                            ampf*ampunc,ib
5                     format(f10.3,2e10.3,i10)
40                 continue
                 end if
                 if(npt.gt.1) then
                   if(itx.lt.2) then
                     call line(x,a,npt)
                   else
                     if(symht.gt.0.) then
                       do 41 j=1,npt
                          call ssymbl(x(j),a(j),symht,1)
41                     continue
                     end if
                   end if
                 else
                   if(symht.gt.0.) call ssymbl(x(1),a(1),symht,1)
                 end if
               else
                 ih=0
                 if(iamout.eq.1) then
                   if(xs(1).ne.xshoth.or.f(1).ne.fidh) then
                     ib=0
                     write(24,5) xs(1),f(1),0.,0
                     xshoth=xs(1)
                     fidh=f(1)
                   end if
                   ib=ib+1
                   ampf=10.**((a(1)-orig)*ascale+ampmin)
                   write(24,5) (x(1)-orig)*xscalt+xmint,ampf,
     +                         ampf*ampunc,ib
                   ampf=10.**((a(2)-orig)*ascale+ampmin)
                   write(24,5) (x(2)-orig)*xscalt+xmint,ampf,
     +                         ampf*ampunc,ib
                 end if
                 call plot(x(1),a(1),3) 
                 call plot(x(2),a(2),2)
                 do 30 j=3,npt
                    if(((x(j)-x(j-1))*(x(j-1)-x(j-2))).ge.0.or.
     +                ih.eq.(j-1)) then
                      if(iamout.eq.1) then
                        ampf=10.**((a(j)-orig)*ascale+ampmin)
                        write(24,5) (x(j)-orig)*xscalt+xmint,
     +                  ampf,ampf*ampunc,ib
                      end if
                      call plot(x(j),a(j),2)
                    else
                      if(iamout.eq.1) then
                        ib=ib+1
                        ampf=10.**((a(j)-orig)*ascale+ampmin)
                        write(24,5) (x(j)-orig)*xscalt+xmint,ampf,
     +                              ampf*ampunc,ib
                      end if
                      call plot(x(j),a(j),3)
                      ih=j 
                    end if
30               continue
                 if(ih.eq.npt.and.symht.gt.0.) 
     +             call ssymbl(x(npt),a(npt),symht,1)
                end if 
                go to 1000
             end if
           end if 
         end if   
10    continue    
      if(iamout.eq.1) write(24,5) 0.,0.,0.,-1
c
      call empty
c
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine pltrms(xrinc,vrmstl,vrmsbl,nrsmth,iaxlab,ivrms,
     +                  iroute)
c                 
c     plot rms velocity versus distance
c                 
      include 'tramp.par'
      real vprms(ppvrms),xprms(ppvrms),vrms(ppvrms),xrms(ppvrms),k
      integer vrmstl,vrmsbl
      include 'tramp.com'
c      
      if(iplots.eq.0) then
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      else
        call aldone
        call erase
      end if
      if(iaxlab.eq.1) then
        call axtick(xmint,xmaxt,xtmint,xtmaxt,ntckxt,ndecxt)
        call axtick(vrmin,vrmax,vrtmin,vrtmax,ntckvr,ndecir)
        call axis(orig,orig,xmint,xmaxt,xmmt,xscalt,0.,1,
     +       xtmint,xtmaxt,ntckxt,ndecxt,'DISTANCE (km)',13,albht)
        call axis(orig,orig,vrmin,vrmax,vrmm,rscale,90.,1,
     +       vrtmin,vrtmax,ntckvr,ndecir,'RMS VELOCITY (km/s)',19,albht)
      end if
      call box(orig,orig,orig+xmmt,orig+vrmm)           
c                 
      nrms=int((xmaxt-xmint)/xrinc+.5)+1
      if(nrms.gt.ppvrms) nrms=ppvrms
      xrinc=(xmaxt-xmint-.02)/float(nrms-1)
      do 100 ii=1,nrms
         xrms(ii)=xmint+.01+float(ii-1)*xrinc
         tsum=0.  
         vsum=0.  
         do 10 i=vrmstl,vrmsbl
            do 20 j=1,nblk(i)
               if(xbnd(i,j,1).le.xrms(ii).and.xbnd(i,j,2).ge.xrms(ii))
     +         then 
                 z1=s(i,j,1)*xrms(ii)+b(i,j,1)
                 z2=s(i,j,2)*xrms(ii)+b(i,j,2)
                 h=z2-z1
                 if(abs(h).lt..001) go to 10
                 vu=(vm(i,j,2)-vm(i,j,1))*(xrms(ii)-xbnd(i,j,1))/
     +              (xbnd(i,j,2)-xbnd(i,j,1))+vm(i,j,1)
                 if(ivrms.lt.0) vu=vu*vsvp(i,j)
                 vl=(vm(i,j,4)-vm(i,j,3))*(xrms(ii)-xbnd(i,j,1))/
     +              (xbnd(i,j,2)-xbnd(i,j,1))+vm(i,j,3)
                 if(ivrms.lt.0) vl=vl*vsvp(i,j)
                 if(vu.le..001.or.vl.le..001) then
                   vsum=0.
                   tsum=1. 
                   go to 110
                 end if 
                 k=abs((vl-vu)/h)
                 if(k.ne.0.) then
                   tsum=tsum+alog(1.+k*h/vu)/k
                 else 
                   tsum=tsum+h/vu
                 end if
                 vsum=vsum+vu*h+0.5*k*h**2
                 go to 10
               end if
20          continue
10       continue 
110      xprms(ii)=(xrms(ii)-xmint)/xscalt+orig
         vrms(ii)=sqrt(vsum/tsum)
         vprms(ii)=(sqrt(vsum/tsum)-vrmin)/rscale+1
100   continue    
      if(nrsmth.gt.0) then
        do 30 i=1,nrsmth
           call smooth(vprms,nrms)
30      continue  
      end if      
c                 
      call line(xprms,vprms,nrms)
      call empty
c                 
      if(idump.eq.1) then
        write(23,5) (xrms(i),i=1,nrms)
        write(23,5) (vrms(i),i=1,nrms)
5       format(10(50f7.2))
      end if      
c                 
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine pltvz(xvz,iaxlab,ivz,ivzp,iroute)
c                 
c     plot velocity versus depth at x=xvz
c                    
      include 'tramp.par'
      real v(player*2),z(player*2),xvz(1)
      include 'tramp.com'
c                 
      if(iplots.eq.0) then
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      else
        call aldone
      end if
c
      iflagp=0
      zscale=zscale/2.
      zmm=2.*zmm
      if(xvz(1).lt.-9999.) xvz(1)=(xmin+xmax)/2.
      do 100 k=1,pnvz
        if(xvz(k).ge.xmin.and.xvz(k).le.xmax) then
          if(idump.eq.1) write(23,5) xvz(k)
5         format(/'x position: ',f7.2,' km'/)
          if(iflagp.eq.0) call erase
          if(iflagp.eq.0.and.iaxlab.eq.1) then
             call axtick(vmin,vmax,vtmin,vtmax,ntickv,ndeciv)
             call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
             call axis(orig,orig+zmm,vmin,vmax,vmm,vscale,0.,-1,
     +            vtmin,vtmax,ntickv,ndeciv,'VELOCITY (km/s)',15,albht)
             call axis(orig,orig,zmin,zmax,zmm,zscale,90.,1,
     +            ztmin,ztmax,ntickz,ndeciz,'DEPTH (km)',10,albht)
          end if
          if(iflagp.eq.0) call box(orig,orig,orig+vmm,orig+zmm)       
c                 
          if(ivz.gt.0) then 
            iw=1  
          else    
            iw=-1 
          end if  
          idash=0  
          io=1    
1000      np=0    
          do 10 i=1,nlayer
             do 20 j=1,nblk(i)
                if(xbnd(i,j,1).le.xvz(k).and.xbnd(i,j,2).ge.xvz(k))
     +          then
                  if(vm(i,j,1).gt.0.) then
                    np=np+1
                    z1=s(i,j,1)*xvz(k)+b(i,j,1)
                    z2=s(i,j,2)*xvz(k)+b(i,j,2)
                    z(np*2-1)=(z1-zmax)/zscale+orig
                    z(np*2)=(z2-zmax)/zscale+orig
                    vu=(vm(i,j,2)-vm(i,j,1))*(xvz(k)-xbnd(i,j,1))/
     +                 (xbnd(i,j,2)-xbnd(i,j,1))+vm(i,j,1)
                    if(iw.eq.-1) vu=vu*vsvp(i,j)
                    vl=(vm(i,j,4)-vm(i,j,3))*(xvz(k)-xbnd(i,j,1))/
     +                 (xbnd(i,j,2)-xbnd(i,j,1))+vm(i,j,3)
                    if(iw.eq.-1) vl=vl*vsvp(i,j) 
                    v(np*2-1)=(vu-vmin)/vscale+orig
                    v(np*2)=(vl-vmin)/vscale+orig
                    if(idump.eq.1.and.io.eq.1) 
     +               write(23,15) i,z1,z2,vu,vl,
     +               vu*vsvp(i,j),vl*vsvp(i,j)
15                   format('layer# ',i2,'   z1=',f7.2,'   z2=',f7.2,
     +                 ' km'/9x,'  vp1=',f7.2,'  vp2=',f7.2,' km/s'/ 
     +                       9x,'  vs1=',f7.2,'  vs2=',f7.2,' km/s') 
                  end if
                  go to 10
                end if
20           continue
10        continue
c
          if(idash.ne.1) then
            call line(v,z,2*np)
          else
            call dashln(v,z,2*np,dash)
          end if
          if(ivzp.eq.1) iflagp=1
          if(iw.eq.-1) go to 110
          if(ivz.eq.2) then 
            iw=-1
            dash=vmm/75.
            idash=1
            io=0  
            go to 1000 
          end if  
          call empty
110       if(iflagp.eq.0) call aldone
        end if
100   continue
c
      return
      end
