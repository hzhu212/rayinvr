c
c     version 1.2  Mar 1992
c
c     Plotting routines for VMODEL
c
c     ----------------------------------------------------------------
c
      subroutine pltmod(ncont,iaxlab,ivel,velht,idash,izort,xtinc,
     +           xsinc,ntsmth,idump,ifrefl,symht,iroute,inode,idcol,
     +           ivcol,ivarz,ivarv,izrefl,ntsmt,npsmt,ndot,dashlf)
c
c     plot the 2-D velocity model versus depth or time
c
      include 'vmodel.par'
      real xcp(ppcntr),zcp(ppcntr),xpts(ppnpts),xppts(ppnpts),
     +     tpts(ppnpts),tppts(ppnpts),xfrefl(ppfref),zfrefl(ppfref),
     +     xcpz(ppnpts),zcpz(ppnpts),xpl(ppnpts),zpl(ppnpts),
     +     xzf(ppzff),zff(ppzff),zpts(ppnpts),vpts(ppnpts)
      integer ivarz(player,ppcntr),ivarv(player,ppvel,2)
      include 'vmodel.com'
      CHARACTER(LEN=50) FMT
c
      call dxmin(ncont)
c
      if(iplots.eq.0) then
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      end if
      call erase
      if(iaxlab.eq.1) then
        call axtick(xmin,xmax,xtmin,xtmax,ntickx,ndecix)
        if(izort.eq.0) then
          call axis(orig,orig+zmm,xmin,xmax,xmm,xscale,0.,-1,
     +         xtmin,xtmax,ntickx,ndecix,'DISTANCE (km)',13,albht)
          call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
          call axis(orig,orig,zmin,zmax,zmm,zscale,90.,1,
     +         ztmin,ztmax,ntickz,ndeciz,'DEPTH (km)',10,albht)
        else
          call axis(orig,orig+tmm,xmin,xmax,xmm,xscale,0.,-1,
     +         xtmin,xtmax,ntickx,ndecix,'DISTANCE (km)',13,albht)
          call axtick(tmin,tmax,ttmin,ttmax,ntickt,ndecit)
          call axis(orig,orig,tmin,tmax,tmm,tscale,90.,1,
     +         ttmin,ttmax,ntickt,ndecit,'TIME (s)',8,albht)
        end if
      end if
      if(izort.eq.0) then
        call box(orig,orig,orig+xmm,orig+zmm)
      else
        call box(orig,orig,orig+xmm,orig+tmm)
      end if
c
      if(izort.eq.0) then
        do 10 i=1,ncont
           nptsc=nzed(i)
           if(nptsc.eq.1) then
             nptsc=2
             xcp(1)=orig
             zcp(1)=(zm(i,1)-zmax)/zscale+orig
             xcp(2)=xmm+orig
             zcp(2)=zcp(1)
           else
             do 20 j=1,nptsc
                xcp(j)=(xm(i,j)-xmin)/xscale+orig
                zcp(j)=(zm(i,j)-zmax)/zscale+orig
20           continue
           end if
           if(idash.ne.1) then
             if(idash.ne.-1) call line(xcp,zcp,nptsc)
           else
             dash=xmm/150.
             call dashln(xcp,zcp,nptsc,dash)
           end if
           if(inode.eq.1.and.i.lt.ncont) then
             call pcolor(idcol)
             if(nzed(i).eq.1) then
               if(ivarz(i,1).eq.1) then
                 xpp=(xmin+xmax)/2.
                 xp=(xpp-xmin)/xscale+orig
                 zp=(zm(i,1)-zmax)/zscale+orig
                 call ssymbl(xp,zp,symht,4)
               end if
             else
               do 3220 j=1,nptsc
                  if(ivarz(i,j).eq.1) then
                    xp=(xm(i,j)-xmin)/xscale+orig
                    zp=(zm(i,j)-zmax)/zscale+orig
                    call ssymbl(xp,zp,symht,4)
                  end if
3220           continue
             end if
             call pcolor(ifcol)
           end if
10      continue
c
        if(ivel.eq.1.or.inode.eq.1) then
          ht=albht*velht
          do 50 i=1,ncont-1
             do 60 j=1,nvel(i,1)
                vu=vf(i,j,1)
                if(nvel(i,1).eq.1) then
                  if(vu.eq.0.) go to 60
                  xp=(xmax-xmin)/2
                else
                  xp=xvel(i,j,1)
                end if
                if(nzed(i).eq.1) then
                  zu=zm(i,1)
                else
                  do 70 k=1,nzed(i)-1
                     if(xp.ge.xm(i,k).and.xp.le.xm(i,k+1)) then
                       zu=(zm(i,k+1)-zm(i,k))*(xp-xm(i,k))/(xm(i,k+1)-
     +                     xm(i,k))+zm(i,k)
                       go to 71
                     end if
70                continue
                end if
71              xpp=(xp-xmin)/xscale+orig
                zpp=(zu-zmax)/zscale+orig
                if(ivel.eq.1)
     +            call number(xpp-1.25*ht,zpp-1.25*ht,ht,vu,0.,2)
                if(inode.eq.1.and.ivarv(i,j,1).eq.1)
     +            call dot(xpp,zpp-0.375*symht,0.75*symht,ivcol)
60           continue
             do 80 j=1,nvel(i,2)
                vl=vf(i,j,2)
                if(nvel(i,2).eq.1) then
                  if(vl.eq.0.) go to 80
                  xp=(xmax-xmin)/2
                else
                  xp=xvel(i,j,2)
                end if
21              if(nzed(i+1).eq.1) then
                  zl=zm(i+1,1)
                else
                  do 90 k=1,nzed(i+1)-1
                     if(xp.ge.xm(i+1,k).and.xp.le.xm(i+1,k+1)) then
                       zl=(zm(i+1,k+1)-zm(i+1,k))*(xp-xm(i+1,k))/
     +                     (xm(i+1,k+1)-xm(i+1,k))+zm(i+1,k)
                       go to 91
                     end if
90                continue
                end if
91              xpp=(xp-xmin)/xscale+orig
                zpp=(zl-zmax)/zscale+orig
                if(ivel.eq.1)
     +            call number(xpp-1.25*ht,zpp+.25*ht,ht,vl,0.,2)
                if(inode.eq.1.and.ivarv(i,j,2).eq.1)
     +            call dot(xpp,zpp+0.375*symht,0.75*symht,ivcol)
80           continue
50        continue
        end if
c
        if(idump.eq.1) then
          if(xsinc.lt.0.) xsinc=(xmax-xmin)/100.
          npts=int((xmax-xmin)/xsinc+.5)+1
          if(npts.gt.ppnpts) npts=ppnpts
          xsinc=(xmax-xmin-.002)/float(npts-1)
          do 160 i=1,npts
             xpts(i)=xmin+.001+float(i-1)*xsinc
160       continue
          write(14,5) (xpts(j),j=1,npts)
5         format(10(50f7.2))
c
c          do 170 i=1,ncont-1
c             do 180 j=1,npts
c                xp=xpts(j)
c                if(nzed(i).eq.1) then
c                  zpts(j)=zm(i,1)
c                else
c                  do 190 k=1,nzed(i)-1
c                     if(xp.ge.xm(i,k).and.xp.le.xm(i,k+1)) then
c                       zpts(j)=(zm(i,k+1)-zm(i,k))*(xp-xm(i,k))/
c     +                    (xm(i,k+1)-xm(i,k))+zm(i,k)
c                       go to 180
c                     end if
c190               continue
c                end if
c180          continue
c             write(14,5) (zpts(j),j=1,npts)
c170       continue
          do 170 i=1,ncont-1
             do 180 j=1,npts
                xp=xpts(j)
                if(nvel(i,1).eq.1) then
                  vpts(j)=vf(i,j,1)
                else
                  do 190 k=1,nvel(i,1)-1
                     if(xp.ge.xvel(i,k,1).and.xp.le.xvel(i,k+1,1))
     +               then
                 vpts(j)=(vf(i,k+1,1)-vf(i,k,1))*(xp-xvel(i,k,1))/
     +                   (xvel(i,k+1,1)-xvel(i,k,1))+vf(i,k,1)
                       go to 180
                     end if
190               continue
                end if
180          continue
             write(14,5) (vpts(j),j=1,npts)
170       continue
        end if
c
        if(ifrefl.eq.1) then
c
          open(unit=30, file='f.in', status='old')
c
790       read(30,745,end=795) nfrefr
745       format(i2)
          WRITE(FMT,'("(3x,", I0, "f7.2)")') nfrefr
          read(30,FMT) (xfrefl(i),i=1,nfrefr)
          read(30,FMT) (zfrefl(i),i=1,nfrefr)
          read(30,765)
C 755       format(3x,<nfrefr>f7.2)
765       format(' ')
c
          ipen=3
          do 7010 ii=1,nfrefr
             xp=xfrefl(ii)
             zp=zfrefl(ii)
             xp=(xp-xmin)/xscale+orig
             zp=(zp-zmax)/zscale+orig
             call plot(xp,zp,ipen)
             ipen=2
c
7010      continue
          go to 790
795       continue
        end if
c
        if(ifrefl.eq.2) then
c
          open(unit=30, file='f.in', status='old')
c
890       read(30,845,end=895) xp,zp
845       format(2f10.3)
c
          xp=(xp-xmin)/xscale+orig
          zp=(zp-zmax)/zscale+orig
          call dot(xp,zp,symht,ifcol)
c         call ssymbl(xp,zp,symht,1)
c
          go to 890
895       continue
        end if
c
        if(izrefl.eq.1) then
c
          open(unit=36, file='zf.in', status='old')
c
395       read(36,365,end=399) np
365       format(i10)
c
          xinc=(xmax-xmin)/float(np-1)
          do 3910 i=1,np
             xcpz(i)=xmin+float(i-1)*xinc
3910      continue
c
          read(36,385) (zcpz(j),j=1,np)
385       format(10f10.3)
c
          if(ntsmt.gt.0) call smooth2(zcpz,np,npsmt,ntsmt)
c
          read(36,3365) npzf
3365      format(i2)
          read(36,3385) (xzf(j),j=1,npzf)
          read(36,3385) (zff(j),j=1,npzf)
3385      format(10f7.2)
c
          do 3950 i=1,np
             if(xzf(1).ge.xcpz(i).and.xzf(1).le.xcpz(i+1)) then
               zp1=(zcpz(i+1)-zcpz(i))/(xcpz(i+1)-xcpz(i))*
     +             (xzf(1)-xcpz(i))+zcpz(i)
             end if
             if(xzf(npzf).ge.xcpz(i).and.xzf(npzf).le.xcpz(i+1)) then
               zp2=(zcpz(i+1)-zcpz(i))/(xcpz(i+1)-xcpz(i))*
     +             (xzf(npzf)-xcpz(i))+zcpz(i)
               go to 3960
             end if
3950      continue
3960      xpl(1)=xzf(1)
          zpl(1)=zp1
          npts=1
c
          do 3930 i=1,np
             if(xcpz(i).lt.xzf(1)) go to 3930
             if(xcpz(i).ge.xzf(1).and.xcpz(i).le.xzf(npzf)) then
               if(abs(xcpz(i)-xzf(1)).gt..01) then
                 npts=npts+1
                 xpl(npts)=xcpz(i)
                 zpl(npts)=zcpz(i)
               end if
             end if
             if(xcpz(i).gt.xzf(npzf)) then
               if(abs(xcpz(i-1)-xzf(npzf)).gt..01) then
                 npts=npts+1
                 xpl(npts)=xzf(npzf)
                 zpl(npts)=zp2
               end if
               go to 3940
             end if
3930      continue
3940      continue
c
          do 3920 i=1,npts
             xcpz(i)=(xpl(i)-xmin)/xscale+orig
             zcpz(i)=(zpl(i)-zmax)/zscale+orig
3920      continue
c
          do 3970 i=1,npts-1
             do 3980 j=1,npzf-1
                if(xpl(i).ge.xzf(j).and.xpl(i).le.xzf(j+1)) ip1=j
                if(xpl(i+1).ge.xzf(j).and.xpl(i+1).le.xzf(j+1)) ip2=j
3980         continue
             zf1=(zff(ip1+1)-zff(ip1))/(xzf(ip1+1)-xzf(ip1))*
     +           (xpl(i)-xzf(ip1))+zff(ip1)
             zf2=(zff(ip2+1)-zff(ip2))/(xzf(ip2+1)-xzf(ip2))*
     +           (xpl(i+1)-xzf(ip2))+zff(ip2)
             zffp=(zf1+zf2)/2.
c
             if(zffp.lt.0.5) then
               nps=ndot
               xinc=(xcpz(i+1)-xcpz(i))/float(nps)
               zinc=(zcpz(i+1)-zcpz(i))/float(nps)
               do 3990 j=1,nps+1
                  xps=xcpz(i)+float(j-1)*xinc
                  zps=zcpz(i)+float(j-1)*zinc
                  call dot(xps,zps,symht,ifcol)
3990           continue
             end if
             if(zffp.ge.0.5.and.zffp.lt.0.8) then
               xcp(1)=xcpz(i)
               xcp(2)=xcpz(i+1)
               zcp(1)=zcpz(i)
               zcp(2)=zcpz(i+1)
               dash=xmm/dashlf
c              call gwicol(-2.,ifcol)
               call dashln(xcp,zcp,2,dash)
             end if
             if(zffp.ge.0.8) then
c              call gwicol(-2.,ifcol)
               call plot(xcpz(i),zcpz(i),3)
               call plot(xcpz(i+1),zcpz(i+1),2)
             end if
c
3970      continue
c
          go to 395
399       continue
c
          call pcolor(ifcol)
c
        end if
c
      else
c
        if(ifrefl.lt.0) go to 1900
c
        if(xtinc.lt.0.) xtinc=(xmax-xmin)/100.
        npts=int((xmax-xmin)/xtinc+.5)+1
        if(npts.gt.ppnpts) npts=ppnpts
        xtinc=(xmax-xmin-.002)/float(npts-1)
        do 100 i=1,npts
           xpts(i)=xmin+.001+float(i-1)*xtinc
           xppts(i)=(xpts(i)-xmin)/xscale+orig
           tpts(i)=0.
100     continue
c
        if(idump.eq.1) write(14,5) (xpts(j),j=1,npts)
c
        do 110 i=1,ncont-1
           do 120 j=1,npts
              xp=xpts(j)
              if(nzed(i).eq.1) then
                zu=zm(i,1)
              else
                do 220 k=1,nzed(i)-1
                   if(xp.ge.xm(i,k).and.xp.le.xm(i,k+1)) then
                     zu=(zm(i,k+1)-zm(i,k))*(xp-xm(i,k))/(xm(i,k+1)-
     +                   xm(i,k))+zm(i,k)
                     go to 221
                   end if
220             continue
              end if
221           if(nzed(i+1).eq.1) then
                zl=zm(i+1,1)
              else
                do 230 k=1,nzed(i+1)-1
                   if(xp.ge.xm(i+1,k).and.xp.le.xm(i+1,k+1)) then
                     zl=(zm(i+1,k+1)-zm(i+1,k))*(xp-xm(i+1,k))/
     +                 (xm(i+1,k+1)-xm(i+1,k))+zm(i+1,k)
                     go to 231
                   end if
230             continue
              end if
231           if(abs(zl-zu).lt..0005) go to 120
              il=i
102           if(nvel(il,1).eq.1) then
                vu=vf(il,1,1)
              else
                do 240 k=1,nvel(il,1)-1
                   if(xp.ge.xvel(il,k,1).and.xp.le.xvel(il,k+1,1))
     +             then
                     vu=(vf(il,k+1,1)-vf(il,k,1))*(xp-xvel(il,k,1))/
     +                  (xvel(il,k+1,1)-xvel(il,k,1))+vf(il,k,1)
                     go to 241
                   end if
240             continue
              end if
241           if(nvel(il,2).eq.1) then
                vl=vf(il,1,2)
              else
                do 250 k=1,nvel(il,2)-1
                   if(xp.ge.xvel(il,k,2).and.xp.le.xvel(il,k+1,2))
     +             then
                     vl=(vf(il,k+1,2)-vf(il,k,2))*(xp-xvel(il,k,2))/
     +                  (xvel(il,k+1,2)-xvel(il,k,2))+vf(il,k,2)
                     go to 251
                   end if
250             continue
              end if
251           if(i.eq.il) then
                vll=vl
                if(vu.gt.0.) go to 252
                il=il-1
                go to 102
              end if
              if(vl.gt.0.) then
                vu=vl
                go to 252
              else
                if(vu.gt.0.) go to 252
                il=il-1
                go to 102
              end if
252           if(vll.eq.0.) then
                vl=vu
              else
                vl=vll
              end if
              if(vu.eq.vl) then
                add=2.*(zl-zu)/vu
              else
                grad=(vl-vu)/(zl-zu)
                add=2.*alog(vl/vu)/grad
              end if
              tpts(j)=tpts(j)+add
120        continue
c
           if(ntsmth.gt.0) then
             do 140 j=1,ntsmth
                call smooth(tpts,npts)
140          continue
           end if
c
           do 130 j=1,npts
              tppts(j)=(tpts(j)-tmax)/tscale+orig
130        continue
c
           if(idash.ne.1) then
             if(idash.ne.-1) call line(xppts,tppts,npts)
           else
             dash=xmm/150.
             call dashln(xppts,tppts,npts,dash)
           end if
           if(idump.eq.1) write(14,5) (tpts(j),j=1,npts)
c
110     continue
c
1900    if(abs(ifrefl).eq.1) then
c
          open(unit=30, file='f.in', status='old')
c
590       read(30,545,end=595) nfrefr
545       format(i2)
          WRITE(FMT,'("(3x,", I0, "f7.2)")') nfrefr
          read(30,FMT) (xfrefl(i),i=1,nfrefr)
          read(30,FMT) (zfrefl(i),i=1,nfrefr)
          read(30,765)
C 555       format(3x,<nfrefr>f7.2)
c
          ipen=3
          do 1010 ii=1,nfrefr
             xp=xfrefl(ii)
             zp=zfrefl(ii)
             tp=0.
             do 1110 i=1,ncont-1
                if(nzed(i).eq.1) then
                  zu=zm(i,1)
                else
                  do 1220 k=1,nzed(i)-1
                     if(xp.ge.xm(i,k).and.xp.le.xm(i,k+1)) then
                       zu=(zm(i,k+1)-zm(i,k))*(xp-xm(i,k))/(xm(i,k+1)-
     +                     xm(i,k))+zm(i,k)
                       go to 1221
                     end if
1220              continue
                end if
1221            if(nzed(i+1).eq.1) then
                  zl=zm(i+1,1)
                else
                  do 1230 k=1,nzed(i+1)-1
                     if(xp.ge.xm(i+1,k).and.xp.le.xm(i+1,k+1)) then
                       zl=(zm(i+1,k+1)-zm(i+1,k))*(xp-xm(i+1,k))/
     +                   (xm(i+1,k+1)-xm(i+1,k))+zm(i+1,k)
                       go to 1231
                     end if
1230              continue
                end if
1231            if(abs(zl-zu).lt..0005) go to 1110
                if(zp.ge.zu.and.zp.lt.zl) then
                  zb=zp
                  isflag=1
                else
                  zb=zl
                  isflag=0
                end if
                il=i
1102            if(nvel(il,1).eq.1) then
                  vu=vf(il,1,1)
                else
                  do 1240 k=1,nvel(il,1)-1
                     if(xp.ge.xvel(il,k,1).and.xp.le.xvel(il,k+1,1))
     +               then
                       vu=(vf(il,k+1,1)-vf(il,k,1))*(xp-xvel(il,k,1))/
     +                    (xvel(il,k+1,1)-xvel(il,k,1))+vf(il,k,1)
                       go to 1241
                     end if
1240              continue
                end if
1241            if(nvel(il,2).eq.1) then
                  vl=vf(il,1,2)
                else
                  do 1250 k=1,nvel(il,2)-1
                     if(xp.ge.xvel(il,k,2).and.xp.le.xvel(il,k+1,2))
     +               then
                       vl=(vf(il,k+1,2)-vf(il,k,2))*(xp-xvel(il,k,2))/
     +                    (xvel(il,k+1,2)-xvel(il,k,2))+vf(il,k,2)
                       go to 1251
                     end if
1250              continue
                end if
1251            if(i.eq.il) then
                  vll=vl
                  if(vu.gt.0.) go to 1252
                  il=il-1
                  go to 1102
                end if
                if(vl.gt.0.) then
                  vu=vl
                  go to 1252
                else
                  if(vu.gt.0.) go to 1252
                  il=il-1
                  go to 1102
                end if
1252            if(vll.eq.0.) then
                  vl=vu
                else
                  vl=vll
                end if
c
                if(vu.eq.vl) then
                  add=2.*(zb-zu)/vu
                else
                  grad=(vl-vu)/(zl-zu)
                  vl=vu+(zb-zu)*grad
                  add=2.*alog(vl/vu)/grad
                end if
                tp=tp+add
                if(isflag.eq.1) go to 1111
1110         continue
c
1111         xp=(xp-xmin)/xscale+orig
             tp=(tp-tmax)/tscale+orig
             call plot(xp,tp,ipen)
             ipen=2
c
1010      continue
          go to 590
595       continue
        end if
c
        if(abs(ifrefl).eq.2) then
c
          open(unit=30, file='f.in', status='old')
c
690       read(30,645,end=695) xp,zp,ifpen
645       format(2f10.3,i10)
          if(xp.lt.xmin.or.xp.gt.xmax) go to 690
c
          tp=0.
          do 2110 i=1,ncont-1
             if(nzed(i).eq.1) then
               zu=zm(i,1)
             else
               do 2220 k=1,nzed(i)-1
                  if(xp.ge.xm(i,k).and.xp.le.xm(i,k+1)) then
                    zu=(zm(i,k+1)-zm(i,k))*(xp-xm(i,k))/(xm(i,k+1)-
     +                  xm(i,k))+zm(i,k)
                    go to 2221
                  end if
2220           continue
             end if
2221         if(nzed(i+1).eq.1) then
               zl=zm(i+1,1)
             else
               do 2230 k=1,nzed(i+1)-1
                  if(xp.ge.xm(i+1,k).and.xp.le.xm(i+1,k+1)) then
                    zl=(zm(i+1,k+1)-zm(i+1,k))*(xp-xm(i+1,k))/
     +                (xm(i+1,k+1)-xm(i+1,k))+zm(i+1,k)
                    go to 2231
                  end if
2230           continue
             end if
2231         if(abs(zl-zu).lt..0005) go to 2110
             if(zp.ge.zu.and.zp.lt.zl) then
               zb=zp
               isflag=1
             else
               zb=zl
               isflag=0
             end if
             il=i
2102         if(nvel(il,1).eq.1) then
               vu=vf(il,1,1)
             else
               do 2240 k=1,nvel(il,1)-1
                  if(xp.ge.xvel(il,k,1).and.xp.le.xvel(il,k+1,1))
     +            then
                    vu=(vf(il,k+1,1)-vf(il,k,1))*(xp-xvel(il,k,1))/
     +                 (xvel(il,k+1,1)-xvel(il,k,1))+vf(il,k,1)
                    go to 2241
                  end if
2240           continue
             end if
2241         if(nvel(il,2).eq.1) then
               vl=vf(il,1,2)
             else
               do 2250 k=1,nvel(il,2)-1
                  if(xp.ge.xvel(il,k,2).and.xp.le.xvel(il,k+1,2))
     +            then
                    vl=(vf(il,k+1,2)-vf(il,k,2))*(xp-xvel(il,k,2))/
     +                 (xvel(il,k+1,2)-xvel(il,k,2))+vf(il,k,2)
                    go to 2251
                  end if
2250           continue
             end if
2251         if(i.eq.il) then
               vll=vl
               if(vu.gt.0.) go to 2252
               il=il-1
               go to 2102
             end if
             if(vl.gt.0.) then
               vu=vl
               go to 2252
             else
               if(vu.gt.0.) go to 2252
               il=il-1
               go to 2102
             end if
2252         if(vll.eq.0.) then
               vl=vu
             else
               vl=vll
             end if
c
             if(vu.eq.vl) then
               add=2.*(zb-zu)/vu
             else
               grad=(vl-vu)/(zl-zu)
               vl=vu+(zb-zu)*grad
               add=2.*alog(vl/vu)/grad
             end if
             tp=tp+add
             if(isflag.eq.1) go to 2111
2110      continue
c
2111      xp=(xp-xmin)/xscale+orig
          tp=(tp-tmax)/tscale+orig
          if(ifpen.eq.0) then
            call dot(xp,tp,symht,ifcol)
          else
            call ssymbl(xp,tp,symht,ifpen)
          end if
c
          go to 690
695       continue
        end if
c
      end if
c
      call empty
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine pltvp(nlayer,xvp,iaxlab,izort,iroute,vtime)
c
c     plot 1-D velocity profile versus depth or time at x=xvp
c
      include 'vmodel.par'
      real v(pinvel),z(pinvel),xvp(1)
      include 'vmodel.com'
c
      call dxmin(nlayer+1)
c
      iaxis=1
      if(xvp(1).lt.-9999.) xvp(1)=(xmin+xmax)/2.
c
      do 100 k=1,pnvp
         if(xvp(k).ge.xmin.and.xvp(k).le.xmax) then
           ifvtime=0
           xp=xvp(k)
           write(13,5) xp
5          format(/'x position: ',f7.2,' km'/)
           if(iplots.eq.0) then
             call plots(xwndow,ywndow,iroute)
             call segmnt(1)
             iplots=1
             iflag1=1
           else
             iflag1=0
           end if
           if(iaxis.eq.1) then
             if(iflag1.ne.1) call aldone
             call erase
             if(iaxlab.eq.1) then
               call axtick(vmin,vmax,vtmin,vtmax,ntickv,ndeciv)
                if(izort.eq.0) then
                  call axis(orig,orig+zmm,vmin,vmax,vmm,vscale,
     +  0.,-1,vtmin,vtmax,ntickv,ndeciv,'VELOCITY (km/s)',15,albht)
                  call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
                  call axis(orig,orig,zmin,zmax,zmm,zscale,90.,1,
     +            ztmin,ztmax,ntickz,ndeciz,'DEPTH (km)',10,albht)
                else
                  call axis(orig,orig+tmm,vmin,vmax,vmm,vscale,
     +  0.,-1,vtmin,vtmax,ntickv,ndeciv,'VELOCITY (km/s)',15,albht)
                  call axtick(tmin,tmax,ttmin,ttmax,ntickt,ndecit)
                  call axis(orig,orig,tmin,tmax,tmm,tscale,90.,1,
     +            ttmin,ttmax,ntickt,ndecit,'TIME (s)',8,albht)
                end if
             end if
             if(izort.eq.0) then
               call box(orig,orig,orig+vmm,orig+zmm)
             else
               call box(orig,orig,orig+vmm,orig+tmm)
             end if
             iaxis=0
           end if
c
1000       np=0
           do 10 i=1,nlayer
              if(nzed(i).eq.1) then
                zu=zm(i,1)
              else
                do 20 j=1,nzed(i)-1
                   if(xp.ge.xm(i,j).and.xp.le.xm(i,j+1)) then
                     zu=(zm(i,j+1)-zm(i,j))*(xp-xm(i,j))/(xm(i,j+1)-
     +                   xm(i,j))+zm(i,j)
                     go to 21
                   end if
20              continue
              end if
21            if(nzed(i+1).eq.1) then
                zl=zm(i+1,1)
              else
                do 30 j=1,nzed(i+1)-1
                   if(xp.ge.xm(i+1,j).and.xp.le.xm(i+1,j+1)) then
                     zl=(zm(i+1,j+1)-zm(i+1,j))*(xp-xm(i+1,j))/
     +                   (xm(i+1,j+1)-xm(i+1,j))+zm(i+1,j)
                     go to 31
                   end if
30              continue
              end if
31            if(abs(zl-zu).lt..0005) go to 10
              np=np+1
              if(nvel(i,1).eq.1) then
                vu=vf(i,1,1)
              else
                do 40 j=1,nvel(i,1)-1
                   if(xp.ge.xvel(i,j,1).and.xp.le.xvel(i,j+1,1))
     +             then
                     vu=(vf(i,j+1,1)-vf(i,j,1))*(xp-xvel(i,j,1))/
     +                  (xvel(i,j+1,1)-xvel(i,j,1))+vf(i,j,1)
                     go to 41
                   end if
40              continue
              end if
41            if(nvel(i,2).eq.1) then
                vl=vf(i,1,2)
              else
                do 50 j=1,nvel(i,2)-1
                   if(xp.ge.xvel(i,j,2).and.xp.le.xvel(i,j+1,2))
     +             then
                     vl=(vf(i,j+1,2)-vf(i,j,2))*(xp-xvel(i,j,2))/
     +                  (xvel(i,j+1,2)-xvel(i,j,2))+vf(i,j,2)
                     go to 51
                   end if
50              continue
              end if
51            if(vu.eq.0.) vu=(v((np-1)*2)-orig)*vscale+vmin
              if(vl.eq.0.) vl=vu
              v(np*2-1)=(vu-vmin)/vscale+orig
              v(np*2)=(vl-vmin)/vscale+orig
              if(izort.eq.0) then
                z(np*2-1)=(zu-zmax)/zscale+orig
                z(np*2)=(zl-zmax)/zscale+orig
                write(13,15) i,zu,zl,vu,vl
15              format('layer# ',i2,'   z1=',f7.2,'   z2=',f7.2,
     +                 ' km'/9x,'  vp1=',f7.2,'  vp2=',f7.2,' km/s')
              else
                if(np.eq.1) then
                  z(1)=0.
                else
                  z(np*2-1)=(z((np-1)*2)-orig)*tscale+tmax
                end if
                if(vu.eq.vl) then
                  add=2.*(zl-zu)/vu
                else
                  grad=(vl-vu)/(zl-zu)
                  add=2.*alog(vl/vu)/grad
                end if
                z(np*2)=z(np*2-1)+add
                write(13,25) i,z(np*2-1),z(np*2),vu,vl
25              format('layer# ',i2,'   t1=',f7.2,'   t2=',f7.2,
     +                 ' s'/9x,'  vp1=',f7.2,'  vp2=',f7.2,' km/s')
                if(ifvtime.eq.0.and.vtime.gt.0..and.
     +            (vtime.le.vu.or.vtime.le.vl)) then
                  if(vtime.le.vu) then
                    timev=z(np*2-1)
                  else
                    if(vu.eq.vl) then
                      timev=z(np*2-1)+2.*(zvtime-zu)/vu
                    else
                      grad=(vl-vu)/(zl-zu)
                      timev=z(np*2-1)+2.*alog(vtime/vu)/grad
                    end if
                  end if
                  write(61,55) xp,timev
55                format(2f10.3)
                  ifvtime=1
                end if
                z(np*2-1)=(z(np*2-1)-tmax)/tscale+orig
                z(np*2)=(z(np*2)-tmax)/tscale+orig
              end if
c
10         continue
c
           call line(v,z,2*np)
c
         end if
100   continue
c
      call empty
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine pltrms(nlayer,xrinc,vrmstl,vrmsbl,nrsmth,iaxlab,
     +                  idump,iroute)
c
c     plot rms velocity versus distance
c
      include 'vmodel.par'
      real vprms(ppnpts),xprms(ppnpts),vrms(ppnpts),xrms(ppnpts),k,
     +     v(pinvel)
      integer vrmstl,vrmsbl
      include 'vmodel.com'
c
      call dxmin(nlayer+1)
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
        call axtick(xmint,xmax,xtmin,xtmax,ntickx,ndecix)
        call axtick(vrmin,vrmax,vrtmin,vrtmax,ntckvr,ndecir)
        call axis(orig,orig,xmin,xmax,xmm,xscale,0.,1,
     +  xtmin,xtmax,ntickx,ndecix,'DISTANCE (km)',13,albht)
        call axis(orig,orig,vrmin,vrmax,vrmm,rscale,90.,1,
     +  vrtmin,vrtmax,ntckvr,ndecir,'RMS VELOCITY (km/s)',19,albht)
      end if
      call box(orig,orig,orig+xmm,orig+vrmm)
c
      nrms=int((xmax-xmin)/xrinc+.5)+1
      if(nrms.gt.ppnpts) nrms=ppnpts
      xrinc=(xmax-xmin-.002)/float(nrms-1)
      do 100 ii=1,nrms
         xrms(ii)=xmin+.001+float(ii-1)*xrinc
         xp=xrms(ii)
         tsum=0.
         vsum=0.
         np=0
         do 10 i=1,nlayer
            if(nzed(i).eq.1) then
              zu=zm(i,1)
            else
              do 20 j=1,nzed(i)-1
                 if(xp.ge.xm(i,j).and.xp.le.xm(i,j+1)) then
                   zu=(zm(i,j+1)-zm(i,j))*(xp-xm(i,j))/(xm(i,j+1)-
     +                 xm(i,j))+zm(i,j)
                   go to 21
                 end if
20            continue
            end if
21          if(nzed(i+1).eq.1) then
              zl=zm(i+1,1)
            else
              do 30 j=1,nzed(i+1)-1
                 if(xp.ge.xm(i+1,j).and.xp.le.xm(i+1,j+1)) then
                   zl=(zm(i+1,j+1)-zm(i+1,j))*(xp-xm(i+1,j))/
     +                 (xm(i+1,j+1)-xm(i+1,j))+zm(i+1,j)
                   go to 31
                 end if
30            continue
            end if
31          if(abs(zl-zu).lt..0005) go to 10
            np=np+1
            if(nvel(i,1).eq.1) then
              vu=vf(i,1,1)
            else
              do 40 j=1,nvel(i,1)-1
                 if(xp.ge.xvel(i,j,1).and.xp.le.xvel(i,j+1,1))
     +           then
                   vu=(vf(i,j+1,1)-vf(i,j,1))*(xp-xvel(i,j,1))/
     +                (xvel(i,j+1,1)-xvel(i,j,1))+vf(i,j,1)
                   go to 41
                 end if
40            continue
            end if
41          if(nvel(i,2).eq.1) then
              vl=vf(i,1,2)
            else
              do 50 j=1,nvel(i,2)-1
                 if(xp.ge.xvel(i,j,2).and.xp.le.xvel(i,j+1,2))
     +           then
                   vl=(vf(i,j+1,2)-vf(i,j,2))*(xp-xvel(i,j,2))/
     +                (xvel(i,j+1,2)-xvel(i,j,2))+vf(i,j,2)
                   go to 51
                 end if
50            continue
            end if
51          if(vu.eq.0.) vu=v((np-1)*2)
            if(vl.eq.0.) vl=vu
            v(np*2-1)=vu
            v(np*2)=vl
            if(i.ge.vrmstl.and.i.le.vrmsbl) then
              h=zl-zu
              k=abs((vl-vu)/h)
              if(abs(k).gt.1.e-6) then
                tsum=tsum+alog(1.+k*h/vu)/k
              else
                tsum=tsum+h/vu
              end if
              vsum=vsum+vu*h+0.5*k*h**2
            end if
10       continue
c
c
         xprms(ii)=(xrms(ii)-xmin)/xscale+orig
         vrms(ii)=sqrt(vsum/tsum)
         vprms(ii)=(vrms(ii)-vrmin)/rscale+orig
100   continue
c
      if(nrsmth.gt.0) then
        do 130 i=1,nrsmth
           call smooth(vprms,nrms)
130     continue
      end if
c
      call line(xprms,vprms,nrms)
c
      call empty
c
      if(idump.eq.1) then
        write(14,5) (xrms(i),i=1,nrms)
        write(14,5) (vrms(i),i=1,nrms)
5       format(10(50f7.2))
      end if
      if(idump.eq.2) then
        write(14,15) (xrms(i),vrms(i),i=1,nrms)
15      format(2f7.2)
      end if
c
      return
      end
