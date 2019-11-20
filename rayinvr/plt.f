c
c     version 1.3  Aug 1992
c
c     Plotting routines for RAYINVR
c
c     ----------------------------------------------------------------
c
      subroutine pltmod(ncont,ibnd,imod,iaxlab,ivel,velht,idash,
     +                  ifrbnd,idata,iroute,i33)
c
c     plot the velocity model
c
      include 'rayinvr.par'
      real xcp(pnsmth),zcp(pnsmth)
      include 'rayinvr.com'
c
      if(iplots.eq.0) then
        if((isep.eq.0.or.isep.eq.2).and.(itx.gt.0.or.idata.ne.0))
     +    sep=sep+orig+tmm
        call plots(xwndow,ywndow,iroute)
        call segmnt(1)
        iplots=1
      end if
      call erase
      ititle=0
      if(iaxlab.eq.1) then
        call axtick(xmin,xmax,xtmin,xtmax,ntickx,ndecix)
        call axtick(zmin,zmax,ztmin,ztmax,ntickz,ndeciz)
c       call axis(orig,sep+zmm,xmin,xmax,xmm,xscale,0.,-1,
c    +       xtmin,xtmax,ntickx,ndecix,'DISTANCE (m)',12,albht)
c       call axis(orig,sep,zmin,zmax,zmm,zscale,90.,1,
c    +       ztmin,ztmax,ntickz,ndeciz,'DEPTH (m)',9,albht)
        call axis(orig,sep+zmm,xmin,xmax,xmm,xscale,0.,-1,
     +       xtmin,xtmax,ntickx,ndecix,'DISTANCE (km)',13,albht)
        call axis(orig,sep,zmin,zmax,zmm,zscale,90.,1,
     +       ztmin,ztmax,ntickz,ndeciz,'DEPTH (km)',10,albht)
      end if
      call box(orig,sep,orig+xmm,sep+zmm)
c
      if(ititle.eq.0.and.title.ne.' ') then
        if(xtitle.lt.-999999.) xtitle=0.
        if(ytitle.lt.-999999.) ytitle=0.
        call symbol(xtitle,ytitle,albht,title,0.,80)
        ititle=1
      end if
c
      if(imod.eq.1) then
c
        call pcolor(mcol(1))
c
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
20           continue
           end if
           if(idash.ne.1) then
             call line(xcp,zcp,nptsc)
             if(i33.eq.1) write(33,335) 0.,0,0,0.,0.
335          format(f12.4,2i4,2f12.4)
           else
             dash=xmm/150.
             call dashln(xcp,zcp,nptsc,dash)
           end if
c
10      continue
c
        call pcolor(mcol(3))
c
        if(ifrbnd.eq.1) then
          do 110 i=1,nfrefl
             do 120 j=1,npfref(i)
                xcp(j)=(xfrefl(i,j)-xmin)/xscale+orig
                zcp(j)=(zfrefl(i,j)-zmax)/zscale+sep
120          continue
             if(idash.ne.1) then
               call line(xcp,zcp,npfref(i))
               if(i33.eq.1) write(33,335) 0.,0,0,0.,0.
             else
               dashfr=xmm/300.
               call dashln(xcp,zcp,npfref(i),dashfr)
             end if
110       continue
        end if
c
        if(ibnd.eq.1) then
c
          call pcolor(mcol(2))
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
                      if(i33.eq.1) write(33,335) 0.,0,0,0.,0.
                    else
                      call dashln(xcp,zcp,2,dash)
                    end if
                  end if
40             continue
             end if
30        continue
        end if
c
        if(ibsmth.eq.2) then
c
          call pcolor(mcol(4))
c
          do 600 i=1,npbnd
             x=xmin+float(i-1)*xsinc
             xcp(i)=(x-xmin)/xscale+orig
600       continue
          do 610 i=1,nlayer+1
             do 620 j=1,npbnd
                zcp(j)=(cosmth(i,j)-zmax)/zscale+sep
620          continue
             call line(xcp,zcp,npbnd)
             if(i33.eq.1) write(33,335) 0.,0,0,0.,0.
610       continue
        end if
      end if
c
      if(ivel.ne.0) then
c
        call pcolor(mcol(5))
c
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
      call pcolor(ifcol)
c
      call empty
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine plttx(ifam,npts,iszero,idata,iaxlab,xshot,idr,nshot,
     +   itxout,ibrka,ivraya,ttunc,itrev,xshotr,fidrr,itxbox,iroute,
     +   iline)
c
c     plot the travel time curves for all rays reaching the surface -
c     for each group of rays a separate curve is drawn for each
c     ray code
c
      include 'rayinvr.par'
      real xh(pnrayf),th(pnrayf),rh(pnrayf),xs(pnrayf),f(pnrayf),
     +     x(pnrayf),t(pnrayf),xshot(1),xsh(pnrayf),fh(pnrayf)
      character tlab*12
      integer npts(1),idr(1),ibrka(1),ivraya(1)
      include 'rayinvr.com'
c
      if(itx.gt.0.or.idata.ne.0) then
        ipflag=1
      else
        ipflag=0
      end if
      xcshot=-999999.
      fc=0.
      dashtl=tmm/100.
      if(vred.ne.0.) then
        rvred=1./vred
      else
        rvred=0.
      end if
c
      if(iaxlab.eq.1) then
        if(vred.eq.0.) then
          nchart=8
          tlab='TIME (s)'
c         nchart=9
c         tlab='TIME (ms)'
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
      end if
c
      if(itrev.ne.1) then
        tadj=tmin
      else
        tadj=tmax
      end if
c
      xshoth=-99999.
      fidh=0.
      iflag1=0
      n=0
      ib=0
      xbmin=orig
      xbmax=orig+xmmt
      tbmin=orig
      tbmax=orig+tmm
c
      do 10 i=1,ifam
c
         nh1=n
c
         if(npts(i).gt.0) then
           do 20 j=1,npts(i)
              xh(j)=range(j+nh1)
              th(j)=tt(j+nh1)
              rh(j)=rayid(j+nh1)
              xsh(j)=xshtar(j+nh1)
              fh(j)=fidarr(j+nh1)
              n=n+1
20         continue
           if(isep.eq.2.and.(abs(xsh(1)-xshotr).gt..001.or.
     +       fh(1).ne.fidrr)) go to 10
           if(iflag1.eq.0) then
             iflaga=1
           else
             if(xsh(1).ne.xshota.and.isep.gt.1) then
               iflaga=1
             else
               iflaga=0
             end if
           end if
           xshota=xsh(1)
           ida=nint(fh(1))
c
           if(iflaga.eq.1) then
             if(ipflag.eq.0) go to 1010
             if(iplots.eq.0) then
               call plots(xwndow,ywndow,iroute)
               call segmnt(1)
               iplots=1
             end if
             if((isep.lt.2.or.isep.eq.3).and.iflag1.ne.0) then
                call empty
                call aldone
             end if
             if((isep.lt.2.or.isep.eq.3).and.(iflag1.ne.0.or.isep.eq.1.
     +       or.isep.eq.3)) then
               call erase
               ititle=0
             end if
c
             if(iaxlab.eq.1) then
               call axtick(xmint,xmaxt,xtmint,xtmaxt,ntckxt,ndecxt)
               call axtick(tmin,tmax,ttmin,ttmax,ntickt,ndecit)
c
c              call axis(orig,orig,xmint,xmaxt,xmmt,xscalt,0.,1,
c    +           xtmint,xtmaxt,ntckxt,ndecxt,'OFFSET (m)',10,albht)
               call axis(orig,orig,xmint,xmaxt,xmmt,xscalt,0.,1,
     +           xtmint,xtmaxt,ntckxt,ndecxt,'DISTANCE (km)',13,albht)
               call axis(orig,orig,tmin,tmax,tmm,tscale,90.,1,
     +           ttmin,ttmax,ntickt,ndecit,tlab,nchart,albht)
             end if
             call box(orig,orig,orig+xmmt,orig+tmm)
c
             if(ititle.eq.0.and.title.ne.' ') then
               if(xtitle.lt.-999999.) xtitle=0.
               if(ytitle.lt.-999999.) ytitle=0.
               call symbol(xtitle,ytitle,albht,title,0.,80)
               ititle=1
             end if
c
             if(itx.eq.4) then
               ntlmax=int(tmax-.001)
               do 310 j=1,ntlmax
                  tp=(float(j)-tadj)/tscale+orig
                  x(1)=orig
                  x(2)=orig+xmmt
                  t(1)=tp
                  t(2)=tp
                  call dashln(x,t,2,dashtl)
310            continue
             end if
c
             if(idata.ne.0)
     +         call pltdat(iszero,idata,xshot,idr,nshot,tadj,xshota,
     +                     xbmin,xbmax,tbmin,tbmax,itxbox,ida)
1010         if((itx.ge.3.or.itxout.eq.3).and.narinv.gt.0) then
               ichold=-999999
               do 120 j=1,narinv
                  if((abs(xscalc(j)-xshota).lt..0001.and.icalc(j)*fh(1).
     +            gt.0..and.isep.eq.2).or.(abs(xscalc(j)-xshota).lt.
     +            .0001.and.isep.eq.3).or.isep.lt.2) then
                    if(itx.ge.3) then
                      xplot=(xcalc(j)-xmint)/xscalt+orig
                      if(itx.eq.3) then
                        tplot=(tcalc(j)-tadj)/tscale+orig
                      else
                        tplot=(tcalc(j)-tobs(j)+float(abs(icalc(j)))-
     +                          tadj)/tscale+orig
                      end if
                      if(itxbox.eq.0.or.(xplot.ge.xbmin.and.xplot.le.
     +                xbmax.and.tplot.ge.tbmin.and.tplot.le.tbmax)) then
                        if(itcol.eq.3) then
                          ipcol=colour(mod(ircalc(j)-1,ncol)+1)
                          call pcolor(ipcol)
                        end if
                        if(iline.eq.0) then
c                         call ssymbl(xplot,tplot,symht,4)
c                         call ssymbl(xplot,tplot,symht,3)
                          call dot(xplot,tplot,symht,ifcol)
                        else
                          if(ircalc(j).eq.ichold) then
                            call plot(xplot,tplot,2)
                          else
                            call plot(xplot,tplot,3)
                            ichold=ircalc(j)
                          end if
                        end if
                      end if
                    end if
                    if(itxout.eq.3) then
                      fcalc=sign(1.,float(icalc(j)))
                      if(xscalc(j).ne.xcshot.or.fcalc.ne.fc) then
                        if(iszero.ne.1) then
                          write(17,5) xscalc(j),fcalc,0.,0
                          xsr=xscalc(j)
                        else
                          write(17,5) 0.,1.,0.,0
                          xsr=0.
                        end if
                      end if
                      xcshot=xscalc(j)
                      fc=fcalc
                      twrite=tcalc(j)+abs(xcalc(j)-xsr)*rvred
c                     write(17,5) xcalc(j),twrite,ttunc,abs(icalc(j))
                      write(17,5) xcalc(j),twrite,uobs(j),abs(icalc(j))
                    end if
                  end if
120            continue
               if(itx.ge.3.and.itcol.eq.3) call pcolor(ifcol)
             end if
           end if
c
           iflag1=1
c
           if(itx.ne.1.and.itx.ne.2.and.itxout.ne.1.and.itxout.ne.2)
     +     go to 10
c
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
             if(ibrka(i).eq.0.or.itx.eq.2.or.npt.lt.3) then
               if(itxout.eq.1.or.itxout.eq.2) then
                 if(xs(1).ne.xshoth.or.f(1).ne.fidh) then
                   ib=0
                   if(iszero.ne.1) then
                     write(17,5) xs(1),f(1),0.,0
                     xsr=xs(1)
                   else
                     write(17,5) 0.,1.,0.,0
                     xsr=0.
                   end if
                   xshoth=xs(1)
                   fidh=f(1)
                 end if
                 if(itxout.eq.1) then
                   ib=ib+1
                 else
                   ib=ivraya(i)
                 end if
                 do 40 j=1,npt
                    xw=(x(j)-orig)*xscalt+xmint
                    tw=(t(j)-orig)*tscale+tadj+abs(xsr-xw)*rvred
                    write(17,5) xw,tw,ttunc,ib
5                   format(3f10.3,i10)
40               continue
               end if
               if(npt.gt.1) then
                 if(itx.ne.2) then
                   if(itx.lt.3.and.itxbox.eq.0) then
                     call line(x,t,npt)
                   else
                     if(x(1).ge.xbmin.and.x(1).le.xbmax.and.t(1).ge.
     +               tbmin.and.t(1).le.tbmax) then
                       if(itx.lt.3) call plot(x(1),t(1),3)
                       iout=0
                     else
                       iout=1
                     end if
                     do 42 j=2,npt
                        ipen=3
                        if(x(j).ge.xbmin.and.x(j).le.xbmax.and.t(j).ge.
     +                  tbmin.and.t(j).le.tbmax) then
                          if(iout.eq.0) ipen=2
                          iout=0
                        else
                          iout=1
                        end if
                        if(itx.lt.3) call plot(x(j),t(j),ipen)
42                   continue
                   end if
                 else
                   if(itx.lt.3.and.symht.gt.0.) then
                     do 41 j=1,npt
                        if(itxbox.eq.0.or.(x(j).ge.xbmin.and.x(j).le.
     +                  xbmax.and.t(j).ge.tbmin.and.t(j).le.tbmax))
     +                  call ssymbl(x(j),t(j),symht,4)
41                   continue
                   end if
                 end if
               else
                 if((itxbox.eq.0.or.(x(1).ge.xbmin.and.x(1).le.xbmax.
     +           and.t(1).ge.tbmin.and.t(1).le.tbmax)).and.symht.gt.0..
     +           and.itx.lt.3) call ssymbl(x(1),t(1),symht,4)
               end if
             else
               ih=0
               if(itxout.eq.1.or.itxout.eq.2) then
                 if(xs(1).ne.xshoth.or.f(1).ne.fidh) then
                   ib=0
                   if(iszero.ne.1) then
                     write(17,5) xs(1),f(1),0.,0
                     xsr=xs(1)
                   else
                     write(17,5) 0.,1.,0.,0
                     xsr=0.
                   end if
                   xshoth=xs(1)
                   fidh=f(1)
                 end if
                 if(itxout.eq.1) then
                   ib=ib+1
                 else
                   ib=ivraya(i)
                 end if
                 xw=(x(1)-orig)*xscalt+xmint
                 tw=(t(1)-orig)*tscale+tadj+abs(xsr-xw)*rvred
                 write(17,5) xw,tw,ttunc,ib
                 xw=(x(2)-orig)*xscalt+xmint
                 tw=(t(2)-orig)*tscale+tadj+abs(xsr-xw)*rvred
                 write(17,5) xw,tw,ttunc,ib
               end if
               if(itxbox.eq.0) then
                 if(itx.lt.3) then
                   call plot(x(1),t(1),3)
                   call plot(x(2),t(2),2)
                 end if
               else
                 if(x(1).ge.xbmin.and.x(1).le.xbmax.and.t(1).ge.tbmin.
     +           and.t(1).le.tbmax.and.x(2).ge.xbmin.and.x(2).le.
     +           xbmax.and.t(2).ge.tbmin.and.t(2).le.tbmax) then
                   if(itx.lt.3) then
                     call plot(x(1),t(1),3)
                     call plot(x(2),t(2),2)
                   end if
                   iout=0
                 else
                   if(x(2).ge.xbmin.and.x(2).le.xbmax.and.t(2).ge.tbmin.
     +             and.t(2).le.tbmax) then
                     iout=0
                   else
                     iout=1
                   end if
                 end if
               end if
               do 30 j=3,npt
                  ipen=3
                  if(((x(j)-x(j-1))*(x(j-1)-x(j-2))).ge.0.or.
     +              ih.eq.(j-1)) then
                    if(itxout.eq.1.or.itxout.eq.2) then
                      xw=(x(j)-orig)*xscalt+xmint
                      tw=(t(j)-orig)*tscale+tadj+abs(xsr-xw)*rvred
                      write(17,5) xw,tw,ttunc,ib
                    end if
                    if(itxbox.eq.0.or.(x(j).ge.xbmin.and.x(j).le.xbmax.
     +              and.t(j).ge.tbmin.and.t(j).le.tbmax.and.iout.eq.0))
     +                ipen=2
                  else
                    if(itxout.eq.1) then
                      ib=ib+1
                    else
                      ib=ivraya(i)
                    end if
                    if(itxout.eq.1.or.itxout.eq.2) then
                      xw=(x(j)-orig)*xscalt+xmint
                      tw=(t(j)-orig)*tscale+tadj+abs(xsr-xw)*rvred
                      write(17,5) xw,tw,ttunc,ib
                    end if
                    ih=j
                  end if
                  if(itx.lt.3) call plot(x(j),t(j),ipen)
                  if(itxbox.ne.0) then
                    if(x(j).ge.xbmin.and.x(j).le.xbmax.and.t(j).ge.
     +              tbmin.and.t(j).le.tbmax) then
                      iout=0
                    else
                      iout=1
                    end if
                  end if
30             continue
               if((itxbox.eq.0.or.(x(npt).ge.xbmin.and.x(npt).le.
     +         xbmax.and.t(npt).ge.tbmin.and.t(npt).le.tbmax)).and.
     +         ih.eq.npt.and.symht.gt.0..and.itx.lt.3)
     +           call ssymbl(x(npt),t(npt),symht,4)
             end if
             go to 1000
           end if
         end if
10    continue
c
      t0=131.
      vrms=.95
      do 101 i=1,19
         xpl=float(i-1)*10.
         tpl=sqrt(t0**2+(xpl/vrms)**2)
c        write(0,*) xpl,tpl
         xppl=(xpl-xmint)/xscalt+orig
         tppl=(tpl-tadj)/tscale+orig
         if(i.eq.1) then
           ipen=3
         else
           ipen=2
         end if
c        call plot(xppl,tppl,ipen)
101   continue
c
      if(ipflag.eq.1) call empty
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine pltdat(iszero,idata,xshot,idr,nshot,tadj,xshota,
     +                  xbmin,xbmax,tbmin,tbmax,itxbox,ida)
c
c     plot observed travel times
c
      include 'rayinvr.par'
      real xshot(1)
      integer idr(1)
      include 'rayinvr.com'
c
      npick=0
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
           xdiff=abs(xshot(i)-xshota)
           if((xdiff.lt..0001.and.ida.eq.idr(i).and.
     +     (imod.gt.0.or.iray.gt.0)).or.
     +     (xdiff.lt..0001.and.imod*iray.eq.0).or.isep.lt.2) then
             xdiff=abs(xshot(i)-xsp)
             if(xdiff.lt..0001.and.idr(i).eq.idp) then
               iplt=1
               npick=isf-nsfc
               isf=isf+1
               nsfc=nsfc+1
               icshot=i
               go to 100
             end if
           end if
10      continue
        iplt=0
        nsfc=nsfc+1
        isf=ilshot(nsfc)
        go to 100
      else
        if(abs(idata).eq.2) then
          npick=npick+1
          iflag=0
          do 20 i=1,narinv
             if(ipinv(i).eq.npick) then
               iflag=1
               go to 30
             end if
20        continue
        else
          iflag=1
        end if
30      if(iplt.eq.1.and.iflag.eq.1) then
          if(iszero.eq.0) then
            xplot=(xp-xmint)/xscalt+orig
          else
            xplot=((xp-xsp)*float(idp)-xmint)/xscalt+orig
          end if
          if(itx.ne.4) then
            tplot=(tp-tadj)/tscale+orig
          else
            tplot=(float(abs(ip))-tadj)/tscale+orig
          end if
          if(itcol.eq.1.or.itcol.eq.2) then
            if(itcol.ne.2) then
              ipcol=colour(mod(ip-1,ncol)+1)
            else
              ipcol=colour(mod(icshot-1,ncol)+1)
            end if
            call pcolor(ipcol)
          end if
          if(itxbox.eq.0.or.(xplot.ge.xbmin.and.xplot.le.xbmax.
     +       and.tplot.ge.tbmin.and.tplot.le.tbmax)) then
            if(idata.gt.0) then
              call plot(xplot,tplot-up/tscale,3)
              call plot(xplot,tplot+up/tscale,2)
            else
c             call ssymbl(xplot,tplot,symht,1)
              call dot(xplot,tplot,symht,ifcol)
            end if
          end if
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
c     plot one ray path
c
      include 'rayinvr.par'
      real x(ppray),z(ppray),xa(ppray),za(ppray),vra(ppray),vpa(ppray)
      character reply*1
      include 'rayinvr.com'
c
      npts=npt-nskip
      if(npts.lt.2) return
c
      do 180 i=1,npts
         x(i)=(xr(i+nskip)-xmin)/xscale+orig
         z(i)=(zr(i+nskip)-zmax)/zscale+sep
         vra(i)=vr(i+nskip,2)
         vpa(i)=vp(i+nskip,2)
180   continue
c
      if(npskp.ne.1) then
        do 190 i=1,nbnd
           if(nbnda(i).gt.nskip+1) then
             nbnds=i
             go to 130
           end if
190     continue
130     nbnd=nbnd-nbnds+1
        do 140 i=1,nbnd
           nbnda(i)=nbnda(i+nbnds-1)-nskip
140     continue
        npts=1
        np1=2
        do 110 j=1,nbnd
           if(nbnda(j).gt.np1.and.nbnda(j)-np1.gt.npskp) then
             do 120 i=np1,nbnda(j)-1,npskp
                npts=npts+1
                x(npts)=x(i)
                z(npts)=z(i)
                vra(npts)=vra(i)
                vpa(npts)=vpa(i)
120          continue
           end if
           npts=npts+1
           x(npts)=x(nbnda(j))
           z(npts)=z(nbnda(j))
           vra(npts)=vra(nbnda(j))
           vpa(npts)=vpa(nbnda(j))
           np1=nbnda(j)+1
110     continue
      end if
c
      if(istep.eq.1) write(6,5) anglew
5     format('take-off angle: ',f10.5$)
c
      if(ircol.ne.0) call pcolor(irrcol)
c
      if(idot.eq.2) go to 100
      if(irayps.ne.1) then
55      format(a3,i10,2f10.3)
        call line(x,z,npts)
      else
        dash=xmm/250.
        na=0
        if(vra(1).eq.vpa(1)) then
          ips=1
        else
          ips=-1
        end if
        do 30 i=1,npts-1
           if(vra(i).eq.vpa(i)) then
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
100   if(idot.gt.0) then
        do 20 i=1,npts
           call ssymbl(x(i),z(i),symht,4)
20      continue
      end if
c
      if(ircol.ne.0) call pcolor(ifcol)
c
      if(istep.eq.1) then
        call empty
        read(5,15) reply
15      format(a1)
        if(reply(1:1).eq.'s') then
          call plotnd(1)
          stop
        end if
        if(reply(1:1).eq.'0') istep=0
      end if
c
      return
      end
