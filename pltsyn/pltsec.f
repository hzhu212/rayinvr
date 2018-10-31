c
c     version 1.2  Mar 1992
c
c     Pltsec routine for PLTSYN
c
c     ----------------------------------------------------------------
c
      subroutine pltsec(vred,xomit,nskip,ipol,tol,nptsw,itrev,idump,
     +                  iamp,twin,imeth,iroute)
c
c     plot synthetic seismic sections using the time, amplitude,
c     and phase of each arrival read in from unit 11
c
      include 'pltsyn.par'
c
      real dist(pseis),s(ppseis),r(ppseis+121),time(ppseis),
     +     sect(pseis,parriv,3),max,wavlet(pwavlt),xamp(parriv*pseis),
     +     tamp(parriv*pseis),
     +     hil(121),xomit(1),sp(ppseis),tp(ppseis),maxamp,minamp
      integer na(pseis),ipamp(parriv*pseis)
      character tlab*12
c
      include 'pltsyn.com'
c
      data wavlet/.01,.072,.121,.12,.044,-.113,-.32,-.512,-.607,
     +           -.538,-.283,.112,.542,.875,1.0,.875,.542,.112,
     +           -.283,-.538,-.607,-.512,-.32,-.113,.044,.12,
     +            .121,.072,.01,71*0./
      nseg=0
      p121=ppseis+121
      pi2=6.283185307
      dens=dens/2.
      fifi=float(ifill)
      tol=tol/xscale
      npts=nint((tmax-tmin)*sps)+1
      nptsp=npts+121
      nomit=0
      if(vred.ne.0.) then
        rvred=1./vred
      else
        rvred=0.
      end if
      nwin=nint(twin*sps)
      do 40 i=1,pseis
         if(xomit(i).lt.-99998.) go to 50
         nomit=nomit+1
40    continue
50    do 80 i=1,npts
         if(itrev.ne.1) then
           time(i)=(float(i-1)/sps)/tscale+orig
         else
           time(i)=(tmin+float(i-1)/sps-tmax)/tscale+orig
         end if
80    continue
      if(ishade.eq.1) then
        ptmm=abs(tscale)*sps
        if(ptmm.lt.dens) then
          nvaplt=nint(dens*float(npts)/ptmm)
        else
          nvaplt=npts
        end if
        if(nvaplt.gt.ppseis) then
          write(6,565)
565       format(/'***  interpolated trace for shading too long  ***'/)
          stop
        end if
        if(nvaplt.gt.npts) call intrtr(time,npts,nvaplt)
      end if
      do 90 i=2,60,2
         fudge=float(i)*.00018604-.00037207
         hil(61+i-1)=.6366198/float(i-1)-fudge
         hil(61-i+1)=-hil(61+i-1)
90    continue
      if(iconv.eq.1) then
        if(iwavlt.lt.2) then
          if(iwavlt.eq.0) then
            nwave=29
            nwave1=30
          else
            nwave=nptsw
            nwave1=nwave+1
            pf=pi2/float(nptsw-1)
            do 430 i=1,nptsw
               wavlet(i)=sin(pf*float(i-1))
430         continue
          end if
        else
          nwave1=1
450       read(20,451,end=440) wavlet(nwave1)
451       format(f10.3)
          nwave1=nwave1+1
          go to 450
440       nwave=nwave1-1
          if(nwave.eq.0) then
            write(6,505)
505         format(/'***  w.in is empty  ***'/)
            stop
          end if
        end if
        wsum=0.
        do 445 i=1,nwave
           wsum=wsum+wavlet(i)
445     continue
        wmean=wsum/float(nwave)
        do 455 i=1,nwave
           wavlet(i)=wavlet(i)-wmean
455     continue
        wmax=0.
        do 465 i=1,nwave
           if(abs(wavlet(i)).gt.wmax) wmax=abs(wavlet(i))
465     continue
        if(wmax.gt..000001) then
          do 475 i=1,nwave
             wavlet(i)=wavlet(i)/wmax
475       continue
        else
          write(6,405)
405       format(/'***  source wavelet is zero  ***'/)
          stop
        end if
        if(abs(wavlet(1)).lt..001) wavlet(1)=sign(.001,wavlet(1))
        if(abs(wavlet(nwave)).lt..001)
     +                wavlet(nwave)=sign(.001,wavlet(nwave))
        if(ipol.eq.-1) then
          do 130 i=1,nwave
             wavlet(i)=-wavlet(i)
130       continue
        end if
        twave=float(nwave-1)/sps
      end if
      if(itrev.ne.1) then
        yorig=orig
        iside=1
      else
        yorig=orig+tmm
        iside=-1
      end if
c
      if(iplot.ge.-1) then
        call plots(xwndow,ywndow,iroute)
        iplots=1
      end if
      read(11,5,end=999) xshot,idum
c
1000  continue
c
      maxamp=0
      ntplt=0
      nsexsp=0
      nhpts=npts
      if(iplot.ge.-1) call erase
      if(iplot.ge.-1.and.iseg.eq.1) call segmnt(1)
      if(iplot.ge.-1) then
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
      call axtick(xmin,xmax,xtmin,xtmax,ntickx,ndecix)
      call axtick(tmin,tmax,ttmin,ttmax,ntickt,ndecit)
      call axis(orig,yorig,xmin,xmax,xmm,xscale,0.,iside,
     +          xtmin,xtmax,ntickx,ndecix,'DISTANCE (km)',13,albht)
      call axis(orig,orig,tmin,tmax,tmm,tscale,90.,1,
     +          ttmin,ttmax,ntickt,ndecit,tlab,nchart,albht)
      call box(orig,orig,orig+xmm,orig+tmm)
      end if
c
      nseis=0
c
180   read(11,5,end=991) dists,nas
5     format(f10.3,i10)
      if(nas.lt.0) then
        xshotn=dists
        go to 99
      end if
      nseis=nseis+1
      dist(nseis)=dists
      na(nseis)=nas
      if(na(nseis).gt.parriv) then
        write(6,555) dist(nseis)
555     format(/
     +  '***  max # of arrivals on seismogram at',f10.3,
     +  ' km exceeded ***'/)
        stop
      end if
      if(na(nseis).gt.0) then
        if(iscale.eq.2) then
          range=abs(dist(nseis)-xshot)
          if(range.eq.0.) range=.001
          rr=range**rcor
        end if
        do 190 j=1,na(nseis)
           read(11,15,end=999) (sect(nseis,j,k),k=1,3)
15         format(3e12.5)
           if(iscale.eq.2) then
             samp=abs(sect(nseis,j,2)*rr)
             if(samp.gt.maxamp) maxamp=samp
           end if
           if(vred.gt.0.) then
             sect(nseis,j,1)=
     +            sect(nseis,j,1)-abs(dist(nseis)-xshot)/vred
           end if
190     continue
      end if
      go to 180
c
991   xshotn=-999999.
99    if(iscale.eq.2) then
        scalef=amp/maxamp
        write(6,25) scalef*xnorm**rcor
25      format(/'***  scalef=',f15.5,'  ***')
      end if
c
      if(iamp.gt.0) then
        namp=0
        rewind(17)
107     read(17,45) xr,tr,ur,ir
45      format(3f10.3,i10)
        if(ir.le.0) then
          if(ir.lt.0) go to 108
          if(abs(xshot-xr).lt..001) then
            iflags=1
          else
            iflags=0
          end if
        else
          if(iflags.eq.1) then
            namp=namp+1
            xamp(namp)=xr
            tamp(namp)=tr-abs(xshot-xr)*rvred
            ipamp(namp)=ir
          end if
        end if
        go to 107
108     if(namp.eq.0) write(*,*, fmt="(/
     +                '***  error in file tx.out  ***'/)")
        if(namp.gt.0) write(18,45) xshot,0.,0.,0
      end if
      idir=1
c
      do 70 i=1,nseis
         npts=nhpts
         if(mod(i,nskip).ne.0) go to 70
         if(dist(i).lt.xmin.or.dist(i).gt.xmax) go to 70
         if(nomit.gt.0) then
           do 75 j=1,nomit
              if(abs(xomit(j)-dist(i)).lt..001) go to 70
75         continue
         end if
         nsexsp=nsexsp+1
         do 170 j=1,nptsp
            r(j)=0.
170      continue
         if(na(i).gt.0) then
           do 270 j=1,na(i)
              if(iconv.eq.1.and.((sect(i,j,1)+twave).lt.tmin.or.
     +           sect(i,j,1).gt.tmax)) go to 270
              nstart=nint((sect(i,j,1)-tmin)*sps)
              sn=sin(sect(i,j,3))*sect(i,j,2)
              kpos=nstart+61
              if(kpos.gt.0.and.kpos.le.nptsp)
     +          r(kpos)=r(kpos)+cos(sect(i,j,3))*sect(i,j,2)
              do 370 k=2,120,2
                 kpos=nstart+k
                 if(kpos.gt.0.and.kpos.le.nptsp)
     +             r(kpos)=r(kpos)-hil(k)*sn
370           continue
270        continue
         end if
         if(nsmth.gt.0) then
           do 280 j=1,nsmth
              call smooth(r,nptsp)
280        continue
         end if
         if(iconv.eq.1) then
           do 470 j=1,npts
              s(j)=0.
              if(j.le.(nwave-60)) then
                kmax=j-59
              else
                kmax=nwave
              end if
              do 570 k=1,kmax
                 s(j)=s(j)+wavlet(k)*r(j-k+60)
570           continue
470        continue
         else
           do 460 j=1,npts
              s(j)=r(j+60)
460        continue
         end if
c
         if(inmo.eq.1) call nmo(s,npts,sps,tmin,dist(i),vrms)
c
         if(idump.eq.1) then
           write(12,35) (s(j),j=1,npts)
35         format(36(100a4))
         end if
c
         if(namp.gt.0) then
           do 109 j=1,namp
              if(abs(xamp(j)-dist(i)).lt..001) then
                n0=nint((tamp(j)-tmin)*sps)+1
                n1=n0+nwin
                if(n0.lt.1) n0=1
                if(n1.lt.1) n1=1
                if(n0.gt.npts) n0=npts
                if(n1.gt.npts) n1=npts
                nwina=n1-n0+1
                minamp=0.
                maxamp=0.
c
                if(n1.gt.n0) then
                  if(imeth.ne.1) then
                    do 360 jj=n0,n1
                       if(s(jj).lt.minamp) minamp=s(jj)
                       if(s(jj).gt.maxamp) maxamp=s(jj)
360                 continue
                    maxamp=(abs(minamp)+abs(maxamp))/2.
                  else
                    sum=0.
                    do 361 jj=n0,n1
                       sum=sum+abs(s(jj))
361                 continue
                    maxamp=sum/float(nwina)
                  end if
                  if(maxamp.gt.0.) maxamp=alog10(maxamp)
                  write(18,45) xamp(j),maxamp,0.,ipamp(j)
                end if
              end if
109        continue
         end if
c
         if(iplot.lt.-1) go to 70
         if(iscale.eq.0) then
           max=0.
           do 480 j=1,npts
              if(abs(s(j)).gt.max) max=abs(s(j))
480        continue
           if(max.gt.0.) then
             f=amp/max
             do 490 j=1,npts
                s(j)=f*s(j)
490          continue
           end if
         end if
         if(iscale.gt.0) then
           range=abs(dist(i)-xshot)
           if(range.eq.0.) range=.001
         end if
         if(iscale.eq.1) then
           srn=scalef*(range/xnorm)**rcor
           do 580 j=1,npts
              s(j)=srn*s(j)
580        continue
         end if
         if(iscale.eq.2) then
           sr=scalef*range**rcor
           do 610 j=1,npts
              s(j)=sr*s(j)
610        continue
         end if
         if(clip.gt.0..and.iscale.eq.1) then
           do 590 j=1,npts
              if(abs(s(j)).gt.clip) s(j)=sign(clip,s(j))
590        continue
         end if
         tmean=(dist(i)-xmin)/xscale+orig
         do 600 j=1,npts
            s(j)=s(j)/xscale+tmean
600      continue
         if(ishade.eq.1) then
           if(nvaplt.gt.npts) call intrtr(s,npts,nvaplt)
           npts=nvaplt
         end if
         if(tol.gt.0..and.npts.gt.2) then
           nplt=1
           sp(nplt)=s(1)
           tp(nplt)=time(1)
           do 700 j=2,npts-1
              if(abs(s(j)-tmean).ge.tol) then
                nplt=nplt+1
                sp(nplt)=s(j)
                tp(nplt)=time(j)
              end if
700        continue
           nplt=nplt+1
           sp(nplt)=s(npts)
           tp(nplt)=time(npts)
         else
           nplt=npts
           do 710 j=1,nplt
              sp(j)=s(j)
              tp(j)=time(j)
710        continue
         end if
         ntplt=ntplt+nplt
         if(idir.eq.1) then
           nl1=1
           nl2=2
           nl3=nplt
           istep=1
         else
           nl1=nplt
           nl2=nplt-1
           nl3=1
           istep=-1
         end if
         call plot(sp(nl1),tp(nl1),3)
         do 620 j=nl2,nl3,istep
            if(ishade.eq.1) then
              if((fifi*(sp(j)-tmean)).ge.0..and.
     +           (fifi*(sp(j-istep)-tmean)).ge.0.)
     +           call plot(tmean,(tp(j)+tp(j-istep))/2.,2)
                 ntplt=ntplt+1
            end if
            call plot(sp(j),tp(j),2)
620      continue
c
         call empty
c
         idir=-idir
70    continue
c
      npts=nhpts
      if(iplot.ge.-1.and.nsexsp.gt.0.and.tol.gt.0.) then
        per=float(ntplt)/float(nsexsp*npts)*100.
        write(6,105) per,tol*xscale
105     format(/'***  ',f6.2,' % of points'
     +         ,' plotted using spmin=',e9.3,'  ***')
      end if
c
      write(6,665) xshot
665   format(/'>>>  shot at ',f10.3,' km completed  <<<')
      if(xshotn.gt.-999998.) then
        xshot=xshotn
        if(iplot.ge.-1.and.iseg.eq.1) call segmnt(0)
        if(iplot.ge.-1.and.iroute.eq.1) call aldone
        go to 1000
      end if
c
      if(namp.gt.0) write(18,45) 0.,0.,0.,-1
c
      return
c
999   write(6,995)
995   format(/'***  insufficient data in sect.out  ***'/)
c
      return
      end
