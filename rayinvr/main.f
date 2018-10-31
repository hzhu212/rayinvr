c
c     version 1.4  Apr 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ********   R A Y I N V R   ********               |
c     |                                                              |
c     |            Two-dimensional Ray Tracing Program               |
c     |           for Traveltime Modeling and Inversion              |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
c     I/O units:
c
c        10 -- input:  program input parameters
c
c        11 -- output: summary ray tracing information for each ray
c                      traced
c
c        12 -- output: details of velocity model and each point
c                      of each ray traced
c
c        17 -- output: calculated traveltime-distance pairs
c
c        18 -- output: partial derivative and traveltime residual
c                      arrays for inversion
c
c        19 -- output: all Calcomp plot calls
c
c        20 -- input:  velocity model
c
c        22 -- output: namelist parameters
c
c        26 -- input:  observed traveltime-distance pairs
c
c        30 -- input:  floating reflectors
c
c        31 -- output: velocity model sampled on a uniform grid
c
c        32 -- output: floating reflectors
c
c        33 -- output: summary info for each ray traced for input
c                      to RAYPLOT to plot minimum traveltime raypaths
c
c        34 -- output: summary info for each observation used for input
c                      to RAYPLOT to plot minimum traveltime raypaths
c
c
c     ----------------------------------------------------------------
c
c
      program main
c
      include 'rayinvr.par'
c
      real amin(prayt),amax(prayt),xshot(pshot),space(prayf),
     +     zshot(pshot),pois(player),poisbl(papois),parunc(3),
     +     zsmth(pnsmth),xshota(pshot2),zshota(pshot2),ta2pt(pr2pt),
     +     ra2pt(pr2pt),tt2pt(pr2pt),xo2pt(pnobsf),ho2pt(pr2pt),
     +     tatan(pitan2),angtan(pitan2)
c
      integer idr(pshot2),nray(prayf),itt(ptrayf),ibrka(ptrayf),
     +        ishot(pshot),ncbnd(prayf),ivraya(ptrayf),ishotr(pshot2),
     +        nrbnd(prayf),rbnd(preflt),nsmax(prayf),ishotw(pshot2),
     +        iturn(prayf),cbnd(pconvt),irayt(prayt),ihead(player),
     +        poisl(papois),poisb(papois),ibreak(prayf),frbnd(prayf),
     +        ifo2pt(pnobsf),ipos(pr2pt),modi(player),nsmin(prayf),
     +        insmth(pncntr)
      character flag*1,title*80
c
      include 'rayinvr.com'
c
      namelist /pltpar/ iplot,imod,ibnd,idash,ivel,iray,irays,
     +                  irayps,idot,itx,idata,iszero,itxout,istep,
     +                  ibreak,idump,isum,symht,velht,nskip,nrskip,
     +                  vred,dvmax,dsmax,isep,iroute,itcol,ircol,mcol,
     +                  itxbox,iseg,xwndow,ywndow,colour,ibcol,ifcol,
     +                  modout,dxmod,dzmod,modi,frz,xmmin,xmmax,
     +                  xmin1d,xmax1d,npskip,ifd,dxzmod,
     +                  title,xtitle,ytitle,iline
c
      namelist /axepar/ iaxlab,itrev,albht,orig,sep,
     +                  xmin,xmax,xtmin,xtmax,xmm,ntickx,ndecix,
     +                  xmint,xmaxt,xtmint,xtmaxt,xmmt,ntckxt,ndecxt,
     +                  zmin,zmax,ztmin,ztmax,zmm,ntickz,ndeciz,
     +                  tmin,tmax,ttmin,ttmax,tmm,ntickt,ndecit
c
      namelist /trapar/ ishot,iraysl,irayt,iturn,isrch,istop,ibsmth,
     +                  imodf,xshot,zshot,ray,nray,space,amin,amax,
     +                  nsmax,aamin,aamax,stol,xsmax,step,smin,smax,
     +                  nrbnd,rbnd,ncbnd,cbnd,pois,poisb,poisl,poisbl,
     +                  crit,hws,npbnd,nbsmth,nhray,frbnd,nsmin,
     +                  i2pt,n2pt,x2pt,ifast,ntan,idiff,xminns,xmaxns,
     +                  insmth,aainc
c
      namelist /invpar/ invr,ivray,ximax,ttunc,velunc,bndunc
c
c     initialize parameters
c
      data amin/prayt*181./,amax/prayt*181./,
     +     space/prayf*0./,zshot/pshot*-9999./,
     +     pois/player*-99./,poisbl/papois*-99./
      data nray/prayf*-1/,itt/ptrayf*0/,modi/player*0/,
     +     ishot/pshot*0/,ncbnd/prayf*-1/,nsmin/prayf*-1/,
     +     nrbnd/prayf*0/,rbnd/preflt*0/,frbnd/prayf*0/,
     +     nsmax/prayf*-1/,iturn/prayf*1/,irayt/prayt*1/,
     +     cbnd/pconvt*-99/,ibreak/prayf*1/,ihead/player*0/,
     +     insmth/pncntr*0/
c
      iline=0
      ititle=0
      xtitle=-9999999.
      ytitle=-9999999.
      title=' '
      ifd=0
      xminns=-999999.
      xmaxns=-999999.
      npskip=1
      idiff=0
      ntan=90
      ifast=1
      xmin1d=-999999.
      xmax1d=-999999.
      istep=0
      iroute=1
      xmmin=-999999.
      xmmax=-999999.
      frz=0.
      dxmod=0.
      dzmod=0.
      dxzmod=0.
      modout=0
      i2pt=0
      n2pt=5
      x2pt=-1.
      itxbox=0
      nhray=-1
      ircol=0
      itcol=0
      isep=0
      velunc=0.1
      bndunc=0.1
      dvmax=1.e10
      dsmax=1.e10
      itrev=0
      xsmax=0.
      ntray=0
      ntpts=0
      isrch=0
      istop=1
      nstepr=15
      narinv=0
      nrskip=1
      itxout=0
      idata=0
      isum=1
      idash=1
      ivel=0
      iaxlab=1
      idot=0
      iszero=0
      ibnd=1
      imod=1
      iray=1
      irayps=0
      nskip=0
      invr=0
      imodf=0
      irays=0
      iraysl=0
      ncont=1
      iflagi=0
      nrayl=0
      ifam=0
      ictbnd=0
      irbnd=0
      nshot=0
      ngroup=0
      ia0=1
      ttunc=.01
      stol=-1.
      velht=.7
      aamin=5.
      aamax=85.
      aainc=.1
      aaimin=1.
      ximax=-1.
      i33=0
c
      open(unit=10, file='r.in', status='old')
c
c     read in program control parameters from unit 10
c
      read(10,pltpar)
      read(10,axepar)
      read(10,trapar)
      read(10,invpar)
c
      if(xmax.lt.-99998) then
        write(6,25)
25      format(/'***  xmax not specified  ***'/)
        stop
      end if
c
      if(imodf.ne.1) then
        iunit=10
        read(10,1)
1       format(' '/' '/' ')
      else
        open(unit=20, file='v.in', status='old')
        iunit=20
      end if
c
      if(mod(ppcntr,10).ne.0.or.mod(ppvel,10).ne.0) then
        write(6,1105)
1105    format(/
     +  '***  array size error for number of model points  ***'/)
        go to 9999
      end if
c
      nrzmax=ppcntr/10
      nrvmax=ppvel/10
      do 20 icont=1,player+1
         nrz=1
         j1=1
         j2=10
11       if(nrz.gt.nrzmax) go to 211
         read(iunit,15,end=999) ilyr,(xm(icont,j),j=j1,j2)
         read(iunit,15,end=999) icnt,(zm(icont,j),j=j1,j2)
         read(iunit,235,end=99) (ivarz(icont,j),j=j1,j2)
15       format(i2,1x,10f7.2)
c235      format(3x,10i7)
235      format(3x,10(5x,i2))
         nrz=nrz+1
         if(icnt.ne.1) go to 211
         j1=j1+10
         j2=j2+10
         go to 11
211      nrv=1
         j1=1
         j2=10
21       if(nrv.gt.nrvmax) go to 311
         read(iunit,15,end=999) ilyr,(xvel(icont,j,1),j=j1,j2)
         read(iunit,15,end=999) icnt,(vf(icont,j,1),j=j1,j2)
         read(iunit,235,end=999) (ivarv(icont,j,1),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 311
         j1=j1+10
         j2=j2+10
         go to 21
311      nrv=1
         j1=1
         j2=10
31       if(nrv.gt.nrvmax) go to 411
         read(iunit,15,end=999) ilyr,(xvel(icont,j,2),j=j1,j2)
         read(iunit,15,end=999) icnt,(vf(icont,j,2),j=j1,j2)
         read(iunit,235,end=999) (ivarv(icont,j,2),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 411
         j1=j1+10
         j2=j2+10
         go to 31
411      ncont=ncont+1
20    continue
c
99    nlayer=ncont-1
c
c     open I/O units
c
      if(iplot.eq.0) iplot=-1
      if(iplot.eq.2) iplot=0
      open(unit=11, file='r1.out')
      if(idump.eq.1) open(unit=12, file='r2.out')
      if(itxout.gt.0) open(unit=17, file='tx.out')
      if(iplot.le.0) open(unit=19, file='p.out')
      if(abs(modout).ne.0) open(unit=31, file='v.out')
      if(i2pt.gt.0.and.iray.eq.2.and.irays.eq.0.and.irayps.eq.0.
     +and.iplot.le.0) then
        open(unit=33, file='ra1.out')
        open(unit=34, file='ra2.out')
        i33=1
      end if
      if(ifd.gt.0) open(unit=35, file='fd.out')

c
c     read in observed data
c
      if(idata.ne.0.or.invr.eq.1.or.i2pt.ne.0) then
c
        open(unit=26, file='tx.in', status='old')
c
        isf=1
        nsf=0
        xshotc=0.
910     read(26,930) xpf(isf),tpf(isf),upf(isf),ipf(isf)
930     format(3f10.3,i10)
        if(ipf(isf).le.0) then
          nsf=nsf+1
          ilshot(nsf)=isf
          xshotc=xpf(isf)
        else
          if(vred.ne.0)
     +      tpf(isf)=tpf(isf)-abs(xshotc-xpf(isf))/vred
        end if
        isf=isf+1
        if(ipf(isf-1).ne.-1) go to 910
      end if
c
c     determine if plot parameters for distance axis of travel
c     time plot are same as model distance parameters
c
      if(xmint.lt.-1000.) xmint=xmin
      if(xmaxt.lt.-1000.) xmaxt=xmax
      if(xtmint.lt.-1000.) xtmint=xtmin
      if(xtmaxt.lt.-1000.) xtmaxt=xtmax
      if(xmmt.lt.-1000.) xmmt=xmm
      if(ndecxt.lt.-1) ndecxt=ndecix
      if(ntckxt.lt.0) ntckxt=ntickx
c
c     calculate scale of each plot axis
c
      xscale=(xmax-xmin)/xmm
      xscalt=(xmaxt-xmint)/xmmt
      zscale=-(zmax-zmin)/zmm
      tscale=(tmax-tmin)/tmm
c
      if(itrev.eq.1) tscale=-tscale
      if(iroute.ne.1) then
c       isep=0
        ibcol=0
      end if
      if(isep.eq.2.and.imod.eq.0.and.iray.eq.0.and.irays.eq.0) isep=3
      if(n2pt.gt.pn2pt) n2pt=pn2pt
      if(i2pt.eq.0) x2pt=0.
      if(x2pt.lt.0.) x2pt=(xmax-xmin)/2000.
c
      ncol=0
      do 1020 i=pcol,1,-1
         if(colour(i).ge.0) then
           ncol=i
           go to 1030
         end if
1020  continue
1030  if(ncol.eq.0) then
        colour(1)=2
        colour(2)=3
        colour(3)=4
        colour(4)=5
        colour(5)=6
        colour(6)=8
        colour(7)=17
        colour(8)=27
        colour(9)=22
        colour(10)=7
        ncol=10
      end if
      do 1040 i=1,5
         if(mcol(i).lt.0) mcol(i)=ifcol
1040  continue
c
      do 40 i=1,prayf
         if(ray(i).lt.1.) go to 50
         ngroup=ngroup+1
40    continue
50    if(ngroup.eq.0) then
        write(6,55)
55      format(/'***  no ray codes specified  ***'/)
        go to 900
      end if
c
      ifrbnd=0
      do 540 i=1,ngroup
         if(frbnd(i).lt.0.or.frbnd(i).gt.pfrefl) then
           write(6,565)
565        format(/'***  error in array frbnd  ***'/)
           stop
         end if
         if(frbnd(i).gt.0) ifrbnd=1
540   continue
c
      if(ifrbnd.eq.1) then
        open(unit=30, file='f.in', status='old')
        nfrefl=0
590     read(30,545,end=595) nfrefr
        if(nfrefr.lt.2.or.nfrefr.gt.ppfref) then
          write(6,585)
585       format(/'***  error in f.in file  ***'/)
          stop
        end if
        nfrefl=nfrefl+1
        npfref(nfrefl)=nfrefr
545     format(i2)
        read(30,555) (xfrefl(nfrefl,i),i=1,npfref(nfrefl))
        read(30,555) (zfrefl(nfrefl,i),i=1,npfref(nfrefl))
        read(30,575) (ivarf(nfrefl,i),i=1,npfref(nfrefl))
555     format(3x,<npfref(nfrefl)>f7.2)
575     format(3x,<npfref(nfrefl)>i7)
        go to 590
595     close(30)
        do 550 i=1,ngroup
           if(frbnd(i).gt.nfrefl) then
             write(6,565)
             stop
           end if
550     continue
      end if
c
c     calculate velocity model parameters
c
      call calmod(ncont,pois,poisb,poisl,poisbl,invr,iflagm,ifrbnd,
     +            xmin1d,xmax1d,insmth,xminns,xmaxns)
c
      if(abs(modout).ne.0.or.ifd.gt.1) then
c
        if(xmmin.lt.-999998) xmmin=xmin
        if(xmmax.lt.-999998) xmmax=xmax
c
        if(dxmod.eq.0.) dxmod=(xmmax-xmmin)/20.
        if(dzmod.eq.0.) dzmod=(zmax-zmin)/10.
        if(dxzmod.eq.0.) dxzmod=(zmax-zmin)/10.
c
        nx=nint((xmmax-xmmin)/dxmod)
        nz=nint((zmax-zmin)/dzmod)
c
        if(abs(modout).ge.2) then
          do 4010 i=1,nz
             do 4020 j=1,nx
                sample(i,j)=0
4020         continue
4010      continue
        end if
c
      end if
c
      if(itx.ge.3.and.invr.eq.0) itx=2
      if(itxout.eq.3.and.invr.eq.0) itxout=2
      if(invr.eq.0.and.abs(idata).eq.2) idata=sign(1,idata)
      if(iflagm.eq.1) go to 9999
c
      if(idvmax.gt.0.or.idsmax.gt.0) then
        write(6,895)
        write(11,895)
        if(idvmax.gt.0) then
          write(6,865) (ldvmax(i),i=1,idvmax)
          write(11,865) (ldvmax(i),i=1,idvmax)
865       format('large velocity gradient in layers: ',100i4)
        end if
        if(idsmax.gt.0) then
          write(6,875) (ldsmax(i),i=1,idsmax)
          write(11,875) (ldsmax(i),i=1,idsmax)
875       format('large slope change in boundaries:  ',100i4)
        end if
        write(6,895)
        write(11,895)
      end if
c
c     assign values to the arrays xshota, zshota and idr
c     and determine nshot
c
      do 290 i=1,pshot
         if(ishot(i).eq.-1.or.ishot(i).eq.2) then
           nshot=nshot+1
           xshota(nshot)=xshot(i)
           zshota(nshot)=zshot(i)
           idr(nshot)=-1
           ishotw(nshot)=-i
           ishotr(nshot)=2*i-1
         end if
         if(ishot(i).eq.1.or.ishot(i).eq.2) then
           nshot=nshot+1
           xshota(nshot)=xshot(i)
           zshota(nshot)=zshot(i)
           idr(nshot)=1
           ishotw(nshot)=i
           ishotr(nshot)=2*i
         end if
290   continue
c
c     calculate the z coordinate of shot points if not specified
c     by the user - assumed to be at the top of the first layer
c
      zshift=abs(zmax-zmin)/10000.
      do 300 i=1,nshot
        if(zshota(i).lt.-1000.) then
         if(xshota(i).lt.xbnd(1,1,1).or.xshota(i).gt.xbnd(1,nblk(1),2))
     +   go to 300
         do 310 j=1,nblk(1)
            if(xshota(i).ge.xbnd(1,j,1).and.xshota(i).le.xbnd(1,j,2))
     +      zshota(i)=s(1,j,1)*xshota(i)+b(1,j,1)+zshift
310      continue
        end if
300   continue
c
c     assign default value to nray if not specified or
c     nray(1) if only it is specified and also ensure that nray<=pnrayf
c
      if(nray(1).lt.0) then
         do 320 i=1,prayf
            nray(i)=10
320      continue
      else
        if(nray(2).lt.0) then
          if(nray(1).gt.pnrayf) nray(1)=pnrayf
          do 330 i=2,prayf
             nray(i)=nray(1)
330       continue
        else
          do 340 i=1,prayf
             if(nray(i).lt.0) nray(i)=10
             if(nray(i).gt.pnrayf) nray(i)=pnrayf
340       continue
        end if
      end if
c
c     assign default value to stol if not specified by the user
c
      if(stol.lt.0.) stol=(xmax-xmin)/3500.
c
c     check array ncbnd for array values greater than pconv
c
      do 470 i=1,prayf
         if(ncbnd(i).gt.pconv) then
           write(6,135)
135        format(/'***  max converting boundaries exceeded  ***/')
           go to 900
         end if
470   continue
c
c     plot velocity model
c
      if((imod.eq.1.or.iray.gt.0.or.irays.eq.1).and.isep.lt.2)
     + call pltmod(ncont,ibnd,imod,iaxlab,ivel,velht,idash,ifrbnd,
     + idata,iroute,i33)
c
c     calculation of smooth layer boundaries
c
      if(ibsmth.gt.0) then
        do 680 i=1,nlayer+1
           zsmth(1)=(cosmth(i,2)-cosmth(i,1))/xsinc
           do 660 j=2,npbnd-1
              zsmth(j)=(cosmth(i,j+1)-cosmth(i,j-1))/(2.*xsinc)
660        continue
           zsmth(npbnd)=(cosmth(i,npbnd)-cosmth(i,npbnd-1))/xsinc
           do 670 j=1,npbnd
              cosmth(i,j)=atan(zsmth(j))
670        continue
680     continue
      end if
      if(isep.gt.1.and.ibsmth.eq.2) ibsmth=1
c
      write(11,35)
35    format('shot  ray i.angle  f.angle   dist     depth',1x,
     +       'red.time  npts code')
      if(idump.eq.1) write(12,45)
45    format(/'gr ray npt   x       z      ang1    ang2    v1 ',1x,
     +       '   v2  lyr bk id iw')
c
      if(nrskip.lt.1) nrskip=1
      if(hws.lt.0.) hws=(xmax-xmin)/25.
      hwsm=hws
      crit=crit/pi18
c
      do 260 i=1,prayf
         if(nrbnd(i).gt.prefl) then
           write(6,125)
125        format(/'***  max reflecting boundaries exceeded  ***'/)
           go to 900
         end if
260   continue
c
      do 710 i=1,preflt
        if(rbnd(i).ge.nlayer) then
         write(6,165)
165   format(/'***  reflect boundary greater than # of layers  ***'/)
        end if
710   continue
c
      if(nsmax(1).lt.0) then
        do 350 i=1,ngroup
           nsmax(i)=10
350     continue
      else
        if(nsmax(2).lt.0.and.ngroup.gt.1) then
          do 360 i=2,ngroup
             nsmax(i)=nsmax(1)
360       continue
        end if
      end if
c
      if(nsmin(1).lt.0) then
        do 351 i=1,ngroup
           nsmin(i)=1000000.
351     continue
      else
        if(nsmin(2).lt.0.and.ngroup.gt.1) then
          do 361 i=2,ngroup
             nsmin(i)=nsmin(1)
361       continue
        end if
      end if
c
c     assign default values to smin and smax if not specified
c
      if(smin.lt.0.) smin=(xmax-xmin)/4500.
      if(smax.lt.0.) smax=(xmax-xmin)/15.
c
      if(ximax.lt.0.) ximax=(xmax-xmin)/20.
      ist=0
      iflagp=0
c
      if(ifast.eq.1) then
        if(ntan.gt.pitan) ntan=pitan
        ntan=2*ntan
        ainc=pi/float(ntan)
        factan=float(ntan)/pi
c
        do 6010 i=1,ntan
           angtan(i)=float(i-1)*pi/float(ntan)
           tatan(i)=tan(angtan(i))
6010    continue
        do 6020 i=1,ntan-1
           mtan(i)=(tatan(i+1)-tatan(i))/ainc
           btan(i)=tatan(i)-mtan(i)*angtan(i)
6020    continue
c
        do 6030 i=2,ntan
           tatan(i)=1./tan(angtan(i))
6030    continue
        do 6040 i=1,ntan-1
           mcotan(i)=(tatan(i+1)-tatan(i))/ainc
           bcotan(i)=tatan(i)-mcotan(i)*angtan(i)
6040    continue
      end if
c
      if(nshot.eq.0) go to 1000
c
      do 60 is=1,nshot
         ist=ist+1
         id=idr(is)
         fid=float(id)
         xshotr=xshota(is)
         zshotr=zshota(is)
         if(ist.eq.1) then
           xsec=xshotr
           zsec=zshotr
           idsec=id
           ics=1
           do 810 i=1,nlayer
              tang(i,1)=999.
              tang(i,2)=999.
              tang(i,3)=999.
              tang(i,4)=999.
810        continue
         else
           if(abs(xshotr-xsec).lt..001.and.abs(zshotr-zsec).lt..001)
     +       then
             if(iflags.eq.1) go to 60
             ics=0
             if(id.ne.idsec) then
               do 830 i=1,nlayer
                  tang(i,1)=999.
                  tang(i,2)=999.
                  tang(i,3)=999.
                  tang(i,4)=999.
830            continue
             end if
           else
             xsec=xshotr
             zsec=zshotr
             idsec=id
             ics=1
             do 820 i=1,nlayer
                tang(i,1)=999.
                tang(i,2)=999.
                tang(i,3)=999.
                tang(i,4)=999.
820          continue
           end if
         end if
         if(ics.eq.1) then
c
           call xzpt(xshotr,zshotr,layer1,iblk1,iflags)
c
           if(iflags.eq.1) then
             write(11,95)
95       format('***  location of shot point outside model  ***')
             go to 60
           end if
c
           if((imod.eq.1.or.iray.gt.0.or.irays.eq.1).and.isep.gt.1) then
             if(iflagp.eq.1) call aldone
             iflagp=1
             call pltmod(ncont,ibnd,imod,iaxlab,ivel,velht,idash,ifrbnd,
     +                   idata,iroute,i33)
           end if
c
         end if
         irbnd=0
         ictbnd=0
c
         do 70 i=1,ngroup
            if(iraysl.eq.1) then
              irpos=(ishotr(is)-1)*ngroup+i
              if(irayt(irpos).eq.0) then
                irbnd=irbnd+nrbnd(i)
                ictbnd=ictbnd+ncbnd(i)
                nrayr=0
                go to 70
              end if
            end if
            id=idr(is)
            fid=float(id)
            fid1=fid
            do 250 j=1,prefl+1
               refll(j)=0
250         continue
            do 251 j=1,pconv+1
               icbnd(j)=-1
251         continue
            ifam=ifam+1
            nrayr=nray(i)
            iflagl=0
            ibrka(ifam)=ibreak(i)
            ivraya(ifam)=ivray(i)
            do 870 j=1,nlayer
               iheadf(j)=ihead(j)
870         continue
            idl=int(ray(i))
            idt=int((ray(i)-float(idl))*10.+.5)
            if(idt.eq.2) then
              if(nrbnd(i).gt.(prefl-1)) then
                write(6,125)
                write(11,65) ishotw(is),ray(i)
                nrayr=0
                go to 69
              end if
              refll(1)=idl
              if(nrbnd(i).gt.0) then
                do 230 j=1,nrbnd(i)
                   irbnd=irbnd+1
                   refll(j+1)=rbnd(irbnd)
230             continue
              end if
            else
              if(nrbnd(i).gt.0) then
                do 240 j=1,nrbnd(i)
                   irbnd=irbnd+1
                   refll(j)=rbnd(irbnd)
240             continue
              end if
            end if
            if(idt.eq.3) then
              iheadf(idl)=1
              ihdwf=1
              ihdwm=1
              hws=hwsm
              tdhw=0.
              if(nhray.lt.0) then
                nrayr=pnrayf
              else
                nrayr=nhray
              end if
            else
              ihdwf=-1
              ihdwm=-1
            end if
            if(ncbnd(i).gt.0) then
              do 630 j=1,ncbnd(i)
                 ictbnd=ictbnd+1
                 icbnd(j)=cbnd(ictbnd)
630           continue
            end if
            idifff=0
c
            if(ircol.eq.1) irrcol=colour(mod(ivray(i)-1,ncol)+1)
            if(ircol.eq.2) irrcol=colour(mod(is-1,ncol)+1)
            if(ircol.eq.3) irrcol=colour(mod(i-1,ncol)+1)
            if(ircol.lt.0) irrcol=-ircol
c
            if(nrayr.le.0) go to 69
c
            if(i2pt.gt.0) then
              iflag2=0
              nsfc=1
              isf=ilshot(nsfc)
1110          xf=xpf(isf)
              tf=tpf(isf)
              uf=upf(isf)
              irayf=ipf(isf)
              if(irayf.lt.0) go to 1200
              if(irayf.eq.0) then
                xshotf=xf
                idf=sign(1.,tf)
                if(abs(xshotr-xshotf).lt..001.and.idr(is).eq.idf) then
                  i2flag=1
                  isf=isf+1
                else
                  i2flag=0
                  nsfc=nsfc+1
                  isf=ilshot(nsfc)
                end if
              else
                if(i2flag.eq.1.and.ivray(i).eq.irayf) then
                  iflag2=1
                  go to 1200
                end if
                isf=isf+1
              end if
              go to 1110
1200          if(iflag2.eq.0) then
                nrayr=0
                go to 69
              end if
            end if
c
            call auto(xshotr,zshotr,i,ifam,idl,idt,aminr,amaxr,
     +         aamin,aamax,layer1,iblk1,aainc,aaimin,nsmax(i),
     +         iflag,iturn(i),amin(ia0),amax(ia0),ia0,stol,
     +         irays,nskip,idot,irayps,xsmax,istep,nsmin(i))
c
            if(iflag.ne.0) then
              write(11,65) ishotw(is),ray(i)
65            format('***  shot#',i4,' ray code',f5.1,
     +               ' no rays traced  ***')
              nrayr=0
              nrayl=nrayl+1
              iflagl=1
              go to 69
            end if
            if(amaxr.eq.aminr.and.ihdwf.ne.1) then
              if(nrayr.gt.1) then
                write(11,665) ishotw(is),ray(i)
665             format('***  shot#',i4,' ray code',f5.1,
     +                 ' 1 ray traced  ***')
                nrayl=nrayl+1
                iflagl=1
              end if
              nrayr=1
            end if
            if(nrayr.gt.1) then
              if(idt.eq.2) then
                if(amaxr.le.aminr) amaxr=90.
              end if
              if(space(i).gt.0.) then
                pinc=space(i)
              else
                if(idt.eq.2) then
                  pinc=2.
                else
                  pinc=1.
                end if
              end if
              ainc=(amaxr-aminr)/float(nrayr-1)**pinc
            else
              pinc=1.
              ainc=0.
            end if
            iend=0
            ninv=0
            ifcbnd=frbnd(i)
            nc2pt=0
c
            if(i2pt.gt.0.and.nrayr.gt.1) then
              ii2pt=i2pt
              ni2pt=1
              no2pt=0
              nco2pt=0
              ic2pt=0
c
              nsfc=1
              isf=ilshot(nsfc)
1100          xf=xpf(isf)
              tf=tpf(isf)
              uf=upf(isf)
              irayf=ipf(isf)
              if(irayf.lt.0) go to 1199
              if(irayf.eq.0) then
                xshotf=xf
                idf=sign(1.,tf)
                if(abs(xshotr-xshotf).lt..001.and.idr(is).eq.idf) then
                  i2flag=1
                  isf=isf+1
                else
                  i2flag=0
                  nsfc=nsfc+1
                  isf=ilshot(nsfc)
                end if
              else
                if(i2flag.eq.1.and.ivray(i).eq.irayf) then
                  no2pt=no2pt+1
                  if(no2pt.gt.min(pnrayf,pnobsf)) then
                    write(6,896)
896                 format('***  pnrayf or pnobsf exceeded  ***')
                    go to 1199
                  end if
                  xo2pt(no2pt)=xf
                  ifo2pt(no2pt)=0
                end if
                isf=isf+1
              end if
              go to 1100
1199          if(no2pt.eq.0) ni2pt=0
            else
              ni2pt=0
              ii2pt=0
            end if
c
91          ir=0
            nrg=0
            ihdwf=ihdwm
            tdhw=0.
            dhw=0.
            i1ray=1
90          ir=ir+1
            if(ir.gt.nrayr.and.ni2pt.le.1) go to 890
               if(i2pt.eq.0.and.iend.eq.1) go to 890
               ircbnd=1
               iccbnd=1
               iwave=1
               ihdw=0
               if(icbnd(1).eq.0) then
                 iwave=-iwave
                 iccbnd=2
               end if
               nptbnd=0
               nbnd=0
               npskp=npskip
               nccbnd=0
               id=idr(is)
               fid=float(id)
               fid1=fid
               iturnt=0
               if(nc2pt.le.1) then
                 angled=aminr+ainc*float(ir-1)**pinc
                 if(amaxr.gt.aminr) then
                   if(angled.gt.amaxr) then
                     angled=amaxr
                     iend=1
                   end if
                 else
                   if(angled.lt.amaxr) then
                     angled=amaxr
                     iend=1
                   end if
                 end if
                 if(ir.eq.nrayr.and.nrayr.gt.1) angled=amaxr
               else
891              nco2pt=nco2pt+1
                 if(nco2pt.gt.no2pt) go to 890
                 if(ifo2pt(nco2pt).ne.0.and.ii2pt.gt.0) go to 891
                 xobs=xo2pt(nco2pt)
                 tt2min=1.e10
                 do 892 j=1,nc2pt-1
                    if((ra2pt(j).ge.xobs.and.ra2pt(j+1).le.xobs).or.
     +              (ra2pt(j).le.xobs.and.ra2pt(j+1).ge.xobs)) then
                      denom=ra2pt(j+1)-ra2pt(j)
                      if(denom.ne.0.) then
                        tpos=(tt2pt(j+1)-tt2pt(j))/denom*
     +                       (xobs-ra2pt(j))+tt2pt(j)
                      else
                        tpos=(tt2pt(j+1)+tt2pt(j))/2.
                      end if
                      if(tpos.lt.tt2min) then
                        tt2min=tpos
                        if(denom.ne.0.) then
                          aort=(ta2pt(j+1)-ta2pt(j))/denom*
     +                           (xobs-ra2pt(j))+ta2pt(j)
                        else
                          aort=(ta2pt(j+1)+ta2pt(j))/2.
                        end if
                        if(ihdwf.ne.1) then
                          angled=aort
                        else
                          tdhw=aort
                          hws=tdhw
                        end if
                        xdiff=min(abs(xobs-ra2pt(j)),
     +                            abs(xobs-ra2pt(j+1)))
                        if(xdiff.lt.x2pt) ifo2pt(nco2pt)=1
                      end if
                    end if
892              continue
                 if(tt2min.gt.1.e9.or.(ifo2pt(nco2pt).ne.0.and.
     +           ii2pt.gt.0)) go to 891
               end if
               angle=fid*(90.-angled)/pi18
               if(ir.gt.1.and.ihdwf.ne.1) then
                 if(angle.eq.am) go to 90
               end if
               am=angle
               if((fid1*angle).lt.0.) then
                 id=-id
                 fid=float(id)
               end if
               layer=layer1
               iblk=iblk1
               npt=1
               xr(1)=xshotr
               zr(1)=zshotr
               ar(1,1)=0.0
               ar(1,2)=angle
               vr(1,1)=0.0
               vp(1,1)=0.0
               vs(1,1)=0.0
               vp(1,2)=vel(xshotr,zshotr)
               vs(1,2)=vp(1,2)*vsvp(layer1,iblk1)
               if(iwave.eq.1) then
                 vr(1,2)=vp(1,2)
               else
                 vr(1,2)=vs(1,2)
               end if
               idray(1)=layer1
               idray(2)=1
c
               if(ii2pt.gt.0) then
                 irs=0
               else
                 irs=ir
               end if
               if(invr.eq.1.and.irs.gt.0) then
                 ninv=ninv+1
                 do 80 j=1,nvar
                    fpart(ninv,j)=0.
80               continue
               end if
               nrg=nrg+1
               nhskip=0
c
               call trace(npt,ifam,irs,iturnt,invr,xsmax,iflag,idl,idt,
     +                    iray,ii2pt,i1ray,modout)
c
               call ttime(ishotw(is),xshotr,npt,irs,angled,ifam,itt,
     +                    iszero,iflag,uf,irayf)
c
c              if(irs.eq.0.and.vr(npt,2).gt.0.) then
               if(irs.eq.0) then
                 ic2pt=ic2pt+1
                 if(ihdwf.ne.1) then
                   ta2pt(ic2pt)=angled
                 else
                   ta2pt(ic2pt)=tdhw-hws
                 end if
                 ra2pt(ic2pt)=xr(npt)
                 tt2pt(ic2pt)=timer
                 if(vr(npt,2).le.0.) tt2pt(ic2pt)=1.e20
               end if
c
               if(((iray.eq.1.or.(iray.eq.2.and.vr(npt,2).gt.0.)).
     +         and.mod(ir-1,nrskip).eq.0.and.irs.gt.0).or.
     +         (irays.eq.1.and.irs.eq.0)) then
                 call pltray(npt,max(nskip,nhskip),idot,irayps,istep,
     +                       angled)
                 if(i33.eq.1) then
                   if(iszero.eq.1) then
                     xwr=abs(xshtar(ntt-1)-xobs)
                   else
                     xwr=xobs
                   end if
                   write(33,335) xshtar(ntt-1),ivraya(ifam),i,xwr,
     +                           tt(ntt-1)
335                format(f12.4,2i4,3f12.4)
                 end if
               end if
c
               if(vr(npt,2).gt.0.and.irs.gt.0.and.abs(modout).ge.2)
     +             call cells(npt,xmmin,dxmod,dzmod)
c
               if(invr.eq.1.and.irs.gt.0) call fxtinv(npt)
c
               if(ihdwf.eq.0) go to 890
               if(ntt.gt.pray) then
                 write(6,995)
995              format(/
     +     '***  max number of rays reaching surface exceeded  ***'/)
                 go to 890
               end if
            go to 90
c
890         continue
c
            if(ii2pt.gt.0) then
              nc2pt=ic2pt
              if(ni2pt.gt.1) then
                call sort3(ta2pt,ipos,nc2pt)
                do 893 j=1,nc2pt
                   ho2pt(j)=ra2pt(j)
893             continue
                do 894 j=1,nc2pt
                   ra2pt(j)=ho2pt(ipos(j))
894             continue
                do 897 j=1,nc2pt
                   ho2pt(j)=tt2pt(j)
897             continue
                do 898 j=1,nc2pt
                   tt2pt(j)=ho2pt(ipos(j))
898             continue
              end if
              ni2pt=ni2pt+1
              if(ni2pt.gt.n2pt) ii2pt=0
              nco2pt=0
              go to 91
            end if
c
            if(ninv.gt.0) call calprt(xshotr,i,ivraya(ifam),
     +                    idr(is),ximax,iflagw,iszero,x2pt)
            if(iflagw.eq.1) iflagi=1
c
            if(iray.gt.0.or.irays.eq.1) call empty
c
            nrayr=nrg
69          if(iflagl.eq.1) then
              flag='*'
            else
              flag=' '
            end if
            write(6,375) ishotw(is),ray(i),nrayr,flag
375         format('shot#',i4,':   ray code',f5.1,':   ',
     +             i3,' rays traced ',a1)
c
            if(ntt.gt.pray) go to 1000
c
70       continue
c
         if(isep.eq.2.and.((itx.gt.0.and.ntt.gt.1).or.idata.ne.0.or.
     +     itxout.gt.0)) call plttx(ifam,itt,iszero,idata,iaxlab,
     +     xshota,idr,nshot,itxout,ibrka,ivraya,ttunc,itrev,xshotr,
     +     float(idr(is)),itxbox,iroute,iline)
c
60    continue
c
1000  if((isep.lt.2.or.isep.eq.3).and.((itx.gt.0.and.ntt.gt.1).or.
     +  idata.ne.0.or.itxout.gt.0)) then
        if(isep.gt.0.and.iplots.eq.1) call aldone
        call plttx(ifam,itt,iszero,idata,iaxlab,xshota,idr,nshot,
     +  itxout,ibrka,ivraya,ttunc,itrev,xshotr,1.,itxbox,iroute,
     +  iline)
      end if
c
      if(itxout.gt.0) write(17,930) 0.,0.,0.,-1
c
      if(irkc.eq.1) then
        write(6,75)
        write(11,75)
75      format(/'***  possible inaccuracies in rngkta  ***')
      end if
c
      if(nrayl.gt.0) then
        write(6,785) nrayl
        write(11,785) nrayl
785     format(/
     +'***  less than nray rays traced for ',i4,' ray groups  ***')
      end if
c
      if(iflagi.eq.1) write(6,775)
775   format(/'***  attempt to interpolate over ximax  ***')
c
      if(i33.eq.1) then
        do 1033 i=1,narinv
           write(34,335) xscalc(i),abs(icalc(i)),ircalc(i),xcalc(i),
     +                   tcalc(i)
1033    continue
      end if
c
      if(abs(modout).ne.0) call
     + modwr(modout,dxmod,dzmod,modi,ifrbnd,frz,xmmin,xmmax)
c
      if(ifd.gt.0) call fd(dxzmod,xmmin,xmmax,ifd)
c
900   ntblk=0
      do 920 i=1,nlayer
         ntblk=ntblk+nblk(i)
920   continue
c
      if(isum.gt.0) write(6,935)
      write(11,935)
935   format(/'|---------------------------------------------------',
     +        '-------------|'/'|',64x,'|')
      if(ntpts.gt.0) then
        if(isum.gt.0) write(6,905) ntray,ntpts
        write(11,905) ntray,ntpts
905     format('| total of ',i6,' rays consisting of ',i8,
     +         ' points were traced |'/'|',64x,'|')
      else
        if(isum.gt.0) write(6,915)
        write(11,915)
915     format('|                   ***  no rays traced  ***',
     +         21x,'|'/'|',64x,'|')
      end if
      if(isum.gt.0) write(6,925) nlayer,ntblk
      write(11,925) nlayer,ntblk
925   format('|           model consists of ',i2,' layers and ',i3,
     +       ' blocks           |'/,'|',64x,'|'/
     +'|----------------------------------------------------------',
     +'------|'/)
c
      if(invr.eq.1) then
c
        open(unit=18, file='i.out')
c
        parunc(1)=bndunc
        parunc(2)=velunc
        parunc(3)=bndunc
        write(18,895)
895     format(' ')
        write(18,835) narinv,nvar
835     format(2i5)
        write(18,895)
        do 801 i=1,nvar
           write(18,805) partyp(i),parorg(i),parunc(partyp(i))
805        format(i5,2f15.5)
801     continue
        write(18,895)
        do 840 i=1,narinv
           write(18,815) (apart(i,j),j=1,nvar)
815        format(5e12.5)
840     continue
        write(18,895)
        write(18,815) (tobs(i)-tcalc(i),i=1,narinv)
        write(18,895)
        write(18,815) (uobs(i),i=1,narinv)
c
        if(narinv.gt.1) then
          sum=0.
          sumx=0.
          do 850 i=1,narinv
             sum=sum+(tobs(i)-tcalc(i))**2
             sumx=sumx+((tobs(i)-tcalc(i))/uobs(i))**2
850       continue
          trms=sqrt(sum/float(narinv))
          chi=sumx/float(narinv-1)
          write(18,895)
          write(18,825) narinv,trms,chi
825       format('Number of data points used: ',i8/
     +           'RMS traveltime residual:    ',f8.3/
     +           'Normalized chi-squared:   ',f10.3)
          write(18,895)
          write(6,895)
          write(6,825) narinv,trms,chi
          write(6,895)
          write(11,895)
          write(11,825) narinv,trms,chi
          write(11,895)
          if(isum.gt.1) then
c
            write(11,968)
968         format(' phase    npts   Trms   chi-squared'/
     +             '-----------------------------------')
            if(isum.eq.3) write(6,968)
            nused=0
            j=0
951         j=j+1
            sum=0.
            sumx=0.
            nars=0
            do 952 i=1,narinv
               if(abs(icalc(i)).eq.j) then
                 sum=sum+(tobs(i)-tcalc(i))**2
                 sumx=sumx+((tobs(i)-tcalc(i))/uobs(i))**2
                 nars=nars+1
                 nused=nused+1
               end if
952         continue
            if(nars.gt.0) then
              trms=sqrt(sum/float(nars))
            end if
            if(nars.gt.1) then
              chi=sumx/float(nars-1)
            else
              chi=sumx
            end if
            if(nars.gt.0) then
              write(11,969) j,nars,trms,chi
969           format(i6,i8,f8.3,f10.3)
              if(isum.eq.3) write(6,969) j,nars,trms,chi
            end if
            if(nused.ge.narinv) go to 953
            go to 951
953         write(11,895)
            if(isum.eq.3) write(6,895)
c
            write(11,868)
868         format(
     +      '     shot  dir   npts   Trms   chi-squared'/
     +      '------------------------------------------')
            if(isum.eq.3) write(6,868)
            xsn=xscalc(1)
            icn=sign(1,icalc(1))
            do 851 j=1,nshot
               iflagn=1
               xsc=xsn
               icc=icn
               sum=0.
               sumx=0.
               nars=0
               do 852 i=1,narinv
                  if(xscalc(i).eq.xsc.and.sign(1,icalc(i)).eq.icc.
     +            and.icalc(i).ne.0) then
                    tdiff=abs(tobs(i)-tcalc(i))
                    sum=sum+tdiff**2
                    sumx=sumx+(tdiff/uobs(i))**2
                    nars=nars+1
                    icalc(i)=0
                  else
                    if(iflagn.eq.1.and.icalc(i).ne.0) then
                      xsn=xscalc(i)
                      icn=sign(1,icalc(i))
                      iflagn=0
                    end if
                  end if
852            continue
               if(nars.gt.0) then
                 trms=sqrt(sum/float(nars))
               end if
               if(nars.gt.1) then
                 chi=sumx/float(nars-1)
               else
                 chi=sumx
               end if
               if(nars.gt.0) then
                 write(11,869) xsc,icc,nars,trms,chi
869              format(f10.3,i3,i8,f8.3,f10.3)
c                write(44,867) xsc,chi
867              format(f10.3,f10.3)
                 if(isum.eq.3) write(6,869) xsc,icc,nars,trms,chi
               end if
851         continue
            write(11,895)
            if(isum.eq.3) write(6,895)
          end if
        end if
      end if
c
      if(iplots.eq.1) call plotnd(1)
c
9999  if(idump.eq.1) then
c
        open(unit=22, file='n.out')
c
        write(22,pltpar)
        write(22,axepar)
        write(22,trapar)
        write(22,invpar)
      end if
c
      stop
c
999   write(6,105)
105   format(/'***  error in velocity model 1  ***'/)
c
      go to 9999
c
      end
