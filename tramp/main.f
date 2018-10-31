c
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ***********  T R A M P  ***********               |
c     |                                                              |
c     |            Two-dimensional Ray Tracing Program               |
c     |           for Traveltime and Amplitude Modelling             |
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
c        13 -- ouptut: summary amplitude information for each ray
c                      traced
c                 
c        14 -- output: details of amplitude calculations for each
c                      ray traced
c                 
c        16 -- input:  receiver locations
c
c        17 -- output: calculated traveltime-distance pairs
c
c        19 -- output: all Calcomp plot calls
c                 
c        20 -- input:  velocity model
c                 
c        21 -- ouptut: calculated traveltime and amplitude of each
c                      arrival of each trace of synthetic sections
c                 
c        22 -- ouptut: namelist parameters
c
c        23 -- output: misc output
c
c        24 -- output: calculated amplitude-distance pairs
c                 
c        26 -- input:  observed traveltime-distance pairs
c                 
c                 
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'tramp.par'
c                
      real amin(prayt),amax(prayt),xshot(pshot),xvz(pnvz),
     +     xshots(pshot2),theta(pnrayf),dist(pnrayf),space(prayf),
     +     zshot(pshot),pois(player),poisbl(papois),qp(player),
     +     qs(player),qpbl(paqpqs),qsbl(paqpqs),zsmth(pnsmth),
     +     xshota(pshot2),zshota(pshot2),tatan(pitan2),angtan(pitan2)
c                 
      integer idr(pshot2),nray(prayf),itt(ptrayf),ishot(pshot),
     +        ncbnd(prayf),nsec(pshot),mir(pnrayf),ishotr(pshot2),
     +        nrbnd(prayf),rbnd(preflt),ncaust(prayf),iturn(prayf),
     +        nsmax(prayf),dbnd(prayf),amprev,iint(prayf),hlayer,
     +        vrmstl,vrmsbl,cbnd(pconvt),ibrka(ptrayf),ishotw(pshot2),
     +        ql(2*paqpqs),qb(2*paqpqs),poisl(papois),poisb(papois),
     +        ibreak(prayf),ihead(player),irayt(prayt),iinta(ptrayf),
     +        mcol(5),modi(player),frbnd(prayf)
c
      include 'tramp.com'                 
c                 
      namelist /pltpar/ iplot,imod,ibnd,idash,ivel,iray,irays,
     +                  irayps,idot,icntr,itx,idata,iszero,itxout,
     +                  iamout,ivz,ivzp,ivrms,igrid,idump,isum,istep,
     +                  symht,velht,nskip,nrskip,vred,iroute,isep,
     +                  ttunc,ampunc,vcntr,xcinc,xclab,ncsmth,
     +                  xvz,vrmstl,vrmsbl,xrinc,nrsmth,xgrid,zgrid,
     +                  itxbox,dvmax,dsmax,iseg,xwndow,ywndow,
     +                  mcol,ircol,itcol,colour,ibcol,ifcol,
     +                  modout,dxmod,dzmod,modi,frz,xmmin,xmmax,
     +                  xmin1d,xmax1d,npskip
c
      namelist /axepar/ iaxlab,itrev,albht,orig,sep,
     +                  xmin,xmax,xtmin,xtmax,xmm,ntickx,ndecix,
     +                  xmint,xmaxt,xtmint,xtmaxt,xmmt,ntckxt,ndecxt,
     +                  zmin,zmax,ztmin,ztmax,zmm,ntickz,ndeciz,
     +                  tmin,tmax,ttmin,ttmax,tmm,ntickt,ndecit,
     +                  ampmin,ampmax,atmin,atmax,amm,nticka,ndecia,
     +                  vmin,vmax,vmm,ntickv,ndeciv,
     +                  vrmin,vrmax,vrmm,ntckvr,ndecir
c                 
      namelist /trapar/ ishot,iraysl,irayt,iturn,isrch,istop,ibsmth,
     +                  imodf,xshot,zshot,ray,nray,space,amin,amax,
     +                  nsmax,aamin,aamax,stol,xsmax,step,smin,smax,
     +                  nrbnd,rbnd,ncbnd,cbnd,pois,poisb,poisl,poisbl,
     +                  crit,hws,angbnd,dbnd,npbnd,nbsmth,nhray,frbnd,
     +                  i2pt,n2pt,x2pt,ifast,ntan
c
      namelist /amppar/ iamp,isect,ibreak,iint,icaust,icomp,icmp,iden,
     +                  surcon,amprev,ampsmt,spamp,denc,denmin,ncaust,
     +                  omega,hedcut,qp,qs,ql,qb,qpbl,qsbl,
     +                  xmins,xmaxs,xincs,izrefl
c
c     initialize parameters
c
      data amin/prayt*181./,amax/prayt*181./,
     +     xvz/pnvz*-99999./,
     +     pois/player*-99./,poisbl/papois*-99./,qp/player*-1./,
     +     qs/player*-1./,qpbl/paqpqs*-1./,qsbl/paqpqs*-1./,
     +     space/prayf*0./,zshot/pshot*-9999./
      data nray/prayf*-1/,itt/ptrayf*0/,
     +     nsec/pshot*0/,nrbnd/prayf*0/,
     +     rbnd/preflt*0/,ncaust/prayf*0/,iturn/prayf*1/,
     +     iint/prayf*0/,nsmax/prayf*-1/,dbnd/prayf*0/,
     +     cbnd/pconvt*-99/,ibreak/prayf*1/,
     +     ishot/pshot*0/,ncbnd/prayf*-1/,irayt/prayt*1/
c     
      istep=0
      ntan=90
      ifast=1
      iroute=1
      itol=0
      ampunc=10.
      itrev=0
      xsmax=0.
      ivzp=1
      nrayl=0
      irays=0
      irayps=0
      iraysl=0
      imodf=0
      ncont=1
      ifam=0
      ictbnd=0
      isec=1
      irbnd=0
      ngroup=0
      nshot=0
      ia0=1
      iden=0   
      icmp=0      
      nrskip=1    
      itxout=0 
      iamout=0
      spamp=1.    
      amprev=0    
      idata=0     
      isum=1      
      stol=-1.    
      ampsmt=0.1  
      ttunc=.01  
      nrsmth=0    
      ivrms=0     
      xrinc=-1.   
      vrmstl=-1   
      vrmsbl=-1   
      icaust=1    
      idash=1     
      ivel=0      
      velht=0.7   
      iaxlab=1    
      splnf=1.    
      sdev=1.     
      aamin=5.    
      aamax=85.   
      aainc=.1    
      aaimin=1.    
      idot=0      
      iszero=0    
      ibnd=1      
      imod=1      
      iray=1      
      iamp=0      
      nskip=0     
      ivz=0       
      isect=0     
c
      open(unit=10, file='r.in', status='old')
c                 
c     read in program control parameters from unit 10
c                 
      read(10,pltpar)
      read(10,axepar)
      read(10,trapar)
      read(10,amppar)
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
         read(iunit,235,end=99) idum
15       format(i2,1x,10f7.2)
235      format(3x,10i7)
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
         read(iunit,235,end=999) idum
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
         read(iunit,235,end=999) idum
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
      if(iplot.eq.0) iplot=-1
      if(iplot.eq.2) iplot=0
c     
c     open I/O units
c
      
      open(unit=11, file='r1.out')
      if(idump.eq.1) open(unit=12, file='r2.out')
      if(iamp.gt.0) open(unit=13, file='a1.out')
      if(iamp.gt.0.and.idump.eq.1) open(unit=14, file='a2.out')
      if(isect.eq.2) open(unit=16, file='rec.in', status='old')
      if(itxout.gt.0) open(unit=17, file='tx.out') 
      if(iplot.le.0) open(unit=19, file='p.out')
      if(isect.gt.0) open(unit=21, file='sect.out') 
      if(((ivrms.ne.0.or.icntr.ne.0.or.ivz.ne.0).and.idump.eq.1).or.
     +  igrid.ne.0.or.iden.eq.1) open(unit=23, file='m.out') 
      if(iamout.eq.1) open(unit=24, file='amp.out') 
c
c     read in observed data
c
      if(idata.ne.0) then
        open(unit=26, file='tx.in', status='old')
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
c     time and amplitude plots are same as model distance parameters
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
      vscale=(vmax-vmin)/vmm
      rscale=(vrmax-vrmin)/vrmm
      ascale=(ampmax-ampmin)/amm
c
      if(itrev.eq.1) tscale=-tscale
      if(iroute.ne.1) then
        isep=0 
        ibcol=0
      end if
      if(isep.eq.2) isep=3
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
c                 
c     calculate velocity model parameters
c                 
      call calmod(ncont,pois,poisb,poisl,poisbl,
     +            qp,qs,qb,ql,qpbl,qsbl,iamp,iden,iflagm)
c
      if(iflagm.eq.1) go to 9999
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
c     assign default values to xrinc, vrmstl, and vrmsbl if not 
c     specified   
c                 
      if(ivrms.ne.0) then
        if(xrinc.lt.0.) xrinc=(xmaxt-xmint)/100.
        if(vrmstl.lt.0) vrmstl=1
        if(vrmsbl.lt.0) vrmsbl=nlayer
      end if      
c                 
c     assign default value to xcinc if not specified
c                 
      if(icntr.ne.0) then
        if(xcinc.lt.0.) xcinc=(xmax-xmin)/100.
      end if      
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
     +  call pltmod(ncont,ibnd,imod,iaxlab,ivel,velht,idash,idata,
     +              iroute)
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
      if(iamp.gt.0) then
        write(13,85)
85      format(' shot  ray#     dist   amplitude     phase  code')
        if(hedcut.le.0.) hedcut=1.
        hedcut=hedcut/100.
        omega=omega*6.283185307
        sdev=sdev/100.
      end if      
c
      if(nrskip.lt.1) nrskip=1
      if(hws.lt.0.) hws=(xmax-xmin)/25.
      crit=crit/pi18 
      angbnd=angbnd/pi18
      do 30 i=1,player
         ihead(i)=0
30    continue    
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
      do 40 i=1,prayf
         if(ray(i).lt.1.) go to 50
         ngroup=ngroup+1
40    continue  
50    if(ngroup.eq.0) then
        write(6,55)
55      format(/'***  no ray codes specified  ***'/)
        go to 900
      end if    
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
      do 3900 k=1,ngroup
         if(izrefl(k).gt.0) then
           open(unit=36, file='zf.in', status='old')
           i=1
395        read(36,365,end=399) npzf(i)
365        format(i2)
           read(36,385) (xzf(i,j),j=1,npzf(i))
           read(36,385) (zff(i,j),j=1,npzf(i))
385        format(10f7.2)
           i=i+1 
           go to 395
399        nzf=i-1
           go to 3910
         end if
3900  continue
3910  continue
c                 
c     assign values to xmins, xmaxs and xincs if not specified
c                 
      if(isect.eq.1) then
        xincss=(xmax-xmin)/50.
        if(xmins.lt.-9999.) xmins=xmin+xincss
        if(xmaxs.lt.-9999.) xmaxs=xmax-xincss
        if(xincs.lt.-9999.) xincs=(xmaxs-xmins)/48.
        iflags=0  
        if(xmins.ge.xmaxs.or.xincs.le.0.) iflags=1
        nseis=0   
        if(xincs.ne.0.) nseis=int((xmaxs-xmins)/xincs+.5)+1
        if(nseis.le.0.or.nseis.gt.pseis) iflags=1
        if(iflags.eq.1) then
         write(6,475) 
475      format(/ 
     +   '***  error in specification of seismogram parameters  ***'/)
         go to 9999
        end if    
      end if      
c                 
c     assign default values to smin and smax if not specified
c                 
      if(smin.lt.0.) smin=(xmax-xmin)/4500.
      if(smax.lt.0.) smax=(xmax-xmin)/15.
c                 
      ist=0       
      iflagp=0
c
      if(nshot.eq.0) go to 1000
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
           if(iszero.eq.0) xshots(isec)=xshotr
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
             isec=isec+1
             xsec=xshotr
             zsec=zshotr
             idsec=id
             if(iszero.eq.0) xshots(isec)=xshotr
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
             call pltmod(ncont,ibnd,imod,iaxlab,ivel,velht,idash,idata,
     +                   iroute)
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
            nsec(isec)=nsec(isec)+1
            nrayr=nray(i)
            caust=float(ncaust(i))
            iintf=iint(i)
            iinta(ifam)=iint(i)
            ibrka(ifam)=ibreak(i)
            difh=dbnd(i)
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
            end if
            if(ncbnd(i).gt.0) then
              do 630 j=1,ncbnd(i)
                 ictbnd=ictbnd+1
                 icbnd(j)=cbnd(ictbnd)
630           continue
            end if
c               
            call auto(xshotr,zshotr,i,ifam,idl,idt,aminr,amaxr,
     +           aamin,aamax,layer1,iblk1,aainc,aaimin,nsmax(i),
     +           iflag,iturn(i),amin(ia0),amax(ia0),ia0,stol,
     +           irays,nskip,idot,irayps,xsmax,istep)
c               
            if(iflag.ne.0) then
              write(11,65) ishotw(is),ray(i)
65            format('***  shot#',i4,' ray code',f5.1,
     +               ' no rays traced  ***')
              nrayr=0
              nrayl=nrayl+1
              go to 69 
            end if 
            if(nrayr.le.0) go to 69
            if(amaxr.eq.aminr.and.ihdwf.ne.1) then
              if(nrayr.gt.1) then
                write(11,665) ishotw(is),ray(i)
665             format('***  shot#',i4,' ray code',f5.1,
     +                 ' 1 ray traced  ***')
                nrayl=nrayl+1
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
            if(iamp.gt.0.and.idt.eq.3) then
              ihamp=1
            else
              ihamp=0
            end if
            namp=0
            iend=0
            nrg=0
c
            do 90 ir=1,nrayr
               if(iend.eq.1) go to 890
               difbnd=difh
               ircbnd=1
               iccbnd=1
               iwave=1
               if(icbnd(1).eq.0) then
                 iwave=-iwave
                 iccbnd=2
               end if
               iep=0
               nbnd=0
               nccbnd=0
               id=idr(is)
               fid=float(id)
               fid1=fid
               iturnt=0
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
               angle=fid*(90.-angled)/pi18
               if(ir.gt.1) then
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
               qr(1)=q(layer1,iblk1,(3-iwave)/2)
               idray(1)=layer1
               idray(2)=1
               nrg=nrg+1
c                 
               call trace(npt,ifam,ir,iturnt,iamp,ihamp,iflag,idl,idt)
c                 
               if(iflag.gt.1) then
                if(ntt.gt.pray) then
                  write(6,995)
                  go to 1000
                end if
                if(iflag.eq.2) then
c                 
                 call hdwave(is,npt,ifam,ir,xshotr,angled,nskip,itt,
     +                hlayer,iszero,idot,iray,iamp,caust,spamp,nrskip,
     +                irayps,xsmax,ishotw(is),ampsmt,istep,angled)
c                 
                 iheadf(hlayer)=0
                else
                 call difrct(is,npt,ifam,ir,xshotr,angled,nskip,itt,
     +              iszero,idot,iray,nrskip,irayps,xsmax,ishotw(is),
     +              istep,anglew)
                end if
               else
c                 
                 call ttime(ishotw(is),xshotr,npt,ir,angled,ifam,itt,
     +                    iszero,iamp,iflag)
c                 
                 if((iray.eq.1.or.(iray.eq.2.and.vr(npt,2).gt.0.)).
     +             and.mod(ir-1,nrskip).eq.0) 
     +             call pltray(npt,nskip,idot,irayps,istep,angled)
c                 
                 if(iamp.gt.0.and.vr(npt,2).ne.0.) 
     +             call calamp(npt,ifam,ir,namp,theta,dist,mir,caust,
     +                xshotr,layer1,iblk1,amprev,spamp,ishotw(is),i)
c                 
               end if 
               if(ntt.gt.pray) then
                 write(6,995)
995              format(/
     +     '***  max number of rays reaching surface exceeded  ***'/)
                 go to 890
               end if 
90          continue
c                 
890         if(iamp.gt.0.and.namp.gt.0) call spread(ifam,theta,dist,
     +        mir,namp,sdev,splnf,ampsmt,iintf,icaust,ibrka(ifam),
     +        xshotr,fid1,idump,ishotw(is))
c                 
            if(iray.gt.0) call empty
c
            nrayr=nrg
69          write(6,375) ishotw(is),ray(i),nrayr
375         format('shot#',i4,':   ray code',f5.1,':   ',
     +             i3,' rays traced')
c
            if(ntt.gt.pray) go to 1000
c
70       continue 
60    continue    
c                 
1000  if((itx.gt.0.and.ntt.gt.1).or.idata.ne.0) then
        if(isep.gt.0.and.iplots.eq.1) call aldone
        call plttx(ifam,itt,iszero,idata,iaxlab,xshota,
     +             idr,nshot,itxout,ibrka,ttunc,itrev,iroute)
      end if
c
      if(irkc.eq.1) write(6,75)
75    format(/'***  possible inaccuracies in rngkta  ***'/)
c                 
      amp1=ampmin
      amp2=ampmax
      if(iamp.eq.2.and.ntt.gt.1) call pltamp(ifam,itt,iaxlab,
     +                       ibrka,iamout,ampunc,iroute,amp1,amp2)
c                 
      if(isect.gt.0.and.ntt.gt.2.and.iamp.gt.0)
     +  call calsec(xshots,itt,nsec,isec,isect,iinta,ibrka,icmp)
c                 
      if(ivrms.ne.0) call pltrms(xrinc,vrmstl,vrmsbl,nrsmth,iaxlab,
     +                           ivrms,iroute)
c                 
      if(ivz.ne.0) call pltvz(xvz,iaxlab,ivz,ivzp,iroute)
c                 
      if(nrayl.gt.0) then
        write(6,785) nrayl
        write(11,785) nrayl
785     format(/
     +'***  less than nray rays traced for ',i4,' ray groups  ***')
      end if
c
900   if(isum.gt.0) then
        ntblk=0     
        do 920 i=1,nlayer
           ntblk=ntblk+nblk(i)
920     continue    
c                 
        write(6,935)
        write(11,935) 
935     format(/'|---------------------------------------------------',
     +          '-------------|'/'|',64x,'|')
        if(ntpts.gt.0) then 
          write(6,905) ntray,ntpts
          write(11,905) ntray,ntpts
905       format('|  total of ',i5,' rays consisting of ',i7,
     +           ' points were traced  |'/'|',64x,'|')
        else        
          write(6,915)
          write(11,915)
915       format('|                   ***  no rays traced  ***',
     +           21x,'|'/'|',64x,'|')
        end if      
        write(6,925) nlayer,ntblk
        write(11,925) nlayer,ntblk
925     format('|           model consists of ',i2,' layers and ',i3,
     +         ' blocks           |'/,'|',64x,'|'/
     +  '|----------------------------------------------------------',
     +  '------|'/) 
c
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
        write(22,amppar)
      end if
c
      stop        
c                 
999   write(6,105)
105   format(/'***  error in velocity model  ***'/)
c                 
      go to 9999  
c
      end         
