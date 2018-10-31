c                 
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |          ***********   P L T S Y N   *************           |
c     |                                                              |
c     |          Synthetic Record Section Plotting Program           |
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
c        10 -- input: program control parameters
c
c        11 -- input: calculated travel time, amplitude, and phase 
c                     of each arrival of each trace of synthetic 
c                     section
c
c        12 -- output: synthetic seismograms
c
c        17 -- input: tx.out file output by RAYINVR
c
c        18 -- output: seismogram amplitudes at times specified in 
c                      tx.out
c
c        20 -- input: source wavelet
c                 
c     ----------------------------------------------------------------
c                 
      program main
c
      include 'pltsyn.par'
c
      real xomit(pseis)
c
      include 'pltsyn.com'
c
      namelist /pltpar/ iplot,iconv,iwavlt,ipol,iscale,ishade,ifill,
     +                  inmo,idump,vred,sps,nsmth,scalef,rcor,xnorm,amp,
     +                  clip,nptsw,nskip,spmin,xomit,dens,vrms,
     +                  iroute,iamp,twin,imeth,iseg,xwndow,ywndow,
     +                  ibcol,ifcol
c
      namelist /axepar/ itrev,albht,orig,
     +                  xmin,xmax,xtmin,xtmax,xmm,ntickx,ndecix,
     +                  tmin,tmax,ttmin,ttmax,tmm,ntickt,ndecit
c
      data xomit/pseis*-99999./
c                 
c     default parameter values
c                 
      iroute=1
      twin=.25
      imeth=1
      iamp=0
      itrev=0     
      idump=0     
      nptsw=19    
      spmin=-1.   
      ipol=0      
      nskip=1     
      vred=8.    
c                 
      open(unit=10, file='s.in', status='old')
c
      read(10,pltpar)
      read(10,axepar)
c                 
      open(unit=11, file='sect.out', status='old')
      if(idump.eq.1) open(unit=12, file='syn.out')
      if(iwavlt.eq.2) open(unit=20, file='w.in', status='old')
      if(iplot.eq.-1) iplot=-2
      if(iplot.eq.0) iplot=-1
      if(iplot.eq.2) iplot=0
      if(iplot.eq.-1.or.iplot.eq.0) open(unit=19, file='p.out')
c
      xscale=(xmax-xmin)/xmm
      tscale=(tmax-tmin)/tmm
      if(itrev.eq.1) tscale=-tscale
c                 
      ipol=-2*ipol+1
      if(amp.lt.0.) amp=(xmax-xmin)/100.
      if(spmin.lt.0.) spmin=(xmax-xmin)/100000.
c
      if(iamp.gt.0) then
        open(17, file='tx.out')
        open(18, file='amp.out')
      end if
c
      if(iroute.ne.1) ibcol=0
c                 
      call pltsec(vred,xomit,nskip,ipol,spmin,nptsw,itrev,idump,
     +            iamp,twin,imeth,iroute)
c                 
      if(iplots.eq.1) call plotnd(1)
c
      stop        
      end         
