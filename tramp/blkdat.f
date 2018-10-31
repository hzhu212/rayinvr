c
c     version 1.3  Aug 1992
c
c     Block data for TRAMP
c                 
c     ----------------------------------------------------------------
c                 
      block data  
c                 
c     assign default and initial values to common block parameters
c                 
      include 'tramp.par'
      include 'tramp.com'
c                 
      data hws,crit,angbnd,omega,hedcut/-1.,1.,20.,5.,25./,
     +  irkc,tol,hdenom,hmin,idump/0,.0005,64.,.01,0/,
     +  itx,vred/0,8./,nzed/pncntr*1/,nvel/pinvel*1/,
     +  step,smin,smax/.05,-1.,-1./,ntt,ray/1,prayf*0./,
     +  xmin,xmax,xmm,ndecix,ntickx/0.,-99999.,250.,-2,-1/,
     +  xmint,xmaxt,xmmt,ndecxt,ntckxt/3*-9999.,-2,-1/,
     +  zmin,zmax,zmm,ndeciz,ntickz/0.,50.,75.,-2,-1/,
     +  tmin,tmax,tmm,ndecit,ntickt/0.,10.,75.,-2,-1/,
     +  vmin,vmax,vmm,ndeciv,ntickv/0.,8.,75.,-2,-1/,
     +  vrmin,vrmax,vrmm,ndecir,ntckvr/2.,7.,75.,-2,-1/,
     +  ampmin,ampmax,amm,ndecia,nticka/-6.,-1.,75.,-2,-1/,
     +  symht,albht/.5,2.5/,xwndow,ywndow/2*0./,colour/pcol*-1/,
     +  xtmin,xtmax,xtmint,xtmaxt,ztmin,ztmax,ttmin,ttmax,
     +  vtmin,vtmax,vrtmin,vrtmax,atmin,atmax/14*-999999./
      data iplot,iplots,orig,sep,iseg,nseg/1,0,12.5,7.5,0,0/,
     +  ibsmth,nbsmth,npbnd/0,10,100/,
     +  surcon/1/,xmins,xmaxs,xincs/3*-99999./,
     +  icntr,xcinc,vcntr,ncsmth,xclab
     +  /0,-1.,pcntr*-999.,0,plcntr*-9999./,
     +  denc,denmin/-.6997,2.2302,-.598,.07036,-.0028311,1.25/,
     +  igrid,xgrid,zgrid/0,-1,-1/,nstepr/15/,ntray,ntpts/0,0/,
     +  icomp/1/,isrch/0/,iq/0/,istop/1/,ibcol,ifcol,sf/0,1,1.2/,
     +  izrefl/prayf*0/
c                 
      end         
