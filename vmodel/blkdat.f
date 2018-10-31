c
c     version 1.2  Mar 1992
c
c     Block data for VMODEL
c                 
c     ----------------------------------------------------------------
c                 
      block data  
c                 
c     assign default and initial values to common block parameters
c                 
      include 'vmodel.par'
      include 'vmodel.com'
c                 
      data 
     +  xmin,xmax,xmm,ndecix,ntickx/2*-99999.,250.,-2,-1/,
     +  zmin,zmax,zmm,ndeciz,ntickz/0.,50.,75.,-2,-1/,
     +  vmin,vmax,vmm,ndeciv,ntickv/0.,8.,75.,-2,-1/,
     +  tmin,tmax,tmm,ndecit,ntickt/0.,10.,75.,-2,-1/,
     +  vrmin,vrmax,vrmm,ndecir,ntckvr/2.,7.,75.,-2,-1/,
     +  xtmin,xtmax,ztmin,ztmax,vtmin,vtmax,ttmin,ttmax,vrtmin,vrtmax
     +  /10*-999999./,
     +  albht,orig/2.5,12.5/,xwndow,ywndow/2*0./,
     +  iplot,iplots,iseg,nseg/1,0,0,0/,sf,ibcol,ifcol/1.2,0,1/
c                 
      end         
