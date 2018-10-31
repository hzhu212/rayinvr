c
c     version 1.3  Aug 1992
c
c     Block data for RAYINVR
c                 
c     ----------------------------------------------------------------
c                 
      block data  
c                 
c     assign default and initial values to common block parameters
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c                 
      data 
     +  irkc,tol,hdenom,hmin,idump/0,.0005,64.,.01,0/,
     +  itx,vred/1,8./,nzed/pncntr*1/,nvel/pinvel*1/,
     +  step,smin,smax/.05,-1.,-1./,ntt,ray/1,prayf*0./,
     +  xmin,xmax,xmm,ndecix,ntickx/0.,-99999.,250.,-2,-1/,
     +  xmint,xmaxt,xmmt,ndecxt,ntckxt/3*-9999.,-2,-1/,
     +  zmin,zmax,zmm,ndeciz,ntickz/0.,50.,75.,-2,-1/,
     +  tmin,tmax,tmm,ndecit,ntickt/0.,10.,75.,-2,-1/,
     +  xtmin,xtmax,xtmint,xtmaxt,ztmin,ztmax,ttmin,ttmax/8*-999999./,
     +  symht,albht/.5,2.5/,ibsmth,nbsmth,npbnd/0,10,100/,
     +  cv/picv*0./,crit,hws/1.,-1./,
     +  iplot,iplots,orig,sep,iseg,nseg/1,0,12.5,7.5,0,0/,
     +  xwndow,ywndow/2*0./,colour/pcol*-1/,mcol/5*-1/,
     +  sf,ibcol,ifcol/1.2,0,1/
c                 
      end         
