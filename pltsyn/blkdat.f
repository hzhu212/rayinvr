c                 
c     version 1.2  Mar 1992
c
c     Block data for PLTSYN
c
c     ----------------------------------------------------------------
c                 
      block data  
c                 
c     assign default and initial values to common block parameters
c                 
      include 'pltsyn.par'
      include 'pltsyn.com'
c
      data xmin,xmax,xmm,ntickx,ndecix/0.,300.,250.,-1,-2/,
     +     tmin,tmax,tmm,ntickt,ndecit/0.,10.,125.,-1,-2/,albht/2.5/,
     +     xtmin,xtmax,ttmin,ttmax/4*-999999./,
     +     rcor,xnorm,scalef,amp,iscale,clip/1.,100.,100.,-1.,0,0./,
     +     nsmth,iconv,sps,iwavlt/0,1,60.,0/,
     +     ishade,ifill,dens/0,1,5./,inmo,vrms/0,2.5/,
     +     iplot,iplots,orig,iseg,nseg/1,0,12.5,0,0/,
     +     xwndow,ywndow/2*0./,ibcol,ifcol/0,1/,sf/1.2/
c                 
      end
