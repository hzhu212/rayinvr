c                 
c     version 1.2  Mar 1992
c                 
c     common blocks for PLTSYN
c
c     ----------------------------------------------------------------
c                 
      common /blk1/ xmin,xmax,xtmin,xtmax,xmm,ntickx,ndecix,xscale,
     +              tmin,tmax,ttmin,ttmax,tmm,ntickt,ndecit,tscale,
     +              albht,iplots,orig
      common /blk2/ rcor,xnorm,scalef,amp,iscale,clip,nsmth,iconv,sps,
     +              iwavlt
      common /blk3/ inmo,vrms
      common /blk4/ ishade,ifill,dens
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
