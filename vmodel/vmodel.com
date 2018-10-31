c      
c     version 1.2  Mar 1992
c
c     common blocks for VMODEL
c           
c     ----------------------------------------------------------------
c
      real xm(pncntr,ppcntr),zm(pncntr,ppcntr),vf(player,ppvel,2),
     +     xvel(player,ppvel,2)
      integer nzed(pncntr),nvel(player,2)
c
      common /blk1/ xm,zm,vf,nzed,nvel,xvel
      common /blk2/ 
     +       xmin,xmax,xtmin,xtmax,xmm,ndecix,xscale,ntickx,
     +       zmin,zmax,ztmin,ztmax,zmm,ndeciz,zscale,ntickz,
     +       vmin,vmax,vtmin,vtmax,vmm,ndeciv,vscale,ntickv,
     +       tmin,tmax,ttmin,ttmax,tmm,ndecit,tscale,ntickt,
     +       vrmin,vrmax,vrtmin,vrtmax,vrmm,ndecir,rscale,ntckvr,
     +       albht,iplots,orig
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
