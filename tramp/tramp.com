c      
c     version 1.3  Aug 1992
c
c     common blocks for TRAMP
c           
c     ----------------------------------------------------------------
c
      integer refll(prefl+1),difbnd,surcon,nblk(player),
     +        iheadf(player),idray(2),irtamp(piray,3),
     +        nptamp(piray),icbnd(pconv+1),ivg(player,ptrap),
     +        nzed(pncntr),nvel(player,2),ipf(prayi+pshot2+1),
     +        ilshot(pshot2+1),colour(pcol),npzf(pzff),izrefl(prayf)
      real*4 c(player,ptrap,11),s(player,ptrap,2),b(player,ptrap,2),
     +       vm(player,ptrap,4),xbnd(player,ptrap,2),
     +       xm(pncntr,ppcntr),zm(pncntr,ppcntr),xvel(player,ppvel,2),
     +       vf(player,ppvel,2),xr(ppray),zr(ppray),range(pray),
     +       ar(ppray,2),vr(ppray,2),tr(ppray),qr(ppray),xshtar(pray),
     +       vp(ppray,2),vs(ppray,2),sh(ppray/2),rayid(pray),tt(pray),
     +       aamp(piray,2),vpamp(piray,2),vsamp(piray,2),fidarr(pray),
     +       ray(prayf),amp(pray),phase(pray),p0(pray),pr(pray),
     +       xclab(plcntr),q(player,ptrap,2),vsvp(player,ptrap),
     +       denc(5),tang(player,4),cosmth(pncntr,pnsmth),vcntr(pcntr),
     +       xpf(prayi+pshot2+1),tpf(prayi+pshot2+1),
     +       upf(prayi+pshot2+1),xzf(pzff,ppzff),zff(pzff,ppzff),
     +       mtan(pitan2),btan(pitan2),mcotan(pitan2),bcotan(pitan2)
c
      common /blk1/ layer,iblk,id,fid,fid1
      common /blk2/ c,ivg
      common /blk3/ s,b,vm
      common /blk4/ xbnd,nblk,nlayer
      common /blk5/ xm,zm,vf,nzed,nvel,xvel
      common /blk6/ xr,zr,ar,vr,tr,qr,vp,vs
      common /blk7/ iheadf,hws,crit,refll,angbnd,sh,omega,hedcut,ircbnd
      common /blk8/ irkc,tol,hdenom,hmin,idump,isrkc,ifast
      common /blk9/ itx,vred,time
      common /blk10/ step,smin,smax
      common /blk11/ range,tt,rayid,xshtar,fidarr,idray,ntt,ray
      common /blk12/ xmin,xmax,xtmin,xtmax,xmm,ndecix,xscale,ntickx,
     +            xmint,xmaxt,xtmint,xtmaxt,xmmt,ndecxt,xscalt,ntckxt,
     +            zmin,zmax,ztmin,ztmax,zmm,ndeciz,zscale,ntickz,
     +            tmin,tmax,ttmin,ttmax,tmm,ndecit,tscale,ntickt,
     +            vmin,vmax,vtmin,vtmax,vmm,ndeciv,vscale,ntickv,
     +            vrmin,vrmax,vrtmin,vrtmax,vrmm,ndecir,rscale,ntckvr,
     +            ampmin,ampmax,atmin,atmax,amm,ndecia,ascale,nticka,
     +            symht,albht,iplots,orig,sep,itcol,colour
      common /blk13/ aamp,irtamp,nptamp,nbnd,vpamp,vsamp
      common /blk14/ amp,phase,p0,pr,surcon
      common /blk15/ xmins,xmaxs,xincs
      common /blk16/ icasel,dstepf,n2,n3,istop
      common /blk17/ difbnd
      common /blk18/ icntr,vcntr,xcinc,ncsmth,xclab
      common /blk19/ iq,q,qa
      common /blk20/ iwave,nccbnd,iccbnd,icbnd,vsvp
      common /blk21/ denc,denmin
      common /blk22/ isrch,tang
      common /blk23/ icomp,iep
      common /blk24/ ntray,ntpts
      common /blk25/ nstepr
      common /blk26/ igrid,xgrid,zgrid
      common /blk27/ ibsmth,nbsmth,npbnd,cosmth,xsinc
      common /blk28/ xpf,tpf,upf,ipf,ilshot
      common /blk36/ izrefl,nzf,npzf,xzf,zff
      common /blktan/ mtan,btan,mcotan,bcotan,factan
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
