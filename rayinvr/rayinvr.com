c      
c     version 1.3  Aug 1992
c
c     common blocks for RAYINVR
c           
c     ----------------------------------------------------------------
c
      integer refll(prefl+1),nblk(player),ivarz(player,ppcntr),
     +        idray(2),partyp(pnvar),iheadf(player),nzed(pncntr),
     +        icbnd(pconv+1),ivg(player,ptrap),izv(player,ptrap,4),
     +        ivv(player,ptrap,4),ivray(prayf),ipinv(prayi),
     +        nvel(player,2),ivarv(player,ppvel,2),ldvmax(player),
     +        ldsmax(pncntr),ipf(prayi+pshot2+1),ilshot(pshot2+1),
     +        icalc(prayi),npfref(pfrefl),colour(pcol),mcol(5),
     +        ivarf(pfrefl,ppfref),ircalc(prayi),
     +        nbnda(piray),sample(pzgrid,pxgrid)
      real*4 c(player,ptrap,11),s(player,ptrap,2),b(player,ptrap,2),
     +       vm(player,ptrap,4),xbnd(player,ptrap,2),xshtar(pray),
     +       xm(pncntr,ppcntr),zm(pncntr,ppcntr),fidarr(pray),
     +       vf(player,ppvel,2),xr(ppray),zr(ppray),parorg(pnvar),
     +       ar(ppray,2),vr(ppray,2),tr(ppray),tfinv(pnrayf),
     +       vp(ppray,2),vs(ppray,2),rayid(pray),xfinv(pnrayf),
     +       ray(prayf),vsvp(player,ptrap),
     +       cv(player,ptrap,4,5),tt(pray),cz(player,ptrap,4,2),
     +       tang(player,4),cosmth(pncntr,pnsmth),range(pray),
     +       apart(prayi,pnvar),fpart(pnrayf,pnvar),
     +       tobs(prayi),uobs(prayi),tcalc(prayi),xvel(player,ppvel,2),
     +       xpf(prayi+pshot2+1),tpf(prayi+pshot2+1),
     +       upf(prayi+pshot2+1),xcalc(prayi),xscalc(prayi),
     +       xfrefl(pfrefl,ppfref),zfrefl(pfrefl,ppfref),
     +       mtan(pitan2),btan(pitan2),mcotan(pitan2),bcotan(pitan2)
      character title*80
c
      common /blk1/ layer,iblk,id,fid,fid1
      common /blk2/ c,ivg
      common /blk3/ s,b,vm
      common /blk4/ xbnd,nblk,nlayer
      common /blk5/ xm,zm,vf,nzed,nvel,xvel
      common /blk6/ xr,zr,ar,vr,tr,vp,vs
      common /blk7/ refll,ircbnd
      common /blk8/ irkc,tol,hdenom,hmin,idump,isrkc,ifast
      common /blk9/ itx,vred,time,timer
      common /blk10/ step,smin,smax
      common /blk11/ range,tt,rayid,xshtar,fidarr,idray,ntt,ray
      common /blk12/ xmin,xmax,xtmin,xtmax,xmm,ndecix,xscale,ntickx,
     +          xmint,xmaxt,xtmint,xtmaxt,xmmt,ndecxt,xscalt,ntckxt,
     +          zmin,zmax,ztmin,ztmax,zmm,ndeciz,zscale,ntickz,
     +          tmin,tmax,ttmin,ttmax,tmm,ndecit,tscale,ntickt,
     +          symht,albht,iplots,orig,sep,
     +          title,ititle,xtitle,ytitle,
     +          ircol,mcol,irrcol,itcol,ncol,colour
      common /blk13/ icasel,dstepf,n2,n3,istop,nstepr
      common /blk14/ iwave,nccbnd,iccbnd,icbnd,vsvp
      common /blk15/ isrch,tang
      common /blk16/ ntray,ntpts
      common /blk17/ ibsmth,nbsmth,npbnd,cosmth,xsinc
      common /blk18/ nptbnd
      common /blk19/ ivarz,ivarv,ivarf,ivv,izv,cz,cv,partyp,nvar,
     +               parorg
      common /blk20/ apart,fpart,tobs,uobs,tcalc,xcalc,xscalc,
     +               ipinv,ivray,narinv,icalc,ircalc
      common /blk21/ ninv,xfinv,tfinv
      common /blk22/ iheadf,hws,crit,dhw,tdhw,ihdwf,ihdw,nhskip,idiff,
     +               idifff
      common /blk23/ dvmax,dsmax,idvmax,idsmax,ldvmax,ldsmax
      common /blk24/ xpf,tpf,upf,ipf,ilshot
      common /blk25/ nfrefl,npfref,xfrefl,zfrefl,ifcbnd
      common /blk26/ sample
      common /blk27/ nbnd,nbnda,npskip,npskp
      common /blktan/ mtan,btan,mcotan,bcotan,factan
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
