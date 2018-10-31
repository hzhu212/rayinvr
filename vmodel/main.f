c                 
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             **********  V M O D E L  **********              |
c     |                                                              |
c     |          Check and edit the velocity model used as           |
c     |           input by the programs RAYINVR and TRAMP            |
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
c        10 -- input:  namelist parameters
c
c        11 -- input:  initial velocity model
c                 
c        12 -- output: final velocity model
c                 
c        13 -- output: warning/error messages and summary info
c
c        14 -- output: layer boundaries versus time and Vrms profile
c
c        19 -- output: all Calcomp plot calls
c
c        30 -- input:  floating reflectors
c                 
c
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'vmodel.par'
c                
      real xmn(ppcntr),zmn(ppcntr),grad(player*ppvel),
     +     azsmth(ppcntr),avsmth(ppvel),thick(player*ppcntr),
     +     xveln(ppvel),vfn(ppvel),xlvz(ppvel),dslope(pncntr),
     +     xps(pncntr,2),dvv(player),xpv(player),dlv(player,2),
     +     xlv(player,2,2),xsortz(ppcntr),zsort(ppcntr),xsortv(ppvel),
     +     vsort(ppvel),xvp(pnvp),xzint(pncntr),xvuint(pncntr),
     +     xvlint(pncntr)
      integer ivarzn(ppcntr),ivarvn(ppvel),
     +        ivarz(player,ppcntr),ivarv(player,ppvel,2),iflagg(player),
     +        mpinch(pncntr,ppcntr,2),igrad(player*ppvel),
     +        izintr(pncntr),ivintr(pinvel),iflagt(pncntr),
     +        nzintr(pncntr),nvintr(pinvel),vrmstl,vrmsbl,
     +        izsmth(pncntr),ivsmth(pinvel),isortz(ppcntr),
     +        nzsmth(pncntr),nvsmth(pinvel),isortv(ppvel),
     +        izint(pncntr),ivuint(player),ivlint(player),
     +        nzint(pncntr),nvuint(player),nvlint(player),
     +        izsmt(pncntr),ivusmt(player),ivlsmt(player),
     +        nzsmt(pncntr),nvusmt(player),nvlsmt(player)
      character charr*1,xfile*13
c                 
      include 'vmodel.com'
c
      namelist /pltpar/ iplot,idash,ivel,velht,imod,idump,xsinc,
     +                  xvp,ivp,izort,xtinc,ntsmth,iroute,symht,
     +                  xrinc,vrmstl,vrmsbl,nrsmth,ivrms,ifrefl,
     +                  ibcol,ifcol,iseg,inode,idcol,ivcol,izrefl,
     +                  ntsmt,npsmt,ndot,dashlf,xwndow,ywndow,iout,
     +                  vtime
      namelist /axepar/ xmin,xmax,xtmin,xtmax,xmm,ndecix,ntickx,
     +                  zmin,zmax,ztmin,ztmax,zmm,ndeciz,ntickz,
     +                  vmin,vmax,vtmin,vtmax,vmm,ndeciv,ntickv,
     +                  tmin,tmax,ttmin,ttmax,tmm,ndecit,ntickt,
     +                  vrmin,vrmax,vrtmin,vrtmax,vrmm,ntckvr,ndecir,
     +                  iaxlab,albht,orig
      namelist /modpar/ izsmt,ivusmt,nzsmt,nvusmt,ivlsmt,nvlsmt,
     +                  izint,ivuint,nzint,nvuint,ivlint,nvlint,
     +                  velmin,velmax,dvvmax,dlvmax,dsmax,
     +                  iorder,iswit,xadd,xzint,xvuint,xvlint
c     
      data izint/pncntr*0/,ivuint/player*0/,ivlint/player*0/,
     +     nzint/pncntr*-1/,nvuint/player*-1/,nvlint/player*-1/,
     +     izsmt/pncntr*0/,ivusmt/player*0/,ivlsmt/player*0/,
     +     nzsmt/pncntr*-1/,nvusmt/player*-1/,nvlsmt/player*-1/,
     +     izintr/pncntr*0/,ivintr/pinvel*0/,
     +     nzintr/pncntr*0/,nvintr/pinvel*0/,
     +     izsmth/pncntr*0/,ivsmth/pinvel*0/,
     +     nzsmth/pncntr*0/,nvsmth/pinvel*0/,
     +     xvp/pnvp*-99999./
c
      open(unit=10, file='vm.in', status='old')
      open(unit=11, file='v.in', status='old')
      open(unit=13, file='vm.out')
c                 
c     default parameter values
c
      vtime=-1.
      iout=0 
      xadd=-1.e20
      ndot=3
      dashlf=200.
      izrefl=0
      ntsmt=0
      npsmt=3
      idcol=-1
      ivcol=-1
      inode=0
      iroute=1
      symht=1.
      ifrefl=0
      idump=0
      izort=0
      imod=0
      ivp=0
      xtinc=-1.
      xsinc=-1.
      ntsmth=0
      xrinc=-1.
      vrmstl=-1
      vrmsbl=-1
      nrsmth=0
      ivrms=0
      iaxlab=1
      idash=1
      ivel=0
      velht=.7
      velmin=0.1
      velmax=8.
      dvvmax=1.0
      dlvmax=0.1
      dsmax=10.
      iorder=0
      iflagw=0
      iswit=1
c
      read(10,pltpar)
      read(10,axepar)
      read(10,modpar)
c
      if(xmax.lt.-99998) then
        write(6,5)
        write(13,5)
5       format(/'***  xmax not specified  ***'/)
        stop
      end if
c
      write(6,555) xmax
      write(13,555) xmax
555   format('xmax = ',f7.2,' km')
c
c     assign values to the arrays nzint and nvint
c
      do 1090 i=1,pncntr
         if(nzint(i).gt.ppcntr) then
           write(6,25) 
           write(13,25) 
25         format(/'***  a value of nzint > ppcntr  ***'/)
           stop     
         end if
1090  continue
c
      do 1110 i=1,player
         if(nvuint(i).gt.ppvel) then
           write(6,35) 
           write(13,35) 
35         format(/'***  a value of nvuint > ppvel  ***'/)
           stop     
         end if
         if(nvlint(i).gt.ppvel) then
           write(6,36) 
           write(13,36) 
36         format(/'***  a value of nvlint > ppvel  ***'/)
           stop     
         end if
1110  continue
c
      if(nzint(1).lt.0) then
        do 1050 i=1,pncntr
           nzint(i)=0
1050    continue 
      else
        if(nzint(2).lt.0) then
          do 1060 i=2,pncntr
             nzint(i)=nzint(1) 
1060      continue
        end if
      end if
c
      if(nvuint(1).lt.0) then
        do 1070 i=1,player
           nvuint(i)=0
1070    continue 
      else
        if(nvuint(2).lt.0) then
          do 1080 i=2,player
             nvuint(i)=nvuint(1) 
1080      continue
        end if
      end if
c
      if(nvlint(1).lt.0) then
        do 1071 i=1,player
           nvlint(i)=0
1071    continue 
      else
        if(nvlint(2).lt.0) then
          do 1081 i=2,player
             nvlint(i)=nvlint(1) 
1081      continue
        end if
      end if
c
c     assign values to the arrays nzsmth and nvsmth
c
      if(nzsmt(1).lt.0) then
        do 2050 i=1,pncntr
           nzsmt(i)=1
2050    continue 
      else
        if(nzsmt(2).lt.0) then
          do 2060 i=2,pncntr
             nzsmt(i)=nzsmt(1) 
2060      continue
        end if
      end if
c
      if(nvusmt(1).lt.0) then
        do 2070 i=1,player
           nvusmt(i)=1
2070    continue 
      else
        if(nvusmt(2).lt.0) then
          do 2080 i=2,player
             nvusmt(i)=nvusmt(1) 
2080      continue
        end if
      end if
c
      if(nvlsmt(1).lt.0) then
        do 2071 i=1,player
           nvlsmt(i)=1
2071    continue 
      else
        if(nvlsmt(2).lt.0) then
          do 2081 i=2,player
             nvlsmt(i)=nvlsmt(1) 
2081      continue
        end if
      end if
c
c     asign values to the arrays izintr and nzintr
c
      do 2710 i=1,pncntr
         if(izint(i).ne.0.and.abs(izint(i)).lt.pncntr) then
           izintr(abs(izint(i)))=sign(1,izint(i))
           nzintr(abs(izint(i)))=nzint(i)
         end if
2710  continue
c
c     asign values to the arrays ivintr and nvintr
c
      do 2720 i=1,player
         if(ivuint(i).ne.0.and.abs(ivuint(i)).lt.player) then
           ivintr(2*abs(ivuint(i))-1)=sign(1,ivuint(i))
           nvintr(2*abs(ivuint(i))-1)=nvuint(i)
         end if
         if(ivlint(i).ne.0.and.abs(ivlint(i)).lt.player) then
           ivintr(2*abs(ivlint(i)))=sign(1,ivlint(i))
           nvintr(2*abs(ivlint(i)))=nvlint(i)
         end if
2720  continue
c
c     asign values to the arrays izsmth and nzsmth
c
      do 2730 i=1,pncntr
         if(izsmt(i).gt.0.and.izsmt(i).lt.pncntr) then
           izsmth(izsmt(i))=1
           nzsmth(izsmt(i))=nzsmt(i)
         end if
2730  continue
c
c     asign values to the arrays ivsmth and nvsmth
c
      do 2740 i=1,player
         if(ivusmt(i).gt.0.and.ivusmt(i).lt.player) then
           ivsmth(2*ivusmt(i)-1)=1
           nvsmth(2*ivusmt(i)-1)=nvusmt(i)
         end if
         if(ivlsmt(i).gt.0.and.ivlsmt(i).lt.player) then
           ivsmth(2*ivlsmt(i))=1
           nvsmth(2*ivlsmt(i))=nvlsmt(i)
         end if
2740  continue
c
c     check for the correct number of lines
c
      nline=0
101   read(11,995,end=990) charr
995   format(a1)
      nline=nline+1
      go to 101
990   if(nline.lt.11.or.mod(nline,3).ne.2) then
        write(6,9995)
        write(13,9995)
9995    format(/'***  v.in has incorrect number of lines  ***'/)
        stop
      end if
c
      rewind(11)
c                 
c     read in velocity model
c
      ncont=1
      nrzmax=ppcntr/10
      nrvmax=ppvel/10
      do 170 icont=1,player+1
         nrz=1
         j1=1
         j2=10
11       if(nrz.gt.nrzmax) go to 211
         read(11,55,end=999) ilyr,(xm(icont,j),j=j1,j2)
         read(11,55,end=999) icnt,(zm(icont,j),j=j1,j2)
         read(11,65,end=99) (ivarz(icont,j),j=j1,j2)
55       format(i2,1x,10f7.2)
c65       format(3x,10i7)
65       format(3x,10(5x,i2))
         nrz=nrz+1
         if(icnt.ne.1) go to 211
         j1=j1+10
         j2=j2+10
         go to 11
211      nrv=1
         j1=1
         j2=10
21       if(nrv.gt.nrvmax) go to 311
         read(11,55,end=999) ilyr,(xvel(icont,j,1),j=j1,j2)
         read(11,55,end=999) icnt,(vf(icont,j,1),j=j1,j2)
         read(11,65,end=999) (ivarv(icont,j,1),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 311
         j1=j1+10
         j2=j2+10
         go to 21
311      nrv=1
         j1=1
         j2=10
31       if(nrv.gt.nrvmax) go to 411
         read(11,55,end=999) ilyr,(xvel(icont,j,2),j=j1,j2)
         read(11,55,end=999) icnt,(vf(icont,j,2),j=j1,j2)
         read(11,65,end=999) (ivarv(icont,j,2),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 411
         j1=j1+10
         j2=j2+10
         go to 31
411      ncont=ncont+1
170   continue    
c               
99    nlayer=ncont-1
c
      write(6,505) nlayer
      write(13,505) nlayer
505   format('number of layers = ',i3)
c
      do 171 i=1,ncont 
         nzed(i)=1
171   continue
      do 172 i=1,nlayer
         nvel(i,1)=1
         nvel(i,2)=1
172   continue
c
      do 180 i=1,ncont
         do 190 j=1,ppcntr 
            if(abs(xm(i,j)-xmax).lt..0001) go to 180
            nzed(i)=nzed(i)+1
190      continue 
180   continue    
c
      do 210 i=1,nlayer  
         do 220 j=1,ppvel  
            if(abs(xvel(i,j,1)-xmax).lt..0001) go to 212
            nvel(i,1)=nvel(i,1)+1
220      continue 
212      if(nvel(i,1).eq.1.and.vf(i,1,1).eq.0.) nvel(i,1)=0
210   continue    
c
      do 240 i=1,nlayer  
         do 250 j=1,ppvel  
            if(abs(xvel(i,j,2)-xmax).lt..0001) go to 260
            nvel(i,2)=nvel(i,2)+1
250      continue 
260      if(nvel(i,2).eq.1.and.vf(i,1,2).eq.0.) nvel(i,2)=0
240   continue    
c
c     check that boundary nodes are listed left to right
c
      iflag=0
      iflags=0
      do 1180 i=1,ncont
         if(nzed(i).eq.ppcntr+1) then
           write(6,1175) i
           write(13,1175) i
1175       format('***  no node specified at xmax for boundary ',
     +            i3,'  ***')
           iflags=iflags+1
           go to 1180
         end if
c
         if(nzed(i).gt.1) then
c
           call dxmin(ncont)
c
           if(abs(xm(i,1)-xmin).gt..0001) then
             write(6,1165) i
             write(13,1165) i
1165         format(
     +         '***  first node not specified at xmin for boundary ',
     +         i3,'  ***')
             iflag=1
           end if
         end if
c
         if(nzed(i).eq.1) then
           if(abs(xm(i,1)-xmax).gt..0001) then
             write(6,1185) i
             write(13,1185) i
1185         format(
     +         '***  single node not specified at xmax for boundary ',
     +         i3,'  ***')
             iflag=1
           end if
         else
           do 1190 j=1,nzed(i)-1
              if(xm(i,j+1).le.xm(i,j)) then
                if(iorder.eq.0) then
                  write(6,1195) xm(i,j+1),i
                  write(13,1195) xm(i,j+1),i
1195              format('***  nodes not specified left to right at ',
     +                   f7.2,' km for boundary ',i3,'  ***') 
                  iflag=1
                else
                  iflagw=1
                  do 1191 k=1,nzed(i)
                     xsortz(k)=xm(i,k)
                     zsort(k)=zm(i,k)
                     isortz(k)=ivarz(i,k)
1191              continue
c
                  call sort(xsortz,zsort,isortz,nzed(i)) 
c
                  do 1192 k=1,nzed(i)
                     xm(i,k)=xsortz(k)
                     zm(i,k)=zsort(k)
                     ivarz(i,k)=isortz(k)
1192              continue
                  go to 1180
                end if
              end if
1190       continue
         end if
1180  continue
c
c     check that velocity nodes are listed left to right
c
      do 2180 i=1,nlayer
         if(nvel(i,1).eq.ppvel+1) then
           write(6,2175) i
           write(13,2175) i
2175       format('***  no upper vel specified at xmax for layer ',
     +            i3,'  ***')
           iflags=iflags+1
           go to 2181
         end if
c
         if(nvel(i,1).gt.1) then
c
           call dxmin(ncont)
c
           if(abs(xvel(i,1,1)-xmin).gt..0001) then
             write(6,2165) i
             write(13,2165) i
2165         format(
     +         '***  first upper vel not specified at xmin for layer ',
     +         i3,'  ***')
             iflag=1
           end if
         end if
c
         if(nvel(i,1).eq.1) then
           if(abs(xvel(i,1,1)-xmax).gt..0001) then
             write(6,2185) i
             write(13,2185) i
2185         format(
     +         '***  single upper vel not specified at xmax for layer',
     +         i3,'  ***')
             iflag=1
           end if
         else
           do 2190 j=1,nvel(i,1)-1
              if(xvel(i,j+1,1).lt.xvel(i,j,1)) then
                if(iorder.eq.0) then
                  write(6,2195) xvel(i,j+1,1),i
                  write(13,2195) xvel(i,j+1,1),i
2195              format(
     +              '***  upper vel not specified left to right at ',
     +              f7.2,' km for layer ',i3,'  ***') 
                  iflag=1
                else
                  iflagw=1
                  do 2191 k=1,nvel(i,1)
                     xsortv(k)=xvel(i,k,1)
                     vsort(k)=vf(i,k,1)
                     isortv(k)=ivarv(i,k,1)
2191              continue
c
                  call sort(xsortv,vsort,isortv,nvel(i,1))
c
                  do 2192 k=1,nvel(i,1)
                     xvel(i,k,1)=xsortv(k)
                     vf(i,k,1)=vsort(k)
                     ivarv(i,k,1)=isortv(k)
2192              continue
                  go to 2181
                end if
              end if
2190       continue
         end if
c
2181     if(nvel(i,2).eq.ppvel+1) then
           write(6,2375) i
           write(13,2375) i
2375       format('***  no lower vel specified at xmax for layer ',
     +            i3,'  ***')
           iflags=iflags+1
           go to 2180
         end if
c
         if(nvel(i,2).gt.1) then
c
           call dxmin(ncont)
c
           if(abs(xvel(i,1,2)-xmin).gt..0001) then
             write(6,2265) i
             write(13,2265) i
2265         format(
     +         '***  first lower vel not specified at xmin for layer ',
     +         i3,'  ***')
             iflag=1
           end if
         end if
c
         if(nvel(i,2).eq.1) then
           if(abs(xvel(i,1,2)-xmax).gt..0001) then
             write(6,2395) i
             write(13,2395) i
2395         format(
     +         '***  single lower vel not specified at xmax for layer',
     +         i3,'  ***')
             iflag=1
           end if
         else
           do 2290 j=1,nvel(i,2)-1
              if(xvel(i,j+1,2).lt.xvel(i,j,2)) then
                if(iorder.eq.0) then
                  write(6,2295) xvel(i,j+1,2),i
                  write(13,2295) xvel(i,j+1,2),i
2295              format(
     +              '***  lower vel not specified left to right at ',
     +              f7.2,' km for layer ',i3,'  ***') 
                  iflag=1
                else
                  iflagw=1
                  do 2291 k=1,nvel(i,2)
                     xsortv(k)=xvel(i,k,2)
                     vsort(k)=vf(i,k,2)
                     isortv(k)=ivarv(i,k,2)
2291              continue
c
                  call sort(xsortv,vsort,isortv,nvel(i,2))
c
                  do 2292 k=1,nvel(i,2)
                     xvel(i,k,2)=xsortv(k)
                     vf(i,k,2)=vsort(k)
                     ivarv(i,k,2)=isortv(k)
2292              continue
                  go to 2180
                end if
              end if
2290       continue
         end if
2180  continue
c
      if(iflags.gt.1) then
        write(6,2500)
        write(13,2500)
2500    format(/'**********************************'/
     +          '***  check for incorrect xmax  ***'/
     +          '**********************************'/)
      end if
c
      if(iflags.gt.0.or.iflag.eq.1) stop
c
c     check for fixed layer thicknesses
c     
      do 471 i=2,ncont
         iflagt(i)=0
         do 481 j=1,nzed(i)
            if(ivarz(i,j).eq.-1) iflagt(i)=1
481      continue
471   continue
c
c     check for fixed velocity gradients
c
      do 371 i=1,nlayer
         iflagg(i)=0
         if(nvel(i,2).gt.0) then
           do 381 j=1,nvel(i,2)
              if(ivarv(i,j,2).eq.-1) iflagg(i)=1
381        continue
         end if
371   continue
c
c     interpolate the layer boundaries
c
      do 3010 i=1,ncont
c
         if(iflagt(i).eq.1) then
           if(abs(izint(i-1)).eq.1) then
             n1=nzint(i-1)
           else
             n1=nzed(i-1)
           end if
           if(abs(izint(i)).eq.1) then
             n2=nzint(i)
           else
             n2=nzed(i)
           end if
           if(n1.ne.n2) then
             write(6,45) i-1
             write(13,45) i-1
45           format(/'***  unequal # of bnd nodes for layer ',i3,
     +               '  ***'/)
             stop
           end if
           do 3110 j=1,nzed(i)
              if(abs(xm(i-1,j)-xm(i,j)).gt..0005) then
                write(6,3105) i-1
                write(13,3105) i-1
3105            format(/'***  different bnd node positions for layer ',
     +                  i3,'  ***'/)
                stop
              end if
3110       continue
         end if
c
         if(izintr(i).ne.0.and.nzintr(i).gt.0.and.
     +   nzed(i).ne.nzintr(i)) then
           iflagw=1
           if(nzintr(i).eq.1) then
             sumz=0.
             do 3020 j=1,nzed(i)
                sumz=sumz+zm(i,j)
3020         continue
             avez=sumz/float(nzed(i))
             nzed(i)=1
             xm(i,1)=xmax
             zm(i,1)=avez
           else
c
             call dxmin(ncont)
c
             xinc=(xmax-xmin)/float(nzintr(i)-1)
             do 3030 j=1,nzintr(i)
                if(izintr(i).lt.0) then
                  xmn(j)=xzint(j)
                else
                  xmn(j)=xmin+xinc*float(j-1)
                end if
                if(nzed(i).eq.1) then
                  zmn(j)=zm(i,1)
                  ivarzn(j)=ivarz(i,1)
                else
                  do 3040 k=1,nzed(i)-1
                     if(xmn(j).ge.xm(i,k).and.xmn(j).le.xm(i,k+1)) then
                       zmn(j)=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-xm(i,k)))
     +                 *(xmn(j)-xm(i,k))+zm(i,k)
                       if(ivarz(i,k).eq.ivarz(i,k+1)) then
                         ivarzn(j)=ivarz(i,k)
                       else
                         if(iswit.eq.0) then
                           ivarzn(j)=0
                         else
                           if(ivarz(i,k).eq.1.or.ivarz(i,k+1).eq.1) 
     +                     then
                             ivarzn(j)=1
                           else
                             ivarzn(j)=-1
                           end if
                         end if
                       end if  
                     end if
3040              continue
                end if
3030         continue
             nzed(i)=nzintr(i)
             do 3050 j=1,nzed(i)
                xm(i,j)=xmn(j)
                zm(i,j)=zmn(j)
                ivarz(i,j)=ivarzn(j)
3050         continue
             xm(i,nzed(i))=xmax
           end if
         end if
3010  continue
c
      do 4010 i=1,nlayer
c
         if(iflagg(i).eq.1) then
           if(i.eq.1.or.nvel(i,1).gt.0) then
             ilyr=i
             ib=1
             ipos=2*i-1
           else 
             do 4110 j=i-1,1,-1
                if(nvel(j,2).gt.0) then
                  ilyr=j
                  ib=2
                  ipos=2*j
                  go to 4100
                end if
                if(nvel(j,1).gt.0) then
                  ilyr=j
                  ib=1
                  ipos=2*j-1
                  go to 4100
                end if
4110         continue                
           end if
4100       if(abs(ivintr(ipos)).eq.1) then
             n1=nvintr(ipos)
           else
             n1=nvel(ilyr,ib)
           end if
           if(abs(ivintr(2*i)).eq.1) then
             n2=nvintr(2*i)
           else
             n2=nvel(i,2)
           end if
           if(n1.ne.n2) then
             write(6,15) i
             write(13,15) i
15           format(/'***  unequal # of vel nodes in layer ',i3,
     +               '  ***'/)
             stop
           end if
           do 4210 j=1,nvel(i,2)
              if(abs(xvel(ilyr,j,ib)-xvel(i,j,2)).gt..0005) then
                write(6,4105) i
                write(13,4105) i
4105            format(/'***  different vel node positions for layer ',
     +                  i3,'  ***'/)
                stop
              end if
4210       continue
         end if 
c
c        interpolate the upper layer velocity values
c
         if(ivintr(2*i-1).ne.0.and.nvintr(2*i-1).gt.0.and.
     +   nvel(i,1).ne.nvintr(2*i-1).and.nvel(i,1).gt.0) then
           iflagw=1
           if(nvintr(2*i-1).eq.1) then
             sumv=0.
             do 4020 j=1,nvel(i,1)
                sumv=sumv+vf(i,j,1)
4020         continue
             avev=sumv/float(nvel(i,1))
             nvel(i,1)=1
             xvel(i,1,1)=xmax
             vf(i,1,1)=avev
           else
c
             call dxmin(ncont)
c
             xinc=(xmax-xmin)/float(nvintr(2*i-1)-1)
             do 4030 j=1,nvintr(2*i-1)
                if(ivintr(2*i-1).lt.0) then
                  xveln(j)=xvuint(j)
                else
                  xveln(j)=xmin+xinc*float(j-1)
                end if
                if(nvel(i,1).eq.1) then
                  vfn(j)=vf(i,1,1)
                  ivarvn(j)=ivarv(i,1,1)
                else
                  do 4040 k=1,nvel(i,1)-1
                     if(xveln(j).ge.xvel(i,k,1).and.xveln(j).le. 
     +               xvel(i,k+1,1)) then
                       vfn(j)=((vf(i,k+1,1)-vf(i,k,1))/
     +                 (xvel(i,k+1,1)-xvel(i,k,1)))*(xveln(j)-
     +                 xvel(i,k,1))+vf(i,k,1)
                       if(ivarv(i,k,1).eq.ivarv(i,k+1,1)) then
                         ivarvn(j)=ivarv(i,k,1)
                       else
                         if(iswit.eq.0) then
                           ivarvn(j)=0
                         else
                           if(ivarv(i,k,1).eq.1.or.ivarv(i,k+1,1).
     +                     eq.1) then
                             ivarvn(j)=1
                           else
                             ivarvn(j)=-1
                           end if 
                         end if
                       end if  
                     end if
4040              continue
                end if
4030         continue
             nvel(i,1)=nvintr(2*i-1)
             do 4050 j=1,nvel(i,1)
                xvel(i,j,1)=xveln(j)
                vf(i,j,1)=vfn(j)
                ivarv(i,j,1)=ivarvn(j)
4050         continue
             xvel(i,nvel(i,1),1)=xmax
           end if
         end if
c
c        interpolate the lower layer velocity values
c
         if(ivintr(2*i).ne.0.and.nvintr(2*i).gt.0.and.
     +   nvel(i,2).ne.nvintr(2*i).and.nvel(i,2).gt.0) then
           iflagw=1
           if(nvintr(2*i).eq.1) then
             sumv=0.
             do 4060 j=1,nvel(i,2)
                sumv=sumv+vf(i,j,2)
4060         continue
             avev=sumv/float(nvel(i,2))
             nvel(i,2)=1
             xvel(i,1,2)=xmax
             vf(i,1,2)=avev
           else
c
             call dxmin(ncont)
c
             xinc=(xmax-xmin)/float(nvintr(2*i)-1)
             do 4070 j=1,nvintr(2*i)
                if(ivintr(2*i).lt.0) then
                  xveln(j)=xvlint(j)
                else
                  xveln(j)=xmin+xinc*float(j-1)
                end if
                if(nvel(i,2).eq.1) then
                  vfn(j)=vf(i,1,2)
                  ivarvn(j)=ivarv(i,1,2)
                else
                  do 4080 k=1,nvel(i,2)-1
                     if(xveln(j).ge.xvel(i,k,2).and.xveln(j).le. 
     +               xvel(i,k+1,2)) then
                       vfn(j)=((vf(i,k+1,2)-vf(i,k,2))/
     +                 (xvel(i,k+1,2)-xvel(i,k,2)))*(xveln(j)-
     +                 xvel(i,k,2))+vf(i,k,2)
                       if(ivarv(i,k,2).eq.ivarv(i,k+1,2)) then
                         ivarvn(j)=ivarv(i,k,2)
                       else
                         if(iswit.eq.0) then
                           ivarvn(j)=0
                         else
                           if(ivarv(i,k,2).eq.1.or.ivarv(i,k+1,2).
     +                     eq.1) then
                             ivarvn(j)=1
                           else
                             ivarvn(j)=-1
                           end if
                         end if
                       end if  
                     end if
4080              continue
                end if
4070         continue
             nvel(i,2)=nvintr(2*i)
             do 4090 j=1,nvel(i,2)
                xvel(i,j,2)=xveln(j)
                vf(i,j,2)=vfn(j)
                ivarv(i,j,2)=ivarvn(j)
4090         continue
             xvel(i,nvel(i,2),2)=xmax
           end if
         end if
4010  continue
c
c     add new nodes
c
      if(xadd.gt.xmin.and.xadd.lt.xmax) then
c
      xfile='xadd   .nodes'
      i1=int(xadd/100)
      i2=int((xadd-i1*100)/10+.05)
      i3=int(xadd-i1*100-i2*10+.05)
      xfile(5:7)=char(i1+48)//char(i2+48)//char(i3+48)
      open(57, file=xfile)
c
c     add a new depth node
c
        do 1170 i=1,ncont
           if(nzed(i).gt.1) then
             do j=1,nzed(i)-1
                if(xadd.ge.xm(i,j).and.xadd.le.xm(i,j+1)) then
                  if(abs(xadd-xm(i,j)).lt..001.or.
     +               abs(xadd-xm(i,j+1)).lt..001) go to 1170
                  iflagw=1
                  zn=(zm(i,j+1)-zm(i,j))/(xm(i,j+1)-xm(i,j))*
     +               (xadd-xm(i,j))+zm(i,j)
                  do k=1,j
                     xmn(k)=xm(i,k)
                     zmn(k)=zm(i,k)
                     ivarzn(k)=ivarz(i,k)
                  end do
                  xmn(j+1)=xadd
                  zmn(j+1)=zn
                  ivarzn(j+1)=1
                  do k=j+1,nzed(i)
                     xmn(k+1)=xm(i,k)
                     zmn(k+1)=zm(i,k)
                     ivarzn(k+1)=ivarz(i,k)
                  end do
                  nzed(i)=nzed(i)+1
                  do k=1,nzed(i)
                     xm(i,k)=xmn(k)
                     zm(i,k)=zmn(k)
                     ivarz(i,k)=ivarzn(k)
                  enddo
                  write(57,1775) 1,i,xadd,zn
1775              format(2i10,4f10.3)
                  go to 1170
                end if
             enddo
           end if
1170    continue
c
c     add new upper velocity nodes
c
        do 1171 i=1,nlayer
           if(nvel(i,1).gt.1) then
             do j=1,nvel(i,1)-1
                if(xadd.ge.xvel(i,j,1).and.xadd.le.xvel(i,j+1,1)) 
     +          then
                  if(abs(xadd-xvel(i,j,1)).lt..001.or.
     +               abs(xadd-xvel(i,j+1,1)).lt..001) go to 1171
                  iflagw=1
                  vn=(vf(i,j+1,1)-vf(i,j,1))/(xvel(i,j+1,1)-
     +                xvel(i,j,1))*(xadd-xvel(i,j,1))+vf(i,j,1)
                  do k=1,j
                     xveln(k)=xvel(i,k,1)
                     vfn(k)=vf(i,k,1)
                     ivarvn(k)=ivarv(i,k,1)
                  end do
                  xveln(j+1)=xadd
                  vfn(j+1)=vn
                  ivarvn(j+1)=1
                  do k=j+1,nvel(i,1)
                     xveln(k+1)=xvel(i,k,1)
                     vfn(k+1)=vf(i,k,1)
                     ivarvn(k+1)=ivarv(i,k,1)
                  end do
                  nvel(i,1)=nvel(i,1)+1
                  do k=1,nvel(i,1)
                     xvel(i,k,1)=xveln(k)
                     vf(i,k,1)=vfn(k)
                     ivarv(i,k,1)=ivarvn(k)
                  enddo
                  write(57,1775) 2,i,xadd,vn
                  go to 1171
                end if
             enddo
           end if
1171    continue
c
c     add new lower velocity nodes
c
        do 1172 i=1,nlayer
           if(nvel(i,2).gt.1) then
             do j=1,nvel(i,2)-1
                if(xadd.ge.xvel(i,j,2).and.xadd.le.xvel(i,j+1,2))
     +          then
                  if(abs(xadd-xvel(i,j,2)).lt..001.or.
     +               abs(xadd-xvel(i,j+1,2)).lt..001) go to 1172
                  iflagw=1
                  vn=(vf(i,j+1,2)-vf(i,j,2))/(xvel(i,j+1,2)-
     +                xvel(i,j,2))*(xadd-xvel(i,j,2))+vf(i,j,2)
                  do k=1,j
                     xveln(k)=xvel(i,k,2)
                     vfn(k)=vf(i,k,2)
                     ivarvn(k)=ivarv(i,k,2)
                  end do
                  xveln(j+1)=xadd
                  vfn(j+1)=vn
                  ivarvn(j+1)=1
                  do k=j+1,nvel(i,2)
                     xveln(k+1)=xvel(i,k,2)
                     vfn(k+1)=vf(i,k,2)  
                     ivarvn(k+1)=ivarv(i,k,2)
                  end do  
                  nvel(i,2)=nvel(i,2)+1
                  do k=1,nvel(i,2)
                     xvel(i,k,2)=xveln(k)
                     vf(i,k,2)=vfn(k)
                     ivarv(i,k,2)=ivarvn(k)
                  enddo
                  write(57,1775) 3,i,xadd,vn
                  go to 1172
                end if
             enddo
           end if
1172    continue
      end if
c
c     check for pinchouts
c
      if(nlayer.gt.1) then
        ipinch=0
        do 281 i=2,nlayer
           do 282 j=1,nzed(i)
              mpinch(i,j,1)=0
              do 283 k=1,i-1
                 do 284 l=1,nzed(k)
                  if(abs(xm(i,j)-xm(k,l)).lt..0005.and.
     +                 abs(zm(i,j)-zm(k,l)).lt..0005) then
                       ipinch=+1
                       mpinch(i,j,1)=k
                       mpinch(i,j,2)=l
                       go to 282
                    end if
284              continue
283           continue
282        continue
281     continue
      end if
c
c     check for fixed layer thicknesses
c     
      nlyr=0
      do 470 i=2,ncont
         do 480 j=1,nzed(i)
            if(ivarz(i,j).eq.-1) then
              nlyr=nlyr+1
              thick(nlyr)=zm(i,j)-zm(i-1,j)
            end if
480      continue
470   continue
c
c     check for fixed velocity gradients
c
      ngrad=0
      do 370 i=1,nlayer
         iflagg(i)=0
         if(nvel(i,2).gt.0) then
           do 380 j=1,nvel(i,2)
              if(ivarv(i,j,2).eq.-1) then
                ngrad=ngrad+1
                xbndc=xvel(i,j,2)
                if(nzed(i).gt.1) then
                  do 390 k=1,nzed(i)-1
                     if(xbndc.ge.xm(i,k).and.xbndc.le.xm(i,k+1))
     +               then
                       zu=(zm(i,k+1)-zm(i,k))*(xbndc-xm(i,k))/
     +                    (xm(i,k+1)-xm(i,k))+zm(i,k)
                       go to 400
                     end if
390               continue
                else
                  zu=zm(i,1)
                end if
400             if(nvel(i,1).gt.0) then
                  vu=vf(i,j,1)
                else
                  do 401 k=i-1,1,-1
                     if(nvel(k,2).gt.0) then
                       vu=vf(k,j,2)
                       go to 402
                     end if
                     if(nvel(k,1).gt.0) then
                       vu=vf(k,j,1)
                       go to 402
                     end if
401               continue
                end if
402             if(nzed(i+1).gt.1) then
                  do 410 k=1,nzed(i+1)-1
                     if(xbndc.ge.xm(i+1,k).and.xbndc.le.xm(i+1,k+1))
     +               then
                       zl=(zm(i+1,k+1)-zm(i+1,k))*(xbndc-xm(i+1,k))/
     +                    (xm(i+1,k+1)-xm(i+1,k))+zm(i+1,k)
                       go to 420
                     end if
410               continue
                else
                  zl=zm(i+1,1)    
                end if
420             vl=vf(i,j,2)
                if(abs(vu-vl).gt..001) then
                  if(abs(zl-zu).gt..001) then
                    grad(ngrad)=(vl-vu)/(zl-zu)
                    igrad(ngrad)=1
                  else
                    grad(ngrad)=vl-vu
                    igrad(ngrad)=0
                    iflagg(i)=1
                  end if
                else
                  grad(ngrad)=0.
                  igrad(ngrad)=1
                end if
              end if 
380        continue
         end if
c
         if(iflagg(i).eq.1) then
           write(6,85) i
           write(13,85) i
85         format('***  check gradient in layer ',i3,'  ***')
         end if
c
370   continue
c
c     smooth the layer boundaries
c
      do 50 i=1,ncont
         if(izsmth(i).gt.0.and.nzsmth(i).gt.0.and.nzed(i).gt.2) then
           iflagw=1
           do 60 j=1,nzed(i)
              azsmth(j)=zm(i,j)
60         continue
           do 70 j=1,nzsmth(i)
              call smooth(azsmth,nzed(i))
70         continue 
           do 80 j=1,nzed(i)
              zm(i,j)=azsmth(j)
80         continue
         end if
50    continue
c
c     smooth the velocity values
c
      do 90 i=1,nlayer
         if(ivsmth(2*i-1).gt.0.and.nvsmth(2*i-1).gt.0.
     +   and.nvel(i,1).gt.2) then
           iflagw=1 
           do 100 j=1,nvel(i,1)
              avsmth(j)=vf(i,j,1)
100        continue
           do 110 j=1,nvsmth(2*i-1)
              call smooth(avsmth,nvel(i,1))
110        continue 
           do 120 j=1,nvel(i,1)
              vf(i,j,1)=avsmth(j)
120        continue
         end if
         if(ivsmth(2*i).gt.0.and.nvsmth(2*i).gt.0.
     +   and.nvel(i,2).gt.2) then
           iflagw=1
           do 130 j=1,nvel(i,2)
              avsmth(j)=vf(i,j,2)
130        continue
           do 140 j=1,nvsmth(2*i)
              call smooth(avsmth,nvel(i,2))
140        continue 
           do 150 j=1,nvel(i,2)
              vf(i,j,2)=avsmth(j)
150        continue
         end if
90    continue
c
c     maintain fixed layer thicknesses
c
      nlyr=0
      do 270 i=1,nlayer
         do 280 j=1,nzed(i)
            if(ipinch.eq.1.and.i.gt.1) then
              if(mpinch(i,j,1).gt.0) then
                zm(i,j)=zm(mpinch(i,j,1),mpinch(i,j,2))
                go to 280
              end if
            end if
            if(ivarz(i,j).eq.-1) then
              nlyr=nlyr+1
              zm(i,j)=zm(i-1,j)+thick(nlyr)
            end if
280      continue
270   continue
c
c     check for nvel(i,1)=0 or nvel(i,2)=0
c
      do 670 i=1,nlayer
         if(nvel(i,1).eq.0) then
           nvel(i,1)=1
           xvel(i,1,1)=xmax    
           vf(i,1,1)=0.
         end if
c
         if(nvel(i,2).eq.0) then
           nvel(i,2)=1
           xvel(i,1,2)=xmax    
           vf(i,1,2)=0.
         end if
670   continue
c
c     maintain fixed velocity gradients
c
      ngrad=0
      do 271 i=1,nlayer
         do 272 j=1,nvel(i,2)
            if(ivarv(i,j,2).eq.-1) then
              ngrad=ngrad+1
              xbndc=xvel(i,j,2)
              if(nzed(i).gt.1) then
                do 440 k=1,nzed(i)-1
                   if(xbndc.ge.xm(i,k).and.xbndc.le.xm(i,k+1))
     +             then
                     zu=(zm(i,k+1)-zm(i,k))*(xbndc-xm(i,k))/
     +                  (xm(i,k+1)-xm(i,k))+zm(i,k)
                     go to 450
                   end if
440             continue
              else
                zu=zm(i,1)
              end if
450           if(nvel(i,1).gt.0.and.vf(i,1,1).gt.0.) then
                vu=vf(i,j,1)
              else
                do 451 k=i-1,1,-1
                   if(nvel(k,2).gt.0.and.vf(k,1,2).gt.0.) then
                     vu=vf(k,j,2)
                     go to 452
                   end if
                   if(nvel(k,1).gt.0.and.vf(k,1,1).gt.0.) then
                     vu=vf(k,j,1)
                     go to 452
                   end if
451             continue
              end if
              if(igrad(ngrad).eq.1) then
452             if(nzed(i+1).gt.1) then
                  do 460 k=1,nzed(i+1)-1
                     if(xbndc.ge.xm(i+1,k).and.xbndc.le.xm(i+1,k+1))
     +               then
                       zl=(zm(i+1,k+1)-zm(i+1,k))*(xbndc-xm(i+1,k))/
     +                    (xm(i+1,k+1)-xm(i+1,k))+zm(i+1,k)
                       vf(i,j,2)=vu+grad(ngrad)*(zl-zu)
                       go to 272
                     end if
460               continue
                else
                  zl=zm(i+1,1)
                  vf(i,j,2)=vu+grad(ngrad)*(zl-zu)
                end if
              else
                vf(i,j,2)=vu+grad(ngrad)
              end if
            end if
272      continue
271   continue
c
c     check for small and large velocities
c
      do 8810 i=1,nlayer
         if(nvel(i,1).eq.1.and.vf(i,1,1).eq.0.) then
           if(i.eq.1) then
             write(6,185)
             write(13,185)
185          format(/'***  zero velocity at top of model  ***'/)
             stop
           end if
           do 8860 j=i-1,1,-1 
              if(nvel(j,2).gt.1.or.vf(j,1,2).ne.0.) then
                ilyr=j
                ib=2
                go to 8861
              end if
              if(nvel(j,1).gt.1.or.vf(j,1,1).ne.0.) then
                ilyr=j
                ib=1
                go to 8861
              end if
8860       continue
         else
           ilyr=i
           ib=1
         end if
c
8861     do 8820 j=1,nvel(ilyr,ib)
            va=vf(ilyr,j,ib)
            if(va.gt.velmin.and.va.lt.velmax) go to 8820
            xa=xvel(ilyr,j,ib)
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 8830 k=1,nzed(i)-1
                 if(xa.ge.xm(i,k).and.xa.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xa-xm(i,k))+zm(i,k)
                   go to 8831
                 end if
8830          continue
            end if
8831        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 8840 k=1,nzed(i+1)-1
                 if(xa.ge.xm(i+1,k).and.xa.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xa-xm(i+1,k))+zm(i+1,k)
                   go to 8832
                 end if
8840          continue
            end if
8832        if(abs(zb-za).gt..0005) then
              if(va.le.velmin) then
                write(6,8805) i,xa,va
                write(13,8805) i,xa,va
8805            format('layer ',i3,' small upper velocity at ',
     +                  f7.2,' km (',f7.2,' km/s)')
              else
                write(6,8815) i,xa,va
                write(13,8815) i,xa,va
8815            format('layer ',i3,' large upper velocity at ',
     +                  f7.2,' km (',f7.2,' km/s)')
              end if
            end if
8820     continue
c
         if(nvel(i,2).eq.1.and.vf(i,1,2).eq.0.) go to 8810
c
         do 9820 j=1,nvel(i,2)
            vb=vf(i,j,2)
            if(vb.gt.velmin.and.vb.lt.velmax) go to 9820
            xb=xvel(i,j,2)
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 9830 k=1,nzed(i)-1
                 if(xb.ge.xm(i,k).and.xb.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xb-xm(i,k))+zm(i,k)
                   go to 9831
                 end if
9830          continue
            end if
9831        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 9840 k=1,nzed(i+1)-1
                 if(xb.ge.xm(i+1,k).and.xb.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xb-xm(i+1,k))+zm(i+1,k)
                   go to 9832
                 end if
9840          continue
            end if
9832        if(abs(zb-za).gt..0005) then
              if(vb.le.velmin) then
                write(6,9805) i,xb,vb
                write(13,9805) i,xb,vb
9805            format('layer ',i3,' small lower velocity at ',
     +                  f7.2,' km (',f7.2,' km/s)')
              else
                write(6,9815) i,xb,vb
                write(13,9815) i,xb,vb
9815            format('layer ',i3,' large lower velocity at ',
     +                  f7.2,' km (',f7.2,' km/s)')
              end if
            end if
9820     continue
c
8810   continue
c
c     check for boundaries crossing
c
      do 710 i=2,ncont
         do 720 j=1,nzed(i)
            xb=xm(i,j)
            zb=zm(i,j)
            if(nzed(i-1).eq.1) then
              za=zm(i-1,1)
            else
              do 730 k=1,nzed(i-1)-1
                 if(xb.ge.xm(i-1,k).and.xb.le.xm(i-1,k+1)) then
                   za=((zm(i-1,k+1)-zm(i-1,k))/(xm(i-1,k+1)-
     +                  xm(i-1,k)))*(xb-xm(i-1,k))+zm(i-1,k)
                   go to 731
                 end if
730           continue
            end if
731         if(za.gt.zb) then
              if(nzed(i).eq.1.and.nzed(i-1).eq.1) then
                write(6,715) i,i-1
                write(13,715) i,i-1
715             format('***  boundary ',i3,' is above boundary ',
     +                  i3,'  ***')
                go to 720
              else
                write(6,705) i,i-1,xb
                write(13,705) i,i-1,xb
705             format('***  boundary ',i3,' crosses boundary ',
     +                  i3,' at ',f7.2,' km  ***')
                go to 720
              end if
            end if
720      continue
710   continue
      do 9110 i=1,ncont-1
         do 9220 j=1,nzed(i)
            xa=xm(i,j)    
            za=zm(i,j)
            if(nzed(i+1).eq.1) then
              zb=zm(i+1,1) 
            else 
              do 9330 k=1,nzed(i+1)-1
                 if(xa.ge.xm(i+1,k).and.xa.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xa-xm(i+1,k))+zm(i+1,k)
                   go to 9331
                 end if    
9330          continue
            end if   
9331        if(za.gt.zb) then
              if(nzed(i).eq.1.and.nzed(i+1).eq.1) then
                write(6,715) i,i+1
                write(13,715) i,i+1
                go to 9220  
              else
                write(6,705) i,i+1,xa
                write(13,705) i,i+1,xa
                go to 9220  
              end if    
            end if 
9220     continue
9110  continue
c
c     check for low-velocity zones 
c
      do 810 i=2,nlayer
         if(nvel(i,1).eq.1.and.vf(i,1,1).eq.0.) go to 810
         if(nvel(i-1,2).eq.1.and.vf(i-1,1,2).eq.0.) then
           if(nvel(i-1,1).gt.1.or.vf(i-1,1,1).ne.0.) then
             ilyr=i-1
             ib=1
           else
             do 860 j=i-2,1,-1 
                if(nvel(j,2).gt.1.or.vf(j,1,2).ne.0.) then
                  ilyr=j
                  ib=2
                  go to 861
                end if
                if(nvel(j,1).gt.1.or.vf(j,1,1).ne.0.) then
                  ilyr=j
                  ib=1
                  go to 861
                end if
860          continue
           end if
         else
           ilyr=i-1
           ib=2
         end if
c
861      do 820 j=1,nvel(ilyr,ib)
            xa=xvel(ilyr,j,ib)
            xlvz(j)=xa
            if(nzed(ilyr).eq.1) then
              za=zm(ilyr,1)
            else
              do 4850 k=1,nzed(ilyr)-1
                 if(xa.ge.xm(ilyr,k).and.xa.le.xm(ilyr,k+1)) then
                   za=((zm(ilyr,k+1)-zm(ilyr,k))/(xm(ilyr,k+1)-
     +                  xm(ilyr,k)))*(xa-xm(ilyr,k))+zm(ilyr,k)
                   go to 4851
                 end if
4850          continue
            end if
4851        if(nzed(ilyr+1).eq.1) then
              zb=zm(ilyr+1,1)
            else
              do 4860 k=1,nzed(ilyr+1)-1
                 if(xa.ge.xm(ilyr+1,k).and.xa.le.xm(ilyr+1,k+1)) then
                   zb=((zm(ilyr+1,k+1)-zm(ilyr+1,k))/(xm(ilyr+1,k+1)-
     +                  xm(ilyr+1,k)))*(xa-xm(ilyr+1,k))+zm(ilyr+1,k)
                   go to 4852
                 end if
4860          continue
            end if
4852        if(abs(zb-za).lt..0005) go to 820
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 4870 k=1,nzed(i)-1
                 if(xa.ge.xm(i,k).and.xa.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xa-xm(i,k))+zm(i,k)
                   go to 4871
                 end if
4870          continue
            end if
4871        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 4880 k=1,nzed(i+1)-1
                 if(xa.ge.xm(i+1,k).and.xa.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xa-xm(i+1,k))+zm(i+1,k)
                   go to 4872
                 end if
4880          continue
            end if
4872        if(abs(zb-za).lt..0005) go to 820
            va=vf(ilyr,j,ib)
            if(nvel(i,1).eq.1) then
              vb=vf(i,1,1)
            else
              do 830 k=1,nvel(i,1)-1
                 if(xa.ge.xvel(i,k,1).and.xa.le.xvel(i,k+1,1)) then
                   vb=((vf(i,k+1,1)-vf(i,k,1))/(xvel(i,k+1,1)-
     +                  xvel(i,k,1)))*(xa-xvel(i,k,1))+vf(i,k,1)
                   go to 831
                 end if
830           continue
            end if
831         if(va.gt.vb) then
              if(nvel(ilyr,ib).eq.1.and.nvel(i,1).eq.1) then
                write(6,815) i,va,vb
                write(13,815) i,va,vb
815             format('layer ',i3,' is low-velocity zone (',
     +                 f7.2,' over ',f7.2,' km/s)')
                go to 820
              else
                write(6,805) i,xa,va,vb
                write(13,805) i,xa,va,vb
805             format('layer ',i3,' is low-velocity zone at ',
     +                  f7.2,' km (',f7.2,' over ',f7.2,' km/s)')
                go to 820
              end if
            end if
820      continue
c
         do 920 j=1,nvel(i,1)
            xb=xvel(i,j,1)
            do 940 k=1,nvel(ilyr,ib)
               if(abs(xb-xlvz(k)).lt..0005) go to 920
940         continue
            if(nzed(ilyr).eq.1) then
              za=zm(ilyr,1)
            else
              do 7850 k=1,nzed(ilyr)-1
                 if(xb.ge.xm(ilyr,k).and.xb.le.xm(ilyr,k+1)) then
                   za=((zm(ilyr,k+1)-zm(ilyr,k))/(xm(ilyr,k+1)-
     +                  xm(ilyr,k)))*(xb-xm(ilyr,k))+zm(ilyr,k)
                   go to 7851
                 end if
7850          continue
            end if
7851        if(nzed(ilyr+1).eq.1) then
              zb=zm(ilyr+1,1)
            else
              do 7860 k=1,nzed(ilyr+1)-1
                 if(xb.ge.xm(ilyr+1,k).and.xb.le.xm(ilyr+1,k+1)) then
                   zb=((zm(ilyr+1,k+1)-zm(ilyr+1,k))/(xm(ilyr+1,k+1)-
     +                  xm(ilyr+1,k)))*(xb-xm(ilyr+1,k))+zm(ilyr+1,k)
                   go to 7852
                 end if
7860          continue
            end if
7852        if(abs(zb-za).lt..0005) go to 920
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 7870 k=1,nzed(i)-1
                 if(xb.ge.xm(i,k).and.xb.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xb-xm(i,k))+zm(i,k)
                   go to 7871
                 end if
7870          continue
            end if
7871        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 7880 k=1,nzed(i+1)-1
                 if(xb.ge.xm(i+1,k).and.xb.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xb-xm(i+1,k))+zm(i+1,k)
                   go to 7872
                 end if
7880          continue
            end if
7872        if(abs(zb-za).lt..0005) go to 920
            vb=vf(i,j,1)
            if(nvel(ilyr,ib).eq.1) then
              va=vf(ilyr,1,ib)
            else
              do 930 k=1,nvel(ilyr,ib)
                 if(xb.ge.xvel(ilyr,k,ib).and.xb.le.xvel(ilyr,k+1,ib)) 
     +           then
                   va=((vf(ilyr,k+1,ib)-vf(ilyr,k,ib))/
     +                (xvel(ilyr,k+1,ib)-xvel(ilyr,k,ib)))*(xb-
     +                xvel(ilyr,k,ib))+vf(ilyr,k,ib)
                   go to 931
                 end if
930           continue
931         end if
            if(va.gt.vb) then
              if(nvel(ilyr,ib).eq.1.and.nvel(i,1).eq.1) then
                write(6,815) i,xb,va,vb
                write(13,815) i,xb,va,vb
                go to 920
              else
                write(6,805) i,xb,va,vb
                write(13,805) i,xb,va,vb
                go to 920
              end if
            end if
920      continue
c
810   continue
c
c     check for negative vertical velocity gradients
c
      do 1810 i=1,nlayer
c
         if(nvel(i,2).eq.1.and.vf(i,1,2).eq.0.) go to 1810
         if(nvel(i,1).eq.1.and.vf(i,1,1).eq.0.) then
           if(i.eq.1) then
             write(6,185)
             write(13,185)
             stop
           end if
           do 1860 j=i-1,1,-1 
              if(nvel(j,2).gt.1.or.vf(j,1,2).ne.0.) then
                ilyr=j
                ib=2
                go to 1861
              end if
              if(nvel(j,1).gt.1.or.vf(j,1,1).ne.0.) then
                ilyr=j
                ib=1
                go to 1861
              end if
1860       continue
         else
           ilyr=i
           ib=1
         end if
c
1861     do 1820 j=1,nvel(ilyr,ib)
            xa=xvel(ilyr,j,ib)
            xlvz(j)=xa
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 6830 k=1,nzed(i)-1
                 if(xa.ge.xm(i,k).and.xa.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xa-xm(i,k))+zm(i,k)
                   go to 6831
                 end if
6830          continue
            end if
6831        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 6840 k=1,nzed(i+1)-1
                 if(xa.ge.xm(i+1,k).and.xa.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xa-xm(i+1,k))+zm(i+1,k)
                   go to 6832
                 end if
6840          continue
            end if
6832        if(abs(zb-za).lt..0005) go to 1820
            va=vf(ilyr,j,ib)
            if(nvel(i,2).eq.1) then
              vb=vf(i,1,2)
            else
              do 1830 k=1,nvel(i,2)-1
                 if(xa.ge.xvel(i,k,2).and.xa.le.xvel(i,k+1,2)) then
                   vb=((vf(i,k+1,2)-vf(i,k,2))/(xvel(i,k+1,2)-
     +                  xvel(i,k,2)))*(xa-xvel(i,k,2))+vf(i,k,2)
                   go to 1831
                 end if
1830          continue
            end if
1831        if(va.gt.vb) then
              if(nvel(ilyr,ib).eq.1.and.nvel(i,2).eq.1) then
                write(6,1815) i,va,vb
                write(13,1815) i,va,vb
1815            format('layer ',i3,' has neg vert vel grad (',
     +                 f7.2,' over ',f7.2,' km/s)')
              else
                write(6,1805) i,xa,va,vb
                write(13,1805) i,xa,va,vb
1805            format('layer ',i3,' has neg vert vel grad at ',
     +                  f7.2,' km (',f7.2,' over ',f7.2,' km/s)')
              end if
            end if
1820     continue
c
         do 1920 j=1,nvel(i,2)
            xb=xvel(i,j,2)
            do 1940 k=1,nvel(ilyr,ib)
               if(abs(xb-xlvz(k)).lt..0005) go to 1920
1940        continue
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 7830 k=1,nzed(i)-1
                 if(xb.ge.xm(i,k).and.xb.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xb-xm(i,k))+zm(i,k)
                   go to 7831
                 end if
7830          continue
            end if
7831        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 7840 k=1,nzed(i+1)-1
                 if(xb.ge.xm(i+1,k).and.xb.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xb-xm(i+1,k))+zm(i+1,k)
                   go to 7832
                 end if
7840          continue
            end if
7832        if(abs(zb-za).lt..0005) go to 1920
            vb=vf(i,j,2)
            if(nvel(ilyr,ib).eq.1) then
              va=vf(ilyr,1,ib)
            else
              do 1930 k=1,nvel(ilyr,ib)-1
                 if(xb.ge.xvel(ilyr,k,ib).and.xb.le.xvel(ilyr,k+1,ib)) 
     +           then
                   va=((vf(ilyr,k+1,ib)-vf(ilyr,k,ib))/
     +                (xvel(ilyr,k+1,ib)-xvel(ilyr,k,ib)))*
     +                (xb-xvel(ilyr,k,ib))+vf(ilyr,k,ib)
                   go to 1931
                 end if
1930          continue
            end if
1931        if(va.gt.vb) then
              if(nvel(ilyr,ib).eq.1.and.nvel(i,2).eq.1) then
                write(6,1815) i,va,vb
                write(13,1815) i,va,vb
              else
                write(6,1805) i,xb,va,vb
                write(13,1805) i,xb,va,vb
              end if
            end if
1920     continue
c
1810  continue
c
c     check for large changes in boundary slope
c
      do 2310 i=1,ncont
         dslope(i)=0.
         xps(i,1)=0.
         xps(i,2)=0.
         if(nzed(i).gt.2) then
           do 2320 j=1,nzed(i)-2
              slope1=(zm(i,j+1)-zm(i,j))/(xm(i,j+1)-xm(i,j))
              slope2=(zm(i,j+2)-zm(i,j+1))/(xm(i,j+2)-xm(i,j+1))
              ds1=atan(slope1)*pi18
              ds2=atan(slope2)*pi18
              ds=abs(ds2-ds1)
              if(ds.gt.dsmax) then
                write(6,2325) i,ds,xm(i,j),xm(i,j+2)
                write(13,2325) i,ds,xm(i,j),xm(i,j+2)
2325            format('bnd ',i3,' has ',f7.2,
     +          ' deg slope change at ',f7.2,' - ',f7.2,' km')
              end if
              if(ds.gt.dslope(i)) then
                dslope(i)=ds
                xps(i,1)=xm(i,j)
                xps(i,2)=xm(i,j+2)
              end if
2320       continue
         end if
2310  continue
c
      write(13,2305)
2305  format(/'boundary   slope change (degrees)       between (km)')
      do 2330 i=1,ncont
         write(13,2315) i,dslope(i),xps(i,1),xps(i,2)
2315     format(i4,f19.4,12x,2f10.3)
2330  continue
c
c     check for large vertical velocity gradients
c
      do 3810 i=1,nlayer
c
         dvv(i)=0.
         xpv(i)=0.
         if(nvel(i,2).eq.1.and.vf(i,1,2).eq.0.) go to 3810
         if(nvel(i,1).eq.1.and.vf(i,1,1).eq.0.) then
           if(i.eq.1) then
             write(6,185)
             write(13,185)
             stop
           end if
           do 3860 j=i-1,1,-1 
              if(nvel(j,2).gt.1.or.vf(j,1,2).ne.0.) then
                ilyr=j
                ib=2
                go to 3861
              end if
              if(nvel(j,1).gt.1.or.vf(j,1,1).ne.0.) then
                ilyr=j
                ib=1
                go to 3861
              end if
3860       continue
         else
           ilyr=i
           ib=1
         end if
c
3861     do 3820 j=1,nvel(ilyr,ib)
            xa=xvel(ilyr,j,ib)
            xlvz(j)=xa
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 4830 k=1,nzed(i)-1
                 if(xa.ge.xm(i,k).and.xa.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xa-xm(i,k))+zm(i,k)
                   go to 4831
                 end if
4830          continue
            end if
4831        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 4840 k=1,nzed(i+1)-1
                 if(xa.ge.xm(i+1,k).and.xa.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xa-xm(i+1,k))+zm(i+1,k)
                   go to 4832
                 end if
4840          continue
            end if
4832        if(abs(zb-za).lt..0005) go to 3820
            va=vf(ilyr,j,ib)
            if(nvel(i,2).eq.1) then
              vb=vf(i,1,2)
            else
              do 3830 k=1,nvel(i,2)-1
                 if(xa.ge.xvel(i,k,2).and.xa.le.xvel(i,k+1,2)) then
                   vb=((vf(i,k+1,2)-vf(i,k,2))/(xvel(i,k+1,2)-
     +                  xvel(i,k,2)))*(xa-xvel(i,k,2))+vf(i,k,2)
                   go to 3831
                 end if
3830          continue
            end if
3831        gradx=(vb-va)/(zb-za)
            if(abs(gradx).gt.dvvmax) then
              if(nvel(ilyr,ib).eq.1.and.nvel(i,2).eq.1) then
                write(6,3815) i,gradx
                write(13,3815) i,gradx
3815            format('layer ',i3,' has large vert vel grad (',
     +                 f10.3,'  km/s/km)')
              else
                write(6,3805) i,xa,gradx
                write(13,3805) i,xa,gradx
3805            format('layer ',i3,' has large vert vel grad at ',
     +                  f7.2,' km (',f10.3,' km/s/km)')
              end if
            end if
            if(abs(gradx).gt.dvv(i)) then
              dvv(i)=abs(gradx)
              xpv(i)=xa
            end if
3820     continue
c
         do 3920 j=1,nvel(i,2)
            xb=xvel(i,j,2)
            do 3940 k=1,nvel(ilyr,ib)
               if(abs(xb-xlvz(k)).lt..0005) go to 3920
3940        continue
            if(nzed(i).eq.1) then
              za=zm(i,1)
            else
              do 9130 k=1,nzed(i)-1
                 if(xb.ge.xm(i,k).and.xb.le.xm(i,k+1)) then
                   za=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xb-xm(i,k))+zm(i,k)
                   go to 9131
                 end if
9130          continue
            end if
9131        if(nzed(i+1).eq.1) then
              zb=zm(i+1,1)
            else
              do 9140 k=1,nzed(i+1)-1
                 if(xb.ge.xm(i+1,k).and.xb.le.xm(i+1,k+1)) then
                   zb=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xb-xm(i+1,k))+zm(i+1,k)
                   go to 9132
                 end if
9140          continue
            end if
9132        if(abs(zb-za).lt..0005) go to 3920
            vb=vf(i,j,2)
            if(nvel(ilyr,ib).eq.1) then
              va=vf(ilyr,1,ib)
            else
              do 3930 k=1,nvel(ilyr,ib)-1
                 if(xb.ge.xvel(ilyr,k,ib).and.xb.le.xvel(ilyr,k+1,ib)) 
     +           then
                   va=((vf(ilyr,k+1,ib)-vf(ilyr,k,ib))/
     +                (xvel(ilyr,k+1,ib)-xvel(ilyr,k,ib)))*
     +                (xb-xvel(ilyr,k,ib))+vf(ilyr,k,ib)
                   go to 3931
                 end if
3930          continue
            end if
3931        gradx=(vb-va)/(zb-za)
            if(abs(gradx).gt.dvvmax) then
              if(nvel(ilyr,ib).eq.1.and.nvel(i,2).eq.1) then
                write(6,3815) i,gradx
                write(13,3815) i,gradx
              else
                write(6,3805) i,xb,gradx
                write(13,3805) i,xb,gradx
              end if
            end if
            if(abs(gradx).gt.dvv(i)) then
              dvv(i)=abs(gradx)
              xpv(i)=xb
            end if
3920     continue
c
3810  continue
c
      write(13,3305)
3305  format(/'layer   vert vel grad (km/s/km)    position (km)')
      do 3330 i=1,nlayer
         write(13,3315) i,dvv(i),xpv(i)
3315     format(i4,f19.4,10x,f10.3)
3330  continue
c
c     check for large lateral velocity gradients
c
      do 5810 i=1,nlayer
         dlv(i,1)=0.
         dlv(i,2)=0.
         xlv(i,1,1)=0.
         xlv(i,1,2)=0.
         xlv(i,2,1)=0.
         xlv(i,2,2)=0.
c
         if(nvel(i,1).eq.1.and.vf(i,1,1).eq.0.) then
           if(i.eq.1) then
             write(6,185)
             write(13,185)
             stop
           end if
           do 5860 j=i-1,1,-1 
              if(nvel(j,2).gt.1.or.vf(j,1,2).ne.0.) then
                ilyr=j
                ib=2
                go to 5861
              end if
              if(nvel(j,1).gt.1.or.vf(j,1,1).ne.0.) then
                ilyr=j
                ib=1
                go to 5861
              end if
5860       continue
         else
           ilyr=i
           ib=1
         end if
c
         if(nvel(ilyr,ib).lt.2) go to 9921
c 
5861     do 5820 j=1,nvel(ilyr,ib)-1
            va1=vf(ilyr,j,ib)
            va2=vf(ilyr,j+1,ib)
            xa1=xvel(ilyr,j,ib)
            xa2=xvel(ilyr,j+1,ib)
            if(nzed(i).eq.1) then
              za1=zm(i,1)
              za2=zm(i,1)
            else
              do 5837 k=1,nzed(i)-1
                 if(xa1.ge.xm(i,k).and.xa1.le.xm(i,k+1)) then
                   za1=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xa1-xm(i,k))+zm(i,k)
                   go to 5833
                 end if
5837          continue
5833          do 5830 k=1,nzed(i)-1
                 if(xa2.ge.xm(i,k).and.xa2.le.xm(i,k+1)) then
                   za2=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xa2-xm(i,k))+zm(i,k)
                   go to 5831
                 end if
5830          continue
            end if
5831        if(nzed(i+1).eq.1) then
              zb1=zm(i+1,1)
              zb2=zm(i+1,1)
            else
              do 5847 k=1,nzed(i+1)-1
                 if(xa1.ge.xm(i+1,k).and.xa1.le.xm(i+1,k+1)) then
                   zb1=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xa1-xm(i+1,k))+zm(i+1,k)
                   go to 5834
                 end if
5847          continue
5834          do 5840 k=1,nzed(i+1)-1
                 if(xa2.ge.xm(i+1,k).and.xa2.le.xm(i+1,k+1)) then
                   zb2=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xa2-xm(i+1,k))+zm(i+1,k)
                   go to 5832
                 end if
5840          continue
            end if
5832        if(abs(zb1-za1).gt..0005.and.abs(zb2-za2).gt..0005) then
              gradl=abs((va2-va1)/((za2-za1)**2+(xa2-xa1)**2)**.5)
              if(gradl.gt.dlvmax) then
                write(6,5805) i,gradl,xa1,xa2
                write(13,5805) i,gradl,xa1,xa2
5805            format('layer ',i3,' large lat upper vel grad of ',
     +                  f7.4,' km/s/km at ',f7.2,' - ',f7.2,' km')
              end if
              if(gradl.gt.dlv(i,1)) then
                dlv(i,1)=gradl
                xlv(i,1,1)=xa1
                xlv(i,1,2)=xa2 
              end if
            end if
5820     continue
c
         if((nvel(i,2).eq.1.and.vf(i,1,2).eq.0.).or.nvel(i,2).lt.2) 
     +   go to 5810
c
9921     do 9920 j=1,nvel(i,2)-1
            vb1=vf(i,j,2)
            vb2=vf(i,j+1,2)
            xb1=xvel(i,j,2)
            xb2=xvel(i,j+1,2)
            if(nzed(i).eq.1) then
              za1=zm(i,1)
              za2=zm(i,1)
            else
              do 9937 k=1,nzed(i)-1
                 if(xb1.ge.xm(i,k).and.xb1.le.xm(i,k+1)) then
                   za1=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xb1-xm(i,k))+zm(i,k)
                   go to 9934
                 end if
9937          continue
9934          do 9930 k=1,nzed(i)-1
                 if(xb2.ge.xm(i,k).and.xb2.le.xm(i,k+1)) then
                   za2=((zm(i,k+1)-zm(i,k))/(xm(i,k+1)-
     +                  xm(i,k)))*(xb2-xm(i,k))+zm(i,k)
                   go to 9931
                 end if
9930          continue
            end if
9931        if(nzed(i+1).eq.1) then
              zb1=zm(i+1,1)
              zb2=zm(i+1,1)
            else
              do 9947 k=1,nzed(i+1)-1
                 if(xb1.ge.xm(i+1,k).and.xb1.le.xm(i+1,k+1)) then
                   zb1=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xb1-xm(i+1,k))+zm(i+1,k)
                   go to 9933
                 end if
9947          continue
9933          do 9940 k=1,nzed(i+1)-1
                 if(xb2.ge.xm(i+1,k).and.xb2.le.xm(i+1,k+1)) then
                   zb2=((zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +                  xm(i+1,k)))*(xb2-xm(i+1,k))+zm(i+1,k)
                   go to 9932
                 end if
9940          continue
            end if
9932        if(abs(zb1-za1).gt..0005.and.abs(zb2-za2).gt..0005) then
              gradl=abs((vb2-vb1)/((zb2-zb1)**2+(xb2-xb1)**2)**.5)
              if(gradl.gt.dlvmax) then
                write(6,5815) i,gradl,xb1,xb2
                write(13,5815) i,gradl,xb1,xb2
5815            format('layer ',i3,' large lat lower vel grad of ',
     +                  f7.4,' km/s/km at ',f7.2,' - ',f7.2,' km')
              end if
              if(gradl.gt.dlv(i,2)) then
                dlv(i,2)=gradl
                xlv(i,2,1)=xb1
                xlv(i,2,2)=xb2 
              end if
            end if
9920     continue
c
5810   continue
c
      write(13,2505)
2505  format(/'layer   upp/low lat vel grad (km/s/km)   between (km)')
      do 2530 i=1,nlayer
         write(13,2515) i,dlv(i,1),xlv(i,1,1),xlv(i,1,2)
         write(13,2515) i,dlv(i,2),xlv(i,2,1),xlv(i,2,2)
2515     format(i4,f19.4,12x,2f10.3)
2530  continue
c
      if(iout.eq.1) then
c
        open(58, file='depth.nodes')
        do i=1,ncont
           do j=1,nzed(i)
              write(58,1155) i,xm(i,j),zm(i,j)
1155          format(i10,3f10.3)
           enddo
        enddo
        close(58)
c
        open(58, file='velu.nodes')
        do i=1,nlayer
           do j=1,nvel(i,1)
              if(nzed(i).gt.1) then
              do k=1,nzed(i)-1
                if(xvel(i,j,1).ge.xm(i,k).and.xvel(i,j,1).
     +            le.xm(i,k+1)) then
                  zv=(zm(i,k+1)-zm(i,k))/(xm(i,k+1)-xm(i,k))*
     +               (xvel(i,j,1)-xm(i,k))+zm(i,k)
                  go to 1890
                end if
              enddo
              else
                zv=zm(i,1)
              end if
1890          write(58,1155) i,xvel(i,j,1),zv,vf(i,j,1)
           enddo
        enddo
        close(58)
c
        open(58, file='vell.nodes')
        do i=1,nlayer
           do j=1,nvel(i,2)
              if(nzed(i+1).gt.1) then
              do k=1,nzed(i+1)-1
                if(xvel(i,j,2).ge.xm(i+1,k).and.xvel(i,j,2).
     +            le.xm(i+1,k+1)) then
                  zv=(zm(i+1,k+1)-zm(i+1,k))/(xm(i+1,k+1)-
     +            xm(i+1,k))*(xvel(i,j,2)-xm(i+1,k))+zm(i+1,k)
                  go to 1891
                end if
              enddo
              else 
                zv=zm(i+1,1)
              end if
1891          write(58,1155) i,xvel(i,j,2),zv,vf(i,j,2)
           enddo
        enddo
        close(58)
c
      end if 
c
c     write out velocity model
c
      if(iflagw.eq.1) then
        open(unit=12, file='v.out')
        do 570 i=1,nlayer
           nstart=1  
  590      j1=nstart
           j2=j1+9 
           if(j2.gt.nzed(i)) j2=nzed(i)
           if(j2.lt.nzed(i)) then
             icnt=1
           else 
             icnt=0
           end if
           write(12,55) i,(xm(i,j),j=j1,j2)
           write(12,55) icnt,(zm(i,j),j=j1,j2)
           write(12,65) (ivarz(i,j),j=j1,j2)
           if(j2.eq.nzed(i)) go to 600 
           nstart=j2+1
           go to 590
600        nstart=1  
620        j1=nstart
           j2=j1+9 
           if(j2.gt.nvel(i,1)) j2=nvel(i,1)
           if(j2.lt.nvel(i,1)) then
             icnt=1
           else 
             icnt=0
           end if
           write(12,55) i,(xvel(i,j,1),j=j1,j2)
           write(12,55) icnt,(vf(i,j,1),j=j1,j2)
           write(12,65) (ivarv(i,j,1),j=j1,j2)
           if(j2.eq.nvel(i,1)) go to 630 
           nstart=j2+1
           go to 620
630        nstart=1  
650        j1=nstart
           j2=j1+9 
           if(j2.gt.nvel(i,2)) j2=nvel(i,2)
           if(j2.lt.nvel(i,2)) then
             icnt=1
           else 
             icnt=0
           end if
           write(12,55) i,(xvel(i,j,2),j=j1,j2)
           write(12,55) icnt,(vf(i,j,2),j=j1,j2)
           write(12,65) (ivarv(i,j,2),j=j1,j2)
           if(j2.eq.nvel(i,2)) go to 570 
           nstart=j2+1
           go to 650
570     continue
c
        write(12,55) ncont,(xm(ncont,j),j=1,nzed(ncont))
        write(12,55) 0,(zm(ncont,j),j=1,nzed(ncont))
      end if
c
c     calculate scale of each plot axis
c
      xscale=(xmax-xmin)/xmm
      zscale=-(zmax-zmin)/zmm
      tscale=-(tmax-tmin)/tmm
      vscale=(vmax-vmin)/vmm
      rscale=(vrmax-vrmin)/vrmm
c
      if(idump.gt.0) open(unit=14, file='m.out')
c
      if(iplot.eq.0) iplot=-1
      if(iplot.eq.2) iplot=0
      if(iplot.le.0) open(unit=19, file='p.out')
      if(iroute.ne.1) ibcol=0
      if(idcol.lt.0) idcol=ifcol
      if(ivcol.lt.0) ivcol=ifcol
c
c     plot 2-D velocity model versus depth or time
c
      if(imod.eq.1) call pltmod(ncont,iaxlab,ivel,velht,idash,izort,
     +  xtinc,xsinc,ntsmth,idump,ifrefl,symht,iroute,inode,
     +  idcol,ivcol,ivarz,ivarv,izrefl,ntsmt,npsmt,ndot,dashlf)
c
c     plot 1-D velocity profiles versus depth or time   
c
      if(ivp.eq.1) call pltvp(nlayer,xvp,iaxlab,izort,iroute,vtime)
c
c     plot RMS velocity versus distance
c
      if(ivrms.eq.1) then
        if(xrinc.lt.0.) xrinc=(xmax-xmin)/100.
        if(vrmstl.lt.0) vrmstl=1
        if(vrmsbl.lt.0) vrmsbl=nlayer
c
        call pltrms(nlayer,xrinc,vrmstl,vrmsbl,nrsmth,iaxlab,idump,
     +              iroute)
c
      end if
c
      if(iplots.eq.1) call plotnd(1)
c
      stop
c
999   write(6,95)
      write(13,95)
95    format(/'***  premature end of file  ***'/)
      stop
c
      end
