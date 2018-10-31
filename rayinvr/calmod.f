c
c     version 1.3  Aug 1992
c
c     Calmod routine for RAYINVR
c
c     ----------------------------------------------------------------
c
      subroutine calmod(ncont,pois,poisb,poisl,poisbl,invr,iflagm,
     +                  ifrbnd,xmin1d,xmax1d,insmth,xminns,xmaxns)
c
c     calculate model parameters now for use later in program
c
c
      include 'rayinvr.par'
      real xa(2*(ppcntr+ppvel)),pois(player),poisbl(papois),
     +     zsmth(pnsmth)
      integer poisb(papois),poisl(papois),insmth(pncntr)
      include 'rayinvr.com'
      iflagm=0
      idvmax=0
      idsmax=0
c
      do 10 i=1,ncont
         do 20 j=1,ppcntr
            if(abs(xm(i,j)-xmax).lt..0001) go to 30
            nzed(i)=nzed(i)+1
20       continue
30       if(nzed(i).gt.1) then
           do 40 j=1,nzed(i)-1
              if(xm(i,j).ge.xm(i,j+1)) write(0,*) 1,i,j
              if(xm(i,j).ge.xm(i,j+1)) go to 999
40         continue
           if(abs(xm(i,1)-xmin).gt..001.or.
     +        abs(xm(i,nzed(i))-xmax).gt..001) write(0,*) 2,i
           if(abs(xm(i,1)-xmin).gt..001.or.
     +        abs(xm(i,nzed(i))-xmax).gt..001) go to 999
         else
            xm(i,1)=xmax
         end if
10    continue
      do 11 i=1,nlayer
         do 21 j=1,ppvel
            if(abs(xvel(i,j,1)-xmax).lt..0001) go to 31
            nvel(i,1)=nvel(i,1)+1
21       continue
31       if(nvel(i,1).gt.1) then
           do 41 j=1,nvel(i,1)-1
              if(xvel(i,j,1).ge.xvel(i,j+1,1)) write(0,*) 3,i,j
              if(xvel(i,j,1).ge.xvel(i,j+1,1)) go to 999
41         continue
           if(abs(xvel(i,1,1)-xmin).gt..001.or.abs(xvel(i,nvel(i,1),1)-
     +        xmax).gt..001) write(0,*) 4,i
           if(abs(xvel(i,1,1)-xmin).gt..001.or.abs(xvel(i,nvel(i,1),1)-
     +        xmax).gt..001) go to 999
         else
           if(vf(i,1,1).gt.0.) then
             xvel(i,1,1)=xmax
           else
             if(i.eq.1) then
              write(0,*) 5
               go to 999
             else
               nvel(i,1)=0
             end if
           end if
         end if
11    continue
      do 12 i=1,nlayer
         do 22 j=1,ppvel
            if(abs(xvel(i,j,2)-xmax).lt..0001) go to 32
            nvel(i,2)=nvel(i,2)+1
22       continue
32       if(nvel(i,2).gt.1) then
           do 42 j=1,nvel(i,2)-1
              if(xvel(i,j,2).ge.xvel(i,j+1,2)) write(0,*) 6,i,j
              if(xvel(i,j,2).ge.xvel(i,j+1,2)) go to 999
42         continue
           if(abs(xvel(i,1,2)-xmin).gt..001.or.abs(xvel(i,nvel(i,2),2)-
     +        xmax).gt..001) write(0,*) 7,i
           if(abs(xvel(i,1,2)-xmin).gt..001.or.abs(xvel(i,nvel(i,2),2)-
     +        xmax).gt..001) go to 999
         else
           if(vf(i,1,2).gt.0.) then
             xvel(i,1,2)=xmax
           else
             nvel(i,2)=0
           end if
         end if
12    continue
c
      do 50 i=1,nlayer
         xa(1)=xmin
         xa(2)=xmax
         ib=2
         ih=ib
         do 60 j=1,nzed(i)
            do 61 k=1,ih
               if(abs(xm(i,j)-xa(k)).lt..005) go to 60
61          continue
            ib=ib+1
            xa(ib)=xm(i,j)
60       continue
         ih=ib
         do 70 j=1,nzed(i+1)
            do 80 k=1,ih
               if(abs(xm(i+1,j)-xa(k)).lt..005) go to 70
80          continue
            ib=ib+1
            xa(ib)=xm(i+1,j)
70       continue
         ih=ib
         if(nvel(i,1).gt.0) then
           il=i
           is=1
         else
           if(nvel(i-1,2).gt.0) then
             il=i-1
             is=2
           else
             il=i-1
             is=1
           end if
         end if
         do 71 j=1,nvel(il,is)
            do 81 k=1,ih
               if(abs(xvel(il,j,is)-xa(k)).lt..005) go to 71
81          continue
            ib=ib+1
            xa(ib)=xvel(il,j,is)
71       continue
         if(nvel(i,2).gt.0) then
           ih=ib
           do 72 j=1,nvel(i,2)
              do 82 k=1,ih
                 if(abs(xvel(i,j,2)-xa(k)).lt..005) go to 72
82            continue
              ib=ib+1
              xa(ib)=xvel(i,j,2)
72         continue
         end if
c
         if(ib.gt.(ptrap+1)) then
           write(6,5) i
5          format(/'***  maximum number of blocks in layer ',
     +       i2,' exceeded  ***'/)
           iflagm=1
           return
         end if
c
         call sort(xa,ib)
c
         nblk(i)=ib-1
         do 90 j=1,nblk(i)
            xbnd(i,j,1)=xa(j)
            xbnd(i,j,2)=xa(j+1)
90       continue
c
50    continue
c
      if(invr.eq.1) then
        do 310 i=1,nlayer
           ivlyr=0 
           do 320 j=1,nzed(i)
              if(ivarz(i,j).ne.1.and.ivarz(i,j).ne.-1) ivarz(i,j)=0
              if(ivarz(i,j).eq.-1) ivlyr=1
320        continue
           if(ivlyr.eq.1) then
             if(i.eq.1) write(0,*) 8
             if(i.eq.1) go to 999
             if(nzed(i).ne.nzed(i-1)) write(0,*) 9,i
             if(nzed(i).ne.nzed(i-1)) go to 999
             do 327 j=1,nzed(i)
                if(xm(i,j).ne.xm(i-1,j)) write(0,*) 10,i,j
                if(xm(i,j).ne.xm(i-1,j)) go to 999
                if(zm(i,j).lt.zm(i-1,j)) write(0,*) 11,i,j
                if(zm(i,j).lt.zm(i-1,j)) go to 999
                if(ivarz(i-1,j).eq.0.and.ivarz(i,j).eq.-1) ivarz(i,j)=0
327          continue
           end if
           do 321 j=1,nvel(i,1)
              if(ivarv(i,j,1).ne.1) ivarv(i,j,1)=0
321        continue
           ivgrad=0
           if(nvel(i,2).gt.0) then
             do 322 j=1,nvel(i,2)
                if(ivarv(i,j,2).ne.1.and.ivarv(i,j,2).ne.-1) 
     +             ivarv(i,j,2)=0
                if(ivarv(i,j,2).eq.-1) ivgrad=1
322          continue
             if(ivgrad.eq.1) then
               iflag=0
               if(nvel(i,1).gt.0) then
                 iflag=1
                 ig=i
                 jg=1
                 go to 326
               end if
               if(i.gt.1.and.nvel(i,1).eq.0) then
                 do 324 j=i-1,1,-1
                   if(nvel(j,2).gt.0) then
                     iflag=1
                     ig=j
                     jg=2
                     go to 326
                   end if
                   if(nvel(j,1).gt.0) then
                     iflag=1
                     ig=j
                     jg=1
                     go to 326
                   end if
324              continue
               end if
326            if(iflag.eq.1.and.nvel(ig,jg).eq.nvel(i,2)) then
                 do 323 j=1,nvel(ig,jg)
          if(xvel(ig,j,jg).ne.xvel(i,j,2)) write(0,*) 12,i,j,ig,jg
                    if(xvel(ig,j,jg).ne.xvel(i,j,2)) go to 999
                    if(ivarv(ig,j,jg).eq.0.and.ivarv(i,j,2).eq.-1)
     +                 ivarv(i,j,2)=0
323              continue
               else
          write(0,*) 13,i,ig,jg
                 go to 999
               end if
             end if
           end if
310     continue
c
        nvar=0
        do 410 i=1,nlayer
           do 420 j=1,nzed(i)
c
c             check for layer pinchouts
c
              iflag=0
              if(i.gt.1) then
                do 427 k=1,i-1
                   do 428 l=1,nzed(k)
                      if(abs(xm(i,j)-xm(k,l)).lt..005.and.
     +                   abs(zm(i,j)-zm(k,l)).lt..005) then
                         iflag=1
                         iv=ivarz(k,l)
                         go to 429
                      end if
428                continue
427             continue
              end if
429           if(ivarz(i,j).eq.1.and.iflag.eq.0) then
                nvar=nvar+1
                ivarz(i,j)=nvar
                partyp(nvar)=1
                parorg(nvar)=zm(i,j)
              else
                if(iflag.eq.1) then
                  ivarz(i,j)=iv
                else
                  if(ivarz(i,j).ne.-1) ivarz(i,j)=0
                end if
              end if
420        continue
c
           do 421 j=1,nvel(i,1)
              if(ivarv(i,j,1).eq.1) then
                nvar=nvar+1
                ivarv(i,j,1)=nvar
                partyp(nvar)=2
                parorg(nvar)=vf(i,j,1)
              else
                ivarv(i,j,1)=0
              end if
421        continue
           if(nvel(i,2).gt.0) then
             do 422 j=1,nvel(i,2)
                if(ivarv(i,j,2).eq.1) then
                  nvar=nvar+1
                  ivarv(i,j,2)=nvar
                  partyp(nvar)=2
                  parorg(nvar)=vf(i,j,2)
                end if
422          continue
           end if
410     continue
c
c       check for inverting floating reflectors
c
        do 426 i=1,nfrefl
           do 431 j=1,npfref(i)
              if(ivarf(i,j).eq.1) then
                nvar=nvar+1
                ivarf(i,j)=nvar
                partyp(nvar)=3
                parorg(nvar)=zfrefl(i,j)
              end if
431        continue
426     continue 
c
        if(nvar.eq.0) then
          write(6,445)
445       format(/'***  no parameters varied for inversion  ***'/)
          invr=0
        end if
c
        if(nvar.gt.pnvar) then
          write(6,455)
455       format(/'***  too many parameters varied for inversion  ***'/)
          iflagm=1
          return
        end if
c
      end if
c
c     calculate slopes and intercepts of each block boundary
c
      do 100 i=1,nlayer
         do 110 j=1,nblk(i)
            xbndc=xbnd(i,j,1)+.001
            if(nzed(i).gt.1) then
              do 120 k=1,nzed(i)-1
                 if(xbndc.ge.xm(i,k).and.xbndc.le.xm(i,k+1)) then
                   dx=xm(i,k+1)-xm(i,k)
                   c1=(xm(i,k+1)-xbnd(i,j,1))/dx
                   c2=(xbnd(i,j,1)-xm(i,k))/dx
                   z1=c1*zm(i,k)+c2*zm(i,k+1)
                   if(ivarz(i,k).ge.0) then
                     iv=ivarz(i,k)
                   else
                     iv=ivarz(i-1,k)
                   end if
                   izv(i,j,1)=iv
                   if(iv.gt.0) then
                     cz(i,j,1,1)=xm(i,k+1)
                     cz(i,j,1,2)=xm(i,k+1)-xm(i,k)
                   end if
                   c1=(xm(i,k+1)-xbnd(i,j,2))/dx
                   c2=(xbnd(i,j,2)-xm(i,k))/dx
                   z2=c1*zm(i,k)+c2*zm(i,k+1)
                   if(ivarz(i,k+1).ge.0) then
                     iv=ivarz(i,k+1)
                   else
                     iv=ivarz(i-1,k+1)
                   end if
                   izv(i,j,2)=iv
                   if(iv.gt.0) then
                     cz(i,j,2,1)=xm(i,k)
                     cz(i,j,2,2)=xm(i,k+1)-xm(i,k)
                   end if
                   go to 130
                 end if
  120         continue
            else
              z1=zm(i,1)
              if(ivarz(i,1).ge.0) then
                iv=ivarz(i,1)
              else
                iv=ivarz(i-1,1)
              end if
              izv(i,j,1)=iv
              if(iv.gt.0) then
                cz(i,j,1,1)=0.
                cz(i,j,1,2)=0.   
              end if
              z2=zm(i,1) 
              izv(i,j,2)=0
            end if
130         s(i,j,1)=(z2-z1)/(xbnd(i,j,2)-xbnd(i,j,1))
            b(i,j,1)=z1-s(i,j,1)*xbnd(i,j,1)
            if(nzed(i+1).gt.1) then
              do 140 k=1,nzed(i+1)-1
                 if(xbndc.ge.xm(i+1,k).and.xbndc.le.
     +             xm(i+1,k+1)) then
                   dx=xm(i+1,k+1)-xm(i+1,k)
                   c1=(xm(i+1,k+1)-xbnd(i,j,1))/dx
                   c2=(xbnd(i,j,1)-xm(i+1,k))/dx
                   z3=c1*zm(i+1,k)+c2*zm(i+1,k+1)
                   if(i.eq.nlayer) then
                     izv(i,j,3)=0
                   else
                     if(ivarz(i+1,k).ge.0) then
                       iv=ivarz(i+1,k)
                     else
                       iv=ivarz(i,k)
                     end if
                     izv(i,j,3)=iv
                     if(iv.gt.0) then
                       cz(i,j,3,1)=xm(i+1,k+1)
                       cz(i,j,3,2)=xm(i+1,k+1)-xm(i+1,k)   
                     end if
                   end if
                   c1=(xm(i+1,k+1)-xbnd(i,j,2))/dx
                   c2=(xbnd(i,j,2)-xm(i+1,k))/dx
                   z4=c1*zm(i+1,k)+c2*zm(i+1,k+1)
                   if(i.eq.nlayer) then
                     izv(i,j,4)=0
                   else
                     if(ivarz(i+1,k+1).ge.0) then
                       iv=ivarz(i+1,k+1)
                     else
                       iv=ivarz(i,k+1)
                     end if
                     izv(i,j,4)=iv
                     if(iv.gt.0) then
                       cz(i,j,4,1)=xm(i+1,k)
                       cz(i,j,4,2)=xm(i+1,k+1)-xm(i+1,k)   
                     end if
                   end if
                   go to 150
                 end if
140           continue
            else
              z3=zm(i+1,1)
              if(i.eq.nlayer) then
                izv(i,j,3)=0
              else
                if(ivarz(i+1,1).ge.0) then
                  iv=ivarz(i+1,1)
                else
                  iv=ivarz(i,1)
                end if
                izv(i,j,3)=iv
                if(iv.gt.0) then
                  cz(i,j,3,1)=0.
                  cz(i,j,3,2)=0.   
                end if
              end if
              z4=zm(i+1,1) 
              izv(i,j,4)=0
            end if
150         s(i,j,2)=(z4-z3)/(xbnd(i,j,2)-xbnd(i,j,1))
            b(i,j,2)=z3-s(i,j,2)*xbnd(i,j,1)
c                 
c           check for layer pinchouts 
c                 
            ivg(i,j)=1
            if(abs(z3-z1).lt..0005) ivg(i,j)=2 
            if(abs(z4-z2).lt..0005) ivg(i,j)=3 
            if(abs(z3-z1).lt..0005.and.abs(z4-z2).lt..0005) ivg(i,j)=-1
c
110      continue 
100   continue    
c
c     assign velocities to each model block
c
      do 160 i=1,nlayer 
c
         if(nvel(i,1).eq.0) then
           do 161 j=i-1,1,-1
              if(nvel(j,2).gt.0) then
                ig=j
                jg=2
                n1g=nvel(j,2)
                go to 162
              end if
              if(nvel(j,1).gt.0) then
                ig=j
                jg=1
                n1g=nvel(j,1)
                go to 162
              end if
161        continue
         else
           ig=i
           jg=1
           n1g=nvel(i,1)
         end if
c
162      if(n1g.gt.1.and.nvel(i,2).gt.1) ivcase=1
         if(n1g.gt.1.and.nvel(i,2).eq.1) ivcase=2
         if(n1g.eq.1.and.nvel(i,2).gt.1) ivcase=3
         if(n1g.eq.1.and.nvel(i,2).eq.1) ivcase=4
         if(n1g.gt.1.and.nvel(i,2).eq.0) ivcase=5
         if(n1g.eq.1.and.nvel(i,2).eq.0) ivcase=6
c
         do 170 j=1,nblk(i)
c
            if(ivg(i,j).eq.-1) go to 170
c
            xbndcl=xbnd(i,j,1)+.001
            xbndcr=xbnd(i,j,2)-.001
c
            go to (1001,1002,1003,1004,1005,1006), ivcase
c
1001        do 180 k=1,n1g-1
               if(xbndcl.ge.xvel(ig,k,jg).and.xbndcl.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(ig,k,jg)
                 vm(i,j,1)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
                 if(invr.eq.1) then
                   iv=ivarv(ig,k,jg)
                   ivv(i,j,1)=iv
                   if(iv.gt.0) then
                     cf=c1/dxx
                     call cvcalc(i,j,1,1,cf)
                     if(ivg(i,j).eq.2) call cvcalc(i,j,1,3,cf)
                   end if
                   if(c2.gt..001) then
                     iv=ivarv(ig,k+1,jg)
                     ivv(i,j,2)=iv
                     if(iv.gt.0) then
                       cf=c2/dxx
                       call cvcalc(i,j,2,1,cf)
                       if(ivg(i,j).eq.2) call cvcalc(i,j,2,3,cf)
                     end if
                   end if
                 end if
                 go to 1811
               end if
180         continue
c
1811        do 1812 k=1,n1g-1
               if(xbndcr.ge.xvel(ig,k,jg).and.xbndcr.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(ig,k,jg)
                 vm(i,j,2)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,2)=vm(i,j,2)+.001
                 if(invr.eq.1) then
                   iv=ivarv(ig,k+1,jg)
                   ivv(i,j,2)=iv
                   if(iv.gt.0) then
                     cf=c2/dxx
                     call cvcalc(i,j,2,2,cf)
                     if(ivg(i,j).eq.3) call cvcalc(i,j,2,4,cf)
                   end if
                   if(c1.gt..001) then
                     iv=ivarv(ig,k,jg)
                     ivv(i,j,1)=iv
                     if(iv.gt.0) then
                       cf=c1/dxx
                       call cvcalc(i,j,1,2,cf)
                       if(ivg(i,j).eq.3) call cvcalc(i,j,1,4,cf)
                     end if
                   end if
                 end if
                 go to 181
               end if
1812        continue
c
181         do 182 k=1,nvel(i,2)-1
               if(xbndcl.ge.xvel(i,k,2).and.xbndcl.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(i,k,2)
                 if(ivg(i,j).ne.2) then
                   vm(i,j,3)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,3)=vm(i,j,1)
                 end if 
                 if(invr.eq.1) then
                   if(ivg(i,j).eq.2) then
                     ivv(i,j,3)=0
                     icorn=0
                   else
                     iv=ivarv(i,k,2)
                     if(iv.gt.0) then
                       ivv(i,j,3)=iv
                       icorn=3
                     else
                       ivv(i,j,3)=0
                       if(iv.lt.0) then
                         icorn=1
                       else 
                         icorn=0
                       end if
                     end if
                   end if
                   if(icorn.gt.0) then
                     cf=c1/dxx
                     call cvcalc(i,j,icorn,3,cf)
                   end if
c
                   if(c2.gt..001) then
                     if(ivg(i,j).eq.2) then
                       ivv(i,j,4)=0
                       icorn=0
                     else
                       iv=ivarv(i,k+1,2)
                       if(iv.gt.0) then
                         ivv(i,j,4)=iv
                         icorn=4
                       else
                         ivv(i,j,4)=0
                         if(iv.lt.0) then
                           icorn=2
                         else
                           icorn=0
                         end if
                       end if
                     end if
                     if(icorn.gt.0) then
                       cf=c2/dxx
                       call cvcalc(i,j,icorn,3,cf)
                     end if
                   end if
                 end if
                 go to 187
               end if
182         continue
c
187         do 1822 k=1,nvel(i,2)-1
               if(xbndcr.ge.xvel(i,k,2).and.xbndcr.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(i,k,2)
                 if(ivg(i,j).ne.3) then
                   vm(i,j,4)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,4)=vm(i,j,2) 
                 end if
                 if(invr.eq.1) then
                   if(ivg(i,j).eq.3) then
                     ivv(i,j,4)=0
                     icorn=0
                   else
                     iv=ivarv(i,k+1,2)
                     if(iv.gt.0) then
                       ivv(i,j,4)=iv
                       icorn=4
                     else
                       ivv(i,j,4)=0
                       if(iv.lt.0) then
                         icorn=2
                       else
                         icorn=0
                       end if
                     end if
                   end if
                   if(icorn.gt.0) then
                     cf=c2/dxx
                     call cvcalc(i,j,icorn,4,cf)
                   end if
c
                   if(c1.gt..001) then
                     if(ivg(i,j).eq.3) then
                       ivv(i,j,3)=0
                       icorn=0
                     else
                       iv=ivarv(i,k,2)
                       if(iv.gt.0) then
                         ivv(i,j,3)=iv
                         icorn=3
                       else
                         ivv(i,j,3)=0
                         if(iv.lt.0) then
                           icorn=1
                         else
                           icorn=0
                         end if
                       end if
                     end if
                     if(icorn.gt.0) then
                       cf=c1/dxx
                       call cvcalc(i,j,icorn,4,cf)
                     end if
                   end if
                 end if
                 go to 171
               end if
1822        continue
c    
1002        do 183 k=1,n1g-1
               if(xbndcl.ge.xvel(ig,k,jg).and.xbndcl.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(ig,k,jg)
                 vm(i,j,1)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
                 if(invr.eq.1) then
                   iv=ivarv(ig,k,jg)
                   ivv(i,j,1)=iv
                   if(iv.gt.0) then
                     cf=c1/dxx
                     call cvcalc(i,j,1,1,cf)
                   end if
                   if(c2.gt..001) then
                     iv=ivarv(ig,k+1,jg)
                     ivv(i,j,2)=iv
                     if(iv.gt.0) then
                       cf=c2/dxx
                       call cvcalc(i,j,2,1,cf)
                     end if
                   end if
                 end if
                 go to 1833
               end if
183         continue
c
1833        do 1832 k=1,n1g-1
               if(xbndcr.ge.xvel(ig,k,jg).and.xbndcr.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(ig,k,jg)
                 vm(i,j,2)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,2)=vm(i,j,2)+.001
                 if(invr.eq.1) then
                   iv=ivarv(ig,k+1,jg)
                   ivv(i,j,2)=iv
                   if(iv.gt.0) then
                     cf=c2/dxx
                     call cvcalc(i,j,2,2,cf)
                   end if
                   if(c1.gt..001) then
                     iv=ivarv(ig,k,jg)
                     ivv(i,j,1)=iv
                     if(iv.gt.0) then
                       cf=c1/dxx
                       call cvcalc(i,j,1,2,cf)
                     end if
                   end if
                 end if
                 go to 184
               end if
1832        continue
c
184         vm(i,j,3)=vf(i,1,2)
            if(ivg(i,j).eq.2) vm(i,j,3)=vm(i,j,1)
            vm(i,j,4)=vf(i,1,2)
            if(ivg(i,j).eq.3) vm(i,j,4)=vm(i,j,2)
            if(invr.eq.1) then
              iv=ivarv(i,1,2)
              ivv(i,j,3)=iv
              ivv(i,j,4)=0 
              if(iv.gt.0) then
                if(ivg(i,j).ne.2) call cvcalc(i,j,3,3,1.)
                if(ivg(i,j).ne.3) call cvcalc(i,j,3,4,1.)
              end if
            end if
            go to 171
c
1003        vm(i,j,1)=vf(ig,1,jg)
            vm(i,j,2)=vf(ig,1,jg)
            if(ig.ne.i) then
              vm(i,j,1)=vm(i,j,1)+.001
              vm(i,j,2)=vm(i,j,2)+.001
            end if
            if(invr.eq.1) then
              iv=ivarv(ig,1,jg)
              ivv(i,j,1)=iv
              ivv(i,j,2)=0 
              if(iv.gt.0) then
                call cvcalc(i,j,1,1,1.)
                call cvcalc(i,j,1,2,1.)
              end if
            end if
c
            do 185 k=1,nvel(i,2)-1
               if(xbndcl.ge.xvel(i,k,2).and.xbndcl.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(i,k,2)
                 if(ivg(i,j).ne.2) then
                   vm(i,j,3)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,3)=vm(i,j,1)
                 end if
                 if(invr.eq.1) then
                   if(ivg(i,j).eq.2) then
                     ivv(i,j,3)=0
                     icorn=1
                   else
                     iv=ivarv(i,k,2)
                     if(iv.gt.0) then
                       ivv(i,j,3)=iv
                       icorn=3
                     else
                       ivv(i,j,3)=0
                       icorn=0
                     end if
                   end if
                   if(icorn.gt.0) then
                     cf=c1/dxx
                     call cvcalc(i,j,icorn,3,cf)
                   end if
c
                   if(c2.gt..001) then
                     if(ivg(i,j).eq.2) then
                       ivv(i,j,3)=0
                       icorn=2
                     else
                       iv=ivarv(i,k+1,2)
                       if(iv.gt.0) then
                         ivv(i,j,4)=iv
                         icorn=4
                       else
                         ivv(i,j,4)=0
                         icorn=0
                       end if
                     end if
                     if(icorn.gt.0) then
                       cf=c2/dxx
                       call cvcalc(i,j,icorn,3,cf)
                     end if
                   end if
                 end if
                 go to 188
               end if
185         continue
c
188         do 1851 k=1,nvel(i,2)-1
               if(xbndcr.ge.xvel(i,k,2).and.xbndcr.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(i,k,2)
                 if(ivg(i,j).ne.3) then
                   vm(i,j,4)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,4)=vm(i,j,2)
                 end if
                 if(invr.eq.1) then
                   if(ivg(i,j).eq.3) then
                     ivv(i,j,4)=0
                     icorn=2   
                   else
                     iv=ivarv(i,k+1,2)
                     if(iv.gt.0) then
                       ivv(i,j,4)=iv
                       icorn=4
                     else
                       ivv(i,j,4)=0
                       icorn=0
                     end if
                   end if
                   if(icorn.gt.0) then
                     cf=c2/dxx
                     call cvcalc(i,j,icorn,4,cf)
                   end if
c
                   if(c1.gt..001) then
                     if(ivg(i,j).eq.3) then
                       ivv(i,j,4)=0
                       icorn=1
                     else
                       iv=ivarv(i,k,2)
                       if(iv.gt.0) then
                         ivv(i,j,3)=iv
                         icorn=3
                       else
                         ivv(i,j,3)=0
                         icorn=0
                       end if
                     end if
                     if(icorn.gt.0) then
                       cf=c1/dxx
                       call cvcalc(i,j,icorn,4,cf)
                     end if
                   end if
                 end if
                 go to 171
               end if
1851        continue
c
1004        vm(i,j,1)=vf(ig,1,jg)
            vm(i,j,2)=vf(ig,1,jg) 
            if(ig.ne.i) then
              vm(i,j,1)=vm(i,j,1)+.001
              vm(i,j,2)=vm(i,j,2)+.001
            end if
            if(invr.eq.1) then
              iv=ivarv(ig,1,jg)
              ivv(i,j,1)=iv
              ivv(i,j,2)=0 
              if(iv.gt.0) then
                call cvcalc(i,j,1,1,1.)
                call cvcalc(i,j,1,2,1.)
              end if
            end if
c
            vm(i,j,3)=vf(i,1,2)
            if(ivg(i,j).eq.2) vm(i,j,3)=vm(i,j,1)
            vm(i,j,4)=vf(i,1,2)
            if(ivg(i,j).eq.3) vm(i,j,4)=vm(i,j,2)
            if(invr.eq.1) then
              iv=ivarv(i,1,2)
              if(iv.gt.0) then
                ivv(i,j,3)=iv
                icorn=3
              else
                ivv(i,j,3)=0
                if(iv.lt.0) then
                  icorn=1
                else
                  icorn=0
                end if
              end if
              ivv(i,j,4)=0
              if(icorn.gt.0) then
                if(ivg(i,j).ne.2) call cvcalc(i,j,icorn,3,1.)
                if(ivg(i,j).ne.3) call cvcalc(i,j,icorn,4,1.)
              end if
            end if
            go to 171
c
1005        do 186 k=1,n1g-1
               if(xbndcl.ge.xvel(ig,k,jg).and.xbndcl.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(ig,k,jg)
                 vm(i,j,1)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
                 vm(i,j,3)=vm(i,j,1)
                 if(invr.eq.1) then
                   iv=ivarv(ig,k,jg)
                   ivv(i,j,1)=iv
                   if(iv.gt.0) then
                     cf=c1/dxx
                     call cvcalc(i,j,1,1,cf)
                     if(ivg(i,j).ne.2) call cvcalc(i,j,1,3,cf)
                   end if
                   if(c2.gt..001) then
                     iv=ivarv(ig,k+1,jg)
                     ivv(i,j,2)=iv
                     if(iv.gt.0) then
                       cf=c2/dxx
                       call cvcalc(i,j,2,1,cf)
                       if(ivg(i,j).ne.2) call cvcalc(i,j,2,3,cf)
                     end if
                   end if
                 end if
                 go to 1861
               end if
186         continue
c
1861        do 1862 k=1,n1g-1
               if(xbndcr.ge.xvel(ig,k,jg).and.xbndcr.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(ig,k,jg)
                 vm(i,j,2)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,2)=vm(i,j,2)+.001
                 vm(i,j,4)=vm(i,j,2)                        
                 if(invr.eq.1) then
                   iv=ivarv(ig,k+1,jg)
                   ivv(i,j,2)=iv
                   if(iv.gt.0) then
                     cf=c2/dxx
                     call cvcalc(i,j,2,2,cf)
                     if(ivg(i,j).ne.3) call cvcalc(i,j,2,4,cf)
                   end if
                   if(c1.gt..001) then
                     iv=ivarv(ig,k,jg)
                     ivv(i,j,1)=iv
                     if(iv.gt.0) then
                       cf=c1/dxx
                       call cvcalc(i,j,1,2,cf)
                       if(ivg(i,j).ne.3) call cvcalc(i,j,1,4,cf)
                     end if
                   end if
                   ivv(i,j,3)=0
                   ivv(i,j,4)=0
                 end if
                 go to 171
               end if
1862        continue
c
1006        vm(i,j,1)=vf(ig,1,jg)
            if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
            vm(i,j,2)=vm(i,j,1) 
            vm(i,j,3)=vm(i,j,1)
            vm(i,j,4)=vm(i,j,1) 
            if(invr.eq.1) then
              iv=ivarv(ig,1,jg)
              ivv(i,j,1)=iv
              ivv(i,j,2)=0 
              ivv(i,j,3)=0 
              ivv(i,j,4)=0 
              if(iv.gt.0) then
                call cvcalc(i,j,1,1,1.)
                call cvcalc(i,j,1,2,1.)
                call cvcalc(i,j,1,3,1.)
                call cvcalc(i,j,1,4,1.)
              end if
            end if
c
c
c           calculate velocity coefficients
c
171         s1=s(i,j,1)
            s2=s(i,j,2)
            b1=b(i,j,1)
            b2=b(i,j,2)
            xb1=xbnd(i,j,1)
            xb2=xbnd(i,j,2)
            if(ivg(i,j).eq.2) then
              z3=s(i,j,2)*xb1+b(i,j,2)+.001
              z4=s(i,j,2)*xb2+b(i,j,2)
              s2=(z4-z3)/(xb2-xb1)
              b2=z3-s2*xb1
            end if
            if(ivg(i,j).eq.3) then
              z3=s(i,j,2)*xb1+b(i,j,2)
              z4=s(i,j,2)*xb2+b(i,j,2)+.001
              s2=(z4-z3)/(xb2-xb1)
              b2=z3-s2*xb1
            end if
            v1=vm(i,j,1)
            v2=vm(i,j,2)
            v3=vm(i,j,3)
            v4=vm(i,j,4)
c
            c(i,j,1)=s2*(xb2*v1-xb1*v2)+b2*(v2-v1)-
     +               s1*(xb2*v3-xb1*v4)-b1*(v4-v3)       
            c(i,j,2)=s2*(v2-v1)-s1*(v4-v3)
            c(i,j,3)=-xb2*v1+xb1*v2+xb2*v3-xb1*v4
            c(i,j,4)=-v2+v1+v4-v3
            c(i,j,5)=b2*(xb2*v1-xb1*v2)-b1*(xb2*v3-xb1*v4)
            c(i,j,6)=(s2-s1)*(xb2-xb1)
            c(i,j,7)=(b2-b1)*(xb2-xb1)
            c(i,j,8)=2.*c(i,j,2)*c(i,j,7)
            c(i,j,9)=c(i,j,2)*c(i,j,6)
            c(i,j,10)=c(i,j,4)*c(i,j,7)-c(i,j,3)*c(i,j,6)
            c(i,j,11)=c(i,j,1)*c(i,j,7)-c(i,j,5)*c(i,j,6)
c
            if(ivg(i,j).eq.-1) then
              vm(i,j,1)=0.
              vm(i,j,2)=0.
              vm(i,j,3)=0.
              vm(i,j,4)=0.
              do 172 k=1,11
                 c(i,j,1)=0.
172           continue
            end if
            if(abs(vm(i,j,1)-vm(i,j,2)).le..001.and.abs(vm(i,j,2)-
     +        vm(i,j,3)).le..001.and.abs(vm(i,j,3)-vm(i,j,4)).le..001. 
     +        and.ivg(i,j).ne.-1) ivg(i,j)=0
c
170      continue
160   continue    
c                 
c     assign values to array vsvp
c                 
      if(pois(1).lt.-10.) then 
        do 190 i=1,nlayer 
           do 200 j=1,nblk(i)
              vsvp(i,j)=0.57735
200        continue
190     continue  
      else        
        if(nlayer.gt.1) then
          if(pois(2).lt.-10.) then
            do 210 j=1,nblk(1)
               vsvp(1,j)=sqrt((1.-2.*pois(1))/(2.*(1.-pois(1))))
210         continue
            do 220 i=2,nlayer
               do 230 j=1,nblk(i)
                  vsvp(i,j)=vsvp(1,1)
230            continue
220         continue
          else    
            do 240 i=1,nlayer
               if(pois(i).lt.-10.) then
                 do 250 j=1,nblk(i)
                    vsvp(i,j)=0.57735
250              continue
               else
                 do 260 j=1,nblk(i)
                    vsvp(i,j)=sqrt((1.-2.*pois(i))/(2.*(1.-pois(i))))
260              continue 
               end if
240         continue 
          end if  
        end if    
      end if      
c                 
c     calculate velocity ratios for specific model blocks specified
c     through the arrays poisbl, poisl and poisb
c                 
      i=1         
270   if(poisbl(i).lt.-10.) go to 400
      vsvp(poisl(i),poisb(i))=
     +  sqrt((1.-2.*poisbl(i))/(2.*(1.-poisbl(i))))
      i=i+1       
      if(i.le.papois) go to 270
c                 
c     calculation of smooth layer boundaries
c
400   if(ibsmth.gt.0) then
        xsinc=(xmax-xmin)/float(npbnd-1)
        do 600 i=1,nlayer+1
           if(i.lt.(nlayer+1)) then
             il=i 
             ib=1 
           else   
             il=i-1
             ib=2 
           end if 
           iblk=1 
           do 610 j=1,npbnd
              x=xmin+float(j-1)*xsinc
              if(x.lt.xmin) x=xmin+.001
              if(x.gt.xmax) x=xmax-.001
620           if(x.ge.xbnd(il,iblk,1).and.x.le.xbnd(il,iblk,2)) then
                cosmth(i,j)=s(il,iblk,ib)*x+b(il,iblk,ib)
                go to 610
              else
                iblk=iblk+1
                go to 620
              end if
610        continue
600     continue  
        n1ns=nint(xminns/xsinc)+1
        n2ns=nint(xmaxns/xsinc)+1
        iflag12=0
        if(n1ns.ge.1.and.n1ns.le.npbnd.and.n2ns.ge.1.and.
     +  n2ns.le.npbnd.and.n1ns.lt.n2ns) iflag12=1
        if(nbsmth.gt.0) then
          do 630 i=1,nlayer+1
             if(xminns.lt.xmin.and.xmaxns.lt.xmin) then
               do 6630 j=1,pncntr
                  if(insmth(j).eq.i) go to 630
6630           continue
             end if
             iflagns=0
             do 6640 j=1,pncntr
                if(insmth(j).eq.i) iflagns=1
6640         continue
             do 640 j=1,npbnd
                zsmth(j)=cosmth(i,j)
640          continue 
             do 650 j=1,nbsmth 
                if(iflag12.eq.1.and.iflagns.eq.1) then
                  call smooth2(zsmth,npbnd,n1ns,n2ns) 
                else
                  call smooth(zsmth,npbnd) 
                end if
650          continue
             do 660 j=1,npbnd
                cosmth(i,j)=zsmth(j)
660          continue 
             if(idump.eq.2) then
               do 670 j=1,npbnd
                  x=xmin+float(j-1)*xsinc
                  write(12,635) x,zsmth(j)
635               format(2f7.2)
670            continue
             end if
630       continue
        end if    
      end if      
c                 
      if(idump.eq.1) then
        write(12,15) nlayer
15      format('***  velocity model:  ***'//'number of layers=',i2)
        do 510 i=1,nlayer
           write(12,25) i,nblk(i)
25         format(/'layer#',i2,'  nblk=',i4,
     +     ' (ivg,x1,x2,z11,z12,z21,z22,s1,b1,s2,b2,vp1,vs1,vp2,vs2,
     +     vp3,vs3,vp4,vs4,c1,c2,...,c11)')
           write(12,35) (ivg(i,j),j=1,nblk(i))
           write(12,45) (xbnd(i,j,1),j=1,nblk(i))
           write(12,45) (xbnd(i,j,2),j=1,nblk(i))
           write(12,45) (s(i,j,1)*xbnd(i,j,1)+b(i,j,1),j=1,nblk(i))
           write(12,45) (s(i,j,1)*xbnd(i,j,2)+b(i,j,1),j=1,nblk(i))
           write(12,45) (s(i,j,2)*xbnd(i,j,1)+b(i,j,2),j=1,nblk(i))
           write(12,45) (s(i,j,2)*xbnd(i,j,2)+b(i,j,2),j=1,nblk(i))
           write(12,45) (s(i,j,1),j=1,nblk(i))
           write(12,45) (b(i,j,1),j=1,nblk(i))
           write(12,45) (s(i,j,2),j=1,nblk(i))
           write(12,45) (b(i,j,2),j=1,nblk(i))
           write(12,45) (vm(i,j,1),j=1,nblk(i))
           write(12,45) (vm(i,j,1)*vsvp(i,j),j=1,nblk(i)) 
           write(12,45) (vm(i,j,2),j=1,nblk(i))
           write(12,45) (vm(i,j,2)*vsvp(i,j),j=1,nblk(i))
           write(12,45) (vm(i,j,3),j=1,nblk(i))
           write(12,45) (vm(i,j,3)*vsvp(i,j),j=1,nblk(i)) 
           write(12,45) (vm(i,j,4),j=1,nblk(i))
           write(12,45) (vm(i,j,4)*vsvp(i,j),j=1,nblk(i))
           write(12,55) (c(i,j,1),j=1,nblk(i))
           write(12,55) (c(i,j,2),j=1,nblk(i))
           write(12,55) (c(i,j,3),j=1,nblk(i))
           write(12,55) (c(i,j,4),j=1,nblk(i))
           write(12,55) (c(i,j,5),j=1,nblk(i))
           write(12,55) (c(i,j,6),j=1,nblk(i))
           write(12,55) (c(i,j,7),j=1,nblk(i))
           write(12,55) (c(i,j,8),j=1,nblk(i))
           write(12,55) (c(i,j,9),j=1,nblk(i))
           write(12,55) (c(i,j,10),j=1,nblk(i))
           write(12,55) (c(i,j,11),j=1,nblk(i))
35         format(100i10)
45         format(100f10.4)
55         format(100e10.3)
510     continue  
c
        xmod=xmax-xmin 
        write(12,65) 
65      format(/'equivalent 1-dimensional velocity model:'/)
        do 520 i=1,nlayer
           z1sum=0.
           z2sum=0.
           vp1=0. 
           vp2=0. 
           vs1=0. 
           vs2=0. 
           vp1sum=0.
           vp2sum=0.
           vs1sum=0.
           vs2sum=0.
           xvmod=0.
           do 530 j=1,nblk(i)
              xblk=xbnd(i,j,2)-xbnd(i,j,1)
              z11=s(i,j,1)*xbnd(i,j,1)+b(i,j,1)
              z12=s(i,j,1)*xbnd(i,j,2)+b(i,j,1)
              z21=s(i,j,2)*xbnd(i,j,1)+b(i,j,2)
              z22=s(i,j,2)*xbnd(i,j,2)+b(i,j,2)
              z1sum=z1sum+xblk*(z11+z12)/2.
              z2sum=z2sum+xblk*(z21+z22)/2.
              if(vm(i,j,1).gt..001) then
                vp1sum=vp1sum+xblk*(vm(i,j,1)+vm(i,j,2))/2.
                vs1sum=vs1sum+xblk*(vm(i,j,1)+vm(i,j,2))*vsvp(i,j)/2.
                vp2sum=vp2sum+xblk*(vm(i,j,3)+vm(i,j,4))/2.
                vs2sum=vs2sum+xblk*(vm(i,j,3)+vm(i,j,4))*vsvp(i,j)/2.
                xvmod=xvmod+xblk
              end if
530        continue
           z1=z1sum/xmod
           z2=z2sum/xmod
           if(xvmod.gt..000001) then
             vp1=vp1sum/xvmod
             vp2=vp2sum/xvmod
             vs1=vs1sum/xvmod
             vs2=vs2sum/xvmod 
           end if 
           write(12,75) i,z1,z2,vp1,vp2,vs1,vs2
75         format('layer# ',i2,'   z1=',f7.2,'   z2=',f7.2,
     +            ' km'/9x,'  vp1=',f7.2,'  vp2=',f7.2,' km/s'/
     +                  9x,'  vs1=',f7.2,'  vs2=',f7.2,' km/s')
520     continue  
c
        if(xmin1d.lt.-999998.) xmin1d=xmin
        if(xmax1d.lt.-999998.) xmax1d=xmax
        write(12,175) xmin1d,xmax1d
175     format(/'1-dimensional P-wave velocity model between ',
     +          f7.2,' and ',f7.2,' km:'/)
        xmod=xmax1d-xmin1d
        do 720 i=1,nlayer
           z1sum=0.
           z2sum=0.
           vp1=0. 
           vp2=0. 
           vs1=0. 
           vs2=0. 
           vp1sum=0.
           vp2sum=0.
           vs1sum=0.
           vs2sum=0.
           xvmod=0.
           do 730 j=1,nblk(i)
              if(xbnd(i,j,1).ge.xmax1d) go to 730
              if(xbnd(i,j,2).le.xmin1d) go to 730
              if(xbnd(i,j,1).lt.xmin1d) then
                xb1=xmin1d
              else
                xb1=xbnd(i,j,1)
              end if
              if(xbnd(i,j,2).gt.xmax1d) then
                xb2=xmax1d
              else
                xb2=xbnd(i,j,2)
              end if
              xblk=xb2-xb1
              z11=s(i,j,1)*xb1+b(i,j,1)
              z12=s(i,j,1)*xb2+b(i,j,1)
              z21=s(i,j,2)*xb1+b(i,j,2)
              z22=s(i,j,2)*xb2+b(i,j,2)
              z1sum=z1sum+xblk*(z11+z12)/2.
              z2sum=z2sum+xblk*(z21+z22)/2.
              if(vm(i,j,1).gt..001) then
                layer=i 
                iblk=j
                v11=vel(xb1,z11)
                v12=vel(xb2,z12)
                v21=vel(xb1,z21)
                v22=vel(xb2,z22)
                vp1sum=vp1sum+xblk*(v11+v12)/2.
                vs1sum=vs1sum+xblk*(v11+v12)*vsvp(i,j)/2.
                vp2sum=vp2sum+xblk*(v21+v22)/2.
                vs2sum=vs2sum+xblk*(v21+v22)*vsvp(i,j)/2.
                xvmod=xvmod+xblk
              end if
730        continue
           z1=z1sum/xmod
           z2=z2sum/xmod
           if(xvmod.gt..000001) then
             vp1=vp1sum/xvmod
             vp2=vp2sum/xvmod
             vs1=vs1sum/xvmod
             vs2=vs2sum/xvmod 
           end if 
           write(12,155) i,xmax,0,z1,0
155        format(i2,1x,f7.2/i2,1x,f7.2/3x,i7)
           write(12,155) i,xmax,0,vp1,0
           write(12,155) i,xmax,0,vp2,0
           if(i.eq.nlayer) write(12,165) i+1,xmax,0,z2
165        format(i2,1x,f7.2/i2,1x,f7.2)
720     continue  
      end if      
c
      if(idump.eq.1) write(12,85) 
85    format(/'layer   max. gradient (km/s/km)   block')
c
      do 910 i=1,nlayer
         delv=0.
         ibd=0
         do 920 j=1,nblk(i)
            if(ivg(i,j).lt.1) go to 920
            if(ivg(i,j).eq.2) then
              delv1=0.
              delv3=0.
            else
              x1=xbnd(i,j,1)
              z1=s(i,j,1)*x1+b(i,j,1)
              denom=c(i,j,6)*x1+c(i,j,7)
              vx=(c(i,j,8)*x1+c(i,j,9)*x1**2+c(i,j,10)*z1+
     +            c(i,j,11))/denom**2
              vz=(c(i,j,3)+c(i,j,4)*x1)/denom
              delv1=(vx**2+vz**2)**.5
              z3=s(i,j,2)*x1+b(i,j,2)
              vx=(c(i,j,8)*x1+c(i,j,9)*x1**2+c(i,j,10)*z3+
     +            c(i,j,11))/denom**2
              vz=(c(i,j,3)+c(i,j,4)*x1)/denom
              delv3=(vx**2+vz**2)**.5
            end if
            if(ivg(i,j).eq.3) then
              delv2=0.
              delv4=0.
            else
              x2=xbnd(i,j,2)
              z2=s(i,j,1)*x2+b(i,j,1)
              denom=c(i,j,6)*x2+c(i,j,7)
              vx=(c(i,j,8)*x2+c(i,j,9)*x2**2+c(i,j,10)*z2+
     +            c(i,j,11))/denom**2
              vz=(c(i,j,3)+c(i,j,4)*x2)/denom
              delv2=(vx**2+vz**2)**.5
              z4=s(i,j,2)*x2+b(i,j,2)
              vx=(c(i,j,8)*x2+c(i,j,9)*x2**2+c(i,j,10)*z4+
     +            c(i,j,11))/denom**2
              vz=(c(i,j,3)+c(i,j,4)*x2)/denom
              delv4=(vx**2+vz**2)**.5
            end if
            delm=amax1(delv1,delv2,delv3,delv4)
            if(delm.gt.delv) then
              delv=delm
              ibd=j
            end if
920      continue
         if(idump.eq.1) write(12,95) i,delv,ibd
95       format(i4,f17.4,i17)
         if(delv.gt.dvmax) then
           idvmax=idvmax+1 
           ldvmax(idvmax)=i
         end if
910   continue
c
      if(idump.eq.1) write(12,105) 
105   format(/'boundary   slope change (degrees)   between points')
c
      do 930 i=1,ncont
         dslope=0.
         ips1=0
         ips2=0
         if(nzed(i).gt.2) then
           do 940 j=1,nzed(i)-2
              slope1=(zm(i,j+1)-zm(i,j))/(xm(i,j+1)-xm(i,j)) 
              slope2=(zm(i,j+2)-zm(i,j+1))/(xm(i,j+2)-xm(i,j+1)) 
              ds1=atan(slope1)*pi18
              ds2=atan(slope2)*pi18
              if(abs(ds2-ds1).gt.dslope) then
                dslope=abs(ds2-ds1)
                ips1=j
                ips2=j+2
              end if
940        continue
         end if
         if(idump.eq.1) write(12,115) i,dslope,ips1,ips2
115      format(i4,f19.4,12x,2i7)
         if(dslope.gt.dsmax) then
           idsmax=idsmax+1
           ldsmax(idsmax)=i
         end if
930   continue
c
      return      
c                 
999   write(6,900)
900   format(/'***  error in velocity model 2 ***'/)
      iflagm=1
      return
      end         
