c
c     version 1.4  Dec 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |          ********   R E C I P R O C I T Y  ******            |   
c     |                                                              |
c     |        Check travel time reciprocity of "tx.in" file         |   
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c     input:  file tx.in used by Rayinvr
c     output: file reciprocity.out
c     output: file r.out
c
      include 'rayinvr.par'
      parameter(pmax=10000)
      real x(pmax),t(pmax),u(pmax),xh(pmax),th(pmax),uh(pmax),
     +     xshot(pshot),x1(pmax),x2(pmax),dpair(pmax),
     +     x1rec(pmax),x2rec(pmax),drec(pmax),trec(pmax),
     +     xrec(pmax),urec(pmax),xirec(pmax),tvect(pshot),
     +     nmat(pshot,pshot),nmati(pshot,pshot),delta(pshot)
      integer ip(pmax),iph(pmax),ipair(pmax),irec(pmax)
      character flag*1
c
      open(unit=10, file='tx.in')
      open(unit=11, file='reciprocity.out')
      open(unit=12, file='r.out')
c
      xmax=5.
      vmax=.25
      tol=99999.
c
      nshot=0
      ipmax=0
c
101   read(10,5,end=999) xx,tt,uu,ii
5     format(3f10.3,i10)
      if(ii.eq.0) then
        if(nshot.gt.0) then
          do 30 i=1,nshot
             if(abs(xshot(i)-xx).lt..001) go to 101
30        continue
        end if
        nshot=nshot+1
        xshot(nshot)=xx
      end if
      if(ii.gt.0) then
        if(ii.gt.ipmax) ipmax=ii
      end if
      go to 101
c
999   rewind(10)
c
      do 460 i=1,nshot
         do 470 j=1,nshot
470         nmat(i,j)=0.
460      tvect(i)=0.
      func=0.
      npairs=0
c
      write(11,115) (i,i=1,nshot)
      write(6,115) (i,i=1,nshot)
115   format(/'shot offsets (km):'/
     +  /'nshot',<nshot>i9/'-----',<nshot>('---------'))
      do 110 i=1,nshot
         write(11,105) i,(xshot(j)-xshot(i),j=1,nshot)
         write(*,105) i,(xshot(j)-xshot(i),j=1,nshot)
105      format(i5,<nshot>f9.2)
110   continue
c
      i=1
100   read(10,5,end=998) x(i),t(i),u(i),ip(i)
      i=i+1
      go to 100
c
998   numf=i-1
c
      iline=1
      npair=0
210   iphase=-2
      nphase=0
200   if(ip(iline).le.0.or.ip(iline).ne.iphase) then
        if(nphase.eq.0) go to 22
        do 10 i=1,nshot
           if(abs(xshot(i)-xshotc).lt..001) go to 10
           dist=99999.
           if(nphase.eq.1) then
             if((xh(1).gt.xshotc.and.xh(1).le.xshot(i)).or.
     +          (xh(1).gt.xshot(i).and.xh(1).le.xshotc)) then      
                distc=abs(xshot(i)-xh(1))
                if(distc.lt.tol) dist=distc
             end if
           else
             do 20 j=1,nphase-1
                if((xh(j).le.xshot(i).and.xh(j+1).ge.xshot(i)).
     +          or.(xh(j).ge.xshot(i).and.xh(j+1).le.xshot(i)))
     +          then
                  dist=0.
                  go to 21
                else
                  if((xh(j).le.xshot(i).and.xh(j).ge.xshotc).
     +            or.(xh(j).ge.xshot(i).and.xh(j).le.xshotc))
     +            then
                    distc=abs(xshot(i)-xh(j))
                    if(distc.lt.dist.and.distc.lt.tol) dist=distc
                  end if
                end if
20           continue
             j=nphase
             if((xh(j).le.xshot(i).and.xh(j).ge.xshotc).or.
     +          (xh(j).ge.xshot(i).and.xh(j).le.xshotc)) then
               distc=abs(xshot(i)-xh(j))
               if(distc.lt.dist.and.distc.lt.tol) dist=distc
             end if
           end if
21         if(dist.lt.99998.) then
             npair=npair+1
             x1(npair)=xshotc
             x2(npair)=xshot(i)
             ipair(npair)=iphase
             dpair(npair)=dist
           end if
10      continue
c
22      if(ip(iline).eq.-1) go to 9998
        if(ip(iline).eq.0) then
          xshotc=x(iline)
          iline=iline+1
          go to 210
        else
          iphase=ip(iline)
          nphase=1
          xh(1)=x(iline)
          th(1)=t(iline)
          uh(1)=u(iline)
          iph(1)=ip(iline)
          iline=iline+1
          go to 200
        end if
      else
        nphase=nphase+1
        xh(nphase)=x(iline)
        th(nphase)=t(iline)
        uh(nphase)=u(iline)
        iph(nphase)=ip(iline)
        iline=iline+1
        go to 200
      end if
c
9998  nrec=0
      iline=1
310   iphase=-2
      nphase=0
300   if(ip(iline).le.0.or.ip(iline).ne.iphase) then
        if(nphase.eq.0) go to 333
        do 50 i=1,nshot
           if(abs(xshot(i)-xshotc).lt..001) go to 50
           d1=-1.
           d2=-1.
           do 60 j=1,npair
              if(abs(x1(j)-xshotc).lt..001.and.
     +           abs(x2(j)-xshot(i)).lt..001.and.
     +           ipair(j).eq.iphase) d1=dpair(j)
              if(abs(x1(j)-xshot(i)).lt..001.and.
     +           abs(x2(j)-xshotc).lt..001.and.
     +           ipair(j).eq.iphase) d2=dpair(j)
60         continue
           if(d1.lt.0..or.d2.lt.0.) go to 50
           if(d1.gt.d2) then
             dist=d1
           else
             dist=d2
           end if
           if(xshot(i).gt.xshotc) then
             xpos=xshot(i)-dist
           else
             xpos=xshot(i)+dist
           end if
           tpick=-99999.
           if(nphase.eq.1) then
             if(abs(xh(1)-xpos).lt..001) then
               tpick=th(1)
               upick=uh(1)
               xipick=0.
             end if
           else
             do 320 j=1,nphase-1
                if((xpos.ge.xh(j).and.xpos.le.xh(j+1)).or.
     +          (xpos.le.xh(j).and.xpos.ge.xh(j+1))) then
                  denom=xh(j+1)-xh(j)
c                 if(denom.ne.0.) then
                    tpick=(th(j+1)-th(j))/denom*(xpos-
     +                     xh(j))+th(j)
                    upick=(uh(j+1)-uh(j))/denom*(xpos-
     +                     xh(j))+uh(j)
c                 else
c                   tpick=(th(j+1)-th(j))/2.
c                   upick=(uh(j+1)-uh(j))/2.
c                 end if
                  xipick=min(abs(xpos-xh(j)),abs(xpos-xh(j+1)))
                  go to 321
                end if
320          continue
           end if
321        if(tpick.gt.-99998) then
             nrec=nrec+1
             x1rec(nrec)=xshotc
             x2rec(nrec)=xshot(i)
             irec(nrec)=iphase
             xrec(nrec)=xpos
             drec(nrec)=dist
             trec(nrec)=tpick
             urec(nrec)=upick
             xirec(nrec)=xipick
           end if
50      continue
c
333     if(ip(iline).eq.-1) go to 9999
        if(ip(iline).eq.0) then
          xshotc=x(iline)
          iline=iline+1
          go to 310
        else
          iphase=ip(iline)
          nphase=1
          xh(1)=x(iline)
          th(1)=t(iline)
          uh(1)=u(iline)
          iph(1)=ip(iline)
          iline=iline+1
          go to 300
        end if
      else
        nphase=nphase+1
        xh(nphase)=x(iline)
        th(nphase)=t(iline)
        uh(nphase)=u(iline)
        iph(nphase)=ip(iline)
        iline=iline+1
        go to 300
      end if
c
9999  continue
      write(11,55)
      nused=0
      do 91 k=1,ipmax
         do 71 i=1,nrec
            if(irec(i).ne.k) go to 71
            do 81 j=1,nrec
               if(abs(x1rec(i)-x2rec(j)).lt..001.and.
     +         abs(x1rec(j)-x2rec(i)).lt..001.and.
     +         irec(i).eq.irec(j)) then
                 if(drec(i).lt..001) then
                   vi=abs(x1rec(i)-x2rec(i))/trec(i)
                   vj=abs(x1rec(i)-x2rec(i))/trec(j)
                   vdiff=abs(vj-vi)
                 else
                   vave=(abs(x1rec(i)-x2rec(i))-drec(i))*
     +                  (1./trec(i)+1./trec(j))/2.
                   if(trec(i).gt.trec(j)) then
                     t1=trec(i)
                   else
                     t1=trec(j)
                   end if
                   vdiff=vave-drec(i)/(t1-(abs(x1rec(i)-x2rec(i))-
     +                   2.*drec(i))/vave)
                 end if
                 if((abs(trec(i)-trec(j)).gt.max(urec(i),urec(j))).and.
     +             (vdiff.gt.vmax.or.drec(i).lt.001)) then
                   flag='*'
                 else
                   flag=' '
                 end if
                 write(11,4) 
     +             x1rec(i),x2rec(i),xirec(i),xirec(j),irec(i),drec(i),
     +             max(urec(i),urec(j)),trec(j)-trec(i),flag
                 nused=nused+1
                 if(nused.ge.nrec) go to 99999
               end if
81          continue
71       continue
91    continue
c
99999 write(11,55)
      write(6,55)
55    format(/'   xshot1   xshot2  xinter1  xinter2 phase',
     +        ' |x2-x1|   pick unc  t2-t1'/
     +        '----------------------------------',
     +        '-----------------------------------')
      write(12,555) xmax
555   format(/'recirpocal pairs with |x2-x1|<',f4.1,' km:'//
     +'    ishot1    ishot2    xshot1    xhsot2     t2-t1')
      nused=0
      nflag=0
      do 90 k=1,ipmax
         do 70 i=1,nrec
            if(irec(i).ne.k) go to 70
            do 80 j=1,nrec
               if(abs(x1rec(i)-x2rec(j)).lt..001.and.
     +         abs(x1rec(j)-x2rec(i)).lt..001.and.
     +         irec(i).eq.irec(j)) then
                 if(drec(i).lt..001) then
                   vi=abs(x1rec(i)-x2rec(i))/trec(i)
                   vj=abs(x1rec(i)-x2rec(i))/trec(j)
                   vdiff=abs(vj-vi)
                 else
                   vave=(abs(x1rec(i)-x2rec(i))-drec(i))*
     +                  (1./trec(i)+1./trec(j))/2.
                   if(trec(i).gt.trec(j)) then
                     t1=trec(i)
                   else
                     t1=trec(j)
                   end if
                   vdiff=vave-drec(i)/(t1-(abs(x1rec(i)-x2rec(i))-
     +                   2.*drec(i))/vave)
                 end if
                 if((abs(trec(i)-trec(j)).gt.max(urec(i),urec(j))).and.
     +             (vdiff.gt.vmax.or.drec(i).lt.001)) then
                   flag='*'
                   nflag=nflag+1
                 else
                   flag=' '
                 end if
                 write(11,4) 
     +             x1rec(i),x2rec(i),xirec(i),xirec(j),irec(i),drec(i),
     +             max(urec(i),urec(j)),trec(j)-trec(i),flag
                 write(6,4) 
     +             x1rec(i),x2rec(i),xirec(i),xirec(j),irec(i),drec(i),
     +             max(urec(i),urec(j)),trec(j)-trec(i),flag
                 if(drec(i).lt.xmax) then 
                   do 430 ii=1,nshot
                      if(abs(xshot(ii)-x1rec(i)).lt..001) is1=ii
                      if(abs(xshot(ii)-x2rec(i)).lt..001) is2=ii
430                continue
                   write(12,45) is1,is2,x1rec(i),x2rec(i),
     +                          trec(j)-trec(i)
45                 format(2i10,3f10.3)
                   nmat(is1,is1)=nmat(is1,is1)+1.
                   nmat(is2,is2)=nmat(is2,is2)+1.
                   nmat(is1,is2)=nmat(is1,is2)-1.
                   nmat(is2,is1)=nmat(is2,is1)-1.
                   tvect(is1)=tvect(is1)+trec(j)-trec(i)
                   tvect(is2)=tvect(is2)+trec(i)-trec(j)
                   func=func+(trec(i)-trec(j))**2
                   npairs=npairs+1
                 end if
4                format(2f9.2,2f9.3,i5,3f9.3,a1,f9.3)
                 irec(i)=-irec(i)
                 irec(j)=-irec(j)
                 nused=nused+2
                 if(nused.ge.nrec) go to 99998
               end if
80          continue
70       continue
90    continue
c
99998 continue
c
      write(11,135) nflag
      write(6,135) nflag
135   format(/'number of possible reciprocity violations: ',i6/
     +'***  did you run ORDER_PICKS on tx.in file before hand?  ***'/)
c
      if(npairs.gt.0) then
        write(12,47) xmax,sqrt(func/float(npairs))
        write(6,47) xmax,sqrt(func/float(npairs))
47      format(
     +'ave time diff for reciprocal pairs having |x2-x1|<',
     + f4.1,' km: ',f8.3,' s'/)
c
        call matinv(nmat,nmati,nshot-1)
c
        do 480 i=1,nshot-1
           delta(i)=0.
           do 480 j=1,nshot-1
480           delta(i)=delta(i)+nmati(i,j)*tvect(j)
        delta(nshot)=0.
c
        deltas=0.
        do 490 i=1,nshot
490        deltas=deltas+delta(i)
        deltas=-deltas/float(nshot)
        do 510 i=1,nshot
510        delta(i)=delta(i)+deltas
        write(12,44) (i,i=1,nshot),(delta(i),i=1,nshot)
        write(6,44) (i,i=1,nshot),(delta(i),i=1,nshot)
44      format('time shifts for each shot:'/
     +  <nshot>i8/<nshot>f8.3/)
c
        nused=0
        func=0.
        write(12,555) xmax
        do 990 k=1,ipmax
           do 970 i=1,nrec
              if(abs(irec(i)).ne.k) go to 970
              do 980 j=1,nrec
                 if(abs(x1rec(i)-x2rec(j)).lt..001.and.
     +           abs(x1rec(j)-x2rec(i)).lt..001.and.
     +           abs(irec(i)).eq.abs(irec(j))) then
                   if(drec(i).lt.xmax) then 
                     do 9430 ii=1,nshot
                        if(abs(xshot(ii)-x1rec(i)).lt..001) is1=ii
                        if(abs(xshot(ii)-x2rec(i)).lt..001) is2=ii
9430                 continue
                     write(12,45) is1,is2,x1rec(i),x2rec(i),
     +               trec(j)+delta(is2)-(trec(i)+delta(is1))
                     func=func+(trec(i)+delta(is1)-
     +                    (trec(j)+delta(is2)))**2
                   end if
                   irec(i)=0
                   irec(j)=0
                   nused=nused+2
                   if(nused.ge.nrec) go to 99988
                 end if
980           continue
970        continue
990     continue
c
99988   continue
        write(12,47) xmax,sqrt(func/float(npairs))
        write(6,47) xmax,sqrt(func/float(npairs))
      end if
c
      stop
      end
c
c     ----------------------------------------------------------------
c
      subroutine matinv(a,y,n)
c
c     invert the nxn matrix a  
c
      include 'rayinvr.par'
c
      real a(pshot,pshot),y(pshot,pshot)
      integer indx(pshot)
c
      do 10 i=1,n
         do 20 j=1,n
            y(i,j)=0.
20       continue
         y(i,i)=1.
10    continue
c
      call ludcmp(a,n,indx,d)
c
      do 30 j=1,n
         call lubksb(a,n,indx,y(1,j))
30    continue
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine ludcmp(a,n,indx,d)
c
c     replace a by its LU decomposition
c
      include 'rayinvr.par'
c
      real a(pshot,pshot),vv(pshot)
      integer indx(n)
c
      tiny=1.0e-20
c
      d=1.
      do 10 i=1,n
         aamax=0.
         do 20 j=1,n
            if(abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
20       continue
         if(aamax.eq.0.) then
           write(6,5)
5          format(/'***  singular matrix  ***'/)
           stop
         end if
         vv(i)=1./aamax
10    continue 
      do 30 j=1,n
         if(j.gt.1) then
           do 40 i=1,j-1
              sum=a(i,j)
              if(i.gt.1) then
                do 50 k=1,i-1
                   sum=sum-a(i,k)*a(k,j)
50              continue
                a(i,j)=sum
              end if
40         continue
         end if
         aamax=0.
         do 60 i=j,n
            sum=a(i,j)
            if(j.gt.1) then
              do 70 k=1,j-1
                 sum=sum-a(i,k)*a(k,j)
70            continue
              a(i,j)=sum
            end if
            dum=vv(i)*abs(sum)
            if(dum.ge.aamax) then
              imax=i
              aamax=dum
            end if
60       continue
         if(j.ne.imax) then
           do 80 k=1,n
              dum=a(imax,k)
              a(imax,k)=a(j,k)
              a(j,k)=dum
80         continue
           d=-d
           vv(imax)=vv(j)
         end if
         indx(j)=imax
         if(j.ne.n) then
           if(a(j,j).eq.0.) a(j,j)=tiny
           dum=1./a(j,j)
           do 90 i=j+1,n
              a(i,j)=a(i,j)*dum
90         continue
         end if
30    continue
      if(a(n,n).eq.0.) a(n,n)=tiny
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine lubksb(a,n,indx,b)
c
c     solve the system of n linear equations ax=b
c
      include 'rayinvr.par'
c
      real a(pshot,pshot),b(n)
      integer indx(n)
c
      ii=0
      do 10 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if(ii.ne.0) then
           do 20 j=ii,i-1
              sum=sum-a(i,j)*b(j)
20         continue
         else if(sum.ne.0.) then
           ii=i
         end if
         b(i)=sum
10    continue
      do 30 i=n,1,-1
         sum=b(i)
         if(i.lt.n) then
           do 40 j=i+1,n
              sum=sum-a(i,j)*b(j)
40         continue
         end if
         b(i)=sum/a(i,i)
30    continue
      return
      end
