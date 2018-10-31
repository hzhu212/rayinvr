c                 
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            *******  D M P L S T S Q R   *******              |
c     |                                                              |
c     |              Damped least-squares solution of                |
c     |                traveltime inversion problem                  |
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
c        10 -- input:  damping parameters
c
c        11 -- input:  matrix of partial derivatives and vector of
c                      traveltime residuals
c                 
c        12 -- output: model parameter adjustments and resolution and
c                      covariance estimates and updated model 
c                 
c        13 -- input/output: initial and updated velocity model
c
c        30 -- input:  initial and updated floating reflectors
c                 
c                 
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'rayinvr.par'
c                
      real apart(prayi,pnvar),tres(prayi),ata(pnvar,pnvar),
     +     att(pnvar),atad(pnvar,pnvar),atadi(pnvar,pnvar),
     +     dx(pnvar),pvar(3),res(pnvar,pnvar),var(pnvar),
     +     parorg(pnvar),xm(pncntr,ppcntr),zm(pncntr,ppcntr),
     +     xvel(player,ppvel,2),vf(player,ppvel,2),grad(player*ppvel),
     +     tunc(prayi),thick(player*ppcntr),parunc(pnvar),
     +     xfrefl(pfrefl,ppfref),zfrefl(pfrefl,ppfref)
      integer partyp(pnvar),nzed(pncntr),nvel(player,2),
     +        ivarz(player,ppcntr),ivarv(player,ppvel,2),
     +        mpinch(pncntr,ppcntr,2),igrad(player*ppvel),
     +        ivarf(pfrefl,ppfref),npfref(pfrefl)
c                 
      namelist /dmppar/ iscrn,dmpfct,velunc,bndunc,xmax
c     
      open(unit=10, file='d.in', status='old')
      open(unit=11, file='i.out', status='old')
      open(unit=12, file='d.out')
      open(unit=13, file='v.in', status='old')
      open(unit=14, file='v.bak')
c                 
c     default parameter values
c
      iscrn=1
      dmpfct=1.0
      velunc=0.1
      bndunc=0.1
      xmax=-99999.
      ifrbnd=0
c
      read(10,dmppar)
c
      if(xmax.lt.-99998) then
        write(6,105)
105     format(/'***  xmax not specified  ***'/)
        stop
      end if
c                 
c     read in matrix of partial derivatives and vector of traveltime
c     residuals
c
      read(11,1)
1     format(' ')
      read(11,115) narinv,nvar
115   format(2i5)
      read(11,1)
      do 6 i=1,nvar
         read(11,5) partyp(i),parorg(i),parunc(i)
5        format(i5,2f15.5)
6     continue
      read(11,1)
      do 10 i=1,narinv
         read(11,15) (apart(i,j),j=1,nvar)
15       format(5e12.5)
10    continue
      read(11,1)
      read(11,15) (tres(i),i=1,narinv)
      read(11,1)
      read(11,15) (tunc(i),i=1,narinv)
c                 
      do 20 i=1,nvar
         do 30 j=1,i
            ata(i,j)=0.
            do 40 k=1,narinv
               ata(i,j)=ata(i,j)+apart(k,i)*apart(k,j)/tunc(k)**2
40          continue 
            if(i.ne.j) ata(j,i)=ata(i,j)
30       continue
20    continue
c
      do 50 i=1,nvar
         att(i)=0.
         do 60 j=1,narinv
            att(i)=att(i)+apart(j,i)*tres(j)/tunc(j)**2
60       continue
50    continue
c
      pvar(1)=bndunc**2
      pvar(2)=velunc**2
      pvar(3)=bndunc**2
c
      do 70 i=1,nvar
         do 80 j=1,nvar
            atad(i,j)=ata(i,j)
            if(i.eq.j) then
              if(partyp(i).eq.1) then
                if(bndunc.le.0.) then
                  parunc(i)=parunc(i)**2
                else
                  parunc(i)=pvar(1)
                end if 
              end if
              if(partyp(i).eq.2) then
                if(velunc.le.0.) then
                  parunc(i)=parunc(i)**2
                else
                  parunc(i)=pvar(2)
                end if 
              end if
              if(partyp(i).eq.3) then
                if(bndunc.le.0.) then
                  parunc(i)=parunc(i)**2
                else
                  parunc(i)=pvar(1)
                end if 
                ifrbnd=1
              end if
              atad(i,j)=atad(i,j)+dmpfct/parunc(i)
            end if
80       continue
70    continue
c
c
      call matinv(atad,atadi,nvar)
c
c
      do 90 i=1,nvar
         dx(i)=0.
         do 100 j=1,nvar
            dx(i)=dx(i)+atadi(i,j)*att(j)
100      continue
90    continue
c
      do 110 i=1,nvar
         do 120 j=1,nvar
            res(i,j)=0.
            do 130 k=1,nvar
               res(i,j)=res(i,j)+atadi(i,k)*ata(k,j)
130         continue
120      continue
110   continue
c
      do 140 i=1,nvar
         var(i)=var(i)+(1.-res(i,i))*parunc(i)
140   continue
c
      write(12,25) dmpfct
25    format(/'overall damping factor: ',f10.5)
      write(12,35)
35    format(/'type  orig. val.  uncert.    adjust.   new val.',
     +        ' resolution std. error')
      write(12,45) (partyp(i),parorg(i),parunc(i)**.5,dx(i),
     +              parorg(i)+dx(i),res(i,i),sqrt(var(i)),i=1,nvar)
45    format(i3,6f11.4)            
c
      if(iscrn.eq.1) then
        write(6,25) dmpfct
        write(6,35)
        write(6,45) (partyp(i),parorg(i),parunc(i)**.5,dx(i),parorg(i)+
     +               dx(i),res(i,i),sqrt(var(i)),i=1,nvar)
      end if
c
      ncont=1
      nrzmax=ppcntr/10
      nrvmax=ppvel/10
      do 170 icont=1,player+1
         nrz=1
         j1=1
         j2=10
11       if(nrz.gt.nrzmax) go to 211
         read(13,55,end=999) ilyr,(xm(icont,j),j=j1,j2)
         read(13,55,end=999) icnt,(zm(icont,j),j=j1,j2)
         read(13,65,end=99) (ivarz(icont,j),j=j1,j2)
55       format(i2,1x,10f7.3)
65       format(3x,10i7)
         nrz=nrz+1
         if(icnt.ne.1) go to 211
         j1=j1+10
         j2=j2+10
         go to 11
211      nrv=1
         j1=1
         j2=10
21       if(nrv.gt.nrvmax) go to 311
         read(13,55,end=999) ilyr,(xvel(icont,j,1),j=j1,j2)
         read(13,55,end=999) icnt,(vf(icont,j,1),j=j1,j2)
         read(13,65,end=999) (ivarv(icont,j,1),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 311
         j1=j1+10
         j2=j2+10
         go to 21
311      nrv=1
         j1=1
         j2=10
31       if(nrv.gt.nrvmax) go to 411
         read(13,55,end=999) ilyr,(xvel(icont,j,2),j=j1,j2)
         read(13,55,end=999) icnt,(vf(icont,j,2),j=j1,j2)
         read(13,65,end=999) (ivarv(icont,j,2),j=j1,j2)
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
c     write out original velocity model
c
      do 1570 i=1,nlayer
         nstart=1  
1590     j1=nstart
         j2=j1+9 
         if(j2.gt.nzed(i)) j2=nzed(i)
         if(j2.lt.nzed(i)) then
           icnt=1
         else 
           icnt=0
         end if
         write(14,55) i,(xm(i,j),j=j1,j2)
         write(14,55) icnt,(zm(i,j),j=j1,j2)
         write(14,65) (ivarz(i,j),j=j1,j2)
         if(j2.eq.nzed(i)) go to 1600 
         nstart=j2+1
         go to 1590
1600     if(nvel(i,1).le.0) then
           write(14,55) i,xmax
           write(14,55) 0,0.
           write(14,65) 0
           go to 1630
         end if
         nstart=1  
1620     j1=nstart
         j2=j1+9 
         if(j2.gt.nvel(i,1)) j2=nvel(i,1)
         if(j2.lt.nvel(i,1)) then
           icnt=1
         else 
           icnt=0
         end if
         write(14,55) i,(xvel(i,j,1),j=j1,j2)
         write(14,55) icnt,(vf(i,j,1),j=j1,j2)
         write(14,65) (ivarv(i,j,1),j=j1,j2)
         if(j2.eq.nvel(i,1)) go to 1630 
         nstart=j2+1
         go to 1620
1630     if(nvel(i,2).le.0) then
           write(14,55) i,xmax
           write(14,55) 0,0.
           write(14,65) 0
           go to 1570
         end if
         nstart=1  
1650     j1=nstart
         j2=j1+9 
         if(j2.gt.nvel(i,2)) j2=nvel(i,2)
         if(j2.lt.nvel(i,2)) then
           icnt=1
         else 
           icnt=0
         end if
         write(14,55) i,(xvel(i,j,2),j=j1,j2)
         write(14,55) icnt,(vf(i,j,2),j=j1,j2)
         write(14,65) (ivarv(i,j,2),j=j1,j2)
         if(j2.eq.nvel(i,2)) go to 1570 
         nstart=j2+1
         go to 1650
1570  continue
c
      write(14,55) ncont,(xm(ncont,j),j=1,nzed(ncont))
      write(14,55) 0,(zm(ncont,j),j=1,nzed(ncont))
c
c     read in floating reflectors
c
      if(ifrbnd.eq.1) then
        open(unit=30, file='f.in', status='old')
        open(unit=31, file='f.bak')
        nfrefl=0
690     read(30,545,end=595) nfrefr
        if(nfrefr.lt.2.or.nfrefr.gt.ppfref) then
          write(6,585)
585       format(/'***  error in f.in file  ***'/)
          stop
        end if
        nfrefl=nfrefl+1
        npfref(nfrefl)=nfrefr
545     format(i2)
        read(30,555) (xfrefl(nfrefl,i),i=1,npfref(nfrefl))
        read(30,555) (zfrefl(nfrefl,i),i=1,npfref(nfrefl))
        read(30,575) (ivarf(nfrefl,i),i=1,npfref(nfrefl))
555     format(3x,<ppfref>f7.3)
575     format(3x,<ppfref>i7)
        go to 690
595     continue
c
c       write out original floating reflectors
c
        if(ifrbnd.eq.1) then
          do 1790 j=1,nfrefl
             write(31,545) npfref(j)
             write(31,565) j,(xfrefl(j,i),i=1,npfref(j))
             write(31,555) (zfrefl(j,i),i=1,npfref(j))
             write(31,575) (ivarf(j,i),i=1,npfref(j))
1790      continue
        end if
      end if
c
c     check for fixed velocity gradients
c
      ngrad=0
      do 370 i=1,nlayer
         iflagg=0
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
                    iflagg=1
                  end if
                else
                  grad(ngrad)=0.
                  igrad(ngrad)=1
                end if
              end if 
380        continue
         end if
c
         if(iflagg.eq.1) then
           write(6,85) i
           write(12,85) i
85         format(/'***  check gradient in layer ',i2,'  ***'/)
         end if
c
370   continue
c
      rewind(13)
c
      write(12,75)
      if(iscrn.eq.1) write(6,75)
75    format(/'velocity model:'/)
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
                    if(abs(xm(i,j)-xm(k,l)).lt..005.and.
     +                 abs(zm(i,j)-zm(k,l)).lt..005) then
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
c     add the parameter adjustments
c
      nvarw=0
      nlyr=0
      do 270 i=1,nlayer
         do 280 j=1,nzed(i)
            if(ipinch.eq.1.and.i.gt.1) then
              if(mpinch(i,j,1).gt.0) then
                zm(i,j)=zm(mpinch(i,j,1),mpinch(i,j,2))
                go to 280
              end if
            end if
            if(ivarz(i,j).gt.0) then
              nvarw=nvarw+1
              zm(i,j)=zm(i,j)+dx(nvarw)
            end if
            if(ivarz(i,j).eq.-1) then
              nlyr=nlyr+1
              zm(i,j)=zm(i-1,j)+thick(nlyr)
            end if
280      continue
c
         if(nvel(i,1).gt.0) then
           do 310 j=1,nvel(i,1)
              if(ivarv(i,j,1).gt.0) then
                nvarw=nvarw+1
                vf(i,j,1)=vf(i,j,1)+dx(nvarw)
              end if
310        continue
         else
           nvel(i,1)=1
           xvel(i,1,1)=xmax
           vf(i,1,1)=0.
         end if
c
         if(nvel(i,2).gt.0) then
           do 430 j=1,nvel(i,2)
              if(ivarv(i,j,2).gt.0) then
                nvarw=nvarw+1
                vf(i,j,2)=vf(i,j,2)+dx(nvarw)
              end if
430        continue
         else
           nvel(i,2)=1
           xvel(i,1,2)=xmax
           vf(i,1,2)=0.
         end if
270   continue
c
      if(ifrbnd.eq.1) then
        do 880 j=1,nfrefl
           do 890 i=1,npfref(j)
              if(ivarf(j,i).eq.1) then
                nvarw=nvarw+1
                zfrefl(j,i)=zfrefl(j,i)+dx(nvarw)
              end if
890        continue
880     continue
      end if
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
c     write out velocity model
c
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
         write(13,55) i,(xm(i,j),j=j1,j2)
         write(13,55) icnt,(zm(i,j),j=j1,j2)
         write(13,65) (ivarz(i,j),j=j1,j2)
         if(iscrn.eq.1) then
           write(6,55) i,(xm(i,j),j=j1,j2)
           write(6,55) icnt,(zm(i,j),j=j1,j2)
           write(6,65) (ivarz(i,j),j=j1,j2)
         end if
         if(j2.eq.nzed(i)) go to 600 
         nstart=j2+1
         go to 590
600      nstart=1  
620      j1=nstart
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
         write(13,55) i,(xvel(i,j,1),j=j1,j2)
         write(13,55) icnt,(vf(i,j,1),j=j1,j2)
         write(13,65) (ivarv(i,j,1),j=j1,j2)
         if(iscrn.eq.1) then
           write(6,55) i,(xvel(i,j,1),j=j1,j2)
           write(6,55) icnt,(vf(i,j,1),j=j1,j2)
           write(6,65) (ivarv(i,j,1),j=j1,j2)
         end if
         if(j2.eq.nvel(i,1)) go to 630 
         nstart=j2+1
         go to 620
630      nstart=1  
650      j1=nstart
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
         write(13,55) i,(xvel(i,j,2),j=j1,j2)
         write(13,55) icnt,(vf(i,j,2),j=j1,j2)
         write(13,65) (ivarv(i,j,2),j=j1,j2)
         if(iscrn.eq.1) then
           write(6,55) i,(xvel(i,j,2),j=j1,j2)
           write(6,55) icnt,(vf(i,j,2),j=j1,j2)
           write(6,65) (ivarv(i,j,2),j=j1,j2)
         end if
         if(j2.eq.nvel(i,2)) go to 570 
         nstart=j2+1
         go to 650
570   continue
c
      write(12,55) ncont,(xm(ncont,j),j=1,nzed(ncont))
      write(12,55) 0,(zm(ncont,j),j=1,nzed(ncont))
      write(13,55) ncont,(xm(ncont,j),j=1,nzed(ncont))
      write(13,55) 0,(zm(ncont,j),j=1,nzed(ncont))
      if(iscrn.eq.1) then
        write(6,55) ncont,(xm(ncont,j),j=1,nzed(ncont))
        write(6,55) 0,(zm(ncont,j),j=1,nzed(ncont))
      end if
c
c     write out floating reflectors
c
      if(ifrbnd.eq.1) then
        write(12,175)
        if(iscrn.eq.1) write(6,175)
175     format(/'floating reflectors:'/)
        rewind(30)
        do 790 j=1,nfrefl
           write(30,545) npfref(j)
           write(12,545) npfref(j)
           write(30,565) j,(xfrefl(j,i),i=1,npfref(j))
           write(12,565) j,(xfrefl(j,i),i=1,npfref(j))
565        format(i2,1x,<ppfref>f7.3)
           write(30,555) (zfrefl(j,i),i=1,npfref(j))
           write(12,555) (zfrefl(j,i),i=1,npfref(j))
           write(30,575) (ivarf(j,i),i=1,npfref(j))
           write(12,575) (ivarf(j,i),i=1,npfref(j))
           if(iscrn.eq.1) then
             write(6,545) npfref(j)
             write(6,565) j,(xfrefl(j,i),i=1,npfref(j))
             write(6,555) (zfrefl(j,i),i=1,npfref(j))
             write(6,575) (ivarf(j,i),i=1,npfref(j))
           end if
790     continue
      end if
c
      stop
c
999   write(6,95)
95    format(/'***  error in velocity model  ***'/)
      stop
c
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
      real a(pnvar,pnvar),y(pnvar,pnvar)
      integer indx(pnvar)
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
      real a(pnvar,pnvar),vv(pnvar)
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
      real a(pnvar,pnvar),b(n)
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
