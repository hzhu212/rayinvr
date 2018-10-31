c                 
c     version 1.2  Nov 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |          **************  V D E P T H  **************         |
c     |                                                              |
c     |    Apply a bulk shift to the z-coordinates of a v.in file    |
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
c        11 -- input:  velocity model
c                 
c        12 -- output: bulk-shifted velocity model
c                 
c
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'rayinvr.par'
c                
      real xm(pncntr,ppcntr),zm(pncntr,ppcntr),
     +     xvel(player,ppvel,2),vf(player,ppvel,2)
      integer nzed(pncntr),nvel(player,2),
     +        ivarz(player,ppcntr),ivarv(player,ppvel,2)
c                 
      open(unit=11, file='v.in', status='old')
      open(unit=12, file='v.out')
c
      write(6,5)
5     format(/'Enter xmax (km)')
      read(5,*) xmax
c
      write(6,6)
6     format(/'Enter bulk shift to apply to z-coordinates (km)')
      read(5,*) zshift
c
      nrzmax=ppcntr/10
      nrvmax=ppvel/10
c
c     read in velocity model
c
      ncont=1
      do 170 icont=1,player+1
         nrz=1
         j1=1
         j2=10
11       if(nrz.gt.nrzmax) go to 211
         read(11,55,end=999) ilyr,(xm(icont,j),j=j1,j2)
         read(11,55,end=999) icnt,(zm(icont,j),j=j1,j2)
         read(11,65,end=99) (ivarz(icont,j),j=j1,j2)
55       format(i2,1x,10f7.2)
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
            if(abs(xvel(i,j,1)-xmax).lt..0001) go to 210
            nvel(i,1)=nvel(i,1)+1
220      continue 
210   continue    
c
      do 240 i=1,nlayer
         do 250 j=1,ppvel  
            if(abs(xvel(i,j,2)-xmax).lt..0001) go to 240
            nvel(i,2)=nvel(i,2)+1
250      continue 
240   continue    
c
      do 1030 i=1,nlayer+1
         do 1040 j=1,nzed(i)
            zm(i,j)=zm(i,j)+zshift
1040     continue
1030  continue
c
c     write out the new velocity model
c
      do 570 i=1,nlayer
         nstart=1  
  590    j1=nstart
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
         if(j2.eq.nvel(i,2)) go to 570 
         nstart=j2+1
         go to 650
570   continue

      write(12,55) ncont,(xm(ncont,j),j=1,nzed(ncont))
      write(12,55) 0,(zm(ncont,j),j=1,nzed(ncont))
c
      stop
c
999   write(6,95)
95    format(/'***  premature end of file  ***'/)
      stop
c
      end
