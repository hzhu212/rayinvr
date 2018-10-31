c                 
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             ***********  V P O I S  ***********              |
c     |                                                              |
c     |         Calculate Poisson's ratio given the P- and           |
c     |       and S-wave velocity models in 2 separate files         |
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
c        11 -- input:  P-wave velocity model
c                 
c        12 -- input:  S-wave velocity model
c                 
c        13 -- output: Poisson's ratio model
c                 
c
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'rayinvr.par'
c                
      real xm1(pncntr,ppcntr),zm1(pncntr,ppcntr),
     +     xm2(pncntr,ppcntr),zm2(pncntr,ppcntr),
     +     xvel1(player,ppvel,2),vf1(player,ppvel,2),
     +     xvel2(player,ppvel,2),vf2(player,ppvel,2)
      integer nzed1(pncntr),nvel1(player,2),
     +        nzed2(pncntr),nvel2(player,2),
     +        ivarz(player,ppcntr),ivarv(player,ppvel,2)
c                 
      open(unit=11, file='vp.in', status='old')
      open(unit=12, file='vs.in', status='old')
      open(unit=13, file='vpois.out')
c
      write(6,5)
5     format(/'Enter xmax (km)')
      read(5,*) xmax
c
      nrzmax=ppcntr/10
      nrvmax=ppvel/10
c
c     read in velocity model #1
c
      ncont1=1
      do 170 icont=1,player+1
         nrz=1
         j1=1
         j2=10
11       if(nrz.gt.nrzmax) go to 211
         read(11,55,end=999) ilyr,(xm1(icont,j),j=j1,j2)
         read(11,55,end=999) icnt,(zm1(icont,j),j=j1,j2)
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
         read(11,55,end=999) ilyr,(xvel1(icont,j,1),j=j1,j2)
         read(11,55,end=999) icnt,(vf1(icont,j,1),j=j1,j2)
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
         read(11,55,end=999) ilyr,(xvel1(icont,j,2),j=j1,j2)
         read(11,55,end=999) icnt,(vf1(icont,j,2),j=j1,j2)
         read(11,65,end=999) (ivarv(icont,j,2),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 411
         j1=j1+10
         j2=j2+10
         go to 31
411      ncont1=ncont1+1
170   continue    
c               
99    nlayer1=ncont1-1
c
c     read in velocity model #2
c
      ncont2=1
      do 1170 icont=1,player+1
         nrz=1
         j1=1
         j2=10
111      if(nrz.gt.nrzmax) go to 1211
         read(12,55,end=999) ilyr,(xm2(icont,j),j=j1,j2)
         read(12,55,end=999) icnt,(zm2(icont,j),j=j1,j2)
         read(12,65,end=199) (ivarz(icont,j),j=j1,j2)
         nrz=nrz+1
         if(icnt.ne.1) go to 1211
         j1=j1+10
         j2=j2+10
         go to 111
1211     nrv=1
         j1=1
         j2=10
121      if(nrv.gt.nrvmax) go to 1311
         read(12,55,end=999) ilyr,(xvel2(icont,j,1),j=j1,j2)
         read(12,55,end=999) icnt,(vf2(icont,j,1),j=j1,j2)
         read(12,65,end=999) (ivarv(icont,j,1),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 1311
         j1=j1+10
         j2=j2+10
         go to 121
1311     nrv=1
         j1=1
         j2=10
131      if(nrv.gt.nrvmax) go to 1411
         read(12,55,end=999) ilyr,(xvel2(icont,j,2),j=j1,j2)
         read(12,55,end=999) icnt,(vf2(icont,j,2),j=j1,j2)
         read(12,65,end=999) (ivarv(icont,j,2),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 1411
         j1=j1+10
         j2=j2+10
         go to 131
1411     ncont2=ncont2+1
1170  continue    
c               
199   nlayer2=ncont2-1
c
      if(nlayer1.ne.nlayer2) go to 998
c
      do 171 i=1,ncont1 
         nzed1(i)=1
         nzed2(i)=1
171   continue
      do 172 i=1,nlayer1
         nvel1(i,1)=1
         nvel1(i,2)=1
         nvel2(i,1)=1
         nvel2(i,2)=1
172   continue
c
      do 180 i=1,ncont1
         do 190 j=1,ppcntr 
            if(abs(xm1(i,j)-xmax).lt..0001) go to 180
            nzed1(i)=nzed1(i)+1
190      continue 
180   continue    
c
      do 181 i=1,ncont1
         do 191 j=1,ppcntr 
            if(abs(xm2(i,j)-xmax).lt..0001) go to 181
            nzed2(i)=nzed2(i)+1
191      continue 
181   continue    
c
      do 192 i=1,ncont1
         if(nzed1(i).ne.nzed2(i)) go to 998
192   continue
c
      do 210 i=1,nlayer1  
         do 220 j=1,ppvel  
            if(abs(xvel1(i,j,1)-xmax).lt..0001) go to 210
            nvel1(i,1)=nvel1(i,1)+1
220      continue 
210   continue    
c
      do 216 i=1,nlayer1  
         do 226 j=1,ppvel  
            if(abs(xvel2(i,j,1)-xmax).lt..0001) go to 216
            nvel2(i,1)=nvel2(i,1)+1
226      continue 
216   continue    
c
      do 212 i=1,nlayer1
         if(nvel1(i,1).ne.nvel2(i,1)) go to 998
212   continue
c
      do 240 i=1,nlayer1  
         do 250 j=1,ppvel  
            if(abs(xvel1(i,j,2)-xmax).lt..0001) go to 240
            nvel1(i,2)=nvel1(i,2)+1
250      continue 
240   continue    
c
      do 241 i=1,nlayer1  
         do 251 j=1,ppvel  
            if(abs(xvel2(i,j,2)-xmax).lt..0001) go to 241
            nvel2(i,2)=nvel2(i,2)+1
251      continue 
241   continue    
c
      do 242 i=1,nlayer1
         if(nvel1(i,2).ne.nvel2(i,2)) go to 998
242   continue
c
c     calculate Poisson's ratio
c
      do 1030 i=1,nlayer1
         do 1040 j=1,nvel1(i,1)
            if(xvel1(i,j,1).ne.xvel2(i,j,1)) go to 998
            ratio=(vf2(i,j,1)/vf1(i,j,1))**2
            vf1(i,j,1)=(0.5-ratio)/(1.0-ratio)
            ivarv(i,j,1)=0
1040     continue
         do 1050 j=1,nvel1(i,2)
            if(xvel1(i,j,2).ne.xvel2(i,j,2)) go to 998
            ratio=(vf2(i,j,2)/vf1(i,j,2))**2
            vf1(i,j,2)=(0.5-ratio)/(1.0-ratio)
            ivarv(i,j,2)=0
1050     continue
1030  continue
c
c     write out the difference of the velocity models
c
      do 570 i=1,nlayer1
         nstart=1  
  590    j1=nstart
         j2=j1+9 
         if(j2.gt.nzed1(i)) j2=nzed1(i)
         if(j2.lt.nzed1(i)) then
           icnt=1
         else 
           icnt=0
         end if
         write(13,55) i,(xm1(i,j),j=j1,j2)
         write(13,55) icnt,(zm1(i,j),j=j1,j2)
         write(13,65) (ivarz(i,j),j=j1,j2)
         if(j2.eq.nzed1(i)) go to 600 
         nstart=j2+1
         go to 590
600      nstart=1  
620      j1=nstart
         j2=j1+9 
         if(j2.gt.nvel1(i,1)) j2=nvel1(i,1)
         if(j2.lt.nvel1(i,1)) then
           icnt=1
         else 
           icnt=0
         end if
         write(13,55) i,(xvel1(i,j,1),j=j1,j2)
         write(13,55) icnt,(vf1(i,j,1),j=j1,j2)
         write(13,65) (ivarv(i,j,1),j=j1,j2)
         if(j2.eq.nvel1(i,1)) go to 630 
         nstart=j2+1
         go to 620
630      nstart=1  
650      j1=nstart
         j2=j1+9 
         if(j2.gt.nvel1(i,2)) j2=nvel1(i,2)
         if(j2.lt.nvel1(i,2)) then
           icnt=1
         else 
           icnt=0
         end if
         write(13,55) i,(xvel1(i,j,2),j=j1,j2)
         write(13,55) icnt,(vf1(i,j,2),j=j1,j2)
         write(13,65) (ivarv(i,j,2),j=j1,j2)
         if(j2.eq.nvel1(i,2)) go to 570 
         nstart=j2+1
         go to 650
570   continue

      write(13,55) ncont1,(xm1(ncont1,j),j=1,nzed1(ncont1))
      write(13,55) 0,(zm1(ncont1,j),j=1,nzed1(ncont1))
c
      stop
c
998   write(6,85)
85    format(/'***  models have different parameterization  ***'/)
      stop
c
999   write(6,95)
95    format(/'***  premature end of file  ***'/)
      stop
c
      end
