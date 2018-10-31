c
c     version 1.2  Oct 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             **************  V E X  **************            |
c     |                                                              |
c     |         Extrapolate the constrained part of a model          |
c     |                   to the edges of the model                  |
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
c        11 -- input:  original velocity model
c
c        12 -- output: extrapolated velocity model
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
      character*72 ofname,ifname
c
      write(*,*, fmt=
     +  "(/'Enter input file name (default is v.in)')")
      read(5,15) ifname
15    format(a72)
      if(ifname.eq.'') ifname='v.in'
      write(*,*, fmt="(/
     +'Enter output file name (default is overwriting input file)')")
      read(5,15) ofname
      if(ofname.eq.'') ofname=ifname
c
      open(unit=11, file=ifname, status='old')
      open(unit=12, file=ofname)

      write(6,5)
5     format(/'Enter xmax (km)')
      read(5,*) xmax
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
      do 1030 i=1,nlayer
         if(nvel(i,1).gt.2) then
           j1=0
           do 1040 j=1,nvel(i,1)/2+1
              if(ivarv(i,j,1).ne.0) then
                j1=j
                go to 1050
              end if
1040       continue
1050       if(j1.gt.1) then
             do 1060 j=1,j1-1
                vf(i,j,1)=vf(i,j1,1)
1060         continue
           end if
           j2=999999
           do 2040 j=nvel(i,1),(nvel(i,1)+1)/2,-1
              if(ivarv(i,j,1).ne.0) then
                j2=j
                go to 2050
              end if
2040       continue
2050       if(j2.lt.nvel(i,1)) then
             do 2060 j=j2+1,nvel(i,1)
                vf(i,j,1)=vf(i,j2,1)
2060         continue
           end if
         end if
         if(nvel(i,2).gt.2) then
           j1=0
           do 1070 j=1,nvel(i,2)/2+1
              if(ivarv(i,j,2).ne.0) then
                j1=j
                go to 1080
              end if
1070       continue
1080       if(j1.gt.1) then
             do 1090 j=1,j1-1
                vf(i,j,2)=vf(i,j1,2)
1090         continue
           end if
           j2=999999
           do 2070 j=nvel(i,2),(nvel(i,2)+1)/2,-1
              if(ivarv(i,j,2).ne.0) then
                j2=j
                go to 2080
              end if
2070       continue
2080       if(j2.lt.nvel(i,2)) then
             do 2090 j=j2+1,nvel(i,2)
                vf(i,j,2)=vf(i,j2,2)
2090         continue
           end if
         end if
         if(nzed(i).gt.2) then
           j1=0
           do 3040 j=1,nzed(i)/2+1
              if(ivarz(i,j).ne.0) then
                j1=j
                go to 3050
              end if
3040       continue
3050       if(j1.gt.1) then
             do 3060 j=1,j1-1
                zm(i,j)=zm(i,j1)
3060         continue
           end if
           j2=999999
           do 4040 j=nzed(i),(nzed(i)+1)/2,-1
              if(ivarz(i,j).ne.0) then
                j2=j
                go to 4050
              end if
4040       continue
4050       if(j2.lt.nzed(i)) then
             do 4060 j=j2+1,nzed(i)
                zm(i,j)=zm(i,j2)
4060         continue
           end if
         end if
1030  continue
c
c     write out the extrapolated velocity model
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
