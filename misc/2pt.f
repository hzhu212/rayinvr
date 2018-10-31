c
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |                   ******  2 P T  ******                      |   
c     |                                                              |
c     |   Produce a bit-mask file using two output files generated   |   
c     |  by RAYINVR to allow a two-point minimum-traveltime raypath  |
c     |              diagram to be produced by RAYPLOT               |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
      parameter(nmax=100000)
      real xshot(nmax),ivray(nmax),ig(nmax),dist(nmax),time(nmax),
     +     xr(nmax),ivr(nmax),igr(nmax),dr(nmax),tr(nmax)
c
      open(unit=11, file='ra1.out', status='old')
      open(unit=12, file='ra2.out', status='old')
      open(unit=13, file='2pt.out')
c
      i=1
100   read(12,5,end=99) xshot(i),ivray(i),ig(i),dist(i),time(i)
5     format(f12.4,2i4,2f12.4)
      i=i+1 
      go to 100
99    n12=i-1
c
      i=1
101   read(11,5,end=999) xr(i),ivr(i),igr(i),dr(i),tr(i)
      i=i+1
      go to 101
999   n11=i-1
c
      do 10 i=1,n11
         if(ivr(i).eq.0) go to 104
         do 20 j=1,n12
            if(abs(xshot(j)-xr(i)).lt..001.and.ivray(j).eq.ivr(i).and.
     +      ig(j).eq.igr(i).and.abs(dist(j)-dr(i)).lt..001) go to 102 
20       continue
         go to 103
c
102      do 30 k=1,n11
            if(k.eq.i) go to 30
            if(abs(xr(i)-xr(k)).lt..001.and.ivr(i).eq.ivr(k).and.
     +      igr(i).eq.igr(k).and.abs(dr(i)-dr(k)).lt..001) then
              if(tr(k).lt.tr(i)) go to 103
            end if
30       continue
c
104      write(13,15) 1
15       format(i1)
         go to 10
c
103      write(13,15) 0
10    continue
c
      stop
      end
