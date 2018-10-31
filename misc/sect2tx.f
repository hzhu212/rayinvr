c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ********   S E C T 2 T X   ********               |   
c     |                                                              |
c     |         Obtain a "tx.in" file from a "sect.out" file         |   
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
      open(11, file='sect.out')
      open(12, file='tx.out')
c
      read(11,*) xshot
      write(12,5) xshot,0.,0.,0
5     format(3f10.3,i10)
100   read(11,*,end=999) dist,na
      if(na.gt.0) then
        tmin=999999.
        do 10 i=1,na
           read(11,*) t,a,p
           if(t.lt.tmin) tmin=t
10      continue
        write(12,5) dist,tmin,0.,1
      end if
      go to 100
999   write(12,5) 0.,0.,0.,-1
      stop
      end 
