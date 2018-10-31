c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |               ********  T X 2 R E C  ********                |   
c     |                                                              |
c     |         Obtain a "rec.in" file for input to TRAMP from       |
c     |                       a "tx.in" file                         |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
      parameter(nmax=10000)
      real xshot(nmax),xrec(nmax)
      open(11, file='rec.in')
      open(12, file='rec.out')
      n=0
c
100   read(11,*,end=999) x,t,u,i
c
      if(i.eq.0) then
        write(12,5) x,1
        xshotc=x
        go to 100
      end if
c
      if(i.gt.0) then
        if(n.gt.0) then
          do 10 i=1,n
             if(abs(xshotc-xshot(i)).lt..001.and.abs(x-xrec(i)).lt..001)
     +       go to 100
10        continue
        end if 
        n=n+1
        xshot(n)=xshotc
        xrec(n)=x
        write(12,5) x,0
        go to 100
      end if
c
      if(i.eq.-1) then
        write(12,5) 0.,-1
        go to 999
      end if
5     format(f10.3,i10)
c
999   write(*,*, fmt="(/'number of receivers: ',i6/)") n
      stop
      end
