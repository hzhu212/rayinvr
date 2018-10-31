c
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ******  S M O O T H_P I C K S  ******             |   
c     |                                                              |
c     |           Smooth a "tx.in" file within each phase            |   
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c     Note: the number of applications of a 3-point averaging filter is
c           the same for each phase if a positive value is input; if a
c           negative value is input the number of applications is equal
c           to the absolute value of the input multipied by the number
c           of points for that particular phase (e.g. if -.5 is input
c           and a phase consists of 100 points, 50 applications are
c           applied)
c
      integer pmax
      parameter(pmax=10000)
      real x(pmax),t(pmax),u(pmax),xh(pmax),th(pmax),uh(pmax),
     +     xw(pmax),tw(pmax),uw(pmax)
      integer ip(pmax),iph(pmax),ipw(pmax)
      character*72 ifname,ofname
c
      write(*,*, fmt="(/'Enter input file name')")
      read(5,15) ifname
15    format(a72)
      write(*,*, fmt=
     +  "(/'Enter output file name (default is input file)')")
      read(5,15) ofname
      if(ofname.eq.'') ofname=ifname
      open(unit=10, file=ifname, status='old')
      open(unit=11, file=ofname)
      write(*,*, fmt="(/
     +'Enter number of smoothing applications (>0) or factor (<0)')")
      read(5,*) rsmth
c
      i=1
100   read(10,5,end=999) x(i),t(i),u(i),ip(i)
5     format(3f10.3,i10)
      i=i+1
      go to 100
c
999   numf=i-1
c
      npicks=0
c
      iline=1
210   iphase=-2
      nphase=0
200   if(ip(iline).le.0.or.ip(iline).ne.iphase) then
        if(nphase.gt.0) then
          if(rsmth.gt.0) then
            nsmth=nint(rsmth)
          else
            nsmth=nint(abs(rsmth)*float(nphase))
            if(nsmth.lt.1) nsmth=1
          end if
c
          do 310 i=1,nsmth
            call smooth(th,nphase)
310       continue
c
          do 10 i=1,nphase
             xw(i)=xh(i)
             tw(i)=th(i)
             uw(i)=uh(i)
             ipw(i)=iph(i)
10        continue
          write(11,5) (xw(i),tw(i),uw(i),ipw(i),i=1,nphase)
          npicks=npicks+nphase
        end if
        if(ip(iline).le.0)
     +    write(11,5) x(iline),t(iline),u(iline),ip(iline)
        if(ip(iline).eq.-1) go to 9999
        if(ip(iline).eq.0) then
          xshot=x(iline)
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
9999  continue
c
      nshot=numf-npicks-1
      write(*,*, fmt="(/'number of lines in file: ',i10)") numf
      write(*,*, fmt="('number of shots:         ',i10)") nshot
      write(*,*, fmt="('number of picks:         ',i10/)") npicks
      stop
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine smooth(x,n) 
c                 
c     three point triangular smoothing filter
c                 
      real x(n) 
      m=n-1       
      a=0.77*x(1)+0.23*x(2) 
      b=0.77*x(n)+0.23*x(m) 
      xx=x(1)     
      xr=x(2)     
      do 10 i=2,m 
         xl=xx    
         xx=xr    
         xr=x(i+1) 
         x(i)=0.54*xx+0.23*(xl+xr) 
 10   continue    
      x(1)=a      
      x(n)=b      
      return      
      end         
