c                 
c     version 1.2  Mar 1992
c                 
c     ----------------------------------------------------------------
c     |                                                              |
c     |            **********  S P R E A D  ***********              |
c     |                                                              |
c     |                Select the receiver locations                 |
c     |              according to a specified geometry               |
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
c        10 -- input:  traveltime-distance pairs for all receivers
c
c        11 -- output:  traveltime-distance pairs for spread geometry
c                 
c                 
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      open(unit=10, file='tx.in', status='old')
      open(unit=11, file='tx.out') 
c
      write(6,15)
15    format(/
     +  'Enter shot gap and minimum offset (km)')
      read(5,*) xgap,xmin
      write(6,25)
25    format(/
     +  'Enter phase number (enter 0 for all phases)')
      read(5,*) iphase
      xgap=xgap-.001
c
100   read(10,5) xf,tf,uf,irayf
5     format(3f10.3,i10)
      if(irayf.le.0) then
        xshot=xf
        write(11,5) xf,tf,uf,irayf
        if(irayf.lt.0) go to 999
      else
        if(iphase.eq.0.or.iphase.eq.irayf) then
          xdiff=abs(xshot-xf)
          if(xdiff.lt.xgap) go to 100
          if(xdiff.lt.xmin) go to 100
        end if
        write(11,5) xf,tf,uf,irayf
      end if
      go to 100
c
999   stop
      end
