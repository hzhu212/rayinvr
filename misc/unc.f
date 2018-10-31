c                 
c     version 1.3  Apr 1993
c                 
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ************  U N C  ***************              |
c     |                                                              |
c     |              Assign an offset-dependent pick                 |
c     |                       uncertainty                            |
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
c        10 -- input:  initial pick (tx.in) file
c
c        11 -- output:  new pick file
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
     +  'Enter minimum and maximum offset (km)'/
     +  'and uncertainites at each (s)')
      read(5,*) xmin,xmax,umin,umax
      write(6,25)
25    format(/
     +  'Enter phase number (enter 0 for all phases)')
      read(5,*) iphase
      slope=(umax-umin)/(xmax-xmin)
      b=umax-slope*xmax
c
100   read(10,5) xf,tf,uf,irayf
5     format(3f10.3,i10)
      if(irayf.le.0) then
        xshot=xf
        write(11,5) xf,tf,uf,irayf
        if(irayf.lt.0) go to 999
      else
        if(iphase.eq.0.or.iphase.eq.irayf) then
          offset=abs(xshot-xf)
          uf=slope*offset+b
        end if
        write(11,5) xf,tf,uf,irayf
      end if
      go to 100
c
999   stop
      end
