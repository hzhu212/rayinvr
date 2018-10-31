c                 
c     version 1.2  May 1998
c                 
c     ----------------------------------------------------------------
c     |                                                              |
c     |             *********  T X Z E R O  **********               |
c     |                                                              |
c     |           Add Gaussian noise to traveltime picks             |
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
c        10 -- input:  original traveltime-distance pairs
c
c        11 -- output:  noisy traveltime-distance pairs
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
100   read(10,35) xr,tr,ur,ir
35    format(3f10.3,i10)
      if(ir.le.0) then
        if(ir.lt.0) write(11,35) xr,tr,ur,ir
        if(ir.lt.0) go to 999
      else
        write(11,35) xr,1.,0.,0
        write(11,35) xr,tr,ur,ir
      end if
      go to 100
c  
999   stop 
      end
