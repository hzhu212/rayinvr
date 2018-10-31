c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |              ********   R E D U C E  ********                |   
c     |                                                              |
c     |            Change the reducing velocity applied              |   
c     |                      to a "tx.in" file                       |   
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c     Note: the old or new reducing velocity may be zero for unreduced
c           data
c
      open(unit=10, file='tx.in')
      open(unit=11, file='tx.out')
c
      write(6,55)
55    format(/'Enter old and new reducing velocities (km/s)')
      read(5,*) vred1,vred2
c
      if(vred1.eq.0.) then
        rvred1=0.
      else
        rvred1=1./vred1
      end if
      if(vred2.eq.0.) then
        rvred2=0.
      else
        rvred2=1./vred2
      end if
c
100   read(10,*,end=999) xf,tf,uf,if
5     format(3f10.3,i10)
c
      if(if.eq.0) xshot=xf
      if(if.gt.0) tf=tf+abs(xf-xshot)*(rvred1-rvred2)
      write(11,5) xf,tf,uf,if
      if(if.eq.-1) go to 999
c
      go to 100
c
999   stop
      end
