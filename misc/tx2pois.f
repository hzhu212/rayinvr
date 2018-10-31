c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ***********  T X 2 P O I S  *********             |
c     |                                                              |
c     |        Calculate Poisson's ratio from P- and S-wave          |
c     |          travel times in two separate tx.in files            |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
      npick=0
      nshot=0
c
      open(11, file='txp.in', status='old')
      open(12, file='txs.in', status='old')
      open(13, file='tx.out')
c
100   read(11,*,end=999) xp,tp,up,ip
      read(12,*,end=999) xs,ts,us,is
c
      if(ip.eq.-1) go to 999
      if(ip.eq.0) then
        write(13,5) xp,tp,up,ip
        xshot=xp
        nshot=nshot+1
        go to 100
      end if
c
      ratio=(tp/ts)**2
      pois=(0.5-ratio)/(1.0-ratio)
      ratio1=((tp-up)/(ts+us))**2
      pois1=(0.5-ratio1)/(1.0-ratio1)
      ratio2=((tp+up)/(ts-us))**2
      pois2=(0.5-ratio2)/(1.0-ratio2)
      upois=abs(pois1-pois2)
      if(pois.ge..2.and.pois.le..4.and.upois.lt..05)  then
        write(13,5) (xp+xshot)/2.,pois,upois,ip
        npick=npick+1
      end if
      go to 100
c
999   write(13,5) xp,tp,up,ip
5     format(3f10.3,i10)
c
      write(*,*, fmt="('number of picks matched:   ',i6)") npick
      write(*,*, fmt="('number of shots:           ',i6/)") nshot
c
      stop
      end
