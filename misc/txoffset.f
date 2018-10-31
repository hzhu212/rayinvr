c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ********  T X O F F S E T  ********               |   
c     |                                                              |
c     |         Convert a "tx.in" file with a shot at 0 km to a      |
c     |              file with shots at any position                 |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
      open(11, file='tx.in', status='old')
      open(12, file='tx.out')
      open(13, file='shots.in', status='old')
c
1000  read(13,*,end=9999) xshot,ishot
c
      if(ishot.eq.-1.or.ishot.eq.2) then
        rewind(11)
c
100     read(11,*,end=999) x,t,u,i
c
        if(i.eq.0) then
          write(12,5) xshot,-1.,0.,0
          go to 100
        end if
c
        if(i.gt.0) then
          write(12,5) xshot-x,t,u,i
          go to 100
        end if
c
      end if
c
      if(ishot.eq.1.or.ishot.eq.2) then
        rewind(11)
c
200     read(11,*,end=999) x,t,u,i
c
        if(i.eq.0) then
          write(12,5) xshot,1.,0.,0
          go to 200
        end if
c
        if(i.gt.0) then
          write(12,5) xshot+x,t,u,i
          go to 200
        end if
c
      end if
c
      go to 1000
c
9999  write(12,5) 0.,0.,0.,-1
5     format(3f10.3,i10)
c
      stop
c
999   write(*,*, fmt="(/'***  end of file reached  ***/')") n
c
      stop
      end
