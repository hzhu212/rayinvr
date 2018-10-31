c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             ********  T X M A T C H   ********               |   
c     |                                                              |
c     |         select picks with the same positions for two         |
c     |           different phases from two "tx.in" files            |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
      npickm=0
      npicku=0
c
      write(6,55)
55    format(/'Enter file 1 and file 2 phases to match')
      read(5,*) ip1,ip2
c
      open(11, file='tx1.in', status='old')
      open(12, file='tx2.in', status='old')
      open(13, file='tx1.out')
      open(14, file='tx2.out')
c
100   read(11,*,end=999) x1,t1,u1,i1
c
      if(i1.eq.-1) go to 999
      if(i1.eq.0) then
        x1shot=x1
        t1shot=t1
        isw=0
        go to 100
      end if
      if(i1.ne.ip1) go to 100
c
101   read(12,*,end=998) x2,t2,u2,i2
      if(i2.eq.-1) go to 998
      if(i2.eq.0) then
        x2shot=x2 
        t2shot=t2 
        u2shot=u2 
        i2shot=i2 
        if(abs(x2shot-x1shot).lt..001.and.abs(t2shot-t1shot).lt..001)
     +  then
          iflag=1
        else
          iflag=0
        end if
        go to 101
      end if
      if(i2.ne.ip2) go to 101
      if(iflag.eq.1) then
        if(abs(x2-x1).lt..001) then
          if(isw.eq.0) then
            write(13,5) x2shot,t2shot,u2shot,i2shot
            write(14,5) x2shot,t2shot,u2shot,i2shot
            isw=1
            nshot=nshot+1
          end if
          write(13,5) x1,t1,u1,i1
          write(14,5) x2,t2,u2,i2
          npickm=npickm+1
          go to 100
        end if
      end if
      go to 101
c
998   npicku=npicku+1
      rewind(12)
      go to 100
999   write(13,5) x1,t1,u1,i1
      write(14,5) x1,t1,u1,i1
5     format(3f10.3,i10)
c
      write(*,*, fmt="('number of picks matched:   ',i6)") npickm
      write(*,*, fmt="('number of picks unmatched: ',i6)") npicku
      write(*,*, fmt="('number of shots:           ',i6/)") nshot
c
      stop
      end
