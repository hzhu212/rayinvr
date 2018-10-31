c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             ********  T X P H A S E  ********                |
c     |                                                              |
c     |         Select particular phases from a "tx.in" file         |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
      parameter(nphase=100)
      integer phase(nphase)
      data phase/nphase*0/
      npick=0
      nshot=0
      np=1
c
      open(11, file='tx.in')
      open(12, file='tx.out')
c
1001  write(6,15)
15    format(/'Enter phase number to select (0 to stop)')
      read(5,*) phase(np)
      if(phase(np).gt.0) then
        np=np+1
        go to 1001
      end if
      np=np-1
c
      if(np.eq.0) then
        write(6,25)
25      format(/'***  no phases selected  ***'/)
        stop
      end if
c
101   write(*,*, fmt="(/'number of phases selected: ',i6)") np
c
100   read(11,*,end=999) x,t,u,i
c
      if(i.eq.0) then
        xshotc=x
        tshotc=t
        ushotc=u
        ishotc=i
        isw=0
        go to 100
      end if
c
      if(i.gt.0) then
        iflag=0
        do 10 j=1,np
           if(i.eq.phase(j)) iflag=1
10      continue
        if(iflag.eq.1) then
          if(isw.eq.0) then
            write(12,5) xshotc,tshotc,ushotc,ishotc
            isw=1
            nshot=nshot+1
          end if
          write(12,5) x,t,u,i
          npick=npick+1
        end if
        go to 100
      end if
c
      if(i.eq.-1) then
        write(12,5) x,t,u,i
        go to 999
      end if
5     format(3f10.3,i10)
c
999   write(*,*, fmt="('number of picks: ',i6)") npick
      write(*,*, fmt="('number of shots: ',i6/)") nshot
c
      stop
      end
