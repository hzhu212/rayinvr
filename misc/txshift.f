c
c     version 1.4  Dec 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             ********  T X S H I F T  ********                |
c     |                                                              |
c     |         Apply a bulk time shift to the picks for             |
c     |              selected shots in a "tx.in" file                |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                   Bullard Laboratories                       |
c     |                  University of Cambridge                     |
c     |                  Cambridge, UK  CB3 0EZ                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
      parameter(nsmax=100)
      real xshot(nsmax),tshift(nsmax)
      data tshift/nsmax*0/
      ns=1
c
      open(11, file='tx.in')
      open(12, file='tx.out')
c
1001  write(6,15)
15    format(/'Enter shot position (km) and time shift (s)'/
     +        '(enter 0,0 to stop)')
      read(5,*) xshot(ns),tshift(ns)
      if(tshift(ns).ne.0.) then
        ns=ns+1
        go to 1001
      end if
      ns=ns-1
c
      if(ns.eq.0) then
        write(6,25)
25      format(/'***  no shots selected  ***'/)
        stop
      end if
c
101   write(*,*, fmt="(/'number of shots selected: ',i6)") ns
c
100   read(11,*,end=999) x,t,u,i
c
      if(i.eq.0) then
        xshotc=x
        nshot=nshot+1
        write(12,5) x,t,u,i
        do 20 i=1,ns
           if(abs(x-xshot(i)).lt..001) then
             tsc=tshift(i)
             go to 100
           end if
20      continue
        tsc=0.
        go to 100
      end if
c
      if(i.gt.0) then
        write(12,5) x,t+tsc,u,i
        npick=npick+1
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
