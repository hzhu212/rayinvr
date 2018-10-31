c
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |               ********   X S H O T  ********                 |   
c     |                                                              |
c     |          Change the shot position and/or direction           |   
c     |                      in a "tx.in" file                       |   
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c
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
c
      write(6,55)
55    format(/'Enter new shot position (km) and direction (+1 or -1)')
      read(5,*) xshotn,idirn
      fidirn=float(idirn)
c
100   read(10,*,end=999) xf,tf,uf,if
5     format(3f10.3,i10)
c
      if(if.eq.0) then
        xw=xshotn
        xshotc=xf
        fidirc=sign(1.,tf)
        tw=fidirn
      end if
      if(if.gt.0) then
        xw=xshotn+fidirn*(xf-xshotc)*fidirc 
        tw=tf
      end if
      if(if.lt.0) then
        xw=xf
        tw=tf
      end if
c
      write(11,5) xw,tw,uf,if
      if(if.eq.-1) go to 999
c
      go to 100
c
999   stop
      end
