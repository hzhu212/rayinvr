c
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ***********  I S E I S  ************              |
c     |                                                              |
c     |               Interpolate traveltime picks at                |
c     |          uniform or specified seismogram locations           |
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
c        11 -- output:  interpolated traveltime-distance pairs
c
c        12 -- input:  observed seismogram locations (uneven spacing)
c
c
c     ----------------------------------------------------------------
c
c
      program main
c
      real seis(10000),x(10000),t(10000),u(10000),
     +     xh(10000),th(10000),uh(10000),uncert(10000)
      integer ir(10000),irhh(10000),ipseis(10000)
      character*72 ifname,ofname
c
      write(*,*, fmt="(/'Enter input file name')")
      read(5,85) ifname
85    format(a72)
      write(*,*, fmt=
     +  "(/'Enter output file name (default is input file)')")
      read(5,85) ofname
      if(ofname.eq.'') ofname=ifname
      open(unit=10, file=ifname, status='old')
      open(unit=11, file=ofname)
c
      nseis=0
      write(6,45)
45    format(/'Enter  0  for uniform or  1  for irregular spacing')
      read(5,*) ispace
      if(ispace.eq.0) then
        write(6,15)
15      format(/'Enter minimum and maximum seismogram location and'
     +         /'station increment (km)')
        read(5,*) xmins,xmaxs,xincs
        nseis=nint((xmaxs-xmins)/xincs)+1
        do 10 i=1,nseis
           seis(i)=xmins+float(i-1)*xincs
10      continue
      else
        write(6,55)
55      format(
     +    'Enter  0  for single or  1  for multiple set of receivers')
        read(5,*) iset
        open(unit=12, file='rec.in', status='old')
        if(iset.eq.0) then
          i=1
101       read(12,25,end=99) seis(i)
25        format(f10.3)
          i=i+1
          go to 101
99        nseis=i-1
        end if
      end if
c
      if(nseis.gt.0) write(6,35) nseis
35    format(/'number of seismograms: ',i4/)
c
      iflag=0
      npts=0
100   read(10,5) xf,tf,uf,irayf
5     format(3f10.3,i10)
      if(irayf.lt.0) then
        iflag=2
        go to 500
      end if
      if(irayf.eq.0) then
500     if(iflag.ne.0.and.npts.gt.1) then
          if(ispace.eq.1.and.iset.eq.1) then
            iflags=0
1001        read(12,5) xs,ts,us,is
            if(is.lt.0) go to 998
            if(is.eq.0) then
              if(iflags.eq.1) go to 998
              if(abs(xs-xshot).lt..005.and.nint(ts).eq.idr) then
                iflags=1
                nseis=0
              end if
            else
              if(iflags.eq.1) then
                nseis=nseis+1
                seis(nseis)=xs
                uncert(nseis)=us
                ipseis(nseis)=is
              end if
            end if
            go to 1001
998         if(nseis.eq.0) then
              write(6,65) xshot,idr
65            format(/'***  no receivers for shot at ',f7.2,
     +                ' km and direction ',i2,'  ***'/)
              rewind(12)
              go to 400
            else
              write(6,75) xshot,idr,nseis
75            format('shot ',f7.2,' km   direction ',i2,
     +               '   number of seismograms: ',i4)
              rewind(12)
            end if
          end if
          irh=ir(1)
          npt=1
          xh(1)=x(1)
          th(1)=t(1)
          uh(1)=u(1)
          irhh(1)=ir(1)
          ipos=2
200       if(ipos.gt.npts) go to 300
          if(ir(ipos).eq.irh) then
            npt=npt+1
            xh(npt)=x(ipos)
            th(npt)=t(ipos)
            uh(npt)=u(ipos)
            irhh(npt)=ir(ipos)
            ipos=ipos+1
            go to 200
          end if
300       if(npt.gt.1) then
            do 30 i=1,nseis
               tmin=999999.
               do 40 j=1,npt-1
                  if((seis(i).ge.xh(j).and.seis(i).le.xh(j+1)).or.
     +            (seis(i).le.xh(j).and.seis(i).ge.xh(j+1))) then
                    if(ispace.eq.1.and.iset.eq.1.and.
     +                ipseis(i).ne.irhh(j)) go to 40
                    if(xh(j+1)-xh(j).ne.0.) then
                      time=(th(j+1)-th(j))*(seis(i)-xh(j))/(xh(j+1)-
     +                      xh(j))+th(j)
                      if(ispace.ne.1.and.iset.ne.1) then
                        unc=(uh(j+1)-uh(j))*(seis(i)-xh(j))/(xh(j+1)-
     +                        xh(j))+uh(j)
                      else
                        unc=uncert(i)
                      end if
                    else
                      time=(th(j+1)+th(j))/2.
                      if(ispace.ne.1.and.iset.ne.1) then
                        unc=(uh(j+1)+uh(j))/2.
                      else
                        unc=uncert(i)
                      end if
                    end if
                    if(time.lt.tmin) then
                      tmin=time
                      umin=unc
                    end if
                  end if
40             continue
               if(tmin.lt.999998.) then
                 write(11,5) seis(i),tmin,umin,irh
               end if
30          continue
          end if
          if(ipos.gt.npts) go to 400
          irh=ir(ipos)
          npt=1
          xh(1)=x(ipos)
          th(1)=t(ipos)
          uh(1)=u(ipos)
          ipos=ipos+1
          go to 200
        end if
400     write(11,5) xf,tf,0.0,irayf
        if(iflag.eq.2) go to 999
        npts=0
        iflag=1
        xshot=xf
        idr=nint(tf)
      else
        npts=npts+1
        x(npts)=xf
        t(npts)=tf
        u(npts)=uf
        ir(npts)=irayf
      end if
      go to 100
c
999   stop
      end
