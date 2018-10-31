c
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            *******  O R D E R_P I C K S  *****               |
c     |                                                              |
c     |             Re-order "tx.in" file in order of                |
c     |          increasing x-coordinate within each phase           |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
      integer pmax,pshot
      parameter(pmax=60000, pshot=400)
      real x(pmax),t(pmax),u(pmax),xh(pmax),th(pmax),uh(pmax)
      integer ip(pmax),iph(pmax),ishot(pshot+1)
      character*72 ifname,ofname
c
      write(*,*, fmt="(/'Enter input file name')")
      read(5,15) ifname
15    format(a72)
      write(*,*, fmt=
     +  "(/'Enter output file name (default is input file)')")
      read(5,15) ofname
      if(ofname.eq.'') ofname=ifname
c
      open(10, file=ifname, status='old')
      open(11, file=ofname)
c
      write(*,*, fmt="(/'Enter nominal receiver spacing (km)')")
      read(5,*) xinc
c
      i=1
100   read(10,5,end=999) x(i),t(i),u(i),ip(i)
5     format(3f10.3,i10)
      i=i+1
      go to 100
c
999   numf=i-1
c
      ipmax=0
      nshot=0
      do 310 i=1,numf
         if(ip(i).gt.ipmax) ipmax=ip(i)
         if(ip(i).le.0) then
           nshot=nshot+1
           ishot(nshot)=i
         end if
310   continue
c
      nshot=nshot-1
      do 320 i=1,nshot
         write(11,5) x(ishot(i)),t(ishot(i)),u(ishot(i)),ip(ishot(i))
         do 330 j=1,ipmax
            do 340 k=ishot(i)+1,ishot(i+1)-1
                   if(ip(k).eq.j) write(11,5) x(k),t(k),u(k),ip(k)
340         continue
330      continue
320   continue
      write(11,5) x(ishot(nshot+1)),t(ishot(nshot+1)),u(ishot(nshot+1)),
     +            ip(ishot(nshot+1))
      rewind(11)
c
      i=1
101   read(11,5,end=998) x(i),t(i),u(i),ip(i)
      i=i+1
      go to 101
c
998   rewind(11)
      numf=i-1
      npicks=0
c
      iline=1
210   iphase=-2
      nphase=0
200   if(ip(iline).le.0.or.ip(iline).ne.iphase) then
        if(nphase.gt.0) then
          if(nphase.eq.1) then
            write(11,5) xh(1),th(1),uh(1),iph(1)
            npicks=npicks+1
          else
            nw=nint(abs(xh(nphase)-xh(1))/xinc)+1
            if(nw.eq.1) then
              xwinc=0.
            else
              xwinc=(xh(nphase)-xh(1))/float(nw-1)
            end if
            xwinc2=xwinc/2.
            do 110 i=1,nw
               xpos=xh(1)+float(i-1)*xwinc
               xw=0.
               tw=0.
               uw=0.
               nwc=0
               do j=1,nphase
                  if(abs(xpos-xh(j)).le.xwinc2) then
                    xw=xw+xh(j)
                    tw=tw+th(j)
                    uw=uw+uh(j)
                    nwc=nwc+1
                    if(iph(j).ne.iph(1)) then
                      write(0,*) 'something is wrong'
                      stop
                    end if
                  end if
               enddo
               if(nwc.gt.0) then
                 write(11,5) xw/nwc,tw/nwc,uw/nwc,iph(1)
                 npicks=npicks+1
               end if
110         continue
          end if
        end if
        if(ip(iline).eq.0)
     +    write(11,5) x(iline),t(iline),u(iline),ip(iline)
        if(ip(iline).eq.-1) go to 9999
        if(ip(iline).eq.0) then
          xshot=x(iline)
          iline=iline+1
          go to 210
        else
          iphase=ip(iline)
          nphase=1
          xh(1)=x(iline)
          th(1)=t(iline)
          uh(1)=u(iline)
          iph(1)=ip(iline)
          iline=iline+1
          go to 200
        end if
      else
        nphase=nphase+1
        xh(nphase)=x(iline)
        th(nphase)=t(iline)
        uh(nphase)=u(iline)
        iph(nphase)=ip(iline)
        iline=iline+1
        go to 200
      end if
c
9999  continue
c
      write(*,*, fmt="(/'number of lines in file: ',i10)") numf
      write(*,*, fmt="('number of picks:         ',i10/)") npicks
      stop
      end
