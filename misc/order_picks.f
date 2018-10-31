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
      real x(pmax),t(pmax),u(pmax),xh(pmax),th(pmax),uh(pmax),
     +     xw(pmax),tw(pmax),uw(pmax)
      integer ip(pmax),iph(pmax),ipw(pmax),ipos(pmax),ishot(pshot+1)
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
          call sort3(xh,ipos,nphase)
          do 10 i=1,nphase
             xw(i)=xh(i)
             tw(i)=th(ipos(i))
             uw(i)=uh(ipos(i))
             ipw(i)=iph(ipos(i))
10        continue
          write(11,5) (xw(i),tw(i),uw(i),ipw(i),i=1,nphase)
          npicks=npicks+nphase
        end if
        if(ip(iline).le.0)
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
      nshot=numf-npicks-1
      write(*,*, fmt="(/'number of lines in file: ',i10)") numf
      write(*,*, fmt="('number of shots:         ',i10)") nshot
      write(*,*, fmt="('number of picks:         ',i10/)") npicks
      stop
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine sort(x,ipos,npts)
c                 
c     sort the elements of array x in order of increasing size using
c     a bubble sort technique
c                 
      parameter(pmax=60000)
      real x(pmax)
      integer ipos(pmax) 
c
      do 30 i=1,npts
         ipos(i)=i
30    continue
c
      do 10 i=1,npts-1
         iflag=0  
         do 20 j=1,npts-1
            if(x(j).gt.x(j+1)) then
              iflag=1
              xh=x(j)
              x(j)=x(j+1)
              x(j+1)=xh
              ih=ipos(j)
              ipos(j)=ipos(j+1)
              ipos(j+1)=ih
            end if
20       continue 
         if(iflag.eq.0) return
10     continue   
      return      
      end         
c
c     ----------------------------------------------------------------
c
      subroutine sort3(ra,rb,n)
c
c     sort the elements of array x in order of increasing size using
c     a heapsort technique
c
      real ra(n)
      integer rb(n)
c
      do 30 i=1,n
         rb(i)=i
30    continue
c
      l=n/2+1
      ir=n
c
10    continue
c
      if(l.gt.1) then
        l=l-1
        rra=ra(l)
        rrb=rb(l)
      else
        rra=ra(ir)
        rrb=rb(ir)
        ra(ir)=ra(1)
        rb(ir)=rb(1)
        ir=ir-1
        if(ir.eq.1) then
          ra(1)=rra
          rb(1)=rrb
          return
        end if
      end if
      i=l
      j=l+l
20    if(j.le.ir) then
        if(j.lt.ir) then
          if(ra(j).lt.ra(j+1)) j=j+1
        end if
        if(rra.lt.ra(j)) then
          ra(i)=ra(j)
          rb(i)=rb(j)
          i=j
          j=j+j
        else 
          j=ir+1
        end if
        go to 20
      end if
      ra(i)=rra
      rb(i)=rrb
      go to 10
c
      end
