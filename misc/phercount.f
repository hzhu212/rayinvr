c
c     version 1.4  Dec 1993
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            *******   P H E R C O U N T  *******              |   
c     |                                                              |
c     |        Calculate number of arrivals for each phase           |   
c     |         and average uncertainites for "tx.in" file           |   
c     |                                                              |
c     |              Original code written by B. Zelt                |
c     |                                                              |
c     |                 Modifications by C. A. Zelt                  |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c     written by T. Thing (aka B. ZELT). (c) 1991 Thing Geophysical.
c     counts # of arrivals for each phase for each shot in tx.in and gives
c     grand totals.  also calculates ave uncertainties for each phase for
c     each shot and gives overall averages.
c     counts of picks are written to file "phcount.out" (and to screen).
c     summary of uncertainties are written to file "ercount.out".
c
      parameter(np=99)
      integer*4 k(np)/np*0/,kg(np)/np*0/
      real*4 ke(np)/np*0/,keg(np)/np*0/,err(np),errg(np),ket,kegt,ketkt
      open(7,file='tx.in',status='old')
      open(8,file='phcount.out',status='unknown')
      open(9,file='ercount.out',status='unknown')
c  
      ipmax=0
100   read(7,*,end=998) d,t,e,iphase
      if(iphase.gt.ipmax) ipmax=iphase
      if(iphase.gt.np) then
        write(6,895)
895     format(/'***  phase number greater than array size  ***'/)
        stop
      end if
      go to 100
998   rewind(7)
c
      if(ipmax.gt.1) then
        write(6,40) (i,i=1,ipmax)
        write(8,40) (i,i=1,ipmax)
        write(9,40) (i,i=1,ipmax)
 40     format(/,'     xshot      ',
     &  i2,<ipmax-1>i5/,
     &  '-------------',<ipmax>('-----'),'-----------')
      else
        write(6,41) (i,i=1,ipmax)
        write(8,41) (i,i=1,ipmax)
        write(9,41) (i,i=1,ipmax)
 41     format(/,'     xshot      ',
     &  i2/,
     &  '-------------',<ipmax>('-----'),'-----------')
      end if
      nshot=0
 10   continue
      read(7,*,end=999) d,t,e,iphase
      if(iphase.lt.0) iphase=0
      if(iphase.eq.0) then
        kt=0.
        ket=0.
        do 50 i=1,ipmax
           kt=kt+k(i)
50         ket=ket+ke(i)
        if(nshot.ge.1) then
          write(6,15) xshot,(k(i),i=1,ipmax),kt
          write(8,15) xshot,(k(i),i=1,ipmax),kt
 15     format(f10.3,5x,<ipmax>(i3,2x),' = ',2x,i4)
          do 888 i=1,ipmax
            if(k(i).eq.0) then
              err(i)=0.
            else
              err(i)=ke(i)/k(i)
            end if
 888      continue
          if(kt.eq.0.) then
            ketkt=0.
          else
            ketkt=ket/kt
          end if
          write(9,115) xshot,(err(i),i=1,ipmax),ketkt
 115      format(f10.3,5x,<ipmax>(f4.3,1x),' = ',2x,f4.3)
        end if
        nshot=nshot+1
        do 20 i=1,ipmax
           k(i)=0
           ke(i)=0.
 20     continue
        xshot=d
        goto 10
      end if
      k(iphase)=k(iphase)+1
      kg(iphase)=kg(iphase)+1
      ke(iphase)=ke(iphase)+e
      keg(iphase)=keg(iphase)+e
      goto 10
 999  continue
      kgt=0.
      kegt=0.
      do 60 i=1,ipmax
         kgt=kgt+kg(i)
60       kegt=kegt+keg(i)
      write(6,27)
      write(8,27)
      write(9,27)
 27   format(13x,
     &<ipmax>('-----'),'-----------')
      if(ipmax.gt.1) then
        write(6,25) (kg(i),i=1,ipmax),kgt
        write(8,25) (kg(i),i=1,ipmax),kgt
 25     format(' totals  ',4x,i5,2x,<ipmax-1>(i4,1x),' = ',2x,i5)
      else
        write(6,26) (kg(i),i=1,ipmax),kgt
        write(8,26) (kg(i),i=1,ipmax),kgt
 26     format(' totals  ',4x,i5,2x,' = ',2x,i5)
      end if
      do 889 i=1,ipmax
        if(kg(i).eq.0) then
          errg(i)=0.
        else
          errg(i)=keg(i)/kg(i)
        end if
 889  continue
      if(ipmax.gt.1) then
        write(9,125) (errg(i),i=1,ipmax),kegt/kgt
 125    format(' totals  ',6x,f4.3,1x,<ipmax-1>(f4.3,1x),' = ',2x,f4.3)
      else
        write(9,126) (errg(i),i=1,ipmax),kegt/kgt
 126    format(' totals  ',6x,f4.3,1x,' = ',2x,f4.3)
      end if
      stop
      end
