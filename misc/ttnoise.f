c                 
c     version 1.2  Mar 1992
c                 
c     ----------------------------------------------------------------
c     |                                                              |
c     |            *********  T T N O I S E  **********              |
c     |                                                              |
c     |           Add Gaussian noise to traveltime picks             |
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
c        11 -- output:  noisy traveltime-distance pairs
c                 
c                 
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      open(unit=10, file='tx.in', status='old')
      open(unit=11, file='tx.out')
c
      write(6,5)
5     format(/'Enter a seed (positive integer)')
      read(5,*) iseed
      rnum=gasdev(-iseed)
c
100   read(10,35) xr,tr,ur,ir
35    format(3f10.3,i10)
      if(ir.le.0) then
        write(11,35) xr,tr,ur,ir
        if(ir.lt.0) go to 999
      else
        rnum=ur*gasdev(iseed)
        tr=tr+rnum
        write(11,35) xr,tr,ur,ir
      end if
      go to 100
c  
999   stop 
      end
c
c     -----------------------------------------------------------------
c
      function gasdev(idum)
c
c     returns a normally distributed number with zero mean and unit 
c     variance
c
      data iset/0/
      if(iset.eq.0) then
1       v1=2.*ran1(idum)-1.
        v2=2.*ran1(idum)-1.
        r=v1**2+v2**2
        if(r.ge.1.) go to 1
        fac=sqrt(-2.*log(r)/r)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      end if
      return
      end
c
c     -----------------------------------------------------------------
c
      function ran1(idum)
c
c     returns a uniform random number between 0.0 and 1.0. Set IDUM to
c     any negative value to initialize or reinitialize the sequence.
c
      dimension r(97)
      parameter (m1=259200, ia1=7141, ic1=54773, rm1=1./m1)
      parameter (m2=134456, ia2=8121, ic2=28411, rm2=1./m2)
      parameter (m3=243000, ia3=4561, ic3=51349)           
      data iff/0/
      if(idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,97
           ix1=mod(ia1*ix1+ic1,m1)
           ix2=mod(ia2*ix2+ic2,m2)
           r(j)=(float(ix1)+float(ix2)*rm2)*rm1 
11      continue
        idum=1
      end if
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j.gt.97.or.j.lt.1) pause
      ran1=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end
