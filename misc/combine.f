c                 
c     version 1.2  Mar 1992
c                 
c     ----------------------------------------------------------------
c     |                                                              |
c     |            *********  C O M B I N E  **********              |
c     |                                                              |
c     |            Combine the partial derivatives and               |
c     |            traveltime residuals from two runs                |
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
c        10 -- input:  partial derivatives/traveltime residuals
c
c        11 -- input:  partial derivatives/traveltime residuals
c
c        12 -- output:  combined partial derivatives/traveltime residuals
c                 
c                 
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'rayinvr.par'
c                
      real tres(prayi),apart1(prayi,pnvar),ures(prayi),
     +     tres1(prayi),apart2(prayi,pnvar),tres2(prayi),parorg(pnvar),
     +     ures1(prayi),ures2(prayi),parunc(pnvar)
      integer partyp(pnvar)
c                 
      open(unit=10, file='i1.in', status='old')
      open(unit=11, file='i2.in', status='old')
      open(unit=12, file='i.out')
c
c     read in matrix of partial derivatives and vector of traveltime
c     residuals 
c
      read(10,1)
1     format(' ')
      read(10,45) narin1,nvar1
45    format(2i5)
c
      read(11,1)
      read(11,45) narin2,nvar2
c
      if(narin1+narin2.gt.prayi) then
        write(6,95) 
95      format(/'***  too many arrivals  ***'/)
        stop 
      end if
c
      if(nvar1.ne.nvar2) then
        write(6,105) 
105     format(/'***  unequal number of parameters varied  ***'/)
        stop 
      end if
c
      nvar=nvar1
      narinv=narin1+narin2
c
      read(10,1)
      do 6 i=1,nvar
         read(10,5) partyp(i),parorg(i),parunc(i)
5        format(i5,2f15.5)
6     continue
      read(10,1)
c
      read(11,1)
      do 7 i=1,nvar
         read(11,5) partyp(i),parorg(i),parunc(i)
7     continue
      read(11,1)
c
      do 10 i=1,narin1
         read(10,15) (apart1(i,j),j=1,nvar)
15       format(5e12.5)
10    continue
c
      do 20 i=1,narin2
         read(11,15) (apart2(i,j),j=1,nvar)
20    continue
c
      read(10,1)
      read(10,15) (tres1(i),i=1,narin1)
c
      read(11,1)
      read(11,15) (tres2(i),i=1,narin2)
c
      read(10,1)
      read(10,15) (ures1(i),i=1,narin1)
c
      read(11,1)
      read(11,15) (ures2(i),i=1,narin2)
c
      do 30 i=1,narin1
         tres(i)=tres1(i)
         ures(i)=ures1(i)
30    continue
      do 40 i=1,narin2
         tres(i+narin1)=tres2(i)
         ures(i+narin1)=ures2(i)
40    continue
c
c     write combined file
c
      write(12,1)
      write(12,45) narinv,nvar
c
      write(12,1)
      do 8 i=1,nvar
         write(12,5) partyp(i),parorg(i),parunc(i)
8     continue
      write(12,1)
c
      do 50 i=1,narin1
         write(12,15) (apart1(i,j),j=1,nvar)
50    continue
c
      do 60 i=1,narin2
         write(12,15) (apart2(i,j),j=1,nvar)
60    continue
c
      write(12,1)
      write(12,15) (tres(i),i=1,narinv)
      write(12,1)
      write(12,15) (ures(i),i=1,narinv)
c                 
      read(10,1)
      read(10,25) narin1,trms1,chi1
25    format(28x,i8/28x,f8.3/26x,f10.3)
c
      read(11,1)
      read(11,25) narin2,trms2,chi2
c
      trms=sqrt((narin1*trms1**2+narin2*trms2**2)/narinv)
      chi=((narin1-1)*chi1+(narin2-1)*chi2)/(narinv-1)
c
      write(12,1)
      write(12,35) narinv,trms,chi
35    format('Number of data points used: ',i8/
     +       'RMS traveltime residual:    ',f8.3/
     +       'Normalized chi-squared:   ',f10.3)
      write(12,1)
c
      stop
      end
