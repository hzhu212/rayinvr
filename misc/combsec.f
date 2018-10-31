c
c     version 1.2  Mar 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ********   C O M B S E C   ********               |   
c     |                                                              |
c     |           Combine two "sect.out" files into one              |   
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c     Note: the second synthetic section file cannot contain sections
c           corresponding to shots not in the first file. Also, the order
c           of sections within each file must be the same for the those
c           sections whose shots are common to both files.
c
c     I/O units:
c
c        10 -- input:  first synthetic section file
c
c        11 -- input:  second synthetic section file
c
c        12 -- output: combined synthetic section file
c
c     ----------------------------------------------------------------
c
      open(unit=10, file='sect1.in', status='old')
      open(unit=11, file='sect2.in', status='old')
      open(unit=12, file='sect.out')
c
      iflag=2
c
1000  read(10,5,end=990) x1,num1
      if(iflag.eq.2) read(11,5,end=991) x2,num2
5     format(f10.3,i10)
c
1001  if(num1.eq.-1) then
        xshot1=x1
        write(12,5) xshot1,-1
        if(iflag.eq.3) go to 1000
        if(num2.eq.-1) then
          xshot2=x2
        else
          write(6,25)
25        format(/'***  files are not consistent  ***'/)
          stop
        end if  
        if(abs(xshot1-xshot2).lt..001) then
          iflag=2
        else
          iflag=1
        end if
        go to 1000
      end if
c
      if(iflag.eq.1.or.iflag.eq.3) then
        write(12,5) x1,num1
      end if
      if(iflag.eq.2) then
        if(abs(x1-x2).gt..001) then
          write(6,25)
          stop
        end if
        write(12,5) x1,num1+num2
      end if
      if(num1.gt.0) then
        do 10 i=1,num1
           read(10,15) s1,s2,s3
15         format(3e12.5)
           write(12,15) s1,s2,s3
10      continue
      end if
      if(iflag.eq.2) then
        if(num2.gt.0) then
          do 20 i=1,num2
             read(11,15) s1,s2,s3
             write(12,15) s1,s2,s3
20        continue
        end if
      end if 
      go to 1000
c
991   iflag=3
      go to 1001 
c
990   stop
      end
