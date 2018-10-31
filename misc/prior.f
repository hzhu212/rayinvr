c
      include 'rayinvr.par'
c
      real apart(prayi,pnvar),tres(prayi),parorg(pnvar),
     +     tunc(prayi),parunc(pnvar)
      integer partyp(pnvar)
c
      real bndp(9),bndu(9),velp(9),velu(9)
c
      namelist /pripar/ bndp,bndu,velp,velu
c
      open(unit=10, file='prior.in', status='old')
      open(unit=11, file='i.out', status='old')
      open(unit=12, file='i.new')
c
      read(10,pripar)
c
c     read in matrix of partial derivatives
c     and vector of traveltime residuals
c
      read(11,1)
      write(12,1)
1     format(' ')
      read(11,115) narinv,nvar
      write(12,115) narinv,nvar
115   format(2i5)
      read(11,1)
      write(12,1)
c
      do 6 i=1,nvar
         read(11,5) partyp(i),parorg(i),parunc(i)
5        format(i5,2f15.5)
c
         if(partyp(i).eq.1) then
           do j=1,9
              if(parorg(i).lt.bndp(j)) then
                parw=bndu(j)
                go to 99
              end if
           enddo
         else
           do j=1,9
              if(parorg(i).lt.velp(j)) then
                parw=velu(j)
                go to 99
              end if
           enddo
         end if
         write(0,*) '*** screw up  ***'
         stop
c
99       write(12,5) partyp(i),parorg(i),parw
6     continue
c
      read(11,1)
      write(12,1)
      do 10 i=1,narinv
         read(11,15) (apart(i,j),j=1,nvar)
         write(12,15) (apart(i,j),j=1,nvar)
15       format(5e12.5)
10    continue
c
      read(11,1)
      write(12,1)
      read(11,15) (tres(i),i=1,narinv)
      write(12,15) (tres(i),i=1,narinv)
      read(11,1)
      write(12,1)
      read(11,15) (tunc(i),i=1,narinv)
      write(12,15) (tunc(i),i=1,narinv)
c
      stop
      end
