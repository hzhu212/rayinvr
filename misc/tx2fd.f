c
      character*10 rfile/'fd  .picks'/
      open(10, file='tx.in', status='old')
      open(11, file='txz.in', status='old')
c
      zrec=-0.17
c
      i=0
      xshotc=-999999.
c
100   read(10,*,end=999) x,t,u,ip
      if(ip.lt.0) go to 999
      if(ip.eq.0) then
        if(abs(x-xshotc).gt.001) then
          xshotc=x
          i=i+1
          id1=i/10
          id2=i-id1*10
          rfile(3:4)=char(id1+48)//char(id2+48)
c
          close(12)
          open(12, file=rfile, form='unformatted')
c
          read(11,*) z
c
          write(12) xshotc,0.,z,0.,0.,-1
        end if
      else
        if(ip.eq.1) write(12) x,0.,zrec,t,u,1
      end if
c
      go to 100
c
999   stop
      end
