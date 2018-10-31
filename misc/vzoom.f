c                 
c     version 1.3  Oct 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |             ***********  V Z O O M  ***********              |
c     |                                                              |
c     |              Change the xmin and/or xmax of a                |
c     |                       velocity model                         |
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
c        11 -- input:  original velocity model
c                 
c        12 -- output: zoomed-in velocity model
c                 
c
c     ----------------------------------------------------------------
c                 
c 
      program main
c
      include 'rayinvr.par'
c                
      real xm1(pncntr,ppcntr),zm1(pncntr,ppcntr),
     +     xm2(pncntr,ppcntr),zm2(pncntr,ppcntr),
     +     xvel1(player,ppvel,2),vf1(player,ppvel,2),
     +     xvel2(player,ppvel,2),vf2(player,ppvel,2)
      integer nzed1(pncntr),nvel1(player,2),
     +        nzed2(pncntr),nvel2(player,2),
     +        ivarz1(player,ppcntr),ivarv1(player,ppvel,2),
     +        ivarz2(player,ppcntr),ivarv2(player,ppvel,2)
c                 
      open(unit=11, file='v.in', status='old')
      open(unit=12, file='v.out')
c
      write(6,5)
5     format(/'Enter original xmax (km)')
      read(5,*) xmax
      write(6,15)
15    format(/'Enter new xmin and xmax (km)')
      read(5,*) xminn,xmaxn
      if(xminn.ge.xmax.or.xminn.ge.xmaxn) then
        write(6,25)
25      format(/
     +  '***  inappropriate value of xmax, xminn or xmaxn  ***'/)
        stop
      end if
c
      nrzmax=ppcntr/10
      nrvmax=ppvel/10
c
c     read in velocity model
c
      ncont=1
      do 170 icont=1,player+1
         nrz=1
         j1=1
         j2=10
11       if(nrz.gt.nrzmax) go to 211
         read(11,55,end=999) ilyr,(xm1(icont,j),j=j1,j2)
         read(11,55,end=999) icnt,(zm1(icont,j),j=j1,j2)
         read(11,65,end=99) (ivarz1(icont,j),j=j1,j2)
55       format(i2,1x,10f7.2)
65       format(3x,10i7)
         nrz=nrz+1
         if(icnt.ne.1) go to 211
         j1=j1+10
         j2=j2+10
         go to 11
211      nrv=1
         j1=1
         j2=10
21       if(nrv.gt.nrvmax) go to 311
         read(11,55,end=999) ilyr,(xvel1(icont,j,1),j=j1,j2)
         read(11,55,end=999) icnt,(vf1(icont,j,1),j=j1,j2)
         read(11,65,end=999) (ivarv1(icont,j,1),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 311
         j1=j1+10
         j2=j2+10
         go to 21
311      nrv=1
         j1=1
         j2=10
31       if(nrv.gt.nrvmax) go to 411
         read(11,55,end=999) ilyr,(xvel1(icont,j,2),j=j1,j2)
         read(11,55,end=999) icnt,(vf1(icont,j,2),j=j1,j2)
         read(11,65,end=999) (ivarv1(icont,j,2),j=j1,j2)
         nrv=nrv+1
         if(icnt.ne.1) go to 411
         j1=j1+10
         j2=j2+10
         go to 31
411      ncont=ncont+1
170   continue    
c               
99    nlayer=ncont-1
c
      do 171 i=1,ncont 
         nzed1(i)=1
171   continue
      do 172 i=1,nlayer
         nvel1(i,1)=1
         nvel1(i,2)=1
172   continue
c
      do 180 i=1,ncont
         do 190 j=1,ppcntr 
            if(abs(xm1(i,j)-xmax).lt..0001) go to 180
            nzed1(i)=nzed1(i)+1
190      continue 
180   continue    
c
      do 210 i=1,nlayer  
         do 220 j=1,ppvel  
            if(abs(xvel1(i,j,1)-xmax).lt..0001) go to 210
            nvel1(i,1)=nvel1(i,1)+1
220      continue 
210   continue    
c
      do 240 i=1,nlayer  
         do 250 j=1,ppvel  
            if(abs(xvel1(i,j,2)-xmax).lt..0001) go to 240
            nvel1(i,2)=nvel1(i,2)+1
250      continue 
240   continue    
c
c     change the xmin values of layer boundaries
c
      do 1010 i=1,ncont
         if(nzed1(i).le.1) go to 1010
         if(xminn.lt.xm1(i,1)) then
           if(nzed1(i).ge.ppcntr) then
             write(6,35) 
35           format(/'***  xminn<xmin and nzed=ppcntr  ***'/)  
             stop
           end if
           xm2(i,1)=xminn
           zm2(i,1)=zm1(i,1)
           ivarz2(i,1)=ivarz1(i,1)
           do 1020 j=1,nzed1(i)
              xm2(i,j+1)=xm1(i,j)
              zm2(i,j+1)=zm1(i,j)
              ivarz2(i,j+1)=ivarz1(i,j)
1020       continue
           nzed2(i)=nzed1(i)+1
           go to 1010
         end if
         if(xminn.gt.xm1(i,1)) then
           do 1030 j=2,nzed1(i)
              if(xminn.le.xm1(i,j)) then
                xm2(i,1)=xminn
                zm2(i,1)=((zm1(i,j)-zm1(i,j-1))/(xm1(i,j)-xm1(i,j-1)))*
     +                   (xminn-xm1(i,j-1))+zm1(i,j-1)
                ivarz2(i,1)=ivarz1(i,j-1)
                kc=1
                do 1040 k=j,nzed1(i)
                   kc=kc+1
                   xm2(i,kc)=xm1(i,k)
                   zm2(i,kc)=zm1(i,k)
                   ivarz2(i,kc)=ivarz1(i,k)
1040            continue
                nzed2(i)=kc
                go to 1010
              end if
1030       continue
         end if
1010  continue
c
c     change the xmax values of layer boundaries
c
      do 2010 i=1,ncont
         if(nzed1(i).le.1) then
           nzed2(i)=1
           xm2(i,1)=xmaxn
           zm2(i,1)=zm1(i,1)
           ivarz2(i,1)=ivarz1(i,1)
           go to 2010
         end if
         if(xmaxn.gt.xm2(i,nzed2(i))) then
           if(nzed2(i).ge.ppcntr) then
             write(6,45) 
45           format(/'***  xmaxn>xmax and nzed=ppcntr  ***'/)  
             stop
           end if
           xm2(i,nzed2(i)+1)=xmaxn
           zm2(i,nzed2(i)+1)=zm2(i,nzed2(i))
           ivarz2(i,nzed2(i)+1)=ivarz2(i,nzed2(i))
           nzed2(i)=nzed2(i)+1
           go to 2010
         end if
         if(xmaxn.lt.xm2(i,nzed2(i))) then
           do 2030 j=nzed2(i)-1,1,-1
              if(xmaxn.ge.xm2(i,j)) then
                zm2(i,j+1)=((zm2(i,j+1)-zm2(i,j))/(xm2(i,j+1)-
     +          xm2(i,j)))*(xmaxn-xm2(i,j))+zm2(i,j)
                xm2(i,j+1)=xmaxn
                nzed2(i)=j+1
                go to 2010
              end if
2030       continue
         end if
2010  continue
c
c     change the xmin values of upper layer velocities
c
      do 3010 i=1,nlayer
         if(nvel1(i,1).le.1) go to 3010
         if(xminn.lt.xvel1(i,1,1)) then
           if(nvel1(i,1).ge.ppvel) then
             write(6,155) 
155          format(/'***  xminn<xmin and nvel1=ppvel  ***'/)  
             stop
           end if
           xvel2(i,1,1)=xminn
           vf2(i,1,1)=vf1(i,1,1)
           ivarv2(i,1,1)=ivarv1(i,1,1)
           do 3020 j=1,nvel1(i,1)
              xvel2(i,j+1,1)=xvel1(i,j,1)
              vf2(i,j+1,1)=vf1(i,j,1)
              ivarv2(i,j+1,1)=ivarv1(i,j,1)
3020       continue
           nvel2(i,1)=nvel1(i,1)+1
           go to 3010
         end if
         if(xminn.gt.xvel1(i,1,1)) then
           do 3030 j=2,nvel1(i,1)
              if(xminn.le.xvel1(i,j,1)) then
                xvel2(i,1,1)=xminn
                vf2(i,1,1)=((vf1(i,j,1)-vf1(i,j-1,1))/(xvel1(i,j,1)-
     +          xvel1(i,j-1,1)))*(xminn-xvel1(i,j-1,1))+vf1(i,j-1,1)
                ivarv2(i,1,1)=ivarv1(i,j-1,1)
                kc=1
                do 3040 k=j,nvel1(i,1)
                   kc=kc+1
                   xvel2(i,kc,1)=xvel1(i,k,1)
                   vf2(i,kc,1)=vf1(i,k,1)
                   ivarv2(i,kc,1)=ivarv1(i,k,1)
3040            continue
                nvel2(i,1)=kc
                go to 3010
              end if
3030       continue
         end if
3010  continue
c
c     change the xmax values of upper layer velocities
c
      do 4010 i=1,ncont
         if(nvel1(i,1).le.1) then
           nvel2(i,1)=1
           xvel2(i,1,1)=xmaxn
           vf2(i,1,1)=vf1(i,1,1)
           ivarv2(i,1,1)=ivarv1(i,1,1)
           go to 4010
         end if
         if(xmaxn.gt.xvel2(i,nvel2(i,1),1)) then
           if(nvel2(i,1).ge.ppvel) then
             write(6,165) 
165          format(/'***  xmaxn>xmax and nvel1=ppvel  ***'/)  
             stop
           end if
           xvel2(i,nvel2(i,1)+1,1)=xmaxn
           vf2(i,nvel2(i,1)+1,1)=vf2(i,nvel2(i,1),1)
           ivarv2(i,nvel2(i,1)+1,1)=ivarv2(i,nvel2(i,1),1)
           nvel2(i,1)=nvel2(i,1)+1
           go to 4010
         end if
         if(xmaxn.lt.xvel2(i,nvel2(i,1),1)) then
           do 4030 j=nvel2(i,1)-1,1,-1
              if(xmaxn.ge.xvel2(i,j,1)) then
                vf2(i,j+1,1)=((vf2(i,j+1,1)-vf2(i,j,1))/(xvel2(i,j+1,1)-
     +          xvel2(i,j,1)))*(xmaxn-xvel2(i,j,1))+vf2(i,j,1)
                xvel2(i,j+1,1)=xmaxn
                nvel2(i,1)=j+1
                go to 4010
              end if
4030       continue
         end if
4010  continue
c
c     change the xmin values of lower layer velocities
c
      do 5010 i=1,nlayer
         if(nvel1(i,2).le.1) go to 5010
         if(xminn.lt.xvel1(i,1,2)) then
           if(nvel1(i,2).ge.ppvel) then
             write(6,75) 
75           format(/'***  xminn<xmin and nvel2=ppvel  ***'/)  
             stop
           end if
           xvel2(i,1,2)=xminn
           vf2(i,1,2)=vf1(i,1,2)
           ivarv2(i,2,2)=ivarv1(i,2,2)
           do 5020 j=1,nvel1(i,2)
              xvel2(i,j+1,2)=xvel1(i,j,2)
              vf2(i,j+1,2)=vf1(i,j,2)
              ivarv2(i,j+1,2)=ivarv1(i,j,2)
5020       continue
           nvel2(i,2)=nvel1(i,2)+1
           go to 5010
         end if
         if(xminn.gt.xvel1(i,1,2)) then
           do 5030 j=2,nvel1(i,2)
              if(xminn.le.xvel1(i,j,2)) then
                xvel2(i,1,2)=xminn
                vf2(i,1,2)=((vf1(i,j,2)-vf1(i,j-1,2))/(xvel1(i,j,2)-
     +          xvel1(i,j-1,2)))*(xminn-xvel1(i,j-1,2))+vf1(i,j-1,2)
                ivarv2(i,2,2)=ivarv1(i,j-1,2)
                kc=1
                do 5040 k=j,nvel1(i,2)
                   kc=kc+1
                   xvel2(i,kc,2)=xvel1(i,k,2)
                   vf2(i,kc,2)=vf1(i,k,2)
                   ivarv2(i,kc,2)=ivarv1(i,k,2)
5040            continue
                nvel2(i,2)=kc
                go to 5010
              end if
5030       continue
         end if
5010  continue
c
c     change the xmax values of lower layer velocities
c
      do 6010 i=1,ncont
         if(nvel1(i,2).le.1) then
           nvel2(i,2)=1
           xvel2(i,1,2)=xmaxn
           vf2(i,1,2)=vf1(i,1,2)
           ivarv2(i,1,2)=ivarv1(i,1,2)
           go to 6010
         end if
         if(xmaxn.gt.xvel2(i,nvel2(i,2),2)) then
           if(nvel2(i,2).ge.ppvel) then
             write(6,85) 
85           format(/'***  xmaxn>xmax and nvel2=ppvel  ***'/)  
             stop
           end if
           xvel2(i,nvel2(i,2)+1,2)=xmaxn
           vf2(i,nvel2(i,2)+1,2)=vf2(i,nvel2(i,2),2)
           ivarv2(i,nvel2(i,2)+1,2)=ivarv2(i,nvel2(i,2),2)
           nvel2(i,2)=nvel2(i,2)+1
           go to 6010
         end if
         if(xmaxn.lt.xvel2(i,nvel2(i,2),2)) then
           do 6030 j=nvel2(i,2)-1,1,-1
              if(xmaxn.ge.xvel2(i,j,2)) then
                vf2(i,j+1,2)=((vf2(i,j+1,2)-vf2(i,j,2))/(xvel2(i,j+1,2)-
     +          xvel2(i,j,2)))*(xmaxn-xvel2(i,j,2))+vf2(i,j,2)
                xvel2(i,j+1,2)=xmaxn
                nvel2(i,2)=j+1
                go to 6010
              end if
6030       continue
         end if
6010  continue
c
c     write out the new velocity model
c
      do 570 i=1,nlayer
         nstart=1  
  590    j1=nstart
         j2=j1+9 
         if(j2.gt.nzed2(i)) j2=nzed2(i)
         if(j2.lt.nzed2(i)) then
           icnt=1
         else 
           icnt=0
         end if
         write(12,55) i,(xm2(i,j),j=j1,j2)
         write(12,55) icnt,(zm2(i,j),j=j1,j2)
         write(12,65) (ivarz2(i,j),j=j1,j2)
         if(j2.eq.nzed2(i)) go to 600 
         nstart=j2+1
         go to 590
600      nstart=1  
620      j1=nstart
         j2=j1+9 
         if(j2.gt.nvel2(i,1)) j2=nvel2(i,1)
         if(j2.lt.nvel2(i,1)) then
           icnt=1
         else 
           icnt=0
         end if
         write(12,55) i,(xvel2(i,j,1),j=j1,j2)
         write(12,55) icnt,(vf2(i,j,1),j=j1,j2)
         write(12,65) (ivarv2(i,j,1),j=j1,j2)
         if(j2.eq.nvel2(i,1)) go to 630 
         nstart=j2+1
         go to 620
630      nstart=1  
650      j1=nstart
         j2=j1+9 
         if(j2.gt.nvel2(i,2)) j2=nvel2(i,2)
         if(j2.lt.nvel2(i,2)) then
           icnt=1
         else 
           icnt=0
         end if
         write(12,55) i,(xvel2(i,j,2),j=j1,j2)
         write(12,55) icnt,(vf2(i,j,2),j=j1,j2)
         write(12,65) (ivarv2(i,j,2),j=j1,j2)
         if(j2.eq.nvel2(i,2)) go to 570 
         nstart=j2+1
         go to 650
570   continue

      write(12,55) ncont,(xm2(ncont,j),j=1,nzed2(ncont))
      write(12,55) 0,(zm2(ncont,j),j=1,nzed2(ncont))
c
      stop
c
999   write(6,95)
95    format(/'***  premature end of file  ***'/)
      stop
c
      end
