c
c     version 1.3  Aug 1992
c
c     Calmod routine for TRAMP
c                 
c     ----------------------------------------------------------------
c                 
      subroutine calmod(ncont,pois,poisb,poisl,poisbl,
     +                  qp,qs,qb,ql,qpbl,qsbl,iamp,iden,iflagm)
c                 
c     calculate model parameters now for use later in program
c                 
c 
      include 'tramp.par'
      real xa(2*(ppcntr+ppvel)),qp(player),qs(player),
     +     pois(player),poisbl(papois),
     +     qpbl(paqpqs),qsbl(paqpqs),vmod(121),zsmth(pnsmth)
      integer poisb(papois),poisl(papois),qb(paqpqs),ql(paqpqs)
      include 'tramp.com'                
      iflagm=0
c
      do 10 i=1,ncont
         do 20 j=1,ppcntr
            if(abs(xm(i,j)-xmax).lt..0001) go to 30
            nzed(i)=nzed(i)+1
20       continue
30       if(nzed(i).gt.1) then
           do 40 j=1,nzed(i)-1
              if(xm(i,j).ge.xm(i,j+1)) go to 999
40         continue
           if(abs(xm(i,1)-xmin).gt..001.or.
     +        abs(xm(i,nzed(i))-xmax).gt..001) go to 999
         else
            xm(i,1)=xmax
         end if
10    continue
      do 11 i=1,nlayer
         do 21 j=1,ppvel
            if(abs(xvel(i,j,1)-xmax).lt..0001) go to 31
            nvel(i,1)=nvel(i,1)+1
21       continue
31       if(nvel(i,1).gt.1) then
           do 41 j=1,nvel(i,1)-1
              if(xvel(i,j,1).ge.xvel(i,j+1,1)) go to 999
41         continue
           if(abs(xvel(i,1,1)-xmin).gt..001.or.abs(xvel(i,nvel(i,1),1)-
     +        xmax).gt..001) go to 999
         else
           if(vf(i,1,1).gt.0.) then
             xvel(i,1,1)=xmax
           else
             if(i.eq.1) then
               go to 999
             else
               nvel(i,1)=0
             end if
           end if
         end if
11    continue
      do 12 i=1,nlayer
         do 22 j=1,ppvel
            if(abs(xvel(i,j,2)-xmax).lt..0001) go to 32
            nvel(i,2)=nvel(i,2)+1
22       continue
32       if(nvel(i,2).gt.1) then
           do 42 j=1,nvel(i,2)-1
              if(xvel(i,j,2).ge.xvel(i,j+1,2)) go to 999
42         continue
           if(abs(xvel(i,1,2)-xmin).gt..001.or.abs(xvel(i,nvel(i,2),2)-
     +        xmax).gt..001) go to 999
         else
           if(vf(i,1,2).gt.0.) then
             xvel(i,1,2)=xmax
           else
             nvel(i,2)=0
           end if
         end if
12    continue
c
      do 50 i=1,nlayer
         xa(1)=xmin
         xa(2)=xmax
         ib=2
         ih=ib
         do 60 j=1,nzed(i)
            do 61 k=1,ih
               if(abs(xm(i,j)-xa(k)).lt..005) go to 60
61          continue
            ib=ib+1
            xa(ib)=xm(i,j)
60       continue
         ih=ib
         do 70 j=1,nzed(i+1)
            do 80 k=1,ih
               if(abs(xm(i+1,j)-xa(k)).lt..005) go to 70
80          continue
            ib=ib+1
            xa(ib)=xm(i+1,j)
70       continue
         ih=ib
         if(nvel(i,1).gt.0) then
           il=i
           is=1
         else
           if(nvel(i-1,2).gt.0) then
             il=i-1
             is=2
           else
             il=i-1
             is=1
           end if
         end if
         do 71 j=1,nvel(il,is)
            do 81 k=1,ih
               if(abs(xvel(il,j,is)-xa(k)).lt..005) go to 71
81          continue
            ib=ib+1
            xa(ib)=xvel(il,j,is)
71       continue
         if(nvel(i,2).gt.0) then
           ih=ib
           do 72 j=1,nvel(i,2)
              do 82 k=1,ih
                 if(abs(xvel(i,j,2)-xa(k)).lt..005) go to 72
82            continue
              ib=ib+1
              xa(ib)=xvel(i,j,2)
72         continue
         end if
c
         if(ib.gt.(ptrap+1)) then
           write(6,5) i
5          format(/'***  maximum number of blocks in layer ',
     +       i2,' exceeded  ***'/)
           iflagm=1
           return
         end if
c
         call sort(xa,ib)
c
         nblk(i)=ib-1
         do 90 j=1,nblk(i)
            xbnd(i,j,1)=xa(j)
            xbnd(i,j,2)=xa(j+1)
90       continue
c
50    continue
c
c     calculate slopes and intercepts of each block boundary
c
      do 100 i=1,nlayer
         do 110 j=1,nblk(i)
            xbndc=xbnd(i,j,1)+.001
            if(nzed(i).gt.1) then
              do 120 k=1,nzed(i)-1
                 if(xbndc.ge.xm(i,k).and.xbndc.le.xm(i,k+1)) then
                   dx=xm(i,k+1)-xm(i,k)
                   c1=(xm(i,k+1)-xbnd(i,j,1))/dx
                   c2=(xbnd(i,j,1)-xm(i,k))/dx
                   z1=c1*zm(i,k)+c2*zm(i,k+1)
                   c1=(xm(i,k+1)-xbnd(i,j,2))/dx
                   c2=(xbnd(i,j,2)-xm(i,k))/dx
                   z2=c1*zm(i,k)+c2*zm(i,k+1)
                   go to 130
                 end if
  120         continue
            else
              z1=zm(i,1)
              z2=zm(i,1) 
            end if
130         s(i,j,1)=(z2-z1)/(xbnd(i,j,2)-xbnd(i,j,1))
            b(i,j,1)=z1-s(i,j,1)*xbnd(i,j,1)
            if(nzed(i+1).gt.1) then
              do 140 k=1,nzed(i+1)-1
                 if(xbndc.ge.xm(i+1,k).and.xbndc.le.
     +             xm(i+1,k+1)) then
                   dx=xm(i+1,k+1)-xm(i+1,k)
                   c1=(xm(i+1,k+1)-xbnd(i,j,1))/dx
                   c2=(xbnd(i,j,1)-xm(i+1,k))/dx
                   z3=c1*zm(i+1,k)+c2*zm(i+1,k+1)
                   c1=(xm(i+1,k+1)-xbnd(i,j,2))/dx
                   c2=(xbnd(i,j,2)-xm(i+1,k))/dx
                   z4=c1*zm(i+1,k)+c2*zm(i+1,k+1)
                   go to 150
                 end if
140           continue
            else
              z3=zm(i+1,1)
              z4=zm(i+1,1) 
            end if
150         s(i,j,2)=(z4-z3)/(xbnd(i,j,2)-xbnd(i,j,1))
            b(i,j,2)=z3-s(i,j,2)*xbnd(i,j,1)
c                 
c           check for layer pinchouts 
c                 
            ivg(i,j)=1
            if(abs(z3-z1).lt..0005) ivg(i,j)=2 
            if(abs(z4-z2).lt..0005) ivg(i,j)=3 
            if(abs(z3-z1).lt..0005.and.abs(z4-z2).lt..0005) ivg(i,j)=-1
c
110      continue 
100   continue    
c
c     assign velocities to each model block
c
      do 160 i=1,nlayer 
c
         if(nvel(i,1).eq.0) then
           do 161 j=i-1,1,-1
              if(nvel(j,2).gt.0) then
                ig=j
                jg=2
                n1g=nvel(j,2)
                go to 162
              end if
              if(nvel(j,1).gt.0) then
                ig=j
                jg=1
                n1g=nvel(j,1)
                go to 162
              end if
161        continue
         else
           ig=i
           jg=1
           n1g=nvel(i,1)
         end if
c
162      if(n1g.gt.1.and.nvel(i,2).gt.1) ivcase=1
         if(n1g.gt.1.and.nvel(i,2).eq.1) ivcase=2
         if(n1g.eq.1.and.nvel(i,2).gt.1) ivcase=3
         if(n1g.eq.1.and.nvel(i,2).eq.1) ivcase=4
         if(n1g.gt.1.and.nvel(i,2).eq.0) ivcase=5
         if(n1g.eq.1.and.nvel(i,2).eq.0) ivcase=6
c
         do 170 j=1,nblk(i)
c
            if(ivg(i,j).eq.-1) go to 170
c
            xbndcl=xbnd(i,j,1)+.001
            xbndcr=xbnd(i,j,2)-.001
c
            go to (1001,1002,1003,1004,1005,1006), ivcase
c
1001        do 180 k=1,n1g-1
               if(xbndcl.ge.xvel(ig,k,jg).and.xbndcl.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(ig,k,jg)
                 vm(i,j,1)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
                 go to 1811
               end if
180         continue
c
1811        do 1812 k=1,n1g-1
               if(xbndcr.ge.xvel(ig,k,jg).and.xbndcr.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(ig,k,jg)
                 vm(i,j,2)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,2)=vm(i,j,2)+.001
                 go to 181
               end if
1812        continue
c
181         do 182 k=1,nvel(i,2)-1
               if(xbndcl.ge.xvel(i,k,2).and.xbndcl.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(i,k,2)
                 if(ivg(i,j).ne.2) then
                   vm(i,j,3)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,3)=vm(i,j,1)
                 end if 
                 go to 187
               end if
182         continue
c
187         do 1822 k=1,nvel(i,2)-1
               if(xbndcr.ge.xvel(i,k,2).and.xbndcr.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(i,k,2)
                 if(ivg(i,j).ne.3) then
                   vm(i,j,4)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,4)=vm(i,j,2) 
                 end if
                 go to 171
               end if
1822        continue
c    
1002        do 183 k=1,n1g-1
               if(xbndcl.ge.xvel(ig,k,jg).and.xbndcl.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(ig,k,jg)
                 vm(i,j,1)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
                 go to 1833
               end if
183         continue
c
1833        do 1832 k=1,n1g-1
               if(xbndcr.ge.xvel(ig,k,jg).and.xbndcr.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(ig,k,jg)
                 vm(i,j,2)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,2)=vm(i,j,2)+.001
                 go to 184
               end if
1832        continue
c
184         vm(i,j,3)=vf(i,1,2)
            if(ivg(i,j).eq.2) vm(i,j,3)=vm(i,j,1)
            vm(i,j,4)=vf(i,1,2)
            if(ivg(i,j).eq.3) vm(i,j,4)=vm(i,j,2)
            go to 171
c
1003        vm(i,j,1)=vf(ig,1,jg)
            vm(i,j,2)=vf(ig,1,jg)
            if(ig.ne.i) then
              vm(i,j,1)=vm(i,j,1)+.001
              vm(i,j,2)=vm(i,j,2)+.001
            end if
c
            do 185 k=1,nvel(i,2)-1
               if(xbndcl.ge.xvel(i,k,2).and.xbndcl.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(i,k,2)
                 if(ivg(i,j).ne.2) then
                   vm(i,j,3)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,3)=vm(i,j,1)
                 end if
                 go to 188
               end if
185         continue
c
188         do 1851 k=1,nvel(i,2)-1
               if(xbndcr.ge.xvel(i,k,2).and.xbndcr.le.xvel(i,k+1,2)) 
     +           then
                 dxx=xvel(i,k+1,2)-xvel(i,k,2)   
                 c1=xvel(i,k+1,2)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(i,k,2)
                 if(ivg(i,j).ne.3) then
                   vm(i,j,4)=(c1*vf(i,k,2)+c2*vf(i,k+1,2))/dxx
                 else
                   vm(i,j,4)=vm(i,j,2)
                 end if
                 go to 171
               end if
1851        continue
c
1004        vm(i,j,1)=vf(ig,1,jg)
            vm(i,j,2)=vf(ig,1,jg) 
            if(ig.ne.i) then
              vm(i,j,1)=vm(i,j,1)+.001
              vm(i,j,2)=vm(i,j,2)+.001
            end if
c
            vm(i,j,3)=vf(i,1,2)
            if(ivg(i,j).eq.2) vm(i,j,3)=vm(i,j,1)
            vm(i,j,4)=vf(i,1,2)
            if(ivg(i,j).eq.3) vm(i,j,4)=vm(i,j,2)
            go to 171
c
1005        do 186 k=1,n1g-1
               if(xbndcl.ge.xvel(ig,k,jg).and.xbndcl.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,1)
                 c2=xbnd(i,j,1)-xvel(ig,k,jg)
                 vm(i,j,1)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
                 vm(i,j,3)=vm(i,j,1)
                 go to 1861
               end if
186         continue
c
1861        do 1862 k=1,n1g-1
               if(xbndcr.ge.xvel(ig,k,jg).and.xbndcr.le.xvel(ig,k+1,jg)) 
     +           then
                 dxx=xvel(ig,k+1,jg)-xvel(ig,k,jg)   
                 c1=xvel(ig,k+1,jg)-xbnd(i,j,2)
                 c2=xbnd(i,j,2)-xvel(ig,k,jg)
                 vm(i,j,2)=(c1*vf(ig,k,jg)+c2*vf(ig,k+1,jg))/dxx
                 if(ig.ne.i) vm(i,j,2)=vm(i,j,2)+.001
                 vm(i,j,4)=vm(i,j,2)                        
                 go to 171
               end if
1862        continue
c
1006        vm(i,j,1)=vf(ig,1,jg)
            if(ig.ne.i) vm(i,j,1)=vm(i,j,1)+.001
            vm(i,j,2)=vm(i,j,1) 
            vm(i,j,3)=vm(i,j,1)
            vm(i,j,4)=vm(i,j,1) 
c
c           calculate velocity coefficients
c
171         s1=s(i,j,1)
            s2=s(i,j,2)
            b1=b(i,j,1)
            b2=b(i,j,2)
            xb1=xbnd(i,j,1)
            xb2=xbnd(i,j,2)
            if(ivg(i,j).eq.2) then
              z3=s(i,j,2)*xb1+b(i,j,2)+.001
              z4=s(i,j,2)*xb2+b(i,j,2)
              s2=(z4-z3)/(xb2-xb1)
              b2=z3-s2*xb1
            end if
            if(ivg(i,j).eq.3) then
              z3=s(i,j,2)*xb1+b(i,j,2)
              z4=s(i,j,2)*xb2+b(i,j,2)+.001
              s2=(z4-z3)/(xb2-xb1)
              b2=z3-s2*xb1
            end if
            v1=vm(i,j,1)
            v2=vm(i,j,2)
            v3=vm(i,j,3)
            v4=vm(i,j,4)
c
            c(i,j,1)=s2*(xb2*v1-xb1*v2)+b2*(v2-v1)-
     +               s1*(xb2*v3-xb1*v4)-b1*(v4-v3)       
            c(i,j,2)=s2*(v2-v1)-s1*(v4-v3)
            c(i,j,3)=-xb2*v1+xb1*v2+xb2*v3-xb1*v4
            c(i,j,4)=-v2+v1+v4-v3
            c(i,j,5)=b2*(xb2*v1-xb1*v2)-b1*(xb2*v3-xb1*v4)
            c(i,j,6)=(s2-s1)*(xb2-xb1)
            c(i,j,7)=(b2-b1)*(xb2-xb1)
            c(i,j,8)=2.*c(i,j,2)*c(i,j,7)
            c(i,j,9)=c(i,j,2)*c(i,j,6)
            c(i,j,10)=c(i,j,4)*c(i,j,7)-c(i,j,3)*c(i,j,6)
            c(i,j,11)=c(i,j,1)*c(i,j,7)-c(i,j,5)*c(i,j,6)
c
            if(ivg(i,j).eq.-1) then
              vm(i,j,1)=0.
              vm(i,j,2)=0.
              vm(i,j,3)=0.
              vm(i,j,4)=0.
              do 172 k=1,11
                 c(i,j,1)=0.
172           continue
            end if
            if(abs(vm(i,j,1)-vm(i,j,2)).le..001.and.abs(vm(i,j,2)-
     +        vm(i,j,3)).le..001.and.abs(vm(i,j,3)-vm(i,j,4)).le..001. 
     +        and.vm(i,j,1).ne.0.) ivg(i,j)=0
c
170      continue
160   continue    
c
c     assign values to array vsvp
c                 
      if(pois(1).lt.-10.) then 
        do 190 i=1,nlayer 
           do 200 j=1,nblk(i)
              vsvp(i,j)=0.57735
200        continue
190     continue  
      else        
        if(nlayer.gt.1) then
          if(pois(2).lt.-10.) then
            do 210 j=1,nblk(1)
               vsvp(1,j)=sqrt((1.-2.*pois(1))/(2.*(1.-pois(1))))
210         continue
            do 220 i=2,nlayer
               do 230 j=1,nblk(i)
                  vsvp(i,j)=vsvp(1,1)
230            continue
220         continue
          else    
            do 240 i=1,nlayer
               if(pois(i).lt.-10.) then
                 do 250 j=1,nblk(i)
                    vsvp(i,j)=0.57735
250              continue
               else
                 do 260 j=1,nblk(i)
                    vsvp(i,j)=sqrt((1.-2.*pois(i))/(2.*(1.-pois(i))))
260              continue 
               end if
240         continue 
          end if  
        end if    
      end if      
c                 
c     calculate velocity ratios for specific model blocks specified
c     through the arrays poisbl, poisl and poisb
c                 
      i=1         
270   if(poisbl(i).lt.-10.) go to 280
      vsvp(poisl(i),poisb(i))=
     +  sqrt((1.-2.*poisbl(i))/(2.*(1.-poisbl(i))))
      i=i+1       
      if(i.le.papois) go to 270
c                 
c     assign switch iq and array q
c                 
280   if(iamp.eq.0) go to 400
      if(qp(1).gt.0..or.qs(1).gt.0..or.qpbl(1).gt.0..or.qsbl(1).gt.0.)
     +  then      
        iq=1      
        if(qp(1).le.0.) then 
          do 290 i=1,nlayer
             qp(i)=1.e10
290       continue 
        else      
          if(nlayer.gt.1) then
            if(qp(2).le.0.) then
              do 300 i=2,nlayer
                 qp(i)=qp(1)
300           continue
            else  
              do 310 i=2,nlayer
                 if(qp(i).le.0.) qp(i)=1.e10
310           continue
            end if
          end if  
        end if    
        if(qs(1).le.0.) then 
          do 320 i=1,nlayer
             qs(i)=1.e10
320       continue 
        else      
          if(nlayer.gt.1) then
            if(qs(2).le.0.) then
              do 330 i=2,nlayer
                 qs(i)=qs(1)
330           continue
            else  
              do 340 i=2,nlayer
                 if(qs(i).le.0.) qs(i)=1.e10
340           continue
            end if
          end if  
        end if    
      end if      
c                 
      do 350 i=1,nlayer
         do 360 j=1,nblk(i)
            q(i,j,1)=qp(i)
            q(i,j,2)=qs(i)
360      continue 
350   continue    
c                 
c     assign q values to specific model blocks specified through
c     the arrays qpbl, qsbl, ql and qb
c                 
      i=1         
370   if(qpbl(i).le.0.) go to 380
      q(ql(i),qb(i),1)=qpbl(i)
      i=i+1       
      if(i.le.paqpqs) go to 370
c                 
380   j=1         
390   if(qsbl(j).le.0.) go to 400
      q(ql(i),qb(i),2)=qsbl(j)
      i=i+1       
      j=j+1       
      if(j.le.paqpqs) go to 390
c                 
400   if(ibsmth.gt.0) then
        xsinc=(xmax-xmin)/float(npbnd-1)
        do 600 i=1,nlayer+1
           if(i.lt.(nlayer+1)) then
             il=i 
             ib=1 
           else   
             il=i-1
             ib=2 
           end if 
           iblk=1 
           do 610 j=1,npbnd
              x=xmin+float(j-1)*xsinc
              if(x.lt.xmin) x=xmin+.001
              if(x.gt.xmax) x=xmax-.001
620           if(x.ge.xbnd(il,iblk,1).and.x.le.xbnd(il,iblk,2)) then
                cosmth(i,j)=s(il,iblk,ib)*x+b(il,iblk,ib)
                go to 610
              else
                iblk=iblk+1
                go to 620
              end if
610        continue
600     continue  
        if(nbsmth.gt.0) then
          do 630 i=1,nlayer+1
             do 640 j=1,npbnd
                zsmth(j)=cosmth(i,j)
640          continue 
             do 650 j=1,nbsmth 
                call smooth(zsmth,npbnd) 
650          continue
             do 660 j=1,npbnd
                cosmth(i,j)=zsmth(j)
660          continue 
630       continue
        end if    
      end if      
c                 
      if(idump.eq.1) then
        write(12,15) nlayer
15      format('***  velocity model:  ***'//'number of layers=',i2)
        do 410 i=1,nlayer
           if(iq.eq.0) then
             write(12,25) i,nblk(i)
25           format(/'layer#',i2,'  nblk=',i2,
     +' (ivg,x1,x2,z11,z12,z21,z22,s1,b1,s2,b2,vp1,vs1,vp2,vs2,
     +vp3,vs3,vp4,vs4,c1,c2,...,c11)')
           else   
             write(12,35) i,nblk(i)
35           format(/'layer#',i2,'  nblk=',i2,
     +' (ivg,x1,x2,z11,z12,z21,z22,s1,b1,s2,b2,vp1,vs1,vp2,vs2,
     +vp3,vs3,vp4,vs4,c1,c2,...,c11,qp,qs)')
           end if 
           write(12,36) (ivg(i,j),j=1,nblk(i))
           write(12,45) (xbnd(i,j,1),j=1,nblk(i))
           write(12,45) (xbnd(i,j,2),j=1,nblk(i))
           write(12,45) (s(i,j,1)*xbnd(i,j,1)+b(i,j,1),j=1,nblk(i))
           write(12,45) (s(i,j,1)*xbnd(i,j,2)+b(i,j,1),j=1,nblk(i))
           write(12,45) (s(i,j,2)*xbnd(i,j,1)+b(i,j,2),j=1,nblk(i))
           write(12,45) (s(i,j,2)*xbnd(i,j,2)+b(i,j,2),j=1,nblk(i))
           write(12,45) (s(i,j,1),j=1,nblk(i))
           write(12,45) (b(i,j,1),j=1,nblk(i))
           write(12,45) (s(i,j,2),j=1,nblk(i))
           write(12,45) (b(i,j,2),j=1,nblk(i))
           write(12,45) (vm(i,j,1),j=1,nblk(i))
           write(12,45) (vm(i,j,1)*vsvp(i,j),j=1,nblk(i)) 
           write(12,45) (vm(i,j,2),j=1,nblk(i))
           write(12,45) (vm(i,j,2)*vsvp(i,j),j=1,nblk(i))
           write(12,45) (vm(i,j,3),j=1,nblk(i))
           write(12,45) (vm(i,j,3)*vsvp(i,j),j=1,nblk(i)) 
           write(12,45) (vm(i,j,4),j=1,nblk(i))
           write(12,45) (vm(i,j,4)*vsvp(i,j),j=1,nblk(i))
           write(12,55) (c(i,j,1),j=1,nblk(i))
           write(12,55) (c(i,j,2),j=1,nblk(i))
           write(12,55) (c(i,j,3),j=1,nblk(i))
           write(12,55) (c(i,j,4),j=1,nblk(i))
           write(12,55) (c(i,j,5),j=1,nblk(i))
           write(12,55) (c(i,j,6),j=1,nblk(i))
           write(12,55) (c(i,j,7),j=1,nblk(i))
           write(12,55) (c(i,j,8),j=1,nblk(i))
           write(12,55) (c(i,j,9),j=1,nblk(i))
           write(12,55) (c(i,j,10),j=1,nblk(i))
           write(12,55) (c(i,j,11),j=1,nblk(i))
           if(iq.eq.1) then
             write(12,55) (q(i,j,1),j=1,nblk(i))
             write(12,55) (q(i,j,2),j=1,nblk(i))
           end if 
36         format(100i10)
45         format(100f10.4)
55         format(100e10.3)
410     continue  
        xmod=xmax-xmin 
        write(12,65) 
65      format(/'equivalent 1-dimensional velocity model:'/)
        do 520 i=1,nlayer
           z1sum=0.
           z2sum=0.
           vp1=0. 
           vp2=0. 
           vs1=0. 
           vs2=0. 
           vp1sum=0.
           vp2sum=0.
           vs1sum=0.
           vs2sum=0.
           xvmod=0.
           do 530 j=1,nblk(i)
              xblk=xbnd(i,j,2)-xbnd(i,j,1)
              z11=s(i,j,1)*xbnd(i,j,1)+b(i,j,1)
              z12=s(i,j,1)*xbnd(i,j,2)+b(i,j,1)
              z21=s(i,j,2)*xbnd(i,j,1)+b(i,j,2)
              z22=s(i,j,2)*xbnd(i,j,2)+b(i,j,2)
              z1sum=z1sum+xblk*(z11+z12)/2.
              z2sum=z2sum+xblk*(z21+z22)/2.
              if(vm(i,j,1).gt..001) then
                vp1sum=vp1sum+xblk*(vm(i,j,1)+vm(i,j,2))/2.
                vs1sum=vs1sum+xblk*(vm(i,j,1)+vm(i,j,2))*vsvp(i,j)/2.
                vp2sum=vp2sum+xblk*(vm(i,j,3)+vm(i,j,4))/2.
                vs2sum=vs2sum+xblk*(vm(i,j,3)+vm(i,j,4))*vsvp(i,j)/2.
                xvmod=xvmod+xblk
              end if
530        continue
           z1=z1sum/xmod
           z2=z2sum/xmod
           if(xvmod.gt..000001) then
             vp1=vp1sum/xvmod
             vp2=vp2sum/xvmod
             vs1=vs1sum/xvmod
             vs2=vs2sum/xvmod 
           end if 
           write(12,75) i,z1,z2,vp1,vp2,vs1,vs2
75         format('layer# ',i2,'   z1=',f7.2,'   z2=',f7.2,
     +            ' km'/9x,'  vp1=',f7.2,'  vp2=',f7.2,' km/s'/
     +                  9x,'  vs1=',f7.2,'  vs2=',f7.2,' km/s')
520     continue  
      end if      
      if(igrid.ne.0) then
        if(xgrid.le.0.) xgrid=(xmax-xmin)/100.
        if(zgrid.le.0.) zgrid=(zmax-zmin)/50.
        nx=nint((xmax-xmin-.002)/xgrid)+1
        xinc=(xmax-xmin-.002)/float(nx-1)
        nz=nint((zmax-zmin-.002)/zgrid)+1
        zinc=(zmax-zmin-.002)/float(nz-1)
        write(23,1760) nx,nz
1760    format(2i10)
        do 1710 i=1,nx
           xp=xmin+.001+float(i-1)*xinc
           do 1720 j=1,nz
              zp=zmin+.001+float(j-1)*zinc 
              call xzpt(xp,zp,la,ib,iflag)
              layer=la
              iblk=ib
              if(iabs(igrid).eq.1) then
                vmod(j)=vel(xp,zp) 
                if(igrid.eq.-1) vmod(j)=vmod(j)*vsvp(la,ib)
              else
                if(iabs(igrid).eq.2) then
                  t1=(c(layer,iblk,8)*xp+c(layer,iblk,9)*xp**2
     +               +c(layer,iblk,10)*zp+c(layer,iblk,11))/
     +               (c(layer,iblk,6)*xp+c(layer,iblk,7))**2
                  t2=(c(layer,iblk,3)+c(layer,iblk,4)*xp)/
     +               (c(layer,iblk,6)*xp+c(layer,iblk,7)) 
                  vmod(j)=sqrt(t1**2+t2**2)
                  if(igrid.eq.-2) vmod(j)=vmod(j)*vsvp(la,ib)
                else 
                  if(igrid.eq.3) vmod(j)=(2.*vsvp(la,ib)**2-1.)/ 
     +                                   (2.*(vsvp(la,ib)**2-1.))
                end if
              end if
1720       continue 
           write(23,1770) (vmod(j),j=1,nz)
c1770       format(121a4)
1770       format(10f10.3)
1710    continue  
      end if      
c                 
      if(iden.eq.1) then
        do 1880 i=1,nlayer 
           if(vm(i,1,1).gt.0.) then
             xl=xbnd(i,1,1)-(xmax-xmin)
             xri=xbnd(i,1,1)
             z11=s(i,1,1)*xri+b(i,1,1)
             z12=z11
             z21=s(i,1,2)*xri+b(i,1,2)
             z22=z21
             vave=(vm(i,1,1)+vm(i,1,2)+vm(i,1,3)+vm(i,1,4))/4.
             vsw=vave*vsvp(i,1)
             den=denc(1)+denc(2)*vave+denc(3)*vave**2+
     +           denc(4)*vave**3+denc(5)*vave**4
             if(den.lt.denmin) den=denmin
             if(vsw.le..001) den=1.0
             write(23,1815) 4,den,den,xl,z11,xri,z12,xri,z22,xl,z21
           end if 
           if(vm(i,nblk(i),1).gt.0.) then
             xl=xbnd(i,nblk(i),2)
             xri=xbnd(i,nblk(i),2)+(xmax-xmin)
             z11=s(i,nblk(i),1)*xl+b(i,nblk(i),1)
             z12=z11
             z21=s(i,nblk(i),2)*xl+b(i,nblk(i),2)
             z22=z21
             vave=(vm(i,nblk(i),1)+vm(i,nblk(i),2)+
     +             vm(i,nblk(i),3)+vm(i,nblk(i),4))/4.
             vsw=vave*vsvp(i,1)
             den=denc(1)+denc(2)*vave+denc(3)*vave**2+
     +           denc(4)*vave**3+denc(5)*vave**4
             if(den.lt.denmin) den=denmin
             if(vsw.le..001) den=1.0
             write(23,1815) 4,den,den,xl,z11,xri,z12,xri,z22,xl,z21
           end if 
1880    continue  
        do 1810 i=1,nlayer
           do 1820 j=1,nblk(i)
              if(vm(i,j,1).gt.0.) then
                xl=xbnd(i,j,1)
                xri=xbnd(i,j,2)
                z11=s(i,j,1)*xl+b(i,j,1)
                z12=s(i,j,1)*xri+b(i,j,1) 
                z21=s(i,j,2)*xl+b(i,j,2)
                z22=s(i,j,2)*xri+b(i,j,2) 
                vave=(vm(i,j,1)+vm(i,j,2)+vm(i,j,3)+vm(i,j,4))/4.
                vsw=vave*vsvp(i,j)
                if(z11.eq.z21) z21=z11+.001
                if(z12.eq.z22) z22=z22+.001
                den=denc(1)+denc(2)*vave+denc(3)*vave**2+
     +              denc(4)*vave**3+denc(5)*vave**4
                if(den.lt.denmin) den=denmin
                if(vsw.le..001) den=1.0
                write(23,1815) 4,den,den,xl,z11,xri,z12,xri,z22,xl,z21
1815            format(i5,2f7.4/10f8.3)
              end if
1820       continue 
1810    continue  
      end if      
c                 
      return      
c                 
999   write(6,900)
900   format(/'***  error in velocity model  ***'/)
      iflagm=1
      return
      end         
