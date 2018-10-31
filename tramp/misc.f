c
c     version 1.3  Aug 1992
c
c     Misc routines for TRAMP
c                 
c     ----------------------------------------------------------------
c                 
      subroutine ttime(is,xshot,npt,nr,a1,ifam,itt,iszero,iamp,iflag)
c                 
c     calculate travel time along a single ray using the equation:
c                 
c                       t=2*h/(v1+v2)
c                 
c     for the travel time between two points a distance h apart
c      
      include 'tramp.par'
      real vave(ppray)
      integer itt(1)
      include 'tramp.com'           
c 
      if(idump.eq.1) write(12,15) ifam,nr,npt,xr(npt),zr(npt),
     +  ar(npt,1)*pi18,ar(npt,2)*pi18,vr(npt,1),vr(npt,2),
     +  layer,iblk,id,iwave
15    format(i2,i3,i4,2f8.3,2f8.2,2f7.2,4i3)
      if(iflag.gt.0) then
        time=0.   
        qa=1.     
      end if      
      do 10 i=1,npt-1
         tr(i)=sqrt((xr(i+1)-xr(i))**2+(zr(i+1)-zr(i))**2)
         vave(i)=(vr(i,2)+vr(i+1,1))/2.
         time=time+tr(i)/vave(i)
10    continue    
      if(iamp.gt.0.and.iq.eq.1) then
        on2=-omega/2.
        do 20 i=1,npt-1
           qf=exp(on2*tr(i)/(vave(i)*qr(i)))
           qa=qa*qf
20      continue  
      end if      
                  
      if(iflag.eq.2) then
        n1=npt    
      else        
        a2=fid1*(90.-fid*ar(npt,1)*pi18)/fid
        if(iflag.ne.-2) then
          nptr=npt
        else      
          nptr=npt+n1
        end if    
        if(vred.eq.0.) then
          timer=time
        else      
          timer=time-abs(xr(npt)-xshot)/vred
        end if    
        rayid(ntt)=float(idray(1))+float(idray(2))/10.
        write(11,5) is,nr,a1,a2,xr(npt),zr(npt),timer,nptr,
     +              rayid(ntt)
5       format(2i4,2f9.3,f9.3,f8.2,f8.3,i6,f6.1)
        if(vr(npt,2).ne.0.) then
          itt(ifam)=itt(ifam)+1
          if(iszero.eq.0) then
            range(ntt)=xr(npt)
          else    
            range(ntt)=abs(xr(npt)-xshot)
          end if  
          tt(ntt)=timer
          xshtar(ntt)=xshot
          fidarr(ntt)=fid1
          p0(ntt)=ar(1,2)
          pr(ntt)=sin(ar(npt,1))/vr(npt,1)
          ntt=ntt+1
        end if    
      end if      
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine sort(x,npts)
c                 
c     sort the elements of array x in order of increasing size using
c     a bubble sort technique
c                 
      real x(1) 
      do 10 i=1,npts-1
         iflag=0  
         do 20 j=1,npts-1
            if(x(j).gt.x(j+1)) then
              iflag=1
              xh=x(j)
              x(j)=x(j+1)
              x(j+1)=xh
            end if
20       continue 
         if(iflag.eq.0) return
10     continue   
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine smooth(x,n) 
c                 
c     three point triangular smoothing filter
c                 
      real x(n) 
      m=n-1       
      a=0.77*x(1)+0.23*x(2) 
      b=0.77*x(n)+0.23*x(m) 
      xx=x(1)     
      xr=x(2)     
      do 10 i=2,m 
         xl=xx    
         xx=xr    
         xr=x(i+1) 
         x(i)=0.54*xx+0.23*(xl+xr) 
 10   continue    
      x(1)=a      
      x(n)=b      
      return      
      end         
