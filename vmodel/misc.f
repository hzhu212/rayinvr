c
c     version 1.2  Mar 1992
c
c     Misc routines for VMODEL
c                 
c     ----------------------------------------------------------------
c                 
      subroutine dxmin(ncont)
c
c     determine the value of xmin
c
      include 'vmodel.par'
      include 'vmodel.com'
c
      if(xmin.gt.-99998.) return
      do 10 i=1,ncont
         if(nzed(i).gt.1) then
           xmin=xm(i,1)
           return
         end if
10    continue
      do 20 i=1,ncont-1
         if(nvel(i,1).gt.1) then
           xmin=xvel(i,1,1)
           return
         end if
         if(nvel(i,2).gt.1) then
           xmin=xvel(i,1,2)
           return
         end if
20    continue
c
      write(6,5)
      write(13,5)
5     format(/'***  xmin not specified  ***'/)
c
      stop
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine smooth(x,n) 
c                 
c     three point triangular smoothing filter
c                 
      real x(n) 
c
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
c                 
c     ----------------------------------------------------------------
c                 
      subroutine sort(x,y,ix,npts)
c                 
c     sort the elements of array x in order of increasing size using
c     a bubble sort technique and sort the arrays y and ix with x
c                 
      real x(1),y(1)
      integer ix(1) 
c
      do 10 i=1,npts-1
         iflag=0  
         do 20 j=1,npts-1
            if(x(j).gt.x(j+1)) then
              iflag=1
              xh=x(j)
              x(j)=x(j+1)
              x(j+1)=xh
              yh=y(j)
              y(j)=y(j+1)
              y(j+1)=yh
              ixh=ix(j)
              ix(j)=ix(j+1)
              ix(j+1)=ixh
            end if
20       continue 
         if(iflag.eq.0) go to 999   
10    continue   
999   return      
      end         
c                 
c     ----------------------------------------------------------------
c
      subroutine smooth2(x,npts,n,nt)
c                 
c     smooth2 applies an n-point moving average operator to the array  
c     x, where n is an odd integer. The procedure is repeated nT times. 
c
      parameter(ixd=1000)
c
      real x(npts),y(ixd)
c        
      ioff=n/2
c
      do 50 k=1,nt
       do 10 i=1,npts
        ndiv=n
        i1=i-ioff
        i2=i+ioff
        if(i1.lt.1) then
         i1=1
        end if
        if(i2.gt.npts) then
         i2=npts
        end if
        y(i)=0.0d0
        ndiv=0
        do 20 j=i1,i2
         y(i)=y(i)+x(j)
         ndiv=ndiv+1
  20    continue
        if(ndiv.eq.0) go to 999
        y(i)=y(i)/dfloat(ndiv)
  10   continue
c
       do 30 i=1,npts
        x(i)=y(i)
  30   continue
c
  50  continue
c
 999  continue
c
      return
      end
