c                 
c                 
c     Misc routines for PLTSYN
c
c     ----------------------------------------------------------------
c                 
      subroutine smooth(x,n) 
c                 
c     three point triangular smoothing filter
c                 
      real x(n) 
      m=n-1       
      a=0.54*x(1)+0.46*x(2) 
      b=0.54*x(n)+0.46*x(m) 
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
      subroutine nmo(x,npts,sps,tmin,dist,vrms)
c                 
c     apply constant vrms velocity nmo to trace x
c                 
      include 'pltsyn.par'
c
      real x(ppseis),y(ppseis)
c                 
      term=(2.*dist/vrms)**2
      do 10 i=1,npts
         t0=tmin+float(i-1)/sps 
         t=sqrt(t0**2+term)
         n1=ifix((t-tmin)*sps+1) 
         n2=n1+1  
         if(n1.lt.1..or.n2.gt.npts) then
           y(i)=0. 
         else     
           t1=tmin+float(n1-1)/sps
           y(i)=(x(n2)-x(n1))*(t-t1)*sps+x(n1)
         end if   
10    continue    
c                 
      do 20 i=1,npts
         x(i)=y(i) 
20    continue    
c                 
      return      
      end         
c                 
c     ----------------------------------------------------------------
c                 
      subroutine intrtr(x,n1,n2)
c                 
c     through linear interpolation extend the array x consisting
c     n1 points so that it consists of n2 points (n2>n1 and n2<=7200)
c                 
      include 'pltsyn.par'
c
      real x(ppseis),y(ppseis),inc,m
c
      scale=float(n1-1)/float(n2-1)
      y(1)=x(1)   
      y(n2)=x(n1) 
c                 
      do 10 i=2,n2-1
         sf=scale*float(i-1) 
         intsf=int(sf)
         inc=sf-float(intsf)
         ip=intsf+1
         m=x(ip+1)-x(ip) 
         y(i)=m*inc+x(ip)
10    continue    
c                 
      do 20 i=1,n2
         x(i)=y(i)
20    continue    
c                 
      return      
      end         
