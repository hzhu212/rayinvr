c
c     version 1.3  Aug 1992
c     
c     Inversion routines for RAYINVR
c                 
c     ----------------------------------------------------------------
c                 
      subroutine cvcalc(i,j,ipos,itype,cf)
c                 
c     calculate coefficients of velocity partial derivatives
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c
      if(itype.eq.1) then
        sign=1.
        sb=s(i,j,2)
        xb=xbnd(i,j,2)
        bb=b(i,j,2)
      end if
      if(itype.eq.2) then
        sign=-1.
        sb=s(i,j,2)
        xb=xbnd(i,j,1)
        bb=b(i,j,2)
      end if
      if(itype.eq.3) then
        sign=-1.
        sb=s(i,j,1)
        xb=xbnd(i,j,2)
        bb=b(i,j,1)
      end if
      if(itype.eq.4) then
        sign=1.
        sb=s(i,j,1)
        xb=xbnd(i,j,1)
        bb=b(i,j,1)
      end if
c
      cv(i,j,ipos,1)=cv(i,j,ipos,1)+cf*sign*(sb*xb-bb)
      cv(i,j,ipos,2)=cv(i,j,ipos,2)-cf*sign*sb
      cv(i,j,ipos,3)=cv(i,j,ipos,3)-cf*sign*xb
      cv(i,j,ipos,4)=cv(i,j,ipos,4)+cf*sign
      cv(i,j,ipos,5)=cv(i,j,ipos,5)+cf*sign*bb*xb
c
      return
      end   
c
c     ----------------------------------------------------------------
c                 
      subroutine velprt(lu,iu,npt)
c                 
c     calculate velocity partial derivatives
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c
      do 10 j=1,4
         if(ivv(lu,iu,j).ne.0) then
           jv=ivv(lu,iu,j)
           sl=sqrt((xr(npt)-xr(npt-1))**2+(zr(npt)-zr(npt-1))**2)
           v1=vr(npt-1,2)
           v2=vr(npt,1)
           dvj1=(cv(lu,iu,j,1)*xr(npt-1)+cv(lu,iu,j,2)*xr(npt-1)**2+
     +          cv(lu,iu,j,3)*zr(npt-1)+cv(lu,iu,j,4)*xr(npt-1)*
     +          zr(npt-1)+cv(lu,iu,j,5))/(c(lu,iu,6)*xr(npt-1)+
     +          c(lu,iu,7))
           dvj2=(cv(lu,iu,j,1)*xr(npt)+cv(lu,iu,j,2)*xr(npt)**2+
     +          cv(lu,iu,j,3)*zr(npt)+cv(lu,iu,j,4)*xr(npt)*
     +          zr(npt)+cv(lu,iu,j,5))/(c(lu,iu,6)*xr(npt)+
     +          c(lu,iu,7))
           dvv2=(dvj1/v1**2+dvj2/v2**2)/2.
           fpart(ninv,jv)=fpart(ninv,jv)-sl*dvv2
         end if
10    continue
c
      return
      end
c                 
c     ----------------------------------------------------------------
c                 
      subroutine bndprt(lu,iu,v1,v2,a1,a2,alpha,npt,itype)
c                 
c     calculate boundary partial derivatives
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c
      if(iabs(itype).eq.2) a2=pi-a2
      if(itype.gt.0) then
        iz1=3
        iz2=4
      else 
        iz1=1
        iz2=2
      end if
c
      do 10 i=iz1,iz2
         if(izv(lu,iu,i).ne.0) then
           jv=izv(lu,iu,i)
           if(cz(lu,iu,i,2).ne.0.) then
             slptrm=cos(alpha)*abs(cz(lu,iu,i,1)-xr(npt))/
     +              cz(lu,iu,i,2)
           else
             slptrm=1.
           end if
           if(itype.lt.0) then
             ibz=0
           else
             ibz=1
           end if
           if(nzed(lu+ibz).eq.1) then
             fudge=2.
           else
             fudge=1.
           end if
           fpart(ninv,jv)=fpart(ninv,jv)+
     +       sign(fudge,1.*itype)*(cos(a1)/v1-cos(a2)/v2)*slptrm
         end if 
10    continue
c
      return
      end 
c                 
c     ----------------------------------------------------------------
c                 
      subroutine frprt(vfr,afr,alpha,npt,ifrpt)
c                 
c     calculate partial derivatives for floating reflectors
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c
      x1=xfrefl(ifcbnd,ifrpt)
      x2=xfrefl(ifcbnd,ifrpt+1)
c
      do 10 i=ifrpt,ifrpt+1
         if(ivarf(ifcbnd,i).gt.0) then
           jv=ivarf(ifcbnd,i)
           if(x2-x1.ne.0.) then
             if(i.eq.ifrpt) then
               ind=ifrpt+1
             else
               ind=ifrpt
             end if
             slptrm=cos(alpha)*abs(xfrefl(ifcbnd,ind)-xr(npt))/
     +              (x2-x1)
           else
             slptrm=1.
           end if
           fpart(ninv,jv)=fpart(ninv,jv)+2.*(cos(afr)/vfr)*slptrm
         end if 
10    continue
c
      return
      end 
c                 
c     ----------------------------------------------------------------
c                 
      subroutine fxtinv(npt)
c
c     keep track of those rays in the current family which reach the 
c     surface of the model
c
      include 'rayinvr.par'
      include 'rayinvr.com'
c
      if(vr(npt,2).ne.0.) then
        xfinv(ninv)=range(ntt-1)
        tfinv(ninv)=tt(ntt-1)
      else
        ninv=ninv-1
      end if
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine calprt(xshotr,ig,iray,idr,ximax,iflagw,iszero,x2pt)
c
c     calculate the partial derivatives at the location of the travel
c     time picks by interpolating across the end points of a single 
c     ray family
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'
c
      npick=0
      iflagw=0
      nsfc=1
      isf=ilshot(nsfc)
c
100   xf=xpf(isf)
      tf=tpf(isf)
      uf=upf(isf)
      irayf=ipf(isf)
c
      if(irayf.lt.0) go to 999
      if(irayf.eq.0) then
        xshotf=xf
        idf=sign(1.,tf)
        if(abs(xshotr-xshotf).lt..001.and.idr.eq.idf) then 
          iflag=1
          npick=isf-nsfc
          isf=isf+1
        else
          iflag=0
          nsfc=nsfc+1
          isf=ilshot(nsfc)
        end if
      else
        npick=npick+1
        if(iflag.eq.1.and.iray.eq.irayf) then
          if(iszero.eq.0) then
            xpick=xf 
          else
            xpick=abs(xshotr-xf)
          end if
          tpick=tf
          upick=uf
          tvmin=99999.
          if(ninv.gt.1) then
            do 10 i=1,ninv-1
               x1=xfinv(i)-x2pt 
               x2=xfinv(i+1)+x2pt 
               x3=xfinv(i)+x2pt 
               x4=xfinv(i+1)-x2pt 
               if((xpick.ge.x1.and.xpick.le.x2).or.
     +           (xpick.le.x3.and.xpick.ge.x4)) then 
                 if(ximax.gt.0..and.min(abs(xfinv(i)-xpick),
     +           abs(xfinv(i+1)-xpick)).gt.ximax) then
                   if(iflagw.eq.0) write(11,15) 
15                 format
     +             ('***  attempt to interpolate over ximax  ***')
                   iflagw=1
                   go to 10
                 end if
                 denom=xfinv(i+1)-xfinv(i)
                 if(denom.ne.0.) then
                   tintr=(xpick-xfinv(i))/denom*
     +                   (tfinv(i+1)-tfinv(i))+tfinv(i)
                 else
                   tintr=(tfinv(i+1)+tfinv(i))/2.
                 end if
                 if(tintr.lt.tvmin) then
                   ifpos=i 
                   tvmin=tintr
                 end if
               end if
10          continue
          else
             x1=xfinv(1)-x2pt
             x2=xfinv(1)+x2pt
             if(xpick.ge.x1.and.xpick.le.x2) then
               if(ximax.gt.0..and.abs(xfinv(1)-xpick).gt.ximax) 
     +         then
                 if(iflagw.eq.0) write(11,15)
                 iflagw=1
                 go to 1010
               end if
               tintr=tfinv(1)
               if(tintr.lt.tvmin) then
                 ifpos=i
                 tvmin=tintr
               end if
             end if
          end if
1010      if(tvmin.gt.9999.) go to 110
          iapos=0
          if(narinv.gt.0) then
            do 20 i=1,narinv
               if(ipinv(i).eq.npick) then
                 iapos=i
                 go to 30
               end if
20          continue
          end if
30        if(iapos.eq.0) then
            narinv=narinv+1
            iapos=narinv
          else  
            if(tcalc(iapos).le.tvmin) go to 110
          end if
          ipinv(iapos)=npick
          tobs(iapos)=tpick
          uobs(iapos)=upick
          tcalc(iapos)=tvmin
          xcalc(iapos)=xpick
          xscalc(iapos)=xshotr
          icalc(iapos)=idr*iray
          ircalc(iapos)=ig
          dx=xfinv(ifpos+1)-xfinv(ifpos)
          if(dx.ne.0.) then
            c1=abs((xpick-xfinv(ifpos))/dx) 
            c2=abs((xpick-xfinv(ifpos+1))/dx) 
          else
            c1=.5
            c2=.5
          end if
          do 40 i=1,nvar
             if(fpart(ifpos,i).ne.0..and.fpart(ifpos+1,i).ne.0.) then
               apart(iapos,i)=c1*fpart(ifpos,i)+c2*fpart(ifpos+1,i)
             else
               apart(iapos,i)=0.
             end if
40        continue
        end if
110     isf=isf+1
      end if
      go to 100
c
999   return
      end
