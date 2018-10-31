c
c     version 1.3  Aug 1992
c
c     Head wave routine for RAYINVR
c                 
c     ----------------------------------------------------------------
c                 
      subroutine hdwave(ifam,ir,n,invr,xsmax,iflag,i1ray,modout)
c                 
c     ray now travels along a layer boundary as a head wave.
c     do not need to use runge kutta routine. after ray has
c     traveled a distance hws further than previous ray, it
c     is directed upward into the model.
c                 
      include 'rayinvr.par'
      include 'rayinvr.com'           
c
      dhw=0.
10    l1=layer+1
      call block(xr(n),l1,ib) 
      if(ivg(l1,ib).ne.-1) go to 20
      layer=layer+1
      if(layer.eq.nlayer) go to 999
      go to 10
20    layer=l1
      iblk=ib
      idray(2)=3
      hwsl=hws
      ilr=(id+3)/2
      irl=(3-id)/2
      alpha=atan(s(layer,iblk,1))
c
1000  if(dhw.ge.(tdhw-.001)) then
        lstart=layer
        istart=iblk
        ihdw=0
        tdhw=tdhw+hws
        if(dhw.eq.0.) then
          vh=vel(xr(n),zr(n))
          if(iwave.eq.-1) vh=vh*vsvp(layer,iblk)
        end if
30      l1=layer-1
        call block(xr(n),l1,ib)
        if(ivg(l1,ib).ne.-1) go to 40
        layer=layer-1
        if(layer.eq.1) go to 999
        go to 30
40      layer=l1
        iblk=ib
        vp(n,2)=vel(xr(n),zr(n))
        vs(n,2)=vp(n,2)*vsvp(layer,iblk)
        nccbnd=nccbnd+1
        if(icbnd(iccbnd).eq.nccbnd) then
          iconvt=1
          iwave=-iwave
          iccbnd=iccbnd+1
        end if
        if(iwave.eq.1) then
          vr(n,2)=vp(n,2)
        else
          vr(n,2)=vs(n,2)
        end if
        if(iwave.eq.-1.and.vr(n,2).le..001) then
          vr(n,2)=0.
          iflag=1
          if(ir.ne.0) write(11,15)
15        format('***  ray stopped - s-wave cannot propagate  ***')
          return
        end if
        if(ibsmth.eq.0) then
          alpha=atan(s(layer,iblk,2))
        else
          npl=ifix((xr(n)-xmin-.001)/xsinc)+1
          npr=npl+1
          xpl=xmin+float(npl-1)*xsinc
          alpha=(cosmth(layer+1,npr)-cosmth(layer+1,npl))*(xr(n)-xpl)
     +          /xsinc+cosmth(layer+1,npl)
        end if
        if(i1ray.eq.1) then
          nhskip=0
          i1ray=0
        else
          nhskip=n-1
        end if
        if(dhw.ne.0.) vh=vr(n,1)
        vratio=vr(n,2)/vh
        if(vratio.ge.1.) then
          if(idiff.eq.0) then
            vr(n,2)=0.
            iflag=1
            return
          else
            a2=.99999*pi2
          end if
        else
          a2=asin(vratio)
        end if
        ar(n,2)=fid*(pi-a2)-alpha
        if(invr.eq.1.and.ir.gt.0) call bndprt
     +    (lstart,istart,vr(n,1),vr(n,2),pi2,a2,alpha,n,-1)
c        if(ir.gt.0.and.modout.eq.1) write(32,155) xr(n),zr(n),4
c155     format(2f10.3,i10)
        iflag=0
        return
      else
c
        if(idump.eq.1) write(12,5) ifam,ir,n,xr(n),zr(n),
     +    ar(n,1)*pi18,ar(n,2)*pi18,vr(n,1),vr(n,2),
     +    layer,iblk,id,iwave
5       format(i2,i3,i4,2f8.3,2f8.2,2f7.2,4i3)
c
      end if
c
500   xrp=xr(n)+fid*hwsl*cos(alpha)
      if((fid*xrp).gt.(fid*xbnd(layer,iblk,ilr))) then
        hl=abs(xbnd(layer,iblk,ilr)-xr(n))/cos(alpha)
        hwsl=hwsl-hl
        if(n.eq.ppray) go to 900
        n=n+1
        nbnd=nbnd+1
        if(nbnd.eq.piray) npskp=1
        if(npskp.ne.1) nbnda(nbnd)=n
        xr(n)=xbnd(layer,iblk,ilr)
        zr(n)=s(layer,iblk,1)*xr(n)+b(layer,iblk,1)
        if(xsmax.gt.0..and.abs(xr(n)-xr(1)).gt.xsmax) then
          vr(n,2)=0.
          iflag=1
          ihdwf=0
          return
        end if
        vp(n,1)=vm(layer,iblk,ilr)
        vs(n,1)=vp(n,1)*vsvp(layer,iblk)
        if(iwave.eq.1) then
          vr(n,1)=vp(n,1)
        else
          vr(n,1)=vs(n,1)
        end if
        ar(n,1)=ar(n-1,2)
        iblk=iblk+id
        if(iblk.lt.1.or.iblk.gt.nblk(layer)) go to 999
50      vp(n,2)=vm(layer,iblk,irl)
        if(vp(n,2).gt.0.) go to 60
        l1=layer+1
        if(l1.gt.nlayer) go to 999
        call block(xr(n)+fid*.001,l1,ib)
        layer=l1
        iblk=ib
        go to 50
60      vs(n,2)=vp(n,1)*vsvp(layer,iblk)
        if(iwave.eq.1) then
          vr(n,2)=vp(n,2)
        else
          vr(n,2)=vs(n,2)
        end if
        alpha=atan(s(layer,iblk,1))
        ar(n,2)=fid*pi2-alpha
        if(invr.eq.1.and.ir.gt.0) call velprt(layer,iblk,n)
c
        if(idump.eq.1) write(12,5) ifam,ir,n,xr(n),zr(n),
     +    ar(n,1)*pi18,ar(n,2)*pi18,vr(n,1),vr(n,2),
     +    layer,iblk,id,iwave
c
        go to 500
      else
        hwsl=hws
        if(n.eq.ppray) go to 900
        n=n+1
        nbnd=nbnd+1
        if(nbnd.eq.piray) npskp=1
        if(npskp.ne.1) nbnda(nbnd)=n
        xr(n)=xrp
        zr(n)=s(layer,iblk,1)*xr(n)+b(layer,iblk,1)
        if(xsmax.gt.0..and.abs(xr(n)-xr(1)).gt.xsmax) then
          vr(n,2)=0.
          iflag=1
          ihdwf=0
          return
        end if
        vp(n,1)=vel(xr(n),zr(n))
        vs(n,1)=vp(n,1)*vsvp(layer,iblk)
        if(iwave.eq.1) then
          vr(n,1)=vp(n,1)
        else
          vr(n,1)=vs(n,1)
        end if
        ar(n,1)=ar(n-1,2)
        vp(n,2)=vp(n,1)
        vs(n,2)=vs(n,1)
        vr(n,2)=vr(n,1)
        ar(n,2)=ar(n,1)
        if(invr.eq.1.and.ir.gt.0) call velprt(layer,iblk,n)
        dhw=dhw+hws
        go to 1000
      end if
      return      
c
900   write(11,25)
25    format('***  ray stopped - consists of too many points  ***')
999   vr(n,2)=0.
      iflag=1
      ihdwf=0
      return
      end         
