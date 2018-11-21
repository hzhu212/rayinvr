c
c     version 1.3  Aug 1992
c
c     Amplitude routines for TRAMP
c
c     ----------------------------------------------------------------
c
      subroutine calamp(npt,ifam,ir,namp,theta,dist,mir,caust,
     +                  xshot,layer1,iblk1,amprev,spamp,is,iray)
c
c     do preliminary calculations for the amplitude of a refracted or
c     reflected ray. the in-plane geometrical spreading is calculated
c     after all rays in the family have been traced in the routine
c     spread.
c
c     complex amplitude, a, is given by the equation
c
c                        a=vy/l
c
c     where  v -- product of ratios of impedances at each boundary
c            y -- product of zoeppritz coefficients at each boundary
c            l -- geometrical spreading
c
      include 'tramp.par'
      integer amprev,mir(1)
      real imped,jo,real,l,theta(1),dist(1)
      complex z,zc,y,cmplx
      include 'tramp.com'
c
      v=imped(vr(1,2),vp(1,2),vs(1,2))/
     +  imped(vr(npt,1),vp(npt,1),vs(npt,1))
      if(qa.eq.0.) then
        write(13,25) is,ir
        if(idump.eq.1) write(14,25) is,ir
25      format(2i5,'  ***  ray attenuated to zero',1x,
     +         'amplitude  ***')
      end if
      y=cmplx(1.,0.)*qa
      l=cos(ar(npt,1))
      namp=namp+1
      theta(namp)=fid1*ar(1,2)
      dist(namp)=fid1*(xr(npt)-xshot)
      mir(namp)=ir
      if(idump.eq.1) write(14,15) ifam,ir
15    format(/'refr. or refl. ray:   ifam=',i3,'  ir=',i5//
     +'irt iw1 iw2  vp1     vs1     vp2     vs2    angle     zc',
     +5x,'phase')
      if(nbnd.gt.1) then
        i1=amprev+1
        i2=2-amprev
        do 10 i=1,nbnd-1
           np=nptamp(i)
           v=v*imped(vr(np,2),vp(np,2),vs(np,2))/
     +         imped(vr(np,1),vp(np,1),vs(np,1))
           if(abs(vpamp(i,2)-vpamp(i,1)).le..002.and.
     +     abs(vsamp(i,2)-vsamp(i,1)).le..002.and.irtamp(i,3).eq.1)
     +     then
             write(13,55) is,ir
             if(idump.eq.1) write(14,55) is,ir
55           format(2i5,'  ***  ray refl. from zero vel.',1x,
     +             'contrast boundary  ***')
             zc=cmplx(0.,0.)
           else
             if(irtamp(i,3).eq.1) then
               call zprtz(1,irtamp(i,i1),irtamp(i,i2),
     +         vpamp(i,1),vpamp(i,2),vsamp(i,1),vsamp(i,2),
     +         aamp(i,i1),zc)
               if(izrefl(iray).gt.0) then
                 ii=izrefl(iray)
                 do 90 j=1,npzf(ii)-1
                    if(xr(np).ge.xzf(ii,j).and.xr(np).le.xzf(ii,j+1))
     +              then
                      zffp=((zff(ii,j+1)-zff(ii,j))/(xzf(ii,j+1)-
     +                     xzf(ii,j)))*(xr(np)-xzf(ii,j))+zff(ii,j)
                      zc=zc*zffp
c                  write(40,*) xr(np),zffp
                      go to 91
                    end if
90               continue
91               continue
               end if
             else
               call zprtz(2,irtamp(i,i1),irtamp(i,i2),
     +         vpamp(i,i1),vpamp(i,i2),vsamp(i,i1),vsamp(i,i2),
     +         aamp(i,i1),zc)
             end if
           end if
           y=y*zc
           if(idump.eq.1) then
             if(aimag(zc).ne.0..or.real(zc).ne.0.) then
               ph=atan2(aimag(zc),real(zc))*pi18
             else
               ph=0.0
             end if
             if(irtamp(i,3).eq.1) then
               iv1=1
               iv2=2
             else
               iv1=i1
               iv2=i2
             end if
             write(14,5) irtamp(i,3),irtamp(i,i1),irtamp(i,i2),
     +        vpamp(i,iv1),vsamp(i,iv1),vpamp(i,iv2),vsamp(i,iv2),
     +        aamp(i,i1)*pi18,cabs(zc),ph
5            format(i2,2i4,6f8.3,f8.2)
           end if
           l=l*cos(aamp(i,1))/cos(aamp(i,2))
10      continue
      end if
      v=sqrt(v)
      if(amprev.eq.1) v=1./v
c
      if(amprev.eq.0) then
        if(icbnd(1).eq.0) then
          iw0=-1
        else
          iw0=1
        end if
        iwend=iwave
        xend=xr(npt)
        vpe=vp(npt,1)
        vse=vs(npt,1)
        lend=layer
        aend=aamp(nbnd,1)
      else
        iw0=iwave
        if(icbnd(1).eq.0) then
          iwend=-1
        else
          iwend=1
        end if
        xend=xr(1)
        vpe=vp(1,2)
        vse=vs(1,2)
        lend=layer1
        alphae=atan(s(layer1,iblk1,1))
        aend=pi-fid1*(ar(1,2)+alphae)
      end if
c
c     calculate the surface conversion coefficient
c
      if(surcon.eq.1) then
        if(abs(ar(1,2)).le.pi2.or.amprev.eq.0) then
          idir=-1
        else
          idir=1
        end if
        la=lend+idir
80      if(la.eq.0.or.la.eq.(nlayer+1)) then
          vpa=.001
          vsa=.001
          go to 60
        end if
        call block(xend,la,iba)
        vpa=(vm(la,iba,4)-vm(la,iba,3))*(xend-xbnd(la,iba,1))/
     +      (xbnd(la,iba,2)-xbnd(la,iba,1))+vm(la,iba,3)
        if(ivg(la,iba).ne.-1) go to 70
        la=la+idir
        go to 80
70      vsa=vpa*vsvp(la,iba)
60      call zprtz(3,iwend,icomp,vpe,vpa,vse,vsa,aend,zc)
        if(idump.eq.1) write(14,5) 3,iwend,0,vpe,vse,vpa,vsa,
     +    aend*pi18,cabs(zc),atan2(aimag(zc),real(zc))*pi18
        y=y*zc
      end if
c
c     calculate the out-of-plane geometrical spreading, jo, by
c     evaluating a single integral of the velocity by the
c     trapezoidal rule
c
      sum=0.
      do 50 i=1,npt-1
         vave=(vr(i,2)+vr(i+1,1))/2.
         sum=sum+vave*tr(i)
50    continue
      jo=sum/vr(1,2)
c
c     the spreading function is given by the product jo*l where l
c     is the product of ratios of cosines at each interface
c
      spread=abs(jo*l)
c
      epfac=(-1)**iep
      if(spread.eq.0.) then
        z=0.
      else
        z=v*y/sqrt(spread)*epfac
      end if
      if(amprev.eq.1) then
        ampcor=vr(npt,1)/vr(1,2)
        z=z*ampcor
      end if
      if(iw0.eq.-1) z=z*spamp
      if(amprev.eq.0) then
        a0=ar(1,2)
        an=ar(npt,1)
      else
        a0=ar(npt,1)
        an=ar(1,2)
      end if
      if(surcon.eq.0) then
        if(icomp.eq.1) then
          if(iwend.eq.1) then
            z=z*abs(cos(an))
          else
            z=z*abs(sin(an))
          end if
        else
          if(iwend.eq.1) then
            z=z*abs(sin(an))
          else
            z=z*abs(cos(an))
          end if
        end if
      end if
      amp(ntt-1)=cabs(z)
      zimag=aimag(z)
      zreal=real(z)
      phase(ntt-1)=caust*pi2
      if(zimag.ne.0..or.zreal.ne.0.)
     +  phase(ntt-1)=phase(ntt-1)+atan2(zimag,zreal)
c
      if(idump.eq.1) then
        write(14,35) nbnd,v,cabs(y),l,jo*sin(a0),qa
35      format(/'nbnd=',i3,' v=',f6.3,' y=',f6.3,' l=',f6.3,1x,
     +  'jo=',f9.3,' q fctr=',f10.4)
        if(amprev.eq.1) write(14,45) ampcor
45      format('amprev fctr=',f10.4)
      end if
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine zprtz(irt,iw1,iw2,vp1,vp2,vs1,vs2,ain,zc)
c
c     calculate the zoeppritz's displacement amplitude coefficient
c
c       irt -- 1 for reflection, 2 for transmission, 3 for surface
c              conversion coefficient
c       iw1 -- incoming wave: +1 for p-wave
c                             -1 for s-wave
c       iw2 -- if irt=1 or 2 - outgoing wave: +1 for p-wave
c                                             -1 for s-wave
c              if irt=3 - seismometer type: 0 for horizontal component
c                                           1 for vertical component
c       vp1 -- incoming p-wave velocity
c       vp2 -- outgoing p-wave velocity
c       vs1 -- incoming s-wave velocity
c       vs2 -- outgoing s-wave velocity
c       ain -- angle of incidence in radians
c       zc  -- zoeppritz coefficient (complex)
c
      real denc(5)
      complex p1,p2,p3,p4,d,zc
      common /blk21/ denc,denmin
c
      a=abs(ain)
      if(vp1.lt..001) vp1=.001
      if(vp2.lt..001) vp2=.001
      if(vs1.lt..001) vs1=.001
      if(vs2.lt..001) vs2=.001
      rho1=denc(1)+denc(2)*vp1+denc(3)*vp1**2+denc(4)*vp1**3+
     +                                        denc(5)*vp1**4
      rho2=denc(1)+denc(2)*vp2+denc(3)*vp2**2+denc(4)*vp2**3+
     +                                        denc(5)*vp2**4
      if(rho1.lt.denmin) rho1=denmin
      if(vs1.le..001) rho1=1.0
      if(rho2.lt.denmin) rho2=denmin
      if(vs2.le..001) rho2=1.0
      if(iw1.eq.1) then
        vel=vp1
      else
        vel=vs1
      end if
      if(irt.eq.-2) then
        theta=1./vel
      else
        theta=sin(a)/vel
      end if
      theta2=theta**2
      qa=2.*(rho2*vs2**2-rho1*vs1**2)
      t1=vp1**2*theta2
      t2=vs1**2*theta2
      t3=vp2**2*theta2
      t4=vs2**2*theta2
      if(t1.gt.1.) then
        a1=0.
        b1=-sqrt(t1-1.)
      else
        a1=sqrt(1.-t1)
        b1=0.
      end if
      if(t2.gt.1.) then
        a2=0.
        b2=-sqrt(t2-1.)
      else
        a2=sqrt(1.-t2)
        b2=0.
      end if
      if(t3.gt.1.) then
        a3=0.
        b3=-sqrt(t3-1.)
      else
        a3=sqrt(1.-t3)
        b3=0.
      end if
      if(t4.gt.1.) then
        a4=0.
        b4=-sqrt(t4-1.)
      else
        a4=sqrt(1.-t4)
        b4=0.
      end if
      x=rho2-qa*theta2
      y=rho1+qa*theta2
      z=rho2-rho1-qa*theta2
      p1=cmplx(a1,b1)
      p2=cmplx(a2,b2)
      p3=cmplx(a3,b3)
      p4=cmplx(a4,b4)
      if(iabs(irt).le.2) then
        d=vp1*vp2*vs1*vs2*theta2*z**2+vp2*vs2*p1*p2*x**2+vp1*vs1*
     +    p3*p4*y**2+rho1*rho2*(vs1*vp2*p1*p4+vp1*vs2*p2*p3)+qa**2*
     +    theta2*p1*p2*p3*p4
        if(iw1.eq.1) then
          if(iabs(irt).eq.1) then
            if(iw2.eq.1) then
              zc=-1.+2.*p1*(vp2*vs2*p2*x**2+vs1*vp2*rho1*rho2*p4+qa**2
     +           *theta2*p2*p3*p4)/d
            else
              zc=-2.*vp1*theta*p1*(qa*p3*p4*y+vp2*vs2*x*z)/d
            end if
          else
            if(iw2.eq.1) then
              if(irt.eq.-2) p1=cmplx(1.,0.)
              zc=2.*vp1*rho1*p1*(vs2*p2*x+vs1*p4*y)/d
            else
              zc=-2.*vp1*rho1*theta*p1*(qa*p2*p3-vs1*vp2*z)/d
            end if
          end if
        else
          if(iabs(irt).eq.1) then
            if(iw2.eq.1) then
              zc=-2.*vs1*theta*p2*(qa*p3*p4*y+vp2*vs2*x*z)/d
            else
              zc=1.-2.*p2*(vp2*vs2*p1*x**2+vp1*vs2*rho1*rho2*p3+
     +           qa**2*theta2*p1*p3*p4)/d
            end if
          else
            if(iw2.eq.1) then
              zc=2.*vs1*rho1*theta*p2*(qa*p1*p4-vp1*vs2*z)/d
            else
              if(irt.eq.-2) p2=cmplx(1.,0.)
              zc=2.*vs1*rho1*p2*(vp1*p3*y+vp2*p1*x)/d
            end if
          end if
        end if
      else
        d=(1.-2.*vs1**2*theta2)**2+4.*p1*p2*theta2*vs1**3/vp1
        if(iw1.eq.1) then
          if(iw2.eq.1) then
            zc=2.*p1*(1.-2.*vs1**2*theta2)/d
          else
            zc=4.*p1*p2*theta*vs1/d
          end if
        else
          if(iw2.eq.1) then
            zc=-4.*p1*p2*vs1**2*theta/(d*vp1)
          else
            zc=2.*p2*(1.-2.*vs1**2*theta2)/d
          end if
        end if
      end if
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine spread(ifam,theta,dist,ir,nray,sdev,splnf,ampsmt,
     +                  iint,icaust,ibreak,xshotr,fidr,idump,is)
c
c     calculate the in-plane geometrical spreading by fitting a cubic
c     spline to the curve defined by the range versus take-off angle
c     for a given ray family (iint=1) or a given ray code within a ray
c     family (iint=0)
c
      include 'tramp.par'
      real theta(1),dist(1),dy(pnrayf),slope(pnrayf),
     +       dum1(pnrayf),dum2(pnrayf),t(pnrayf),d(pnrayf),
     +       as(pnrayf),rayid(pray),ray(prayf),amp(pray),
     +       phase(pray),p0(pray),pr(pray),range(pray),tt(pray),
     +       xshtar(pray),fidarr(pray),xzf(pzff,ppzff),zff(pzff,ppzff)
      double precision w(11*pnrayf+14)
      integer ir(1),icid(pnrayf),isid(pnrayf),idray(2),surcon,
     +        npzf(pzff),izrefl(prayf)
      common /blk11/ range,tt,rayid,xshtar,fidarr,idray,ntt,ray
      common /blk14/ amp,phase,p0,pr,surcon
      common /blk36/ izrefl,nzf,npzf,xzf,zff
c
      nstart=ntt-1-nray
      nh=0
      if(nray.eq.1) then
        amp(nstart+1)=-1.
        return
      end if
      isid(1)=1
      ng=1
      do 40 i=2,nray
         if(rayid(nstart+i).eq.rayid(nstart+i-1).or.iint.eq.1) then
           if(ng.eq.1.or.ibreak.eq.0) then
             ng=ng+1
             isid(i)=isid(i-1)
           else
             if(((dist(i)-dist(i-1))*(dist(i-1)-dist(i-2))).ge.0.)
     +         then
               ng=ng+1
               isid(i)=isid(i-1)
             else
               ng=1
               isid(i)=isid(i-1)+1
             end if
           end if
         else
           ng=1
           isid(i)=isid(i-1)+1
         end if
40    continue
1000  i=1
      t(i)=theta(i+nh)
      d(i)=dist(i+nh)
      isidm=isid(i+nh)
      if((i+nh).eq.nray) go to 65
900   if(isid(i+nh+1).ne.isidm) go to 65
      i=i+1
      t(i)=theta(i+nh)
      d(i)=dist(i+nh)
      if((i+nh).eq.nray) go to 65
      go to 900
65    npts=i
      if(icaust.eq.1) then
        do 70 i=1,npts
           icid(i)=int((rayid(nstart+nh+i)-float(int(rayid
     +             (nstart+nh+i))))*10.+.5)
70      continue
      end if
      iflag=0
      do 10 i=1,npts
         if(amp(nstart+nh+i).gt.0.) iflag=1
10    continue
      if(iflag.eq.0) go to 2000
      if(npts.eq.1) then
        amp(nstart+nh+1)=-1.
        go to 2000
      end if
      if(npts.eq.2) then
        slope(1)=(d(2)-d(1))/(t(2)-t(1))
        slope(2)=slope(1)
        if(idump.eq.1) then
          write(14,75) ifam
75        format(/'straight line:'//'ifam=',i3,'  nray=   2'//
     +      '       ray #    theta      dist      slope')
          write(14,85) (ir(i+nh),t(i),d(i),slope(i),i=1,2)
85        format(i10,2f12.5,e12.5)
        end if
        amp(nstart+nh+1)=amp(nstart+nh+1)/sqrt(abs(slope(1)))
        amp(nstart+nh+2)=amp(nstart+nh+2)/sqrt(abs(slope(1)))
        if(icaust.eq.1.and.slope(1).gt.0.) then
          if(icid(1).eq.1) phase(nstart+nh+1)=phase(nstart+nh+1)+pi2
          if(icid(2).eq.1) phase(nstart+nh+2)=phase(nstart+nh+2)+pi2
        end if
        d1=xshotr+fidr*d(1)
        d2=xshotr+fidr*d(2)
        if(amp(nstart+nh+1).gt.0.) write(13,25) is,ir(nh+1),d1,
     +  amp(nstart+nh+1),phase(nstart+nh+1)*pi18,rayid(nstart+nh+1)
        if(amp(nstart+nh+2).gt.0.) write(13,25) is,ir(nh+2),d2,
     +  amp(nstart+nh+2),phase(nstart+nh+2)*pi18,rayid(nstart+nh+2)
        go to 2000
      end if
      sn=splnf*float(npts)
      do 80 i=1,npts
         dy(i)=d(i)*sdev
80    continue
      if(t(2).lt.t(1)) then
        fics=-1.0
        do 30 i=1,npts
           t(i)=pi2-t(i)
30      continue
      else
        fics=1.0
      end if
c
      if(idump.eq.1) write(14,45) ifam,npts,splnf,sdev*100.
45    format(/'cubic spline:'//
     +      'ifam=',i3,'  nray=',i4,'  splnf=',f10.4,'  sdev=',f10.4//
     +      '       ray #    theta      dist      spln dist   slope')
c
      isflag=0
      call cspln(t,d,dy,sn,npts,w,dum1,slope,dum2,isflag)
      if(isflag.ne.0) go to 999
c
      do 20 i=1,npts
         amp(nstart+nh+i)=amp(nstart+nh+i)/sqrt(abs(slope(i)))
         if(icaust.eq.1.and.(slope(i)*fics).gt.0.and.icid(i).eq.1)
     +     phase(nstart+nh+i)=phase(nstart+nh+i)+pi2
20    continue
c
      nsmth=nint(ampsmt*(npts-2))
      if(nsmth.gt.0) then
        do 110 i=1,npts
           if(amp(nstart+nh+i).gt.0.) then
             as(i)=amp(nstart+nh+i)
           else
             as(i)=0.
           end if
110     continue
        do 120 i=1,nsmth
           call smooth(as,npts)
120     continue
        do 130 i=1,npts
           amp(nstart+nh+i)=as(i)
130     continue
      end if
c
      do 140 i=1,npts
         di=xshotr+fidr*d(i)
         if(amp(nstart+nh+i).gt.0.) write(13,25) is,ir(i+nh),di,
     +   amp(nstart+nh+i),phase(nstart+nh+i)*pi18,rayid(nstart+nh+i)
25       format(2i5,f10.3,e12.3,f10.2,f6.1)
140   continue
c
      if(idump.eq.1) write(14,55) (ir(i+nh),t(i),d(i),dum1(i),
     +                             slope(i),i=1,npts)
55    format(i10,3f12.5,e12.3)
      go to 2000
c
999   write(13,100) ifam
100   format('ray group=',i5,'  ***  problem with cubic spline  ***')
      do 200 i=1,npts
         amp(nstart+nh+i)=-1.
200   continue
2000  nh=nh+npts
      if(nh.lt.nray) go to 1000
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine cspln(x,y,dy,s,n,w,yy,yy1,yy2,isflag)
c
c     routine to fit a cubic spline to the n values of y versus x
c
      implicit double precision (a-h,o-z)
      real*4 x,y,dy,s,yy,yy1,yy2
      common /spln$/ ncount
      dimension x(n),y(n),dy(n),w(1)
      dimension yy(1),yy1(1),yy2(1)
      integer a,b,c,d,r,r1,r2,t,t1,u,v
c
      if (n.lt.3) then
        isflag=1
        return
      end if
      ncount=1
      ss=s
      ds=ss
      if(s.le.0.0) ss=1.0d-8
      ee=ds*0.5d-6
      if (s.le.0.0) ee=0.5d-6
      nm=n-1
      n2=n+2
      nb=n
      nc=nb+n
      nd=nc+n
      nr=nd+n
      nr1=nr+n2
      nr2=nr1+n2
      nt=nr2+n2
      nt1=nt+n2
      nu=nt1+n2
      nv=nu+n2
      w(nr+1)=0.0d0
      w(nr+2)=0.0d0
      w(nr2-1)=0.0d0
      w(nr2)=0.0d0
      w(nt-1)=0.0d0
      w(nt)=0.0d0
      w(nu+1)=0.0d0
      w(nu+2)=0.0d0
      w(nv-1)=0.0d0
      w(nv)=0.0d0
      p=0.0d0
c
c     set a to qt*y
c     set r,r1,r2, to d*q
c
      h=x(2)-x(1)
      f=(y(2)-y(1))/h
      do 1 i=2,nm
         if(h.le.0.0d0) then
           isflag=2
           return
         end if
         a=i
         g=h
         h=x(i+1)-x(i)
         e=f
         f=(y(i+1)-y(i))/h
         w(i)=f-e
         t=nt+i
         w(t+1)=.6666667d0*(g+h)
         t1=nt1+i
         w(t1+1)=.3333333d0*h
         r2=nr2+i
         w(r2+1)=dy(i-1)/g
         r=nr+i
         w(r+1)=dy(i+1)/h
         r1=nr1+i
1        w(r1+1)=-dy(i)/g-dy(i)/h
c
c     set b,c and d to qt*d*d*q
c
      do 2 i=2,nm
         b=nb+i
         r=nr+i
         r1=nr1+i
         r2=nr2+i
         w(b)=w(r+1)*w(r+1)+w(r1+1)*w(r1+1)+w(r2+1)*w(r2+1)
         c=nc+i
         w(c)=w(r+1)*w(r1+2)+w(r1+1)*w(r2+2)
         d=nd+i
2        w(d)=w(r+1)*w(r2+3)
3     continue
c
c     do a ldu decomposition of qt*d*d*q+p*t
c
      do 4 i=2,nm
         a=i
         r1=nr1+i
         r=nr+i
         r2=nr2+i
         u=nu+i
         t1=nt1+i
         b=nb+i
         c=nc+i
         d=nd+i
         t=nt+i
         w(r1)=f*w(r)
         w(r2-1)=g*w(r-1)
         w(r+1)=1.0d0/(w(b)+p*w(t+1)-f*w(r1)-g*w(r2-1))
         w(u+1)=w(a)-w(r1)*w(u)-w(r2-1)*w(u-1)
         f=w(c)+p*w(t1+1)-h*w(r1)
         g=h
4        h=w(d)
c
c     do back substitution
c
      do 5 i=2,nm
         ii=n2-i
         u=nu+ii
         r=nr+ii
         r1=nr1+ii
         r2=nr2+ii
5        w(u)=w(r)*w(u)-w(r1)*w(u+1)-w(r2)*w(u+2)
c
      e=0.0d0
      h=0.0d0
c
c     the v is really d*d*q*u
c
      do 6 i=1,nm
         u=nu+i
         v=nv+i
         g=h
         h=(w(u+2)-w(u+1))/(x(i+1)-x(i))
         w(v+1)=(h-g)*dy(i)*dy(i)
6        e=e+w(v+1)*(h-g)
c
      g=-h*dy(n)*dy(n)
      v=nv+n
      w(v+1)=g
      e=e-g*h
      if (e.le.ds.or.dabs(e-ds).le.ee) go to 8
c
c     calculate f and g
c     for g consider ... a*x=i ; a*(x*y)=i*y=y
c
      f=0.0d0
      g=0.0d0
      do 7 i=2,nm
         r=nr+i
         r1=nr1+i
         r2=nr2+i
         u=nu+i
         t=nt+i
         t1=nt1+i
         h=w(u)*w(t1)+w(u+1)*w(t+1)+w(u+2)*w(t1+1)
         f=f+w(u+1)*h
         h=h-w(r1)*w(r)-w(r2-1)*w(r-1)
         g=g+h*w(r+1)*h
7        w(r+1)=h
c
      h=f-p*g
      if (h.le.0.0d0) go to 8
      ncount=ncount+1
      if (ncount.gt.100) then
        isflag=3
        return
      end if
c
      p=p+dsqrt(e/ss)*(e-dsqrt(ds*e))/h
c
      go to 3
8     continue
c
      do 9 i=1,n
         c=nc+i
         v=nv+i
         u=nu+i
         a=i
         w(a)=y(i)-w(v+1)
9        w(c)=p*w(u+1)
c
      do 10 i=1,nm
         d=nd+i
         a=i
         c=nc+i
         b=nb+i
         h=x(i+1)-x(i)
         w(d)=(w(c+1)-w(c))/(3.0d0*h)
10       w(b)=(w(a+1)-w(a))/h-(h*w(d)+w(c))*h
c
c     calculate the derivates of the cubic spline
c
      if (n.le.0) return
      j=1
      a=1
      b=a+n
      c=b+n
      d=c+n
      do 15 i=1,n
11       if (x(i).lt.x(j)) go to 14
         if (x(i).lt.x(j+1)) go to 13
         if (j.lt.nm) go to 12
         if (x(i).eq.x(j+1)) go to 13
         isflag=4
         return
12       j=j+1
         a=j
         b=a+n
         c=b+n
         d=c+n
         go to 11
13       diff=x(i)-x(j)
         yy(i)=w(a)+diff*(w(b)+diff*(w(c)+diff*w(d)))
         yy1(i)=w(b)+diff*(2.0d0*w(c)+3.0d0*diff*w(d))
         yy2(i)=2.0d0*w(c)+6.0d0*w(d)*diff
         go to 15
14       if(j.eq.1) then
           isflag=5
           return
         end if
         j=1
         a=j
         b=a+n
         c=b+n
         d=c+n
         go to 11
15    continue
c
      return
      end
c
c     ----------------------------------------------------------------
c
      function imped(vw,vpw,vsw)
c
c     calculate the impedance imped=p*v where
c
c        p is given by fourth order polynomial   if vs > 0.0
c        p = 1.0                                 if vs = 0.0
c
      real imped,denc(5)
      common /blk21/ denc,denmin
c
      p=denc(1)+denc(2)*vpw+denc(3)*vpw**2+denc(4)*vpw**3+
     +  denc(5)*vpw**4
      if(p.lt.denmin) p=denmin
      if(vsw.le..001) p=1.0
      imped=vw*p
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine calsec(xshots,itt,nsec,isec,isect,iint,ibreak,icmp)
c
c     calculate a synthetic seismic section (one for each shot) using
c     the calculated distance, time, amplitude and phase of each ray
c     reaching the surface. the time, amplitude, and phase at a
c     particular seismogram is determined by linear interpolation
c     between two calculated values of the same ray family and ray
c     code (iint=0) or two values of the same ray family (iint=1)
c
      include 'tramp.par'
      real x(pnrayf),t(pnrayf),a(pnrayf),p(pnrayf),
     +     dist(pseis),sect(pseis,parriv,4),xshots(1),
     +     xh(pnrayf),th(pnrayf),ah(pnrayf),
     +     ph(pnrayf),rh(pnrayf),disth(pseis),p0h(pnrayf),
     +     prh(pnrayf),p1(pnrayf),p2(pnrayf)
      integer itt(1),na(pseis),nsec(1),iint(1),
     +        iflags(pseis),ibreak(1),isint(pnrayf)
      include 'tramp.com'
      data iflags/pseis*0/
c
      n=0
      ifam=0
      if(icmp.ne.1) then
        if(isect.eq.1) then
          nseis=int((xmaxs-xmins)/xincs+.5)+1
          do 30 i=1,nseis
             dist(i)=xmins+float(i-1)*xincs
30        continue
        end if
      else
        nseis=int((xmaxs-xmins)/xincs+.5)+1
        do 31 i=1,nseis
           disth(i)=xmins+float(i-1)*xincs
31      continue
      end if
c
      do 2000 ii=1,isec
         if(nsec(ii).eq.0) go to 2000
         if(isect.eq.2) then
           rewind(16)
           nseis=0
300        read(16,*) dists,iflagr
           if(iflagr.eq.-1) go to 399
           if(iflagr.gt.0) then
             if(abs(dists-xshots(ii)).lt..001.or.iflagr.eq.2) then
               irec=1
             else
               irec=0
             end if
           else
             if(irec.eq.1) then
               nseis=nseis+1
               dist(nseis)=dists
             end if
           end if
           go to 300
399        if(nseis.eq.0.or.nseis.gt.pseis) then
            write(6,475)
475         format(/'***  error in file rec.in  ***'/)
            go to 2000
           end if
         else
           if(icmp.eq.1) then
             do 32 i=1,nseis
                dist(i)=disth(i)+xshots(ii)
32           continue
           end if
         end if
         do 40 i=1,nseis
            na(i)=0
40       continue
         do 10 i=1,nsec(ii)
            ifam=ifam+1
            nh1=n
            if(itt(ifam).gt.0) then
c
c             do not interpolate across rays of different ray code
c
              nptsa=0
              if(iint(ifam).eq.0) then
                do 20 j=1,itt(ifam)
                   if(amp(j+nh1).gt.0.) then
                     nptsa=nptsa+1
                     xh(nptsa)=range(j+nh1)
                     th(nptsa)=tt(j+nh1)
                     ah(nptsa)=amp(j+nh1)
                     ph(nptsa)=phase(j+nh1)
                     rh(nptsa)=rayid(j+nh1)
                     p0h(nptsa)=p0(j+nh1)
                     prh(nptsa)=pr(j+nh1)
                   end if
                   n=n+1
20              continue
                if(nptsa.gt.1) then
                  nh2=0
1000              if(nh2.lt.nptsa) then
                    k=0
100                 k=k+1
                    if((k+nh2).le.nptsa) then
                      x(k)=xh(k+nh2)
                      t(k)=th(k+nh2)
                      a(k)=ah(k+nh2)
                      p(k)=ph(k+nh2)
                      p1(k)=p0h(k+nh2)
                      p2(k)=prh(k+nh2)
                      if(k.eq.1) then
                        rhc=rh(k+nh2)
                        go to 100
                      else
                        if(rh(k+nh2).eq.rhc) go to 100
                      end if
                    end if
                    npt=k-1
                    nh2=nh2+npt
                    if(npt.gt.1) then
                     isint(1)=1
                     if(npt.gt.2) then
                       isflag=0
                       do 150 l=2,npt-1
                          if(((x(l+1)-x(l))*(x(l)-x(l-1))).ge.0..or.
     +                    isflag.eq.1.or.ibreak(ifam).eq.0) then
                            isint(l)=1
                            isflag=0
                          else
                            isint(l)=0
                            isflag=1
                          end if
150                    continue
                     end if
                     do 50 j=1,nseis
                      if(iflags(j).eq.1) go to 50
                      itwice=0
                      do 60 l=1,npt-1
                       if(itwice.eq.0) then
                        if((dist(j).ge.x(l).and.dist(j).le.x(l+1)).or.
     +                   (dist(j).le.x(l).and.dist(j).ge.x(l+1))) then
                          if(na(j).eq.parriv) then
                            write(13,35) dist(j)
35            format('***  too many arrivals for seismogram at '
     +               ,f7.2,' km  ***')
                            iflags(j)=1
                            go to 50
                          end if
                          if(isint(l).eq.0) go to 60
                          if(dist(j).eq.x(l+1)) then
                            itwice=1
                          else
                            itwice=0
                          end if
                          na(j)=na(j)+1
                          denom=x(l+1)-x(l)
                          sect(j,na(j),1)=((t(l+1)-t(l))*dist(j)+
     +                      t(l)*x(l+1)-t(l+1)*x(l))/denom
                          sect(j,na(j),2)=((a(l+1)-a(l))*dist(j)+
     +                      a(l)*x(l+1)-a(l+1)*x(l))/denom
                          sect(j,na(j),3)=((p(l+1)-p(l))*dist(j)+
     +                      p(l)*x(l+1)-p(l+1)*x(l))/denom
                          sect(j,na(j),4)=rhc
                        end if
                       else
                        itwice=0
                       end if
60                    continue
50                   continue
                    end if
                    go to 1000
                  end if
                end if
              else
c
c               interpolate across rays of the same family but not
c               necessarily the same ray code
c
                do 70 j=1,itt(ifam)
                   if(amp(j+nh1).gt.0.) then
                     nptsa=nptsa+1
                     x(nptsa)=range(j+nh1)
                     t(nptsa)=tt(j+nh1)
                     a(nptsa)=amp(j+nh1)
                     p(nptsa)=phase(j+nh1)
                     p1(nptsa)=p0(j+nh1)
                     p2(nptsa)=pr(j+nh1)
                   end if
                   n=n+1
70              continue
                if(nptsa.gt.1) then
                  isint(1)=1
                  if(nptsa.gt.2) then
                    isflag=0
                    do 160 l=2,nptsa-1
                       if(((x(l+1)-x(l))*(x(l)-x(l-1))).ge.0..or.
     +                 isflag.eq.1.or.ibreak(ifam).eq.0) then
                         isint(l)=1
                         isflag=0
                       else
                         isint(l)=0
                         isflag=1
                       end if
160                 continue
                  end if
                  do 80 j=1,nseis
                     if(iflags(j).eq.1) go to 80
                     itwice=0
                     do 90 l=1,nptsa-1
                       if(itwice.eq.0) then
                        if((dist(j).ge.x(l).and.dist(j).le.x(l+1)).or.
     +                  (dist(j).le.x(l).and.dist(j).ge.x(l+1))) then
                          if(na(j).eq.parriv) then
                            write(13,35) dist(j)
                            iflags(j)=1
                            go to 80
                          end if
                          if(isint(l).eq.0) go to 90
                          if(dist(j).eq.x(l+1)) then
                            itwice=1
                          else
                            itwice=0
                          end if
                          na(j)=na(j)+1
                          denom=x(l+1)-x(l)
                          sect(j,na(j),1)=((t(l+1)-t(l))*dist(j)+
     +                      t(l)*x(l+1)-t(l+1)*x(l))/denom
                          sect(j,na(j),2)=((a(l+1)-a(l))*dist(j)+
     +                      a(l)*x(l+1)-a(l+1)*x(l))/denom
                          sect(j,na(j),3)=((p(l+1)-p(l))*dist(j)+
     +                      p(l)*x(l+1)-p(l+1)*x(l))/denom
                          sect(j,na(j),4)=rhc
                        end if
                      else
                        itwice=0
                      end if
90                   continue
80                continue
                end if
              end if
            end if
10       continue
         write(21,15) xshots(ii),-1
         write(210,15) xshots(ii),-1
         do 180 i=1,nseis
            write(21,15) dist(i),na(i)
            write(210,15) dist(i),na(i)
15          format(f10.3,i10)
25          format(2f10.3,i10)
            if(na(i).gt.0) then
              do 190 j=1,na(i)
                 if(vred.eq.0.) then
                   sectt=sect(i,j,1)
                 else
                   sectt=sect(i,j,1)+abs(dist(i)-xshots(ii))/vred
                 end if
                 write(21,45) sectt,(sect(i,j,k),k=2,3)
                 write(210,455) sectt,(sect(i,j,k),k=2,4)
45               format(3e12.5)
455              format(3e12.5,f4.1)
190           continue
            end if
180      continue
2000  continue
      return
      end
