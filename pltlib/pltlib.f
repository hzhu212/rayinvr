c
c     version 1.3  Aug 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |           ***********  P L T L I B  ***********              |
c     |                                                              |
c     |         Plot library to convert Calcomp-like calls           |
c     |          to Uniras or other local graphics system            |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c     The purpose of each routine is described below the subroutine
c     statement.
c
c     The routine pcolor is required only if the local graphics allows
c     colour plotting. The routine segmnt creates a Uniras device-
c     independent plot metafile and can be ignored if Uniras is not
c     the local graphics system.
c
c     The following is a list and description of all Uniras plot calls:
c
c     groute - select graphics device
c     gshmes - set Uniras error-message display mode
c     gopen  - open uniras
c     grpsiz - inquire plot area size (for background colour only)
c     rrect  - plot rectangle (for background colour only)
c     gclear - erase screen
c     gempty - empty graphics buffer
c     gwicol - set polyline colour
c     gvect  - plot polyline
c     rtxcol - set text colour
c     rtxhei - set text height
c     rtxang - set text angle
c     rtxn   - plot floating point number
c     rtx    - plot text
c     gdot   - plot dot
c     gsegcr - create segment file (device-independent plot metafile)
c     gsegcl - close segment file (device-independent plot metafile)
c     gclose - terminate Uniras
c
c     ----------------------------------------------------------------
c
      subroutine plots(x,y,iroute)
c
c     initialize the plot; this must be the first plot call
c
c     iroute = 1 plots to the screen
c              2 creates a postscript file
c              3 creates a uniras file for the VERSATEC plotter
c              4 creates a postscript file for the colour plotter
c              5 creates a legal size postscript file
c              6 creates an A3 postscript file
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      character*50 proute
c
      x=x*sf
      y=y*sf
c
      if(iplot.ge.0) then
        if(iroute.lt.2.or.iroute.gt.7) then
          if(x.le.0..or.y.le.0.) then
            proute='select mx11;exit'
          else
            if(x.gt.350.) x=350.
            if(y.gt.265.) y=265.
            ix=nint(x)
            iy=nint(y)
            ixd1=ix/100
            ixd2=(ix-ixd1*100)/10
            ixd3=ix-ixd1*100-ixd2*10
            iyd1=iy/100
            iyd2=(iy-iyd1*100)/10
            iyd3=iy-iyd1*100-iyd2*10
            proute='select mx11;area    ,   ;exit'
            proute(18:20)=char(ixd1+48)//char(ixd2+48)//char(ixd3+48)
            proute(22:24)=char(iyd1+48)//char(iyd2+48)//char(iyd3+48)
          end if
        else
          if(iroute.eq.2) proute='select mpost;exit'
          if(iroute.eq.3) proute='select gv7236;exit'
          if(iroute.eq.4) proute='select hcposta;exit'
          if(iroute.eq.5) proute='select hpostl;exit'
          if(iroute.eq.6) proute='select hposta3;exit'
          if(iroute.eq.7) proute='select hposta4;exit'
        end if
c
c       call groute(proute)
c       call gshmes('SUP','SUP')
c       call gopen
c
        if(ibcol.ne.0) then
c         call grpsiz(xwndow,ywndow)
c         call rrect(0.,0.,xwndow,ywndow,ibcol,0.)
        end if
        if(ifcol.ne.1) then
c         call gwicol(-1.,ifcol)
c         call gwicol(.001,ifcol)
c         call rtxcol(ifcol,ifcol)
        end if
c
      end if
c
      if(iplot.le.0) write(19,5) -1,ibcol,ifcol
5     format(i2/2i10)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine erase
c
c     erase the screen
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      if(iplot.ge.0) then
        if(ibcol.ne.0) then
c         call rrect(0.,0.,xwndow,ywndow,ibcol,0.)
        else
c         call gclear
        end if
      end if

c
      if(iplot.le.0) write(19,5) -2
5     format(i2)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine plot(x,y,ipen)
c
c     move to the point (x,y); pen up if ipen=3, pen down if ipen=2
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
c     if(iplot.ge.0) call gvect(x*sf,y*sf,3-ipen)
c
      if(iplot.le.0) write(19,5) 1,x,y,ipen
5     format(i2/2e15.5,i10)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine number(x,y,ht,xnum,ang,ndeci)
c
c     plot the floating point number
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      if(iplot.ge.0) then
c       call rtxhei(ht)
c       call rtxang(ang)
c       call rtxn(xnum,ndeci,x*sf,y*sf)
      end if
c
      if (iplot.le.0) write(19,5) 2,x,y,ht,xnum,ang,ndeci
5     format(i2/5e15.5,i10)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine symbol(x,y,ht,label,ang,nchar)
c
c     plot the character string label; special symbols can be plotted
c     with the routine ssymbol
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      character label(nchar)
c
      if(iplot.ge.0) then
c       call rtxhei(ht)
c       call rtxang(ang)
c       call rtx(nchar,label,x*sf,y*sf)
      end if
c
      if(iplot.le.0) write(19,5)
     +  3,nchar,x,y,ht,ang,(label(i),i=1,nchar)
5     format(i2/i10/4e15.5,<nchar>a1)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine line(x,y,npts)
c
c     connect the (x,y) points with straight line segments
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      real x(npts),y(npts)
c
      if(iplot.ge.0) then
        do 10 i=1,npts
           x(i)=x(i)*sf
10         y(i)=y(i)*sf
c
c       call gvect(x,y,npts)
c
        do 20 i=1,npts
          x(i)=x(i)/sf
20        y(i)=y(i)/sf
      end if
c
      if(iplot.le.0) then
        write(19,5) 9,npts
        write(19,15) (x(i),y(i),i=1,npts)
5       format(i2/i10)
15      format(2e15.5)
      end if
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine empty
c
c     flush the graphics buffer
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
c     if(iplot.ge.0) call gempty
c
      if(iplot.le.0) write(19,5) 4
5     format(i2)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine dot(x,y,size,icol)
c
c     plot a dot centred at (x,y) of size isize pixels or mm and colour
c     icol
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      if(iplot.ge.0) then
c       call gwicol(size,icol)
c       call gdot(sf*x,sf*y,1)
c       call gwicol(-1.,1)
c       call gwicol(.001,ifcol)
      end if
c
      if(iplot.le.0) write(19,5) 10,x,y,size,icol
5     format(i2/3e15.5,i10)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine pcolor(icol)
c
c     set the colour for polylines
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
c     if(iplot.ge.0) call gwicol(-1.,icol)
c     if(iplot.ge.0) call gwicol(.001,icol)
c
      if(iplot.le.0) write(19,5) 7,icol
5     format(i2/i10)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine segmnt(iflag)
c
c     open and close Uniras segments
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      if(iplot.ge.0) then
        if(iseg.eq.0) return
c       if(nseg.gt.0) call gsegcl(nseg)
        if(iflag.eq.1) then
          nseg=nseg+1
c         call gsegcr(nseg)
        end if
      end if
c
      if(iplot.le.0) write(19,5) 8,iflag
5     format(i2/i10)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine aldone
c
c     wait unitl the user is ready for the next plot
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      character reply*1
c
      if(iplot.ge.0) then
        write(6,15)
15      format(/'Enter <CR> to continue')
        read(5,25) reply
25      format(a1)
        if(reply(1:1).eq.'s') then
          call plotnd
          stop
        end if
        if(reply(1:1).eq.'0') isep=0
        if(reply(1:1).eq.'1') isep=1
        if(reply(1:1).eq.'2') isep=2
        if(reply(1:1).eq.'3') isep=3
        call segmnt(1)
      end if
c
      if(iplot.le.0) write(19,5) 5
5     format(i2)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine plotnd
c
c     terminate all Uniras plotting
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      if(iplot.ge.0) then
        call segmnt(0)
c       call gclose
      end if
c
      if(iplot.le.0) write(19,5) 6
5     format(i2)
c
      return
      end
