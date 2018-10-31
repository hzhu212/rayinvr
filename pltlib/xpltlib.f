c
c     version 1.0  June 1992
c
c     ----------------------------------------------------------------
c     |                                                              |
c     |           ***********  X P L T L I B  ***********            |  
c     |                                                              |
c     |         Plot library to convert Calcomp-like calls           |   
c     |                     to xbuplot calls                         |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |           Modified from pltlib by T. J. Cote (GSC)           |
c     |     Postscript added by D. Demanet (Liege Univ, Belgium)     |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c
c     The routine pcolor is required only if the local graphics allows
c     colour plotting. 
c
c     ----------------------------------------------------------------
c
      subroutine plots (x, y, iroute)
c
c     initialize the plot; this must be the first plot call
c
c     iroute = 1 plots to the screen
c              2 creates a postscript file
c              3 creates a uniras file for the VERSATEC plotter
c              4 creates a postscript file for the colour plotter
c              5 creates a legal size postscript file
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
      common /cps/ ipscol,xcurr,ycurr
c
      integer fill
      real x0, y0
c
      if (iplot .ge. 0) then
        call xbuinit (x, y)
        xwndow = x
        ywndow = y
ccc     if (ibcol .ne. 0) then
          call xbucolour (ibcol)
          x0 = 0.0
          y0 = 0.0
ccc       fill = 1
          fill = ibcol
          call xburect (x0, y0, x, y, fill)
ccc     endif
        call xbucolour (ifcol)
      end if
c
      if (iplot .le. 0) then
c       if (iroute .eq. 4) then
	  ipscol=1
c       else
c         ipscol=0
c       endif
	write(19,1) '%! xpltlib'
1       format(a10)
ccz	write(19,*) 'initgraphics'
	write(19,*) '/initxplt'
	write(19,*) '{90 rotate'
	write(19,*) '50 -600 translate'
	write(19,*) '2.834 2.834 scale'
	write(19,*) '0.2 setlinewidth'
	write(19,*) '/Helvetica findfont 12 scalefont setfont'
	write(19,*) '} def'
	write(19,*) 'initxplt'
      endif
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
      integer fill
      real x0, y0
c 
      if (iplot .ge. 0) then
        call xbuerase
c       if background colour is not black then draw colour rectangle
ccc     if (ibcol .ne. 0) then
          call xbucolour (ibcol)
          x0 = 0.0
          y0 = 0.0
ccc       fill = 1
          fill = ibcol
          call xburect (x0, y0, xwndow, ywndow, fill)
ccc     end if
        call xbucolour (ifcol)
      end if
      
c
      if (iplot .le. 0) then
	write(19,*) 'erasepage'
      endif
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine plot (x, y, ipen)
c
c     move to the point (x,y); pen up if ipen=3, pen down if ipen=2
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
      common /cps/ ipscol,xcurr,ycurr
c
      integer pen
c
      pen = 3 - ipen 
      if (iplot .ge. 0) call xbuplot (x, y, pen)
c
      if (iplot .le. 0) then
	if (ipen .eq. 3) then
          write (19,*)  x, y, ' moveto'
	else
          write (19,*)  xcurr, ycurr, ' moveto'
          write (19,*)  x, y, ' lineto'
          write (19,*)  'stroke'
	endif
      xcurr=x
      ycurr=y
      endif
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine number (x, y, ht, xnum, ang, ndeci)
c
c     plot the floating point number 
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
      character*30 buffor
      integer lon
c 
      if (iplot .ge. 0) then
        call xbunumber (x, y, ht, xnum, ang, ndeci)
      end if
c
      if (iplot .le. 0) then
        write (19,*) x, y, ' moveto'
        write (19,*) 'gsave'
        write (19,*) ang, ' rotate'
        write (19,*) '/Helvetica findfont ',ht,' scalefont setfont'
	if (xnum .eq. 0.) then
	  lon=1
        else
	  if (xnum .gt. 0.) then
	    lon=int(log10(abs(xnum)))+1
          else
	    lon=int(log10(abs(xnum)))+2
           endif
        endif
	if (ndeci .lt. 0) then
	  write (buffor,6) '(a1,i',lon,',a6)'
  6       format(a5,i2,a4)
	  write(19,buffor) '(',int(xnum),') show'
	else
	  lon=lon+ndeci+2
	  write (buffor,5) '(a1,f',lon,'.',ndeci,',a6)'
  5       format(a5,i2,a1,i2,a4)
	  write(19,buffor) '(',xnum,') show'
        endif
        write (19,*) 'grestore'
      endif
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine symbol (x, y, ht, label, ang, nchar)
c
c     plot the character string label
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      character label(nchar)
      character*30 buffor
c
      if (iplot .ge. 0) then
        call xbusymbol (x, y, ht, label, ang, nchar)
      end if
c
      if (iplot .le. 0) then
        write (19,*) x, y, ' moveto'
        write (19,*) 'gsave'
        write (19,*) ang, ' rotate'
        write (19,*) '/Helvetica findfont ',ht,' scalefont setfont'
	write (buffor,5) '(a1,',nchar,'a1,a6)'
  5     format(a4,i3,a6)
	write(19,buffor) '(',(label(i),i=1,nchar),') show'
        write (19,*) 'grestore'
      endif
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine line (x, y, npts)
c
c     connect the (x,y) points with straight line segments
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      real x (npts),y (npts)
c
      if (iplot .ge. 0) call xbulines (x, y, npts)
c
      if (iplot .le. 0) then
        write (19,*) 'newpath'
        write (19,*) x(1), y(1), ' moveto'
	do 5 i=2,npts
  5       write (19,*) x(i), y(i), ' lineto'
        write (19,*) 'stroke'
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
      if (iplot .ge. 0) call xbuflush
c
      if (iplot .le. 0) then
        write (19,*) '% flushbuffer'
      end if
c
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine dot (x, y, size, icol)
c
c     plot a dot centred at (x,y) of size isize pixels or mm and colour
c     icol
c
      integer fill_flag
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      if (iplot.ge.0) then
ccc
c       call xbucolour (icol)
c       call ssymbl(x,y,size/2.,1)
ccc
        call pcolor (icol)
        fill_flag = 1
        call xbucircle (x, y, size, fill_flag)
      end if
c
      if (iplot .le. 0) then
	write(19,*) 'gsave'
	call pcolor(icol)
	write(19,*) 'newpath'
	write(19,*) x, y, size/2, ' 0 360 arc fill'
	write(19,*) 'grestore'
      endif
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine pcolor (icol)
c
c     set the colour for polylines
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
      common /cps/ ipscol,xcurr,ycurr
c 
      if (iplot .ge. 0) call xbucolour (icol)
c
      if (iplot .le. 0) then
	if (ipscol .eq. 1) then
	  if (icol.eq.0) then
	    write(19,*) '1 setgray'
          elseif (icol .eq. 1) then
	    write(19,*) '0 setgray'
          elseif (icol .eq. 2) then
	    write(19,*) '1 0 0 setrgbcolor'
          elseif (icol .eq. 3) then
	    write(19,*) '0 1 0 setrgbcolor'
          elseif (icol .eq. 4) then
	    write(19,*) '0 0 1 setrgbcolor'
          elseif (icol .eq. 5) then
	    write(19,*) '1 1 0 setrgbcolor'
          elseif (icol .eq. 6) then
	    write(19,*) '1 0 1 setrgbcolor'
          elseif (icol .eq. 7) then
	    write(19,*) '0 1 1 setrgbcolor'
          elseif (icol .eq. 8) then
	    write(19,*) '1 0.5 0 setrgbcolor'
          elseif (icol .eq. 9) then
	    write(19,*) '1 0 0.5 setrgbcolor'
          else
	    write(19,*) '0 setgray'
	  endif
	else
	  if (icol.eq.0) then
	    write(19,*) '1 setgray'
          elseif (icol .eq. 1) then
	    write(19,*) '0 setgray'
          else
	    write(19,*) real(icol)/10., ' setgray'
	  endif
	endif
      endif
5     format (i2/i10)
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine segmnt (nsegd)
c
c     open and close Uniras segments
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c 
      if (iplot .ge. 0) then
      end if
c
      if (iplot .le. 0) then
	write(19,*) '% segment'
      endif
c
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
      if (iplot .ge. 0) then
        write (6,15)
15      format (/'Enter <CR> to continue')
        read (5,25) reply
25      format (a1)
        if (reply (1:1) .eq. 's') then
          call plotnd(1)
          stop   
        end if
        if (reply (1:1) .eq. '0') isep = 0
        if (reply (1:1) .eq. '1') isep = 1
        if (reply (1:1) .eq. '2') isep = 2
        if (reply (1:1) .eq. '3') isep = 3
        call segmnt (1)
      end if
c
      if (iplot .le. 0) then
	write(19,*) 'showpage'
	write(19,*) 'initxplt'
      endif
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine plotnd(iflag)
c
c     terminate all Uniras plotting
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      character reply*1
c
      if (iplot .ge. 0) then
        if(iflag.gt.0) then
          write (6,15)
15        format (/'Enter  q  to quit')
101       read (5,25) reply
25        format (a1)
          if (reply .ne. 'q' .and. reply .ne. 'Q') go to 101
        end if
c 
        call xbuclose
      end if
c
      if (iplot .le. 0) then
	write(19,*) 'showpage'
      endif
c
      return
      end
c
c     ----------------------------------------------------------------
c
      subroutine get_event (x, y, button, key)
c
c     get a key or a button press event
c     returns x, y of cursor. button and key pressed
c
ccc
      integer button
      character*1 key
c
      common /cplot/ iplot,isep,iseg,nseg,xwndow,ywndow,ibcol,ifcol,sf
c
      if (iplot.ge.0) then
        call xbuevent (x, y, button, key)
      end if
c
      return
      end
