c
c     version 1.2  Mar 1992
c                 
c     ----------------------------------------------------------------
c     |                                                              |
c     |            ***********  F T E S T  ************              |
c     |                                                              |
c     |            Compute the significance of two given             |
c     |             chi-squared values uisng the F-test              |
c     |                                                              |
c     |                   Written by C. A. Zelt                      |
c     |                                                              |
c     |                Geological Survey of Canada                   |   
c     |                  Ottawa, Canada K1A 0Y3                      |
c     |                                                              |
c     ----------------------------------------------------------------
c                 
c                 
      program main
c
      data pmin/.05/
c
      write(6,*) 'Enter chi-squared and number of points'
      read(5,*) chi1,n1 
      write(6,*) 'Enter chi-squared and number of points'
      read(5,*) chi2,n2 
c     
      if(chi1.gt.chi2) then
        f=chi1/chi2
        df1=n1-1
        df2=n2-1
      else
        f=chi2/chi1
        df1=n2-1
        df2=n1-1
      end if
c
      prob=betai(0.5*df2,0.5*df1,df2/(df2+df1*f))+
     +     (1.-betai(0.5*df1,0.5*df2,df1/(df1+df2/f)))
c
      write(6,5) f,(1.-prob)*100. 
5     format(/'F = ',f8.3//'There is a ',
     + f6.2,'% chance of significantly different variances'//)
c
      if(prob.lt.pmin) then
        write(6,15)
15      format('===>  models are significantly different'//)
      else
        write(6,25)
25      format('===>  models are NOT significantly different'//)
      end if
c
      stop
      end
c
c     ---------------------------------------------------------------------
c
      function betai(a,b,x)
c
c     Returns the incomplete beta function
c
      if(x.lt.0..or.x.gt.1.) then
        write(6,*) '***  bad argument x in betai  ***'
        stop
      end if
c
      if(x.eq.0..or.x.eq.1.) then
        bt=0.
      else
        bt=exp(gammln(a+b)-gammln(a)-gammln(b)+
     +     a*alog(x)+b*alog(1.-x))
      end if
c
      if(x.lt.(a+1.)/(a+b+2.)) then
        betai=bt*betacf(a,b,x)/a
      else    
        betai=1.-bt*betacf(b,a,1.-x)/b
      end if
c
      return
      end
c
c     ---------------------------------------------------------------------
c
      function betacf(a,b,x)
c
c     continued fraction for incomplete beta function used by betai
c
      parameter(itmax=100, eps=3.e-7) 
c
      am=1.
      bm=1.
      az=1.
      qab=a+b
      qap=a+1. 
      qam=a-1.
      bz=1.-qab*x/qap
c
      do 10 m=1,itmax
         em=m
         tem=em+em
         d=em*(b-m)*x/((qam+tem)*(a+tem))
         ap=az+d*am
         bp=bz+d*bm
         d=-(a+em)*(qab+em)*x/((a+tem)*(qap+tem))
         app=ap+d*az
         bpp=bp+d*bz
         aold=az
         am=ap/bpp
         bm=bp/bpp
         az=app/bpp
         bz=1.
         if(abs(az-aold).lt.eps*abs(az)) go to 1
10    continue
c
      write(6,*) '***  a or b too big or itmax too small  ***'
      stop
c
1     betacf=az
c
      return
      end
c
c     ---------------------------------------------------------------------
c
      function gammln(xx)
c
c     Returns the value of ln(gamma function) for xx>1
c
      real*8 cof(6),stp,half,one,fpf,x,tmp,ser
      data cof,stp/76.18009173d0,-86.50532033d0,24.01409822d0,
     +    -1.231739516d0,.120858003d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/
c
      if(xx.lt.1.) then
        write(6,*) '***  argument too small for gammln  ***'
        stop
      end if
c 
      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp 
      ser=one
c
      do 10 j=1,6
         x=x+one
         ser=ser+cof(j)/x
10    continue
c
      gammln=tmp+log(stp*ser)
c
      return
      end
