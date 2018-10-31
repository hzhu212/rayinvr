#
# makefile for RAYINVR
#
EXEDIR=../build/
PLTLIB=../pltlib/
#
# compile flags
#
FFLAGS =  -O -fdollar-ok
#CFLAGS =  -O -I/usr/openwin/include
CFLAGS =  -O -I/usr/openwin/share/include
#
#rules for compiling .f files:
#
.f.o:
	$(COMPILE.f) $(FFLAGS) $< $(INCS)
.c.o:
	$(COMPILE.c) $(CFLAGS) $< $(INCS)
#---------------------------------------------------------------------------
#
# Source objects
#
RAYINVR_OBJS= main.o trc.o atrc.o adjpt.o rngkta.o plt.o misc.o inv.o hdw.o calmod.o blkdat.o $(PLTLIB)/pltsub.o $(PLTLIB)/pltlib.o
XRAYINVR_OBJS= main.o trc.o atrc.o adjpt.o rngkta.o plt.o misc.o inv.o hdw.o calmod.o blkdat.o $(PLTLIB)/pltsub.o $(PLTLIB)/xpltlib.o $(PLTLIB)/xbuplot.o
#
# Executable name
#
RAYINVR_EXEC=$(EXEDIR)rayinvr
XRAYINVR_EXEC=$(EXEDIR)xrayinvr
#
rayinvr: $(RAYINVR_EXEC)
$(RAYINVR_EXEC): ${RAYINVR_OBJS}
	f77 -o main ${RAYINVR_OBJS}
	mv main ${RAYINVR_EXEC}
#	unilink ${RAYINVR_OBJS}
#
xrayinvr: $(XRAYINVR_EXEC)
$(XRAYINVR_EXEC): ${XRAYINVR_OBJS}
	f77 -o main ${XRAYINVR_OBJS} -L/usr/openwin/lib -lX11
	mv main ${XRAYINVR_EXEC}
#
# dependencies
#
main.o     :  rayinvr.par rayinvr.com
trc.o      :  rayinvr.par rayinvr.com
atrc.o     :  rayinvr.par rayinvr.com
adjpt.o    :  rayinvr.par rayinvr.com
plt.o      :  rayinvr.par rayinvr.com
misc.o     :  rayinvr.par rayinvr.com
inv.o      :  rayinvr.par rayinvr.com
blkdat.o   :  rayinvr.par rayinvr.com
hdw.o      :  rayinvr.par rayinvr.com
calmod.o   :  rayinvr.par rayinvr.com
rngkta.o   :  rayinvr.par rayinvr.com

$(PLTLIB)/pltsub.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:

$(PLTLIB)/pltlib.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:

$(PLTLIB)/xpltlib.o: FORCE
	cd $(@D); $(MAKE) $(@F)
FORCE:


$(PLTLIB)/xbuplot.o: FORCE
	cd $(@D);cc -O -I/usr/openwin/share/include -c xbuplot.c
FORCE:

