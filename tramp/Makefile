#
# makefile for TRAMP
#
EXEDIR=../build/
PLTLIB=../pltlib/
#
# compile flags
#
FFLAGS =  -O -fdollar-ok
CFLAGS =  -O -I/usr/openwin/include
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
TRAMP_OBJS= main.o trc.o atrc.o adjpt.o rngkta.o plt.o misc.o amp.o hdw.o calmod.o blkdat.o $(PLTLIB)/pltsub.o $(PLTLIB)/pltlib.o
XTRAMP_OBJS= main.o trc.o atrc.o adjpt.o rngkta.o plt.o misc.o amp.o hdw.o calmod.o blkdat.o $(PLTLIB)/pltsub.o $(PLTLIB)/xpltlib.o $(PLTLIB)/xbuplot.o
#
# Executable name
#
TRAMP_EXEC=$(EXEDIR)tramp
XTRAMP_EXEC=$(EXEDIR)xtramp
#
tramp: $(TRAMP_EXEC)
$(TRAMP_EXEC): ${TRAMP_OBJS}
	# unilink ${TRAMP_OBJS}
	f77 -o main ${TRAMP_OBJS}
	mv main ${TRAMP_EXEC}
#
xtramp: $(XTRAMP_EXEC)
$(XTRAMP_EXEC): ${XTRAMP_OBJS}
	f77 -o main ${XTRAMP_OBJS} -L/usr/openwin/lib -lX11
	mv main ${XTRAMP_EXEC}
#
# dependencies
#
main.o     :  tramp.par tramp.com
trc.o      :  tramp.par tramp.com
atrc.o     :  tramp.par tramp.com
adjpt.o    :  tramp.par tramp.com
plt.o      :  tramp.par tramp.com
misc.o     :  tramp.par tramp.com
amp.o      :  tramp.par tramp.com
blkdat.o   :  tramp.par tramp.com
hdw.o      :  tramp.par tramp.com
calmod.o   :  tramp.par tramp.com
rngkta.o   :  tramp.par tramp.com

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
	cd $(@D); $(MAKE) $(@F)
FORCE:

