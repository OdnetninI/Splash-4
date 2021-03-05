TARGET = WATER-NSQUARED
OBJS = bndry.o cnstnt.o cshift.o initia.o interf.o intraf.o kineti.o mdmain.o poteng.o predcor.o syscons.o water.o

include ../../Makefile.config


bndry.o: bndry.c.in split.h mdvar.h parameters.h mddata.h global.h
cnstnt.o: cnstnt.c.in water.h wwpot.h cnst.h frcnst.h fileio.h parameters.h global.h
cshift.o: cshift.c.in water.h global.h
initia.o: initia.c.in split.h mdvar.h water.h cnst.h fileio.h parameters.h mddata.h global.h
interf.o: interf.c.in mdvar.h water.h wwpot.h cnst.h global.h split.h parameters.h mddata.h
intraf.o: intraf.c.in frcnst.h mdvar.h water.h wwpot.h parameters.h mddata.h global.h split.h
kineti.o: kineti.c.in global.h split.h mdvar.h parameters.h mddata.h
mdmain.o: mdmain.c.in  parameters.h mdvar.h water.h wwpot.h cnst.h mddata.h fileio.h frcnst.h randno.h global.h split.h
poteng.o: poteng.c.in global.h split.h mdvar.h frcnst.h water.h wwpot.h parameters.h mddata.h
predcor.o: predcor.c.in split.h mdvar.h parameters.h mddata.h global.h
syscons.o: syscons.c.in parameters.h mdvar.h water.h wwpot.h cnst.h mddata.h global.h
water.o: water.c.in parameters.h mdvar.h water.h wwpot.h cnst.h mddata.h fileio.h frcnst.h randno.h global.h split.h
