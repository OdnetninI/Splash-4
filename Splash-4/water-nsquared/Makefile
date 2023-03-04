TARGET = WATER-NSQUARED
OBJS = bndry.o cnstnt.o cshift.o initia.o interf.o intraf.o kineti.o mdmain.o poteng.o predcor.o syscons.o water.o

include ../../Makefile.config


bndry.o: split.h mdvar.h parameters.h mddata.h global.h
cnstnt.o: water.h wwpot.h cnst.h frcnst.h fileio.h parameters.h global.h
cshift.o: water.h global.h
initia.o: split.h mdvar.h water.h cnst.h fileio.h parameters.h mddata.h global.h
interf.o: mdvar.h water.h wwpot.h cnst.h global.h split.h parameters.h mddata.h
intraf.o: frcnst.h mdvar.h water.h wwpot.h parameters.h mddata.h global.h split.h
kineti.o: global.h split.h mdvar.h parameters.h mddata.h
mdmain.o: parameters.h mdvar.h water.h wwpot.h cnst.h mddata.h fileio.h frcnst.h randno.h global.h split.h
poteng.o: global.h split.h mdvar.h frcnst.h water.h wwpot.h parameters.h mddata.h
predcor.o: split.h mdvar.h parameters.h mddata.h global.h
syscons.o: parameters.h mdvar.h water.h wwpot.h cnst.h mddata.h global.h
water.o: parameters.h mdvar.h water.h wwpot.h cnst.h mddata.h fileio.h frcnst.h randno.h global.h split.h
