TARGET = BARNES
OBJS = code.o code_io.o load.o grav.o getparam.o util.o

include ../../Makefile.config


stdinc.h: code.h defs.h util.h vectmath.h load.h code_io.h grav.h getparam.h 
code.o: stdinc.h
code_io.o: stdinc.h
getparam.o: stdinc.h
grav.o: stdinc.h
load.o: stdinc.h
util.o: stdinc.h

