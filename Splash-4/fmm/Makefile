TARGET = FMM
OBJS = box.o construct_grid.o cost_zones.o defs.o fmm.o interactions.o memory.o particle.o partition_grid.o

include ../../Makefile.config


defs.c: defs.h
box.c: box.h
construct_grid.c: construct_grid.h
cost_zones.c: cost_zones.h
interactions.c: interactions.h
memory.c: memory.h
particle.c: particle.h
partition_grid.c: partition_grid.h

fmm.o: fmm.c.in defs.h memory.h particle.h box.h partition_grid.h cost_zones.h construct_grid.h interactions.h
interactions.o: interactions.c.in defs.h memory.h particle.h box.h partition_grid.h interactions.h
shell.o: shell.c.in defs.h memory.h particle.h box.h partition_grid.h interactions.h
construct_grid.o: construct_grid.c.in defs.h memory.h particle.h box.h partition_grid.h construct_grid.h
cost_zones.o: cost_zones.c.in defs.h memory.h box.h partition_grid.h cost_zones.h
#orb.o: orb.c.in defs.h memory.h box.h partition_grid.h orb.h
partition_grid.o: partition_grid.c.in defs.h memory.h box.h partition_grid.h
box.o: box.c.in defs.h memory.h particle.h box.h
particle.o: particle.c.in defs.h memory.h particle.h
memory.o: memory.c.in defs.h memory.h
defs.o: defs.c.in defs.h memory.h
memory.h: defs.h particle.h box.h
