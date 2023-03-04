TARGET = VOLREND
OBJS = adaptive.o file.o main.o map.o normal.o octree.o opacity.o option.o raytrace.o render.o view.o libtiff/libtiff.a

include ../../Makefile.config

CFLAGS := $(CFLAGS) -I./libtiff
LDFLAGS := -static -L./libtiff $(LDFLAGS) -ltiff

main.c:	incl.h anl.h
file.c:	incl.h
option.c: incl.h
map.c: incl.h
normal.c: incl.h anl.h address.h
opacity.c: incl.h anl.h
octree.c: incl.h anl.h
view.c:	incl.h
render.c: incl.h
adaptive.c: incl.h anl.h
raytrace.c: incl.h address.h

incl.h:	user_options.h const.h my_types.h global.h macros.h address.h

libtiff/libtiff.a:
	make -C libtiff

clean:
	$(MAKE) -C libtiff clean
	$(RM) $(OBJS) VOLREND
