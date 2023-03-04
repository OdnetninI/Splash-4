TARGET = RADIOSITY

OBJS = display.o elemman.o modelman.o patchman.o rad_main.o rad_tools.o room_model.o smallobj.o taskman.o visible.o
HHS = model.h.in parallel.h.in patch.h.in radiosity.h.in task.h.in structs.h.in
HS = model.h parallel.h patch.h radiosity.h task.h structs.h

include ../../Makefile.config


CFLAGS := $(CFLAGS) -I./glibdumb -I./glibps

#CCOPTS = -I -float -I/usr/include $(CFLAGS)
#LDFLAGS = -lm -lpthread
#CFLAGS = $(CINCLUDE) -O4 -march=i686 -mcpu=i686 -malign-double -funroll-loops -finline-functions -fomit-frame-pointer -D_REENTRANT
#CINCLUDE = -I$(GLIBDIR) -I./glibps

LOADLIBES = glibdumb/glib.a glibps/glibps.a

$(TARGET): $(OBJS) $(LOADLIBES)
	$(CC) $(OBJS) $(CFLAGS) $(LOADLIBES) -o $(TARGET) $(LDFLAGS)

radiosity.h: patch.h model.h parallel.h task.h
rad_main.c: $(HS)
patchman.c: $(HS)
elemman.c: $(HS)
taskman.c: $(HS)
modelman.c: $(HS)
smallobj.c: $(HS)
display.c: $(HS)
visible.c: $(HS)
rad_tools.c: $(HS)
room_model.c: $(HS)

glibdumb/glib.a:
	$(MAKE) -C glibdumb

glibps/glibps.a:
	$(MAKE) -C glibps

clean:
	$(MAKE) -C glibdumb clean
	$(MAKE) -C glibps clean
	$(RM) $(OBJS) RADIOSITY
