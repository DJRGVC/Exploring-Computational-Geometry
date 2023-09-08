PLATFORM = $(shell uname)


## Compilation flags
##comment out one or the other 
##debugging 
CFLAGS = -g 
##release
#CFLAGS = -O3 -DNDEBUG
LDFLAGS=

CFLAGS+= -Wall

ifeq ($(PLATFORM),Darwin)
## Mac OS X
CFLAGS += -m64  -Wno-deprecated
INCLUDEPATH=-I/system/usr/local/include 
LDFLAGS+= -m64 -lc -framework AGL -framework OpenGL -framework GLUT -framework Foundation
else
## Linux
CFLAGS += -m64
INCLUDEPATH  = -I/usr/include/GL/ 
LIBPATH = -L/usr/lib64 -L/usr/X11R6/lib
LDFLAGS+=  -lGL -lglut -lrt -lGLU -lX11 -lm  -lXmu -lXext -lXi
endif


CC = g++ -O0 -Wall $(INCLUDEPATH)


PROGS = vis

default: $(PROGS)

vis: vis.o geom.o
	$(CC) -o $@ vis.o geom.o  $(LDFLAGS)

vis.o: vis.cpp  geom.h 
	$(CC) -c $(CFLAGS)   vis.cpp  -o $@

geom.o: geom.cpp geom.h 
	$(CC) -c $(CFLAGS)  geom.cpp -o $@


clean:
	rm *.o
	rm vis


