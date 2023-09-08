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


CC = g++ -O3 -Wall $(INCLUDEPATH)


PROGS = thePaintMan

default: $(PROGS)

thePaintMan: thePaintMan.o mondrian.o 
	$(CC) -o $@ thePaintMan.o mondrian.o $(LDFLAGS)

thePaintMan.o: thePaintMan.cpp mondrian.h
	$(CC) -c $(CFLAGS)   thePaintMan.cpp  -o $@

mondrian.o: mondrian.cpp mondrian.h
	$(CC) -c $(CFLAGS)   mondrian.cpp  -o $@

clean:
	rm *.o
	rm thePaintMan
