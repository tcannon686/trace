
OBJECTS = trace.o matrix.o lodepng.o shade.o hashtable.o linkedlist.o texture.o cmd.o xwin.o halton.o
PROGRAM = trace.exe

DEBUG_OPTIONS = -Wall -g
RELEASE_OPTIONS = 

OPTIONS = $(DEBUG_OPTIONS)
LIBS = -lm -lpthread

CC = gcc

ifeq ($(INCLUDE_GUI), true) 
LIBS += -lX11
OPTIONS += -DINCLUDE_GUI
endif

all : $(PROGRAM)
$(PROGRAM) : $(OBJECTS)
	$(CC) $(OBJECTS) -o $(PROGRAM) $(LIBS)

trace.o : trace.c trace.h
	$(CC) -c trace.c -o trace.o $(OPTIONS)

shade.o : shade.c shade.h
	$(CC) -c shade.c -o shade.o $(OPTIONS)

matrix.o : matrix.c matrix.h
	$(CC) -c matrix.c -o matrix.o $(OPTIONS)

hashtable.o : hashtable.c hashtable.h
	$(CC) -c hashtable.c -o hashtable.o $(OPTIONS)

xwin.o : xwin.c xwin.h
	$(CC) -c xwin.c -o xwin.o $(OPTIONS)

linkedlist.o : linkedlist.c linkedlist.h
	$(CC) -c linkedlist.c -o linkedlist.o $(OPTIONS)

texture.o : texture.c texture.h
	$(CC) -c texture.c -o texture.o $(OPTIONS)

cmd.o : cmd.c cmd.h
	$(CC) -c cmd.c -o cmd.o $(OPTIONS)

halton.o : halton.c halton.h
	$(CC) -c halton.c -o halton.o $(OPTIONS)

lodepng.o : lodepng.c lodepng.h
	$(CC) -c lodepng.c -o lodepng.o

clean :
	rm $(OBJECTS) $(PROGRAM)

