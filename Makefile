
OBJECTS = trace.o matrix.o lodepng.o shade.o hashtable.o
PROGRAM = trace.exe

DEBUG_OPTIONS = -Wall -g
RELEASE_OPTIONS = 

OPTIONS = $(DEBUG_OPTIONS)
LIBS = -lm -lpthread

$(PROGRAM) : $(OBJECTS)
	gcc $(OBJECTS) -o $(PROGRAM) $(LIBS)

trace.o : trace.c trace.h
	gcc -c trace.c -o trace.o $(OPTIONS)

shade.o : shade.c shade.h
	gcc -c shade.c -o shade.o $(OPTIONS)

matrix.o : matrix.c matrix.h
	gcc -c matrix.c -o matrix.o $(OPTIONS)

hashtable.o : hashtable.c hashtable.h
	gcc -c hashtable.c -o hashtable.o $(OPTIONS)

lodepng.o : lodepng.c lodepng.h
	gcc -c lodepng.c -o lodepng.o

clean :
	rm $(OBJECTS) $(PROGRAM)


