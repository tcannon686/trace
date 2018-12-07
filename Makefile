
OBJECTS = trace.o matrix.o lodepng.o shade.o hashtable.o
PROGRAM = trace.exe

DEBUG_OPTIONS = -Wall -g
RELEASE_OPTIONS = 
OPTIONS = $(DEBUG_OPTIONS)

$(PROGRAM) : $(OBJECTS)
	gcc $(OBJECTS) -o $(PROGRAM)

trace.o :
	gcc -c trace.c -o trace.o $(OPTIONS)

shade.o :
	gcc -c shade.c -o shade.o $(OPTIONS)

matrix.o : 
	gcc -c matrix.c -o matrix.o $(OPTIONS)

hashtable.o :
	gcc -c hashtable.c -o hashtable.o $(OPTIONS)

lodepng.o :
	gcc -c lodepng.c -o lodepng.o

clean :
	rm $(OBJECTS) $(PROGRAM)


