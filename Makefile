
trace.exe : trace.o
	gcc trace.o matrix.o lodepng.o shade.o -o trace.exe

trace.o : matrix.o lodepng.o shade.o
	gcc -c trace.c -o trace.o -Wall -g

shade.o :
	gcc -c shade.c -o shade.o -Wall -g

matrix.o : 
	gcc -c matrix.c -o matrix.o -Wall -g

lodepng.o :
	gcc -c lodepng.c -o lodepng.o

clean :
	rm matrix.o trace.o shade.o trace.exe


