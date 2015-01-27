OBJFILE = objectives.cpp
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

nsga2 : objectives objectives.h nsga2.h nsga2.cpp
	$(CC) $(CFLAGS) nsga2.cpp

objectives : objectives.h
	$(CC) $(CFLAGS) ${OBJFILE} -o objectives.o

clean :
	-rm -f *.o *~ nsga2 &> /dev/null
	-rm -f *.o *~ objectives &> /dev/null

