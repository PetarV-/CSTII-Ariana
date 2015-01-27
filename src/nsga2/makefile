OBJFILE = objectives.cpp
CC = g++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

nsga2 : objectives nsga2.cpp
	$(CC) $(CFLAGS) nsga2.cpp -o ../../build/nsga2.o

objectives : 
	$(CC) $(CFLAGS) ${OBJFILE} -o ../../build/objectives.o

clean :
	-rm -f ../../build/*.o *~ nsga2 &> /dev/null
	-rm -f ../../build/*.o *~ objectives &> /dev/null

