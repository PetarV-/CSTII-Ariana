OBJS = nsga2_tester.o ../../src/nsga2/nsga2.o ../../src/nsga2/objectives.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
OBJFILE = objectives.cpp

nsga2_tester : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o nsga2_tester

nsga2_tester.o : nsga2_tester.cpp
	$(CC) $(CFLAGS) nsga2_tester.cpp

../../src/nsga2/nsga2.o ../../src/nsga2/objectives.o :
	$(MAKE) -C ../../src/nsga2 OBJFILE=${OBJFILE}

clean :
	-rm -f *.o *~ nsga2_tester &> /dev/null
	$(MAKE) -C ../../src/nsga2 clean

