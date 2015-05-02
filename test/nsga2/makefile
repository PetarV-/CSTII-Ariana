OBJS = ../../build/nsga2_tester.o ../../build/nsga2.o ../../build/objectives.o
CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
OBJFILE = objectives.cpp

nsga2_tester : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o ../../bin/nsga2_tester

../../build/nsga2_tester.o : nsga2_tester.cpp
	$(CC) $(CFLAGS) nsga2_tester.cpp -o ../../build/nsga2_tester.o

../../build/nsga2.o : 
	$(MAKE) -C ../../src/nsga2

../../build/objectives.o :
	$(CC) $(CFLAGS) $(OBJFILE) -o ../../build/objectives.o

clean :
	-rm -f ../../build/*.o *~ ../../bin/nsga2_tester &> /dev/null
	$(MAKE) -C ../../src/nsga2 clean
