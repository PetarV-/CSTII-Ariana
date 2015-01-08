OBJS = nsga2.o objectives.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

nsga2 : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o nsga2

nsga2.o : objectives.h nsga2.cpp
	$(CC) $(CFLAGS) nsga2.cpp

objectives.o : objectives.h
	$(CC) $(CFLAGS) ${OBJFILE} -o objectives.o

clean :
	\rm *.o *~ nsga2 &> /dev/null
	\rm *.o *~ objectives &> /dev/null

