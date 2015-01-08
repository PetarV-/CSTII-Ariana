OBJS = multiplex.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

multiplex : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o multiplex

multiplex.o : multiplex.h multiplex.cpp
	$(CC) $(CFLAGS) multiplex.cpp

clean :
	\rm *.o *~ multiplex &> /dev/null

