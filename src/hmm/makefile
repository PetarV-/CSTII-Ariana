OBJS = hmm.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -Wall --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

hmm : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o hmm

hmm.o : hmm.h hmm.cpp
	$(CC) $(CFLAGS) hmm.cpp

clean :
	\rm *.o *~ hmm &> /dev/null

