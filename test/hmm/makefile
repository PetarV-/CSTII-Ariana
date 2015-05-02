OBJS = ../../build/hmm_tester.o ../../build/hmm.o
CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
OBJFILE = objectives.cpp

hmm_tester : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o ../../bin/hmm_tester

../../build/hmm_tester.o : hmm_tester.cpp
	$(CC) $(CFLAGS) hmm_tester.cpp -o ../../build/hmm_tester.o

../../build/hmm.o : 
	$(MAKE) -C ../../src/hmm hmm

clean :
	-rm -f ../../build/*.o *~ ../../bin/hmm_tester &> /dev/null
	$(MAKE) -C ../../src/hmm clean
