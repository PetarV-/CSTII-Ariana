OBJS_SINGLE_CHAIN = hmm/simple_chain_gmhmm.o ariana_single_chain.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

single_chain : $(OBJS_SINGLE_CHAIN)
	$(CC) $(LFLAGS) $(OBJS_SINGLE_CHAIN) -o ariana_single_chain

ariana_single_chain.o : hmm/hmm.h ariana_single_chain.cpp
	$(CC) $(CFLAGS) ariana_single_chain.cpp

hmm/simple_chain_gmhmm.o : 
	$(MAKE) -C hmm/ simple_chain_gmhmm

clean :
	-rm -f *.o *~ ariana_single_chain &> /dev/null
	$(MAKE) -C hmm/ clean
