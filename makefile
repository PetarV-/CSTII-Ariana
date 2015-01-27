OBJS_SINGLE_CHAIN = build/simple_chain_gmhmm.o build/ariana_single_chain.o
CC = g++ -Iinclude/
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

single_chain : $(OBJS_SINGLE_CHAIN)
	$(CC) $(LFLAGS) $(OBJS_SINGLE_CHAIN) -o bin/ariana_single_chain

build/ariana_single_chain.o : src/ariana_single_chain.cpp
	$(CC) $(CFLAGS) src/ariana_single_chain.cpp -o build/ariana_single_chain.o

build/simple_chain_gmhmm.o : 
	$(MAKE) -C src/hmm/ simple_chain_gmhmm

clean :
	-rm -f build/*.o src/*~ bin/ariana_single_chain &> /dev/null
	$(MAKE) -C src/hmm/ clean
