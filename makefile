OBJS_SINGLE_CHAIN = build/simple_chain_gmhmm.o build/ariana_single_chain.o
OBJS_MULTIPLEX_CHAIN = build/simple_chain_gmhmm.o build/ariana_multiplex_chain.o build/hmm_chain_multiplex.o build/nsga2.o build/vector_cmp.o
CC = g++ -Iinclude/
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

single_chain : $(OBJS_SINGLE_CHAIN)
	$(CC) $(LFLAGS) $(OBJS_SINGLE_CHAIN) -o bin/ariana_single_chain

multiplex_chain : $(OBJS_MULTIPLEX_CHAIN)
	$(CC) $(LFLAGS) $(OBJS_MULTIPLEX_CHAIN) -o bin/ariana_multiplex_chain

build/ariana_single_chain.o : src/ariana_single_chain.cpp
	$(CC) $(CFLAGS) src/ariana_single_chain.cpp -o build/ariana_single_chain.o

build/ariana_multiplex_chain.o : src/ariana_multiplex_chain.cpp
	$(CC) $(CFLAGS) src/ariana_multiplex_chain.cpp -o build/ariana_multiplex_chain.o

build/simple_chain_gmhmm.o : 
	$(MAKE) -C src/hmm/ simple_chain_gmhmm

build/hmm_chain_multiplex.o : 
	$(MAKE) -C src/multiplex/ hmm_chain_multiplex

build/nsga2.o : 
	$(MAKE) -C src/nsga2/ nsga2

build/vector_cmp.o : 
	$(MAKE) -C src/vector_cmp/ vector_cmp

clean :
	-rm -f build/*.o src/*~ bin/* &> /dev/null
	$(MAKE) -C src/hmm/ clean
	$(MAKE) -C src/layers/ clean
	$(MAKE) -C src/matrix_lib/ clean
	$(MAKE) -C src/multiplex/ clean
	$(MAKE) -C src/nsga2/ clean
	$(MAKE) -C src/vector_cmp/ clean
