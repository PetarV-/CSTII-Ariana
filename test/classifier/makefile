OBJS = ../../build/classifier_evaluator.o ../../build/classifier_single_chain.o ../../build/classifier_multiplex_chain.o ../../build/classifier_gmhmm.o ../../build/classifier_multiplex_gmhmm.o ../../build/classifier_generic_single_layer.o ../../build/classifier_generic_multiplex.o ../../build/simple_chain_gmhmm.o ../../build/hmm_chain_multiplex.o ../../build/gmhmm.o ../../build/multiplex_gmhmm.o ../../build/layers.o ../../build/multiplex.o ../../build/nsga2.o ../../build/phi.o ../../build/vector_cmp.o ../../build/matrix_lib.o ../../build/simple_tester.o
CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
LFLAGS = -lpthread -Wall $(DEBUG)

simple_tester : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o ../../bin/simple_tester

../../build/classifier_evaluator.o : classifier_evaluator

../../build/classifier_single_chain.o :
	$(MAKE) -C ../../src/classifier/ classifier_single_chain

../../build/classifier_multiplex_chain.o :
	$(MAKE) -C ../../src/classifier/ classifier_multiplex_chain

../../build/classifier_gmhmm.o :
	$(MAKE) -C ../../src/classifier/ classifier_gmhmm

../../build/classifier_multiplex_gmhmm.o :
	$(MAKE) -C ../../src/classifier/ classifier_multiplex_gmhmm

../../build/classifier_generic_single_layer.o :
	$(MAKE) -C ../../src/classifier/ classifier_generic_single_layer

../../build/classifier_generic_multiplex.o :
	$(MAKE) -C ../../src/classifier/ classifier_generic_multiplex

../../build/simple_chain_gmhmm.o :
	$(MAKE) -C ../../src/hmm/ simple_chain_gmhmm

../../build/hmm_chain_multiplex.o :
	$(MAKE) -C ../../src/multiplex/ hmm_chain_multiplex

../../build/gmhmm.o :
	$(MAKE) -C ../../src/hmm/ gmhmm

../../build/multiplex_gmhmm.o :
	$(MAKE) -C ../../src/multiplex/ multiplex_gmhmm

../../build/layers.o :
	$(MAKE) -C ../../src/layers/ layers

../../build/multiplex.o :
	$(MAKE) -C ../../src/multiplex/ multiplex

../../build/nsga2.o :
	$(MAKE) -C ../../src/nsga2/ nsga2

../../build/phi.o : 
	$(MAKE) -C ../../src/hmm/ phi

../../build/vector_cmp.o : 
	$(MAKE) -C ../../src/vector_cmp/ vector_cmp

../../build/matrix_lib.o :
	$(MAKE) -C ../../src/matrix_lib/ matrix_lib

../../build/simple_tester.o :
	$(CC) $(CFLAGS) simple_tester.cpp -o ../../build/simple_tester.o

classifier_evaluator : classifier_evaluator.cpp
	$(CC) $(CFLAGS) classifier_evaluator.cpp -o ../../build/classifier_evaluator.o

clean :
	-rm -f ../../build/*.o *~ classifier_evaluator &> /dev/null
	-rm -f ../../bin/simple_tester *~ simple_tester &> /dev/null
