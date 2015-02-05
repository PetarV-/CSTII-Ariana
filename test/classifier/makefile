OBJS = ../../build/classifier_evaluator.o ../../build/classifier_single_chain.o ../../build/classifier_multiplex_chain.o ../../build/simple_chain_gmhmm.o ../../build/hmm_chain_multiplex.o ../../build/nsga2.o ../../build/vector_cmp.o ../../build/simple_tester.o
CC = g++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

simple_tester : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o ../../bin/simple_tester

../../build/classifier_evaluator.o : classifier_evaluator

../../build/classifier_single_chain.o :
	$(MAKE) -C ../../src/classifier/ classifier_single_chain

../../build/classifier_multiplex_chain.o :
	$(MAKE) -C ../../src/classifier/ classifier_multiplex_chain

../../build/simple_chain_gmhmm.o :
	$(MAKE) -C ../../src/hmm/ simple_chain_gmhmm

../../build/hmm_chain_multiplex.o :
	$(MAKE) -C ../../src/multiplex/ hmm_chain_multiplex

../../build/nsga2.o :
	$(MAKE) -C ../../src/nsga2/ nsga2

../../build/vector_cmp.o : 
	$(MAKE) -C ../../src/vector_cmp/ vector_cmp

../../build/simple_tester.o :
	$(CC) $(CFLAGS) simple_tester.cpp -o ../../build/simple_tester.o

classifier_evaluator : classifier_evaluator.cpp
	$(CC) $(CFLAGS) classifier_evaluator.cpp -o ../../build/classifier_evaluator.o

clean :
	-rm -f ../../build/*.o *~ classifier_evaluator &> /dev/null
	-rm -f ../../bin/simple_tester *~ simple_tester &> /dev/null
