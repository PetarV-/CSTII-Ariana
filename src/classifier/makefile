CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

classifier_single_chain : classifier_single_chain.cpp
	$(CC) $(CFLAGS) classifier_single_chain.cpp -o ../../build/classifier_single_chain.o

classifier_multiplex_chain : classifier_multiplex_chain.cpp
	$(CC) $(CFLAGS) classifier_multiplex_chain.cpp -o ../../build/classifier_multiplex_chain.o

classifier_gmhmm : classifier_gmhmm.cpp
	$(CC) $(CFLAGS) classifier_gmhmm.cpp -o ../../build/classifier_gmhmm.o

classifier_multiplex_gmhmm : classifier_multiplex_gmhmm.cpp
	$(CC) $(CFLAGS) classifier_multiplex_gmhmm.cpp -o ../../build/classifier_multiplex_gmhmm.o

classifier_generic_single_layer : classifier_generic_single_layer.cpp
	$(CC) $(CFLAGS) classifier_generic_single_layer.cpp -o ../../build/classifier_generic_single_layer.o

classifier_generic_multiplex : classifier_generic_multiplex.cpp
	$(CC) $(CFLAGS) classifier_generic_multiplex.cpp -o ../../build/classifier_generic_multiplex.o

clean :
	-rm -f ../../build/*.o *~ classifier &> /dev/null
