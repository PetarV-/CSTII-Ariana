CC = g++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

classifier_single_chain : classifier_single_chain.cpp
	$(CC) $(CFLAGS) classifier_single_chain.cpp -o ../../build/classifier_single_chain.o

classifier_multiplex_chain : classifier_multiplex_chain.cpp
	$(CC) $(CFLAGS) classifier_multiplex_chain.cpp -o ../../build/classifier_multiplex_chain.o

clean :
	-rm -f ../../build/*.o *~ classifier &> /dev/null
