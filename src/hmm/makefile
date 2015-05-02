CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

hmm : hmm.cpp
	$(CC) $(CFLAGS) hmm.cpp -o ../../build/hmm.o

simple_chain_gmhmm : simple_chain_gmhmm.cpp
	$(CC) $(CFLAGS) simple_chain_gmhmm.cpp -o ../../build/simple_chain_gmhmm.o

gmhmm : gmhmm.cpp
	$(CC) $(CFLAGS) gmhmm.cpp -o ../../build/gmhmm.o

phi : phi.cpp
	$(CC) $(CFLAGS) phi.cpp -o ../../build/phi.o

clean :
	-rm -f ../../build/*.o *~ hmm &> /dev/null
