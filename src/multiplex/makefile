CC = g++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

multiplex : multiplex.cpp
	$(CC) $(CFLAGS) multiplex.cpp -o ../../build/multiplex.o

hmm_chain_multiplex : hmm_chain_multiplex.cpp
	$(CC) $(CFLAGS) hmm_chain_multiplex.cpp -o ../../build/hmm_chain_multiplex.o

clean :
	-rm -f ../../build/*.o *~ multiplex &> /dev/null
