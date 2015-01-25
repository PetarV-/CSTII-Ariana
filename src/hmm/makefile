CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

hmm : hmm.h hmm.cpp
	$(CC) $(CFLAGS) hmm.cpp

simple_chain_gmhmm : hmm.h simple_chain_gmhmm.cpp
	$(CC) $(CFLAGS) simple_chain_gmhmm.cpp

tdhmm : hmm.h tdhmm.cpp
	$(CC) $(CFLAGS) tdhmm.cpp

clean :
	-rm -f *.o *~ hmm &> /dev/null
