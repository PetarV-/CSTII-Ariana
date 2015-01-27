CC = g++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

multiplex : multiplex.cpp
	$(CC) $(CFLAGS) multiplex.cpp -o ../../build/multiplex.o

clean :
	-rm -f ../../build/*.o *~ multiplex &> /dev/null


