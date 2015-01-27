OBJS = multiplex.o ../layers/layers.o ../matrix_lib/matrix_lib.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

multiplex : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o multiplex

multiplex.o : multiplex.h multiplex.cpp
	$(CC) $(CFLAGS) multiplex.cpp

../layers/layers.o :
	$(MAKE) -C ../layers

../matrix_lib/matrix_lib.o :
	$(MAKE) -C ../matrix_lib

clean :
	-rm -f *.o *~ multiplex &> /dev/null
	$(MAKE) -C ../layers clean
	$(MAKE) -C ../matrix_lib clean

