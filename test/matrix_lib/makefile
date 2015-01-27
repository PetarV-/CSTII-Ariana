OBJS = exp_tester.o ../../src/matrix_lib/matrix_lib.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

exp_tester : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o exp_tester

exp_tester.o : ../../src/matrix_lib/matrix_lib.h exp_tester.cpp
	$(CC) $(CFLAGS) exp_tester.cpp

../../src/matrix_lib/matrix_lib.o :
	$(MAKE) -C ../../src/matrix_lib

clean :
	-rm -f *.o *~ exp_tester &> /dev/null
	$(MAKE) -C ../../src/matrix_lib clean

