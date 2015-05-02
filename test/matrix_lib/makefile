OBJS = ../../build/exp_tester.o ../../build/matrix_lib.o
CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

exp_tester : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o ../../bin/exp_tester

../../build/exp_tester.o : exp_tester.cpp
	$(CC) $(CFLAGS) exp_tester.cpp -o ../../build/exp_tester.o

../../build/matrix_lib.o :
	$(MAKE) -C ../../src/matrix_lib

clean :
	-rm -f ../../build/*.o *~ ../../bin/exp_tester &> /dev/null
	$(MAKE) -C ../../src/matrix_lib clean

