CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

matrix_lib : matrix_lib.h matrix_lib.cpp
	$(CC) $(CFLAGS) matrix_lib.cpp

clean :
	-rm -f *.o *~ matrix_lib &> /dev/null

