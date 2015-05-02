CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

matrix_lib : matrix_lib.cpp
	$(CC) $(CFLAGS) matrix_lib.cpp -o ../../build/matrix_lib.o

clean :
	-rm -f ../../build/*.o *~ matrix_lib &> /dev/null

