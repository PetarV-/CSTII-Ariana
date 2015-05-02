CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

vector_cmp : vector_cmp.cpp
	$(CC) $(CFLAGS) vector_cmp.cpp -o ../../build/vector_cmp.o

clean :
	-rm -f ../../build/*.o *~ vector_cmp &> /dev/null
