CC = clang++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

layers : layers.cpp
	$(CC) $(CFLAGS) layers.cpp -o ../../build/layers.o

clean :
	-rm -f ../../build/*.o *~ layers &> /dev/null

