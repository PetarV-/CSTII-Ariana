CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

layers : layers.h layers.cpp
	$(CC) $(CFLAGS) layers.cpp

clean :
	-rm -f *.o *~ layers &> /dev/null

