OBJS = syn_gen.o
CC = g++
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

syn_gen : $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o syn_gen

syn_gen.o : syn_gen.cpp
	$(CC) $(CFLAGS) syn_gen.cpp

clean :
	-rm -f *.o *~ syn_gen &> /dev/null

