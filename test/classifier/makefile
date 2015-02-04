CC = g++ -I../../include
DEBUG = -g
CFLAGS = -std=c++11 -O3 -Wall -Wextra -Werror -Weffc++ -Wstrict-aliasing --pedantic -c $(DEBUG)

classifier_evaluator : classifier_evaluator.cpp
	$(CC) $(CFLAGS) classifier_evaluator.cpp -o ../../build/classifier_evaluator.o

clean :
	-rm -f ../../build/*.o *~ classifier_evaluator &> /dev/null
