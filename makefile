CC = g++
CFLAGS = -I./boost -I. -O7 -g

all :
	$(CC) Driver.cpp NetworkSimplifier.cpp Equation.cpp Solution.cpp -o Program $(CFLAGS)