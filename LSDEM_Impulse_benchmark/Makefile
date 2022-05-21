
CC = mpicxx
CFLAGS = -std=c++17

CC_INCLUDE = -I/global/home/users/tanpeng/EIGEN_LIB/eigen3

TARGETS = Main

all: $(TARGETS)
Main: MainTest.cpp $(wildcard *.hpp)
	$(CC) $< -O3 -o $@ -Wall $(CC_INCLUDE) $(CFLAGS) -fopenmp

clean:
	rm -f $(TARGETS)

again: clean $(TARGETS)
