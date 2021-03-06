CXX = g++
OBJECTS = Complex.o
FLAGS = -Wextra -Wall -std=c++17
LC_F = --leak-check=full
SPL_Y = --show-possibly-lost=yes
SR_Y = --show-reachable=yes
UVE_Y = --undef-value-errors=yes
TARFILES = timeChecker.cpp Matrix.hpp README Makefile
ARG1 = 500
ARG2 = 50

all: timeChecker
	./timeChecker $(ARG1)
	./timeChecker $(ARG2)

timeChecker: $(OBJECTS)
	$(CXX) $(FLAGS) -c timeChecker.cpp -o TimeChecker.o
	$(CXX) $(FLAGS) $(OBJECTS) TimeChecker.o -o timeChecker

Matrix: Matrix.hpp.gch

Matrix.hpp.gch: Matrix.hpp
	$(CXX) $(FLAGS) Matrix.hpp -o Matrix.hpp.gch

TimeChecker.o: TimeChecker.cpp
	$(CXX) $(FLAGS) -c TimeChecker.cpp

Complex.o: Complex.cpp Complex.h
	$(CXX) $(FLAGS) -c Complex.cpp

clean:
	rm -f *.o timeChecker Matrix Matrix.hpp.gch

tar:
	tar cvf ex3.tar $(TARFILES)

valdbg: timeChecker
	valgrind $(LC_F) $(SPL_Y) $(SR_Y) $(UVE_Y) ./timeChecker $(ARG)
