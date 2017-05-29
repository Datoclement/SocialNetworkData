CC=mpic++

# CC=g++

relation.o: relation.cpp relation.hpp makefile
	$(CC) -c $< -o $@

relation_test.o: relation_test.cpp
	$(CC) -c $< -o $@

test: relation.o relation_test.o
	$(CC) $^ -o $@

run: test
	mpirun -np 2 ./test

clear:
	rm *.o
