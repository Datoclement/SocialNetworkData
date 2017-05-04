relation.o: relation.cpp relation.hpp
	g++ -c relation.cpp

relation_test.o: relation_test.cpp relation.hpp
	g++ relation.hpp -c relation_test.cpp

relation_test: relation_test.o relation.o
	g++ relation.o relation_test.o -o relation_test
