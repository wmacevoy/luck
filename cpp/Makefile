CXXFLAGS=-O3 -std=c++20
# CXXFLAGS=-lthreads

all : btluck binomial_luck_test poisson_luck_test geometric_luck_test

btluck : btluck.cpp
	$(CXX) $(CXXFLAGS) -o btluck btluck.cpp $(CXXLIBS)

binomial_luck_test : binomial_luck_test.cpp binomial_luck.h luck.h
	$(CXX) $(CXXFLAGS) -I. -o binomial_luck_test binomial_luck_test.cpp $(CXXLIBS)

poisson_luck_test : poisson_luck_test.cpp poisson_luck.h luck.h
	$(CXX) $(CXXFLAGS) -I. -o poisson_luck_test poisson_luck_test.cpp $(CXXLIBS)

geometric_luck_test : geometric_luck_test.cpp geometric_luck.h luck.h
	$(CXX) $(CXXFLAGS) -I. -o geometric_luck_test geometric_luck_test.cpp $(CXXLIBS)

test : binomial_luck_test poisson_luck_test geometric_luck_test
	./binomial_luck_test
	./poisson_luck_test
	./geometric_luck_test

luck_test : luck_test.cpp luck.h
	$(CXX) $(CXXFLAGS) -I. -o luck_test luck_test.cpp $(CXXLIBS)

clean :
	rm -f btluck binomial_luck_test poisson_luck_test geometric_luck_test
