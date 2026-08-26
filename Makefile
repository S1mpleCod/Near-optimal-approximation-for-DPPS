# dpps_grid -- single translation unit (dpps_grid.cpp includes all headers)
CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++23
HDRS     := dpss_common.hpp layered_graph.hpp grid.hpp peeler.hpp config.hpp coreduce.hpp driver.hpp

dpps_grid: dpps_grid.cpp $(HDRS)
	$(CXX) $(CXXFLAGS) -o $@ dpps_grid.cpp

debug: CXXFLAGS = -O0 -g -std=c++23 -fsanitize=address,undefined
debug: dpps_grid

test: dpps_grid
	./dpps_grid test_queries/graph_1k.txt 0.2 --quiet
	./dpps_grid test_queries/graph_1k.txt 0.2 --demo --quiet

clean:
	rm -f dpps_grid

.PHONY: debug test clean
