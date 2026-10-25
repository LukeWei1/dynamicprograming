CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

all: matrix_test matrix_benchmark smith_waterman_test align_salmonella

matrix_test: matrix_test.cpp matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) matrix_test.cpp matrix.cpp -o $@

matrix_benchmark: matrix_benchmark.cpp matrix.cpp matrix.h
	$(CXX) $(CXXFLAGS) matrix_benchmark.cpp matrix.cpp -o $@

smith_waterman_test: smith_waterman_test.cpp smith_waterman.cpp smith_waterman.h
	$(CXX) $(CXXFLAGS) smith_waterman_test.cpp smith_waterman.cpp -o $@

align_salmonella: align_salmonella.cpp smith_waterman.cpp dna_data.cpp smith_waterman.h dna_data.h
	$(CXX) $(CXXFLAGS) align_salmonella.cpp smith_waterman.cpp dna_data.cpp -o $@

test: matrix_test smith_waterman_test
	./matrix_test && ./smith_waterman_test

bench: matrix_benchmark
	./matrix_benchmark

align: align_salmonella
	./align_salmonella

clean:
	rm -f matrix_test matrix_benchmark smith_waterman_test align_salmonella
