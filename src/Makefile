CC=gcc
CXX=g++ -std=c++11
CFLAGS=-g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -Wno-literal-suffix -DVERBOSE
BIN=../bin

all: repkiller

repkiller: repkiller.cpp
	$(CXX) $(CFLAGS) alignment_functions.cpp -lm commonFunctions.cpp -lm comparisonFunctions.cpp evolutionaryEventsFunctions.cpp class_structs.cpp repkiller.cpp -lpthread -o $(BIN)/repkiller

clean:
	rm -rf $(BIN)/repkiller
