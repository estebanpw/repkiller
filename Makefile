CC=gcc
CXX=g++ -std=c++11
CFLAGS=-g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -Wno-literal-suffix -DVERBOSE
SRC=src
BIN=bin

all: repkiller

repkiller:
	$(CXX) $(CFLAGS) $(SRC)/alignment_functions.cpp -lm $(SRC)/commonFunctions.cpp -lm $(SRC)/comparisonFunctions.cpp $(SRC)/evolutionaryEventsFunctions.cpp $(SRC)/class_structs.cpp $(SRC)/repkiller.cpp -lpthread -o $(BIN)/repkiller

clean:
	rm -rf $(BIN)/repkiller
