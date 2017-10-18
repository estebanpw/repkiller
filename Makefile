CC=gcc
CXX=g++ -std=c++1y
CFLAGS=-g -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall -Wno-literal-suffix
SRC=src
BIN=bin

all: repkiller

repkiller:
	$(CXX) $(CFLAGS) $(SRC)/commonFunctions.cpp -lm $(SRC)/comparisonFunctions.cpp $(SRC)/class_structs.cpp $(SRC)/repkiller.cpp -lpthread -o $(BIN)/repkiller

clean:
	rm -rf $(BIN)/repkiller
