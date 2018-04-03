CC=gcc
CXX=g++ -std=c++14
CFLAGS=-g -O3 -D_FILE_OFFSET_BITS=64 -Wall
SRC=src
BIN=bin

all: repkiller

repkiller:
	$(CXX) $(CFLAGS) $(SRC)/SaverQueue.cpp $(SRC)/FragmentsDatabase.cpp $(SRC)/SequenceOcupationList.cpp $(SRC)/commonFunctions.cpp -lm $(SRC)/class_structs.cpp $(SRC)/repkiller.cpp -lpthread -o $(BIN)/repkiller

clean:
	rm -rf $(BIN)/repkiller
