#pragma once
#pragma pack(1)

#include <inttypes.h>
#include <iostream>
#include <list>
#include <vector>
#include <set>
#include <forward_list>

using namespace std;

#define SEQ_REALLOC 5000000
//#define MAX_READ_SIZE 5000 //Maximum length of read to have a portion of the table allocated
#define MAX_WINDOW_SIZE 70 //Maximum window length to explore NW table
#define SYN_TABLE_REALLOC 10000
#define INIT_SEQS 20
#define INIT_ANNOTS 5000
#define INIT_BLOCKS_TO_ADD 50
#define INIT_REARRANGEMENT_CAPACITY 100
#define WRITE_BUFFER_CAPACITY 100000 //100 KB
#define READLINE 2000
#define READBUF 50000000 //50MB
#define INIT_TRIM_FRAGS 10000
#define TABLE_RATE 100 //hash table lengh divisor
#define INIT_CANDIDATES_ALIGN 100 //For wordbucket list
#define PRINT_RATE 70
#define MAX_LINE 2048

#define MAX_MEM_POOLS 256
#define POOL_SIZE 1024*1024*128 //128 MB

#define NOFRAG 0
#define OPENFRAG 1
#define COVERFRAG 2
#define CLOSEFRAG 3

#define FORWARD 1
#define REVERSE 2
#define MIXED 3

#define INSERTION 1
#define DELETION 2
#define NOTHING 3

#define POINT 4

#define UINT64_T_MAX 0xFFFFFFFFFFFFFFFF
#define SEQUENCE_INDELS_LEN 1000*1000*2 // 1 M

#define IGAP -24
#define EGAP -8

#define EVENTS_WRITE_TIME 5

#define PRINT_BLOCKS 1
#define PRINT_BLOCKS_AND_FRAGS 2

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

// To globally measure RAM usage
extern uint64_t total_bytes_in_use;

//Class prototypes
class memory_pool;
class sequence_manager;

//Struct for FragHits, af2png and leeFrag programs
struct FragFile {
    //Diagonal where the frag is located
    //This value is calculated as:
    //posX - posY
    int64_t diag;
    //Start position in sequence X
    uint64_t xStart;
    //Start position in Sequence Y
    uint64_t yStart;
    //End position in Sequence X
    uint64_t xEnd;
    //End position in Sequence Y
    uint64_t yEnd;
    //Fragment Length
    //For ungaped aligment is:
    //xEnd-xStart+1
    uint64_t length;
    //Number of identities in the
    //fragment
    uint64_t ident;
    //Score of the fragment. This
    //depends on the score matrix
    //used
    uint64_t score;
    //Percentage of similarity. This
    //is calculated as score/scoreMax
    //Where score max is the maximum
    //score possible
    float similarity;
    //sequence number in the 'X' file
    uint64_t seqX;
    //sequence number in the 'Y' file
    uint64_t seqY;
    //synteny block id
    int64_t block;
    //'f' for the forward strain and 'r' for the reverse
    char strand;
    //E-value of fragment
    long double evalue;
};

//Sequence descriptor
typedef struct sequence{
    uint64_t id;    //Label of the sequence
    uint64_t len;   //Length in nucleotides of the sequence
    uint64_t acum;  //Accumulated length from the sequences found before in the file (if any)
    uint32_t coverage; //The percentage of bases covered by fragments of a minimum length
    uint64_t n_frags; //Number of fragments that the sequence had
    char * seq; //DNA sequence
} Sequence;

// A collection of fragments
typedef struct frags_list{
    struct FragFile * f;
    struct frags_list * next;
} Frags_list;

typedef struct linked_list_pos{
    uint64_t pos;
    struct linked_list_pos * next;
} llpos;



//A block that belongs to a genome and that has some synteny level (conserved block)
typedef struct block{
    uint64_t start;     //Starting coordinate
    uint64_t end;       //Ending coordinate
    uint64_t order;     //Order of block according to the genome
    Frags_list * f_list;    //List of fragments that compose it
    Sequence * genome;      //A pointer to the genome to which it belongs
    struct synteny_list * present_in_synteny;   //To tell whether it has already been used in a synteny block
    unsigned char strand_in_synteny;    //The strand that it has at the synteny block
    struct block * prev;
    struct block * next;
    uint64_t id;
} Block;

//A struct to hold arrays of blocks
typedef struct holder{
    Block ** cons_AB;
    Block ** cons_BC;
    uint64_t n_sequences;
} Holder;

//Word struct that identifies a kmer in a sequence
typedef struct word{
    uint64_t hash;
    uint64_t pos;
    char strand;
    Block * b;
} Word;

typedef struct point{
    uint64_t x;
    uint64_t y;
} Point;

//A synteny block is a collection of blocks
typedef struct synteny_block{
    Block * b;
    struct synteny_block * next;
} Synteny_block;

//A synteny list is a collection of synteny blocks
typedef struct synteny_list{
    Synteny_block * sb;
    uint64_t synteny_level;
    struct synteny_list * next;
    struct synteny_list * prev;
    uint64_t id;
} Synteny_list;

/*
class FragsGroup {
private:
  set<FragFile*> groups;
  set<uint64_t> xstarts, ystarts;
public:
  FragsGroup();
  bool sim(FragsGroup & fg, uint64_t proxim);
  void insert(FragFile * f);
  void merge(const FragsGroup & fg);
  void clear() { groups.clear(); xstarts.clear(); ystarts.clear(); };
  set<FragFile*>::iterator begin() { return groups.begin(); }
  set<FragFile*>::iterator end()   { return groups.end(); }
};*/

typedef forward_list<FragFile*> FragsGroup;
typedef forward_list<FragsGroup*> FGList;

//Class for allocating memory only once and requesting particular amounts of bytes
class memory_pool{

private:
    std::vector<char *> * mem_pool;
    std::vector<uint64_t> * base;
    uint64_t current_pool;
    uint64_t max_pools;
    uint64_t pool_size;

public:
    memory_pool(uint64_t pool_size);
    void * request_bytes(uint64_t bytes);
    void reset_n_bytes(uint64_t bytes);
    void reset_to(uint64_t pool, uint64_t position){ this->current_pool = 0; this->base->at(this->current_pool) = 0;}
    void full_reset();
    uint64_t get_pools_used(){ return current_pool; }
    ~memory_pool();
};

//Bucket structs to hold blocks in hash table (1) and words in dictionary (2)
typedef struct bucket {
	Block b;
	struct bucket * next;
} Bucket;

typedef struct wordbucket{
    Word w;
    struct wordbucket * next;
} Wordbucket;

// Sequence class management
class sequence_manager
{
private:
    Sequence * sequences;   //A pointer to the sequences
    uint64_t n_sequences;   //Number of sequences
public:
    sequence_manager();
    uint64_t load_sequences_descriptors(FILE * lengths_file);
    Sequence * get_sequence_by_label(uint64_t label) const;
    uint64_t get_maximum_length() const;
    uint64_t get_number_of_sequences() const { return n_sequences; }
    void print_sequences_data() const;
    ~sequence_manager();
};

// Bucket of linked list for heuristic sorting
struct hs_bucket{
  uint64_t id;
  uint64_t value;
  struct hs_bucket * next;
};

class heuristic_sorted_list{
  private:
    hs_bucket * head = NULL;
    void insert_new(uint64_t id, uint64_t value, hs_bucket * prev);
    void insert_at(uint64_t id, uint64_t value, hs_bucket * prev);
  public:
    heuristic_sorted_list();
    void insert(uint64_t id, uint64_t value);
    uint64_t get_first() { return this->head->id; };
    bool contains(uint64_t id);
    void print();
    void clear();
    ~heuristic_sorted_list();
};

class FragmentsDatabase {
private:
  FragFile * loaded_frags;
  uint64_t num_frags = 0;
public:
  FragmentsDatabase(FILE * frags_file, FILE * lengths_file, sequence_manager & seq_manager);
  FragFile * getFragAt(size_t index) const { return &loaded_frags[index]; }
  uint64_t getTotalFrags() const { return num_frags; }
  ~FragmentsDatabase();
};

class SequenceOcupationList {
private:
  double len_pos_ratio, threshold;
  struct Ocupation {
    uint64_t center;
    uint64_t length;
    FragsGroup * group;
    Ocupation(uint64_t center, uint64_t length, FragsGroup * group) : center(center), length(length), group(group) {};
  };
  std::forward_list<Ocupation> ocupations;
public:
  SequenceOcupationList(double len_pos_ratio, double threshold);
  FragsGroup * get_associated_group(uint64_t center, uint64_t length) const;
  void insert(uint64_t center, uint64_t length, FragsGroup * group) { ocupations.push_front(Ocupation(center, length, group)); };
};
