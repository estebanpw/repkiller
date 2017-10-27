#pragma once
#pragma pack(1)

#include <inttypes.h>
#include <iostream>
#include <vector>

using namespace std;

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
typedef struct sequence {
        uint64_t id; //Label of the sequence
        uint64_t len; //Length in nucleotides of the sequence
        uint64_t acum; //Accumulated length from the sequences found before in the file (if any)
        uint32_t coverage; //The percentage of bases covered by fragments of a minimum length
        uint64_t n_frags; //Number of fragments that the sequence had
} Sequence;

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

typedef vector<const FragFile*> FragsGroup;
typedef vector<FragsGroup*> FGList;

// Sequence class management
class sequence_manager
{
  private:
    Sequence * sequences;       //A pointer to the sequences
    uint64_t n_sequences;       //Number of sequences

  public:
    sequence_manager();
    uint64_t load_sequences_descriptors(FILE * lengths_file);
    Sequence * get_sequence_by_label(uint64_t label) const;
    uint64_t get_maximum_length() const;
    uint64_t get_number_of_sequences() const { return n_sequences; }
    void print_sequences_data() const;
    ~sequence_manager();
};
