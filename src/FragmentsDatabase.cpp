#include "FragmentsDatabase.h"

ostream & operator<<(ostream & os, const FragFile & f) {
  os << "FRAG::(" << f.xStart << ", " << f.yStart << ")";
  os <<  " to (" << f.xEnd << ", " << f.yEnd << ") : ";
  os << "SCORE[" << f.score << "], SIM[" << f.similarity << "], ";
  os << "STRAND[" << f.strand << "], LEN[" << f.length <<  "]";
  return os;
}


void endianessConversion(char * data, size_t n) {
  reverse(data, data + n);
}


bool readFragment(FragFile * frag, ifstream & f) {
  static_assert(__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__ or __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__,
      "Undefined endianness");
  f.read((char*) &frag->diag,       sizeof(int64_t));
  f.read((char*) &frag->xStart,     sizeof(uint64_t));
  f.read((char*) &frag->yStart,     sizeof(uint64_t));
  f.read((char*) &frag->xEnd,       sizeof(uint64_t));
  f.read((char*) &frag->yEnd,       sizeof(uint64_t));
  f.read((char*) &frag->length,     sizeof(uint64_t));
  f.read((char*) &frag->ident,      sizeof(uint64_t));
  f.read((char*) &frag->score,      sizeof(uint64_t));
  f.read((char*) &frag->similarity, sizeof(float));
  f.read((char*) &frag->seqX,       sizeof(uint64_t));
  f.read((char*) &frag->seqY,       sizeof(uint64_t));  frag->seqY = 1;
  f.read((char*) &frag->block,      sizeof(int64_t));
  f.read((char*) &frag->strand,     sizeof(char));
  f.read((char*) &frag->evalue,     sizeof(long double));
  if (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__) {
    endianessConversion((char*) (&frag->diag),       sizeof(int64_t));
    endianessConversion((char*) (&frag->xStart),     sizeof(uint64_t));
    endianessConversion((char*) (&frag->yStart),     sizeof(uint64_t));
    endianessConversion((char*) (&frag->xEnd),       sizeof(uint64_t));
    endianessConversion((char*) (&frag->yEnd),       sizeof(uint64_t));
    endianessConversion((char*) (&frag->length),     sizeof(uint64_t));
    endianessConversion((char*) (&frag->ident),      sizeof(uint64_t));
    endianessConversion((char*) (&frag->score),      sizeof(uint64_t));
    endianessConversion((char*) (&frag->similarity), sizeof(float));
    endianessConversion((char*) (&frag->seqX),       sizeof(uint64_t));
    //endianessConversion((char*) (&frag->seqY),       sizeof(uint64_t)); // Already overwritten
    endianessConversion((char*) (&frag->block),      sizeof(int64_t));
    endianessConversion((char*) (&frag->evalue),     sizeof(long double));
  }
  return f.good();
}


FragmentsDatabase::FragmentsDatabase(ifstream & frags_file, ifstream & lengths_file, sequence_manager & seq_manager, uint64_t threshold) {
  //Compute number of fragments in file
  frags_file.seekg(0, ios::end);
  uint64_t total_frags = frags_file.tellg();
  frags_file.seekg(2 * sizeof(uint64_t), ios::beg);

  //Divide by size of frag to get the number of fragments
  //Plus one because it might have padding, thus rounding up to bottom and missing 1 struct
  total_frags = 1 + total_frags / sizeof(FragFile);
  nc = 1 + seq_manager.get_sequence_by_label(0).len / DIV;
  nr = 1 + seq_manager.get_sequence_by_label(1).len / DIV;

  loaded_frags = new vector<FragFile>*[nr];
  if (!loaded_frags) throw runtime_error("Could not allocate memory for fragments!");
  for (size_t i = 0; i < nc; i++) {
    loaded_frags[i] = new vector<FragFile>[nr];
    if (!loaded_frags[i]) throw runtime_error("Could not allocate memory for fragments!");
  }

  //To keep track of current frag
  frags_count = 0;

  FragFile temp_frag;
  while(!frags_file.eof()) {
    if (!readFragment(&temp_frag, frags_file)) break;
    if (temp_frag.seqX != 0 or temp_frag.strand == 'r') continue;
    this->add(temp_frag);
    if(frags_count > total_frags) throw runtime_error("Unexpected number of fragments");
  }
  // Sort each cell of table by length of fragments
  for (size_t i = 0; i < nc; i++) {
    for (size_t j = 0; j < nr; j++) {
      auto comp = [](const FragFile & f1, const FragFile & f2){ return f1.length > f2.length; };
      std::sort(loaded_frags[i][j].begin(), loaded_frags[i][j].end(), comp);
    }
  }
}

FragmentsDatabase::~FragmentsDatabase() {
  //for (size_t i = 0; i < nc; i++) delete[] loaded_frags[i];
  //delete[] loaded_frags;
}
