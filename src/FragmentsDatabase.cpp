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


bool readFragment(FragFile * frag, const string & line) {
  try {
    vector<string> vect;
    istringstream s(line);
    string buf;

    for (size_t i = 0; i < 14; ++i) {
        getline(s, buf, ',');
        if (buf.empty()) return false;
        vect.push_back(buf);
    }

    if (vect[0] != "Frag") return false;
    frag->xStart = atoll(vect[1].c_str());
    frag->yStart = atoll(vect[2].c_str());
    frag->diag = (int64_t)frag->xStart - (int64_t)frag->yStart;
    frag->xEnd = atoll(vect[3].c_str());
    frag->yEnd = atoll(vect[4].c_str());
    frag->strand = vect[5][0];
    frag->block = atoll(vect[6].c_str());
    frag->length = atoll(vect[7].c_str());
    frag->score = atoll(vect[8].c_str());
    frag->ident = (uint64_t) stof(vect[10]);;
    frag->similarity = stof(vect[10]);
    frag->seqX = 0;
    frag->seqY = 1;
    frag->evalue = 0;

    return true;
  } catch (...) {
    return false;
  }

}

string final_read_header;

FragmentsDatabase::FragmentsDatabase(ifstream & frags_file, sequence_manager & seq_manager) {
  string line;

  for (size_t i = 0; i < 7; ++i) {
    getline(frags_file, line);
    header.append(line).append("\n");
  }

  seq_manager.sequences.emplace_back(0, atoll(line.substr(line.find(':') + 1).c_str()) + 1);
  getline(frags_file, line);
  header.append(line).append("\n");
  seq_manager.sequences.emplace_back(1, atoll(line.substr(line.find(':') + 1).c_str()) + 1);

  for (size_t i = 0; i < 5; ++i) {
    getline(frags_file, line);
    header.append(line).append("\n");
  }

  uint64_t total_frags = atoll(line.substr(line.find(':') + 1).c_str());

  for (size_t i = 0; i < 4; ++i) {
    getline(frags_file, line);
    header.append(line).append("\n");
  }

  seq_manager.read_header(header);


  //Divide by size of frag to get the number of fragments
  //Plus one because it might have padding, thus rounding up to bottom and missing 1 struct
  vsize = 1 + seq_manager.get_sequence_by_label(0).len / 10;
  loaded_frags = unique_ptr<vector<FragFile>[]>(new vector<FragFile>[vsize]);
  if (!loaded_frags) throw runtime_error("Could not allocate memory for fragments!");

  //To keep track of current frag
  frags_count = 0;

  FragFile temp_frag;
  while(!frags_file.eof()) {
    getline(frags_file, line);
    if (!readFragment(&temp_frag, line)) continue;
    if (temp_frag.seqX != 0) continue;
    size_t index = temp_frag.xStart / 10;
    loaded_frags[index].push_back(temp_frag);
    ++frags_count;
    if(frags_count > total_frags) throw runtime_error("Unexpected number of fragments");
  }
}
