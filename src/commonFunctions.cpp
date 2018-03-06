#include "commonFunctions.h"

void print_help(){
        cout << "Repkiller v0.9.b\n";
        cout << "Usage: ./repkiller <input_file_path> <output_file_path> <length_ratio> <position_ratio>\n";
        cout << flush;
}

void init_args(const vector<string> & args, ifstream & multifrags, string & out_file_base_path,
               string & path_frags, queue<pair<double, double>> & params){
        out_file_base_path.clear();
        path_frags.clear();

        if (args.size() < 4) throw invalid_argument("Invalid number of arguments.");

        path_frags = args.at(1);
        multifrags.open(path_frags, ifstream::in | ifstream::binary);
        if (!multifrags) throw runtime_error("Could not open input file " + path_frags + ".");

        out_file_base_path = args.at(2);
        if (out_file_base_path.empty()) throw runtime_error("Output file name is missing");

        for (size_t i = 3; i < args.size(); i += 2) {
          double len_ratio = stod(args.at(i),     nullptr);
          double pos_ratio = stod(args.at(i + 1), nullptr);
          if (len_ratio <= 0) throw invalid_argument("Ratio between length and position must be greater than zero");
          if (pos_ratio <= 0) throw invalid_argument("Position proximity must be greater than zero");
          params.push(make_pair(len_ratio, pos_ratio));
        }
}

void terror(const char *s) {
  fprintf(stderr, "[ERROR] %s \n", s);
  exit(-1);
}

void printFragment(const FragFile & f){
  fprintf(stdout, "FRAG::(%" PRIu64 ", %" PRIu64 ") to (%" PRIu64 ", %" PRIu64 ") @%" PRId64 ": [%" PRIu64 "]-[%" PRIu64 "] %c LEN[%" PRIu64 "]\n", f.xStart, f.yStart, f.xEnd, f.yEnd, f.diag, f.seqX, f.seqY, f.strand, f.length);
}

size_t generate_fragment_groups(const FragmentsDatabase & frags_db, FGList & efrags_groups, const sequence_manager & seq_manager, double lensim, double possim) {
  // Create sequence ocupation lists
  const auto seqxlen = seq_manager.get_sequence_by_label(0).len;
  const auto seqylen = seq_manager.get_sequence_by_label(1).len;
  SequenceOcupationList solxf(lensim, possim, seqxlen);
  SequenceOcupationList solyf(lensim, possim, seqylen);
  SequenceOcupationList solxr(lensim, possim, seqxlen);
  SequenceOcupationList solyr(lensim, possim, seqylen);

  // Iterate over all fragments database and map fragments to groups
  for (const auto & fl : frags_db) for (const auto & f : fl) {
    auto & solx = f.strand == 'f' ? solxf : solxr;
    auto & soly = f.strand == 'f' ? solyf : solyr;

    auto agx = solx.get_associated_group(f.xStart + f.length / 2, f.length);
    if (agx != nullptr) {
      // Add to agx and update soly
      agx->push_back(&f);
      soly.insert(f.yStart + f.length / 2, f.length, agx);
      continue;
    }

    auto agy = soly.get_associated_group(f.yStart + f.length / 2, f.length);
    if (agy != nullptr) {
      // Add to agy and update solx
      agy->push_back(&f);
      solx.insert(f.xStart + f.length / 2, f.length, agy);
      continue;
    }

    // Create new group with fragment f, add to solx and soly
    FragsGroup * ngroup = new FragsGroup();
    ngroup->push_back(&f);
    efrags_groups.push_back(ngroup);
    solx.insert(f.xStart + f.length / 2, f.length, ngroup);
    soly.insert(f.yStart + f.length / 2, f.length, ngroup);
  }
  cout << flush;
  return efrags_groups.size();
}

void write_header(ofstream & f, uint64_t sx_len, uint64_t sy_len) {
  f << "CSV file\n";
  f << "[Jul.15 -- < bitlab - Departamento de Arquitectura de Computadores >\n";
  f << "SeqX filename	: DATA1.dat\n";
  f << "SeqY filename	: DATA2.dat\n";
  f << "SeqX name	: S1\n";
  f << "SeqY name	: S2\n";
  f << "SeqX length	: " << sx_len << "\n";
  f << "SeqY length	: " << sy_len << "\n";
  f << "Min.fragment.length	: 0\n";
  f << "Min.Identity	: 0.0\n";
  f << "Total hits	: 0\n";
  f << "Total hits (used)	: 0\n";
  f << "Total fragments	: 0\n";
  f << "Total CSBs:	: 0\n";
  f << "Frag/CSB,xStart,yStart,xEnd,yEnd,strand,block,length,score,ident,similarity,identity,geneX,geneY\n";
  f << "========================================================\n";
}

inline static void store_frag(ofstream & out_file, const FragFile * f, uint64_t gid, unsigned repval) {
  out_file << "Frag," << f->xStart << "," << f->yStart << "," << f->xEnd << "," << f->yEnd << "," << f->strand << "," << gid;
  out_file << "," << f->length << "," << f->score << "," << f->ident << "," << f->similarity << "," << ((float)f->ident * 100 / (float)f->length) << ",0," << repval << "\n";
}

void save_frags_from_group(ofstream & out_file, FragsGroup & fg, uint64_t gid) {
  if (fg.size() == 1) {
    store_frag(out_file, fg.front(), gid, 0);
  } else {
    store_frag(out_file, fg.front(), gid, 1);
    for (auto f_ptr = fg.begin() + 1; f_ptr != fg.end(); f_ptr++) {
      store_frag(out_file, *f_ptr, gid, 2);
    }
  }
}

void save_frag_pair(ofstream & out_file, uint64_t seq1_label, uint64_t seq2_label, const sequence_manager & seq_mngr, const FGList &fgl) {
  const Sequence & seq1 = seq_mngr.get_sequence_by_label(seq1_label), & seq2 = seq_mngr.get_sequence_by_label(seq2_label);

  uint64_t gid = 0;
  //int repetitions;

  //write_header(out_file, seq1.len, seq2.len);
  out_file << read_header;

  for (auto fg : fgl) {
    save_frags_from_group(out_file, *fg, gid++);
  }
}

void save_all_frag_pairs(const string & out_file_base_path, const sequence_manager & seq_manager, const FGList & fgl){
  // Iterators
  uint64_t i, j;
  // Number of sequences involved
  uint64_t n_seq;
  ofstream out_file;
  n_seq = seq_manager.get_number_of_sequences();
  // For each pair of sequences
  for(i=0; i<n_seq; i++) for(j=i+1; j<n_seq; j++) {
    string out_file_name = out_file_base_path;
    out_file.open(out_file_name, ofstream::out);
    if (!out_file) throw runtime_error("Could not open output directory " + out_file_name);
    save_frag_pair(out_file, i, j, seq_manager, fgl);
    out_file.close();
  }
}

void sort_groups(FGList & fgl, const size_t * diag_func) {
  auto comp = [&diag_func](const FragFile * a, const FragFile * b){
                                   size_t ha, hb;
                                   size_t dx;
                                   dx = diag_func[a->xStart / 10];
                                   ha = a->yStart > dx ? a->yStart - dx : dx - a->yStart;
                                   dx = diag_func[b->xStart / 10];
                                   hb = b->yStart > dx ? b->yStart - dx : dx - b->yStart;
                                   return ha < hb;
                           };
  for (auto fg : fgl) if (fg->size() > 1) sort(fg->begin(), fg->end(), comp);
}

void generate_diagonal_func(const FragmentsDatabase & fdb, size_t * diag_func) {
  size_t i = 0;
  for (const auto & fl : fdb) {
    double oh = numeric_limits<double>::infinity();
    if (fl.empty()) {
      diag_func[i] = i == 0 ? 0 : diag_func[i - 1];
      i++;
      continue;
    }
    double nh;
    for (const auto & f : fl) {
      nh = f.length;
      if (nh < oh) diag_func[i] = f.yStart;
    }
    i++;
  }
}
