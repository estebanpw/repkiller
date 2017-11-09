#include "commonFunctions.h"

void print_help(){
        cout << "Repkiller v0.9.a by Carles Bordas\n";
        cout << "Usage: ./repkiller <input_file_path> <output_file_path>.csv <length_ratio> <position_ratio>\n";
        cout << flush;
}

void init_args(const vector<string> & args, ifstream & multifrags, ifstream & lengths_file, ifstream & inf_file,
              ifstream & fasta_file_0, ifstream & fasta_file_1, string & out_file_base_path, string & path_frags,
              queue<pair<double, double>> & params){
        out_file_base_path.clear();
        path_frags.clear();

        if (args.size() < 6) throw invalid_argument("Invalid number of arguments.");

        path_frags = args.at(1);
        multifrags.open(path_frags, ifstream::in | ifstream::binary);
        if (!multifrags) throw runtime_error("Could not open input file " + path_frags + ".");

        auto path_lengths = path_frags + ".lengths";
        lengths_file.open(path_lengths, ifstream::in | ifstream::binary);
        if (!lengths_file) throw runtime_error("Could not find lengths file " + path_lengths + ".");

        auto path_inf = path_frags + ".INF";
        inf_file.open(path_inf, ifstream::in);
        if (!inf_file) throw runtime_error("Could not find information file " + path_inf + ".");

        auto path_fasta_0 = args.at(2);
        fasta_file_0.open(path_fasta_0, ifstream::in);
        if (!fasta_file_0) throw runtime_error("Could not find fasta 1 file " + path_fasta_0 + ".");

        auto path_fasta_1 = args.at(3);
        fasta_file_1.open(path_fasta_1, ifstream::in);
        if (!fasta_file_1) throw runtime_error("Could not find fasta 1 file " + path_fasta_1 + ".");

        out_file_base_path = args.at(4);
        if (out_file_base_path.empty()) throw runtime_error("Output file name is missing");

        for (size_t i = 5; i < args.size(); i += 2) {
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
  fprintf(stdout, "FRAG::(%" PRIu64 ", %" PRIu64 ") to (%" PRIu64 ", %" PRIu64 "): [%" PRIu64 "]-[%" PRIu64 "] %c LEN[%" PRIu64 "]\n", f.xStart, f.yStart, f.xEnd, f.yEnd, f.seqX, f.seqY, f.strand, f.length);
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
  vector<FragFile> ** loaded_frags = frags_db.loaded_frags;
  for (size_t i = 0; i < frags_db.nc; i++) for (size_t j = 0; j < frags_db.nr; j++) for (const auto & f : loaded_frags[i][j]) {
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
  uint64_t gid = 0;
  //int repetitions;

  const auto & seq1 = seq_mngr.get_sequence_by_label(seq1_label);
  const auto & seq2 = seq_mngr.get_sequence_by_label(seq2_label);

  write_header(out_file, seq1.len, seq2.len);
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

/*
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
      if (nh < oh) {
        diag_func[i] = f.yStart;
      }
    }
    i++;
  }
}
*/

forward_list<FragFile> * generateSubFragsFromCell(size_t i, size_t j, const FragmentsDatabase & fdb, uint64_t threshold, const sequence_manager & seq_mngr) {
  cout << i << ", " << j << endl;
  forward_list<FragFile> * nfrags = new forward_list<FragFile>;
  vector<FragFile> & ff = fdb.loaded_frags[i][j];
  size_t nc = fdb.nc;
  size_t nr = fdb.nr;
  cout << "A" << endl;
  for (const auto & f : ff) {
    cout << "F:" << f.xStart << endl;
    // COLUMN WISE |---|
    for (size_t c = f.xStart / DIV; c <= f.xEnd / DIV && c < nc; c++) {
      cout << "c" << endl;
      for (size_t r = 0; r < nr; r++) {
        for (const auto & f2 : fdb.loaded_frags[c][r]) {
          if (f2.xStart >= f.xStart && (f2.xEnd <= f.xEnd || (f.xEnd >= f2.xStart && f.xEnd - f2.xStart > threshold))) {
            uint64_t l = min(f.xEnd, f2.xEnd) - f2.xStart + 1;
            if (l == f.length) continue;
            FragFile nf1;
            nf1.length = l;
            nf1.xStart = f2.xStart;
            nf1.yStart = f.yStart + (f2.xStart - f.xStart);
            nf1.xEnd = nf1.xStart + l - 1;
            nf1.yEnd = nf1.yStart + l - 1;
            nf1.seqX = 0;
            nf1.seqY = 1;
            nf1.strand = 'f';

            /*printf("Intersection of X:\n");
            printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", f.xStart, f.xEnd, f.yStart, f.yEnd, f.length);
            auto pair = seq_mngr.getSequencePair(f);
            std::cout << "\t X:" << pair.first << "\n" << std::flush;
            std::cout << "\t Y:" << pair.second << "\n" << std::flush;
            printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", f2.xStart, f2.xEnd, f2.yStart, f2.yEnd, f2.length);
            pair = seq_mngr.getSequencePair(f2);
            std::cout << "\t X:" << pair.first << "\n" << std::flush;
            std::cout << "\t Y:" << pair.second << "\n" << std::flush;
            printf("\t LEAD TO\n");
            printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", nf1.xStart, nf1.xEnd, nf1.yStart, nf1.yEnd, l);
            pair = seq_mngr.getSequencePair(nf1);
            std::cout << "\t X:" << pair.first << "\n" << std::flush;
            std::cout << "\t Y:" << pair.second << "\n" << std::flush;
            */

            nfrags->push_front(nf1);
            if (f2.xEnd > f.xEnd) {
              FragFile nf2;
              nf2.length = l;
              nf2.xStart = f2.xStart;
              nf2.yStart = f.yStart + (f2.xStart - f.xStart);
              nf2.xEnd = nf2.xStart + l - 1;
              nf2.yEnd = nf2.yStart + l - 1;
              nf2.seqX = 0;
              nf2.seqY = 1;
              nf2.strand = 'f';
              nfrags->push_front(nf2);
              /*printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", nf2.xStart, nf2.xEnd, nf2.yStart, nf2.yEnd, l);
              pair = seq_mngr.getSequencePair(nf2);
              std::cout << "\t X:" << pair.first << "\n" << std::flush;
              std::cout << "\t Y:" << pair.second << "\n" << std::flush;
              */
            }
            //printf("\n");
            //getchar();
          }
        }
      }
    }

    // ROW WISE =
    for (size_t r = f.yStart / DIV; r <= f.yEnd / DIV && r < fdb.nr; r++) {
      for (size_t c = 0; c < fdb.nc; c++) {
        for (const auto & f2 : fdb.loaded_frags[c][r]) {
          if (f2.yStart >= f.yStart && (f2.yEnd <= f.yEnd || (f.yEnd >= f2.yStart && f.yEnd - f2.yStart > threshold))) {
            uint64_t l = min(f.yEnd, f2.yEnd) - f2.yStart + 1;
            if (l == f.length) continue;
            FragFile nf1;
            nf1.length = l;
            nf1.yStart = f2.yStart;
            nf1.xStart = f.xStart + (f2.yStart - f.yStart);
            nf1.yEnd = nf1.yStart + l - 1;
            nf1.xEnd = nf1.xStart + l - 1;
            nf1.seqX = 0;
            nf1.seqY = 1;
            nf1.strand = 'f';

            /*printf("Intersection of X:\n");
            printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", f.xStart, f.xEnd, f.yStart, f.yEnd, f.length);
            auto pair = seq_mngr.getSequencePair(f);
            std::cout << "\t X:" << pair.first << "\n" << std::flush;
            std::cout << "\t Y:" << pair.second << "\n" << std::flush;
            printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", f2.xStart, f2.xEnd, f2.yStart, f2.yEnd, f2.length);
            pair = seq_mngr.getSequencePair(f2);
            std::cout << "\t X:" << pair.first << "\n" << std::flush;
            std::cout << "\t Y:" << pair.second << "\n" << std::flush;
            printf("\t LEAD TO\n");
            printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", nf1.xStart, nf1.xEnd, nf1.yStart, nf1.yEnd, l);
            pair = seq_mngr.getSequencePair(nf1);
            std::cout << "\t X:" << pair.first << "\n" << std::flush;
            std::cout << "\t Y:" << pair.second << "\n" << std::flush;
*/

            nfrags->push_front(nf1);
            if (f2.yEnd > f.yEnd) {
              FragFile nf2;
              nf2.length = l;
              nf2.yStart = f2.yStart;
              nf2.xStart = f.xStart + (f2.yStart - f.yStart);
              nf2.yEnd = nf2.yStart + l - 1;
              nf2.xEnd = nf2.xStart + l - 1;
              nf2.seqX = 0;
              nf2.seqY = 1;
              nf2.strand = 'f';
              nfrags->push_front(nf2);
              /*printf("\tX: (%lu, %lu), Y: (%lu, %lu) L: %lu\n", nf2.xStart, nf2.xEnd, nf2.yStart, nf2.yEnd, l);
              pair = seq_mngr.getSequencePair(nf2);
              std::cout << "\t X:" << pair.first << "\n" << std::flush;
              std::cout << "\t Y:" << pair.second << "\n" << std::flush;
              */
            }
            //printf("\n");
            //getchar();
          }
        }
      }
    }
  }
  return nfrags;
}

void generateSubfragments(FragmentsDatabase & fdb, uint64_t threshold, const sequence_manager & seq_mngr) {
  future<forward_list<FragFile>*> ** futures = new future<forward_list<FragFile>*>*[fdb.nc];
  for (size_t i = 0; i < fdb.nc; i++) {
    futures[i] = new future<forward_list<FragFile>*>[fdb.nr];
  }
  for (size_t i = 0; i < fdb.nc; i++) {
    for (size_t j = 0; j < fdb.nr; j++) {
      futures[i][j] = std::async(std::launch::async | std::launch::deferred, generateSubFragsFromCell, i, j, fdb, threshold, seq_mngr);
      futures[i][j].get();
    }
  }
  for (size_t i = 0; i < fdb.nc; i++) {
    for (size_t j = 0; j < fdb.nr; j++) {
      auto l = futures[i][j].get();
      for (const auto & f : *l) {
        cout << "ADDING" << endl;
        fdb.add(f);
      }
      delete l;
    }
    delete[] futures[i];
  }
  delete[] futures;
}

void divideGroups(FGList & fgl, double threshold) {



}
