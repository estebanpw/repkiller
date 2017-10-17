
#include "commonFunctions.h"

void terror(const char *s) {
    fprintf(stderr, "[ERROR] %s \n", s);
    exit(-1);
}

void generate_fragment_groups(const FragmentsDatabase & frags_db, FGList * efrags_groups, double lensim, double possim) {
  // Create sequence ocupation lists
  SequenceOcupationList solx(lensim, possim), soly(lensim, possim);
  uint64_t n_frags = 0;

  // Iterate over all fragments database and map fragments to groups
  for (uint64_t i = 0; i < frags_db.getTotalFrags(); i++) {
    FragFile * f = frags_db.getFragAt(i);
    ++n_frags;
    auto agx = solx.get_associated_group(f->xStart + f->length / 2, f->length);
    if (agx != nullptr) {
      // Add to agx and update soly
      agx->push_front(f);
      soly.insert(f->yStart + f->length / 2, f->length, agx);
    } else {
      auto agy = soly.get_associated_group(f->yStart + f->length / 2, f->length);
      if (agy != nullptr) {
        // Add to agy and update solx
        agy->push_front(f);
        solx.insert(f->xStart + f->length / 2, f->length, agy);
      } else {
        // Create new group with fragment f, add to solx and soly
        FragsGroup * ngroup = new FragsGroup();
        ngroup->push_front(f);
        efrags_groups->push_front(ngroup);
        soly.insert(f->yStart + f->length / 2, f->length, ngroup);
        solx.insert(f->xStart + f->length / 2, f->length, ngroup);
      }
    }
    print_load(100.0 * n_frags / frags_db.getTotalFrags());
  }
}

void print_memory_usage(){
    //fprintf(stdout, "[INFO] Current RAM usage: %"PRIu64" (MB)\n", total_bytes_in_use/(1024*1024));
    fprintf(stdout, "[INFO] Current RAM usage: %"PRIu64" bytes, which are %"PRIu64" Megabytes. \n", total_bytes_in_use, total_bytes_in_use/(1024*1024));
}

void write_header(FILE * f, uint64_t sx_len, uint64_t sy_len){
    fprintf(f, "CSV file\n");
    fprintf(f, "[Jul.15 -- < bitlab - Departamento de Arquitectura de Computadores >\n");
    fprintf(f, "SeqX filename	: DATA1.dat\n");
    fprintf(f, "SeqY filename	: DATA2.dat\n");
    fprintf(f, "SeqX name	: S1\n");
    fprintf(f, "SeqY name	: S2\n");
    fprintf(f, "SeqX length	: %"PRIu64"\n", sx_len);
    fprintf(f, "SeqY length	: %"PRIu64"\n", sy_len);
    fprintf(f, "Min.fragment.length	: 0\n");
    fprintf(f, "Min.Identity	: 0.0\n");
    fprintf(f, "Total hits	: 0\n");
    fprintf(f, "Total hits (used)	: 0\n");
    fprintf(f, "Total fragments	: 0\n");
    fprintf(f, "Total CSBs:	: 0\n");
    fprintf(f, "Frag/CSB,xStart,yStart,xEnd,yEnd,strand,block,length,score,ident,similarity,identity,geneX,geneY\n");
    fprintf(f, "========================================================\n");
}

void save_frags_from_group(FILE * out_file, FragsGroup & fg, heuristic_sorted_list * hsl, uint64_t gid) {
  uint64_t i, t;
  int v;
  uint64_t heuristic_val;

  i = 0;
  for (auto f : fg) {
    heuristic_val = f->xStart >= f->yStart ? f->xStart - f->yStart : f->yStart - f->xStart;
    hsl->insert(i++, heuristic_val);
  }

  t = hsl->get_first();
  v = i == 1 ? 0 : 1;
  i = 0;
  for (auto f : fg){
    if (v != 0)
      v = i == t ? 1 : 2;
    fprintf(out_file, "Frag,%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%c,", f->xStart, f->yStart, f->xEnd, f->yEnd, f->strand);
    fprintf(out_file, "%"PRIu64",%"PRIu64",%"PRIu64",%"PRIu64",%.2f,%.2f,0,%d\n", gid, f->length, f->score, f->ident, f->similarity, ((float)f->ident * 100 / (float)f->length), v);
    i++;
  }
  hsl->clear();
}

void save_frag_pair(FILE * out_file, uint64_t seq1_label, uint64_t seq2_label, const sequence_manager & seq_mngr, FGList &fgl) {
  Sequence * seq1, * seq2;
  heuristic_sorted_list hsl;
  uint64_t gid = 0;
  //int repetitions;

  seq1 = seq_mngr.get_sequence_by_label(seq1_label);
  seq2 = seq_mngr.get_sequence_by_label(seq2_label);

  write_header(out_file, seq1->len, seq2->len);
  for (auto fg : fgl){
    save_frags_from_group(out_file, *fg, &hsl, gid);
    gid++;
  }
}

void save_all_frag_pairs(char * out_file_base_path, const sequence_manager & seq_manager, FGList & fgl){
  // Iterators
  uint64_t i, j;
  // Path to svc
  char out_file_name[512];
  // Number of sequences involved
  uint64_t n_seq;
  FILE * out_file;

  n_seq = seq_manager.get_number_of_sequences();
  // For each pair of sequences
  for(i=0;i<n_seq;i++) for(j=i+1;j<n_seq;j++){
    sprintf(out_file_name, "%s-%"PRIu64"-%"PRIu64".csv", out_file_base_path, i, j);
    out_file = fopen64(out_file_name, "wt");
    if (out_file == NULL) terror("Could not load output file");
    save_frag_pair(out_file, i, j, seq_manager, fgl);
    fclose(out_file);
  }
}

void print_load(double percentage) {
  size_t i;
  printf("[");
  for (i = 0; i < percentage / 4; i++) printf("|");
  for (i = percentage / 4; i < 25; i++) printf(".");
  printf("] %.2lf%%\r", percentage);
  fflush(stdout);
}
