#include "FragmentsDatabase.h"

ostream & operator<<(ostream & os, const FragFile & f) {
        os << "FRAG::(" << f.xStart << ", " << f.yStart << ")";
        os <<  " to (" << f.xEnd << ", " << f.yEnd << ") : ";
        os << "STRAND[" << f.strand << "], LEN[" << f.length <<  "]";
        return os;
}


void endianessConversion(char *source, char *target, int numberOfBytes) {
        int i, j;
        for (i = numberOfBytes - 1; i >= 0; i--) {
                j = numberOfBytes - 1 - i;
                target[j] = source[i];
        }
}


bool FragmentsDatabase::readFragment(FragFile * frag, FILE *f) {
        if (htons(1) == 1) {
                //big endian
                if (fread(&frag->diag, sizeof(int64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->xStart, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->yStart, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->xEnd, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->yEnd, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->length, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->ident, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->score, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->similarity, sizeof(float), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->seqX, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                if (fread(&frag->seqY, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                frag->seqY = 1;
                if (fread(&frag->block, sizeof(int64_t), 1, f) != 1) {
                        return false;
                }
                frag->strand = fgetc(f);
                if (fread(&frag->evalue, sizeof(long double), 1, f) != 1) {
                        return false;
                }
        } else {
                char tmpArray[sizeof(long double)];
                //little endian
                if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->diag), sizeof(int64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->xStart), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->yStart), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->xEnd), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->yEnd), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->length), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->ident), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->score), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(float), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->similarity), sizeof(float));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->seqX), sizeof(uint64_t));
                if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->seqY), sizeof(uint64_t));
                frag->seqY = 1;
                if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->block), sizeof(int64_t));
                frag->strand = fgetc(f);
                if (fread(tmpArray, sizeof(long double), 1, f) != 1) {
                        return false;
                }
                endianessConversion(tmpArray, (char *) (&frag->evalue), sizeof(long double));
        }
        return true;
}


FragmentsDatabase::FragmentsDatabase(FILE * frags_file, FILE * lengths_file, sequence_manager & seq_manager) {
        uint64_t total_frags;

        (void) seq_manager.load_sequences_descriptors(lengths_file);

        //Compute number of fragments in file
        fseeko(frags_file, 0L, SEEK_END);
        total_frags = ftello(frags_file) - 2*(sizeof(uint64_t)); //Remove the headers
        fseeko(frags_file, 0L, SEEK_SET);

        //Divide by size of frag to get the number of fragments
        //Plus one because it might have padding, thus rounding up to bottom and missing 1 struct
        total_frags = 1 + total_frags / sizeof(FragFile);
        //Allocate memory for all frags
        loaded_frags = (FragFile *)malloc(total_frags * sizeof(FragFile));
        if(loaded_frags == nullptr) throw runtime_error("Coult not allocate memory for all fragments");

        //Skip headers
        fseeko(frags_file, 16, SEEK_SET);
        //To keep track of current frag
        num_frags = 0;

        while(!feof(frags_file)) {
                if (!readFragment(&loaded_frags[num_frags], frags_file)) break;
                if (loaded_frags[num_frags].seqX != 0) continue;
                ++num_frags;
                if(num_frags > total_frags) throw runtime_error("Unexpected number of fragments");
        }
}

FragmentsDatabase::~FragmentsDatabase() {
        free(loaded_frags);
}
