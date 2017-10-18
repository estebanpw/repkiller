#define __STDC_FORMAT_MACROS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"



void endianessConversion(char *source, char *target, int numberOfBytes) {
    int i, j;
    for (i = numberOfBytes - 1; i >= 0; i--) {
        j = numberOfBytes - 1 - i;
        target[j] = source[i];
    }
}


/**
 * Function to read a fragment from the specified file
 */
bool readFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];

    if (htons(1) == 1) {
        //big endian
        if (fread(&frag->diag, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return false;
            terror("Error reading the HSP diagonal");
        }
        if (fread(&frag->xStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        if (fread(&frag->yStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        if (fread(&frag->xEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        if (fread(&frag->yEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        if (fread(&frag->length, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        if (fread(&frag->ident, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identities");
        }
        if (fread(&frag->score, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        if (fread(&frag->similarity, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP similarity");
        }
        if (fread(&frag->seqX, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        if (fread(&frag->seqY, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        frag->seqY = 1;
        if (fread(&frag->block, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        frag->strand = fgetc(f);
        if (fread(&frag->evalue, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
    } else {
        //little endian
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return false;
            terror("Error reading the HSP diagonal");
        }
        endianessConversion(tmpArray, (char *) (&frag->diag), sizeof(int64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->xStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->yStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->xEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->yEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        endianessConversion(tmpArray, (char *) (&frag->length), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identity");
        }
        endianessConversion(tmpArray, (char *) (&frag->ident), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        endianessConversion(tmpArray, (char *) (&frag->score), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP float");
        }
        endianessConversion(tmpArray, (char *) (&frag->similarity), sizeof(float));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqX), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqY), sizeof(uint64_t));
        frag->seqY = 1;
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        endianessConversion(tmpArray, (char *) (&frag->block), sizeof(int64_t));
        frag->strand = fgetc(f);
        if (fread(tmpArray, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
        endianessConversion(tmpArray, (char *) (&frag->evalue), sizeof(long double));
    }
    return true;
}

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];
    if (htons(1) == 1) {
        //Big endian
        fwrite(&frag->diag, sizeof(int64_t), 1, f);
        fwrite(&frag->xStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->yStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->xEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->yEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->length, sizeof(uint64_t), 1, f);
        fwrite(&frag->ident, sizeof(uint64_t), 1, f);
        fwrite(&frag->score, sizeof(uint64_t), 1, f);
        fwrite(&frag->similarity, sizeof(float), 1, f);
        fwrite(&frag->seqX, sizeof(uint64_t), 1, f);
        fwrite(&frag->seqY, sizeof(uint64_t), 1, f);
        fwrite(&frag->block, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        fwrite(&frag->evalue, sizeof(long double), 1, f);
    } else {
        //Little endian
        endianessConversion((char *) (&frag->diag), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        endianessConversion((char *) (&frag->xStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->xEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->length), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->ident), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->score), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->similarity), tmpArray, sizeof(float));
        fwrite(tmpArray, sizeof(float), 1, f);
        endianessConversion((char *) (&frag->seqX), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->seqY), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->block), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        endianessConversion((char *) (&frag->evalue), tmpArray, sizeof(long double));
        fwrite(tmpArray, sizeof(long double), 1, f);
    }
}

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f) {
    char tmpArray[8];
    if (htons(1) == 1) {
        //big endian
        if (fread(length, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading sequence length");
        }
    } else {
        //little endian
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading sequence length");
        }
        endianessConversion(tmpArray, (char *) length, sizeof(uint64_t));
    }
}

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f) {
    char tmpArray[8];
    if (htons(1) == 1) {
        //big endian
        fwrite(length, sizeof(uint64_t), 1, f);
    } else {
        //little endian
        endianessConversion((char *) length, tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
    }
}

void printFragment(struct FragFile * f){

    fprintf(stdout, "FRAG::(%"PRIu64", %"PRIu64") to (%"PRIu64", %"PRIu64"): [%"PRIu64"]-[%"PRIu64"] %c LEN[%"PRIu64"]\n", f->xStart, f->yStart, f->xEnd, f->yEnd, f->seqX, f->seqY, f->strand, f->length);
}
