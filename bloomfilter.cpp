#include <cstdint>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include "hash.h"


using namespace std;

struct Parameters {
    char * filename;
    uint16_t k;  // Size of kmers
    uint64_t n;  // Size of bloom filter
    uint16_t nf; // Number of hashing functions
    uint32_t r;  // Number of requests
};

void help() {
    printf("Usage:\n./bloomfilter <filepath> <kmerSize> <bitvecSize> <hashes> <requests>\nExample:cargo run data/ecoli.fasta 31 456637 3 10000");
    exit(0);
}

/** Returns the next nucleotide in the DNA sequence,
skipping over headers and N chars */
char next_nucl(fstream &f) {
    char c;

    c = f.get();

    if (c == '>') {
        while (f.get() != '\n')
            continue;
        c = f.get();
    }

    while (c == 'N') {
        f.get();
    }

    return (c);
}

/** 
*   A = 1000 00 1 = 0
*   C = 1000 01 1 = 1
*   T = 1010 10 0 = 2
*   G = 1000 11 1 = 3
*/ 

/** Encodes A,C,G,T characters into 2,3,5,4 respectively. */
u_int16_t nucltoi(char n) {
    return ((u_int16_t) ((n >> 1) & 0b11) + 2);
}

uint16_t * next_kmer(uint16_t * currkmer, fstream &f) {
    return (currkmer);
}

uint32_t kmertoi(uint16_t * currkmer, uint16_t kmersize) {
    uint32_t kmerint = 1;

    for (int i = 0; i < kmersize; i++) {
        if (currkmer[i] == 4) {
            kmerint *= 7;
        }
        else {
            kmerint *= currkmer[i];
        }
    }

    return kmerint;
}

int main(int argc, char ** argv) {
    if (argc != 6) help();

    Parameters params = {};
    params.filename = argv[1];
    params.k = atoi(argv[2]);
    params.n = atoi(argv[3]);
    params.nf = atoi(argv[4]);
    params.r = atoi(argv[5]);

    fstream fasta_stream(params.filename);
    uint16_t kmer[params.k];

    for (int i = 0; i < params.k; i++) {
        kmer[i] = nucltoi(next_nucl(fasta_stream));
    }

    printf("%d ", kmertoi(kmer, params.k));
}