#include <cstdint>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <cmath>
#include <string.h>
#include "hash.h"
#include <deque>
#include <iterator>
#include "bloomfilter.cpp"

using namespace std;


struct Parameters {
    char * filename;
    uint16_t k;  // Size of kmers
    uint64_t n;  // Size of bloom filter
    uint64_t nf; // Number of hashing functions
    uint32_t r;  // Number of requests
};

void help() {
    printf("Usage:\n\t./bloomfilter <filepath> <kmerSize> <bitvecSize> <hashes> <requests>\n\tExample : ./bloomfilter data/ecoli.fasta 31 456637 3 10000");
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
*   A = 01000001 & 00011111 = 1
*   C = 01000011 & 00011111 = 3
*   G = 01000111 & 00011111 = 7
*   T = 01010100 & 00011111 = 20
*/ 

/** Encodes A,C,G,T ASCII characters into 1,3,7,20 respectively.
    works the same for a,c,g,t */
u_int8_t nucltoi(char n) {
    return ((u_int8_t) (n & 0b11111));
}

deque<u_int8_t> next_kmer(deque<u_int8_t> currkmer, fstream &f) {
    currkmer.pop_front();
    currkmer.push_back(nucltoi(next_nucl(f)));

    return(currkmer);
}

/** Rolling hash function */
uint64_t hash_kmer(deque<u_int8_t> currkmer, uint16_t kmersize) {
    uint64_t kmerint = 0;
    deque<u_int8_t>::iterator it = currkmer.begin();
    int i = 0;

    while (it != currkmer.end()) {
        kmerint += *it++ * pow(4, i);
        i++;
    }

    return (kmerint);
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

    deque<u_int8_t> kmer;
    deque<u_int8_t>::iterator kmer_it;

    for (int i = 0; i < params.k; i++) {
        kmer.push_back(nucltoi(next_nucl(fasta_stream)));
    }
    
    kmer_it = kmer.begin();

    while (kmer_it != kmer.end()) {
       printf("%d ", *kmer_it++);
    }

    Bloomfilter bf = Bloomfilter(params.n, params.nf);

    u_int64_t hashedkmer = hash_kmer(kmer, params.k);

    printf("%d\n", bf.is_present(hashedkmer));

    bf.add_value(hashedkmer);

    printf("%d", bf.is_present(hashedkmer));
   
}