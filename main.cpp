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
#include <map>
#include <random>

using namespace std;


struct Parameters {
    char * filename;
    uint16_t k;  // Size of kmers
    uint64_t n;  // Size of bloom filter
    uint64_t nf; // Number of hashing functions
    uint32_t r;  // Number of requests
};

void help() {
    printf("Usage:\n\t./bloomfilter <filepath> <kmerSize> <bitvecSize> <hashes> <requests>\n\t");
    printf("Example : ./bloomfilter data/ecoli.fasta 31 456637 3 10000");
    exit(0);
}

void printkmer(uint64_t kmer, uint16_t kmersize) {
  std::map<int,char> first;

  first[0b00]='A';
  first[0b01]='C';
  first[0b11]='G';
  first[0b10]='T';

  int i = kmersize * 2 - 2;
  while (i > -1) {
      uint64_t c = kmer;
      printf("%c", first[c>>i & 0b11]);
      i -= 2;
  }
}


uint64_t random_kmer(uint16_t kmersize) {
    uint64_t kmer = 0;

    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<uint64_t> dis;

    return (dis(gen));
}

/** Returns the next nucleotide in the DNA sequence,
skipping over headers and N chars */
char next_nucl(ifstream &f) {
    char c;

    c = f.get();

    if (c == '>') {
        while (f.get() != '\n')
            continue;
        c = f.get();
    }

    if (c == 'N') {
        while (f.get() == 'N') {
            continue;
        }
    }

    return (c);
}

/** 
*   A = 1000 00 1 
*   C = 1000 01 1
*   G = 1000 11 1 
*   T = 1010 10 0
*/

/** Encodes A,C,G,T ASCII characters into 0,1,3,2 respectively.
    works the same for a,c,g,t */
u_int8_t nucltoi(char n) {
    return ((n >> 1) & 0b11);
}

uint64_t next_kmer(uint64_t currkmer, ifstream &f) {
    currkmer <<= 2;
    return(currkmer | nucltoi(next_nucl(f)));
}

/** Compares a kmer and it's reverse complement,
returns wether the kmer is smaller lexicographically
than it's reverse complement */
bool comp_kmer(uint64_t kmer, uint64_t rev, uint16_t kmersize) {
    uint8_t i = (kmersize * 2) - 2;
    uint64_t kmercopy = kmer;
    uint64_t revcopy = rev;

    while (i > 0) {
        kmercopy = kmer;
        revcopy = rev;

        kmercopy = (kmercopy >> i & 0b11);
        revcopy = (revcopy >> i & 0b11);

        if (kmercopy == revcopy) {
            i -= 2;
        }
        else {
            //if kmerG and revT return kmer
            if ((kmercopy) == 0b11 && (revcopy) == 0b10) {
                return (true); //kmer
            }
            //if kmerT and revG return rev
            if ((revcopy) == 0b11 && (kmercopy) == 0b10) {
                return (false); //rev
            }
            else {
                return (kmercopy < revcopy);
            }
        }
    }

    return (true);
}

/*
*  000000 00 = A   000000 11 = G 
* -2              -2
*  111111 10 = T   000000 01 = C
*
*  000000 01 = C   000000 10 = T
* -2              -2
*  111111 11 = G   000000 00 = A 
*/ 

// test for performance: https://www.biostars.org/p/113640/

uint64_t encode_kmer(uint64_t currkmer, uint16_t kmersize) {
    uint64_t rev_kmer = 0;
    uint64_t c = currkmer;
    
    for (int i = kmersize * 2 - 2; i > -1; i -= 2) {
        c = (c >> i) & 0b11;
        c = (c - 2) & 0b11;
        c <<= 62;

        rev_kmer |= c;
        rev_kmer >>= 2;

        c = currkmer;
    }

    rev_kmer >>= 2 * (31 - kmersize);
    
    return (comp_kmer(currkmer, rev_kmer, kmersize) ? currkmer : rev_kmer);
}

int main(int argc, char ** argv) {
    if (argc != 6) help();

    Parameters params = {};
    params.filename = argv[1];
    params.k = atoi(argv[2]);
    params.n = atoi(argv[3]);
    params.nf = atoi(argv[4]);
    params.r = atoi(argv[5]);

    ifstream fasta_stream(params.filename);
    
    if (!fasta_stream.is_open()) { 
        perror("Error opening file"); 
        exit(1);
    }

    uint64_t kmer = 0;
    
    Bloomfilter bf = Bloomfilter(params.n, params.nf);
    
    for (int i = 0; i < params.k; i++) {
        kmer = next_kmer(kmer, fasta_stream);
    }

    while (fasta_stream.peek() != EOF) {
        bf.add_value(encode_kmer(kmer, params.k));

        kmer = next_kmer(kmer, fasta_stream);
    }

    for (int i = 0; i < params.r; i++) {
        uint64_t randkmer = random_kmer(params.k);

        printf("Testing presence for: ");
        
        printkmer(randkmer, params.k);

        uint64_t x = encode_kmer(randkmer, params.k);

        printf("\nIs present: %d\n\n", bf.is_present(x));
    }

    fasta_stream.close();

    return(0);
}