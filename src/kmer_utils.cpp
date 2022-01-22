#include "kmer_utils.hpp"

/*
*  000000 00 = A   000000 11 = G 
* -2              -2
*  111111 10 = T   000000 01 = C
*
*  000000 01 = C   000000 10 = T
* -2              -2
*  111111 11 = G   000000 00 = A 
*/ 

/** Compares a kmer and it's reverse complement,
returns wether the kmer is smaller lexicographically
than it's reverse complement */
bool comp_kmer(uint64_t kmer, uint64_t rev, uint8_t kmersize) {
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

// test for performance: https://www.biostars.org/p/113640/
uint64_t choose_kmer_or_rev(uint64_t currkmer, uint8_t kmersize) {
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


void print_kmer(uint64_t kmer, uint8_t kmersize) {
  map<int,char> first;

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

uint64_t random_kmer(uint8_t kmersize) {
    uint64_t kmer = 0;
    int shift = 64 - kmersize * 2;
    
    random_device rd;
    mt19937_64 gen(rd());
    uniform_int_distribution<uint64_t> dis;

    kmer = dis(gen) << shift;
    kmer >>= shift;
    return (kmer);
}