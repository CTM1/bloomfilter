#include "fasta_utils.hpp"

/** Returns the next nucleotide in the DNA sequence,
skipping over headers and N chars */
char next_nucl(ifstream &f) {
    char c;

    c = f.get();

    if (c == 'N') {
        while (c == 'N')
            c = f.get();
    }

    return (c);
}

void skip_line(ifstream &f) {
    char c = f.get();

    while (c != '\n')
        c = f.get();
}

/** 
*   A = 1000 00 1 
*   C = 1000 01 1
*   G = 1000 11 1 
*   T = 1010 10 0
*/

/** Encodes A,C,G,T ASCII characters into 0,1,3,2 respectively.
    works the same for a,c,g,t */
uint8_t nucltoi(char n) {
    return ((n >> 1) & 0b11);
}

uint64_t next_kmer(uint64_t currkmer, ifstream &f, uint8_t kmersize) {
    int shift = (kmersize - 1) * 2;
    currkmer <<= 64 - shift;
    currkmer >>= 62 - shift;
    return(currkmer | nucltoi(next_nucl(f)));
}