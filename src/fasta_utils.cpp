#include "fasta_utils.hpp"
#include <cstdint>

uint8_t encoding[4] = {0b00, 0b01, 0b11, 0b10};

uint8_t next_nucl(ifstream &f) {
    char c;

    c = f.get();

    if (c == 'N') {
        while (c == 'N')
            c = f.get();
    }

    return (encoding[((c >> 1) & 0b11)]);
}

void skip_line(ifstream &f) {
    char c = f.get();

    while (c != '\n')
        c = f.get();
}

uint64_t next_kmer(uint64_t currkmer, ifstream &f, uint8_t kmersize) {
    int shift = (kmersize - 1) * 2;
    currkmer <<= 64 - shift;
    currkmer >>= 62 - shift;
    return(currkmer | next_nucl(f));
}