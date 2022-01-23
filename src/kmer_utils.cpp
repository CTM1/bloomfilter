#include "kmer_utils.hpp"

char char_encoding[4] = {'A', 'C', 'G', 'T'};
random_device rd;
mt19937_64 gen(rd());
uniform_int_distribution<uint64_t> dis;

// This solution inspired from https://www.biostars.org/p/113640/ provides
// an overall performance boost from the previous implementation.
uint64_t choose_kmer_or_rev(uint64_t kmer, uint8_t kmersize) {
    uint64_t res = ~kmer;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    res >>= (2 * (32 - kmersize));

    return (kmer < res ? kmer : res);
}

void print_kmer(uint64_t kmer, uint8_t kmersize) {
  int i = (kmersize - 1) * 2;
  while (i > -1) {
      uint64_t c = kmer;
      printf("%c", char_encoding[c>>i & 0b11]);
      i -= 2;
  }
  printf("\n");
}

uint64_t random_kmer(uint8_t kmersize) {
    uint64_t kmer = 0;
    int shift = 64 - kmersize * 2;

    kmer = dis(gen) << shift;
    kmer >>= shift;
    return (kmer);
}