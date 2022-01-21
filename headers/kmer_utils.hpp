#ifndef KMERUTILS_HPP
#define KMERUTILS_HPP

#include <cstdint>
#include <map>
#include <random>
using namespace std;

uint64_t choose_kmer_or_rev(uint64_t currkmer, uint8_t kmersize);
uint64_t random_kmer(uint8_t kmersize);
void print_kmer(uint64_t kmer, uint8_t kmersize);

#endif