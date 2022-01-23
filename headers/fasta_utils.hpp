#ifndef FASTAUTILS_HPP
#define FASTAUTILS_HPP

#include <cstdint>
#include <iostream>
#include <cstdio>
#include <fstream>

using namespace std;

uint8_t next_nucl(ifstream &f);
uint64_t next_kmer(uint64_t currkmer, ifstream &f, uint8_t kmersize);
void skip_line(ifstream &f);

#endif