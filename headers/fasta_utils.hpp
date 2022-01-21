#ifndef FASTAUTILS_HPP
#define FASTAUTILS_HPP

#include <cstdint>
#include <iostream>
#include <cstdio>
#include <fstream>

using namespace std;

char next_nucl(ifstream &f);
uint8_t nucltoi(char n);
uint64_t next_kmer(uint64_t currkmer, ifstream &f);

#endif