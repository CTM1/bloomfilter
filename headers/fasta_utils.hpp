#ifndef FASTAUTILS_HPP
#define FASTAUTILS_HPP

#include <cstdint>
#include <iostream>
#include <cstdio>
#include <fstream>

using namespace std;

/**
 * @brief Reads next charactes from file stream and encodes A,C,G,T ASCII characters 
 *        into 0,1,3,2 respectively. Skips over N characters.
 * 
 * @param f File stream to read characters from.
 * @return uint8_t Encoded nucleotide
 */
uint8_t next_nucl(ifstream &f);

/**
 * @brief Returns the next kmer of kmersize in the filestream f
 * Removes two rightmost bits from currkmer witha left shift
 * and adds the next_nucl from the file stream to the two rightmost
 * bits
 * @param currkmer The current kmer representation.
 * @param f File stream.
 * @param kmersize Size of kmer.
 */
uint64_t next_kmer(uint64_t currkmer, ifstream &f, uint8_t kmersize);

/**
 * @brief Skips a line
 * Reads characters from file until a '\n' character is reached.
 * @param f File stream.
 */
void skip_line(ifstream &f);

#endif