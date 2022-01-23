/**
 * @file kmer_utils.hpp
 * @brief Utilities used to deal with kmers represented as uint64_t's.
 * 
 */

#ifndef KMERUTILS_HPP
#define KMERUTILS_HPP

#include <cstdint>
#include <map>
#include <random>
using namespace std;

/**
 * @brief Returns a kmer or it's reverse complement, whichever is smaller
 *        lexicographically.
 * Creates the reverse complement of kmer and compares both,
 * supposing a {0:A, 1:C, 2:G, 3:T} encoding is used, 
 * it will return whichever one is smaller lexicographically.
 * @param kmer kmer representation.
 * @param kmersize 
 */
uint64_t choose_kmer_or_rev(uint64_t currkmer, uint8_t kmersize);

/**
 * @brief Returns a random uint64_t representing a kmer.
 * 
 * @param kmersize Size of kmer that the uint64_t represents.
 */
uint64_t random_kmer(uint8_t kmersize);

/**
 * @brief Prints kmer as characters to stdout.
 * 
 * @param kmer Unsigned integer representing kmer.
 * @param kmersize Size of kmer.
 */
void print_kmer(uint64_t kmer, uint8_t kmersize);

#endif