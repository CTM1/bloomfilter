#ifndef KMERUTILS_HPP
#define KMERUTILS_HPP

#include <cstdint>
#include <map>
#include <random>
using namespace std;

/**
 * @brief Returns a kmer or it's reverse complement, whichever is smaller
 *        lexicographically.
 * Creates the reverse complement of kmer and compares both using
 * comp_kmer, returns whichever is smaller lexicographically.
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
 * @brief Tests if a kmer is smaller than it's reverse complement
 * Returns a boolean indicating wether the kmer is smaller 
 * lexicographically than it's reverse complement.
 * @param kmer Unsigned integer representing a kmer.
 * @param rev Unsigned integer representing it's reverse complement.
 * @param kmersize Size of kmer.
 * @return true Returned if kmer is lexicographically smaller than rev.
 * @return false Returned if rev is lexicographically smaller than kmer.
 */
bool comp_kmer(uint64_t kmer, uint64_t rev, uint8_t kmersize)

/**
 * @brief Prints kmer as characters to stdout
 * 
 * @param kmer Unsigned integer representing kmer
 * @param kmersize Size of kmer
 */
void print_kmer(uint64_t kmer, uint8_t kmersize);

#endif