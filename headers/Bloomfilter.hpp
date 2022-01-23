/**
 * @file Bloomfilter.hpp
 * @brief Bloomfilter class.
 * 
 */
#ifndef BLOOMFILTER_HPP
#define BLOOMFILTER_HPP

#include <vector>
#include <cmath>
#include "hash.hpp"

using namespace std;
/**
 * @brief Bloomfilter class
 * 
 */
class Bloomfilter {
private:
    uint64_t size;
    uint32_t n_hashes;
    vector<bool> bit_set;
    
    /**
     * @brief Estimates the required size of bit vector in Bloomfilter
     *        according to it's false positive rate and expected items.
     *
     * @param f_pos_rate False positive rate.
     * @param items Expected number of items.
     * @return uint64_t Optimal size of bit vector.
     */
    uint64_t needed_bits(double f_pos_rate, uint64_t items);
    /**
     * @brief Estimates n_hashes in Bloomfilter
     *        according to it's false positive rate and expected items.
     *
     * @param bits False positive rate.
     * @param items Expected number items.
     * @return uint64_t Optimal number of hashes.
     */
    uint64_t optimal_nhashes(uint64_t bits, uint64_t items);
public:
    /**
     * @brief Construct a new empty Bloomfilter object.
     * 
     */
    Bloomfilter();

    /**
     * @brief Construct a new Bloomfilter with it's bit vector size.
     * 
     * @param n The size of the bloomfilter's bit vector in bits.
     * @param hashes The number of hashes for each kmer.
     */
    Bloomfilter(uint64_t n, uint8_t hashes);
    
    /**
     * @brief Construct a new Bloomfilter with it's false positive rate
     *        and expected number of kmers to store in it.
     * @param rate 
     * @param expected_items 
     */
    Bloomfilter(double rate, uint64_t expected_items);

    /**
     * @brief Hashes a kmer n_hashes times and sets corresponding bits
     *        to 1.
     * 
     * @param x kmer represented as an uint64_t.
     */
    void add_value(uint64_t x);

    /**
     * @brief Hashes a kmer n_hashes times and check if its bits in the 
     *        Bloomfilter are set to 1.
     * 
     * @param x kmer represented as an uint64_t.
     */
    bool is_present(uint64_t x);
};

#endif