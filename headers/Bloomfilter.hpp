#ifndef BLOOMFILTER_HPP
#define BLOOMFILTER_HPP

#include <vector>
#include <cmath>
#include "hash.hpp"

using namespace std;

class Bloomfilter {
private:
    uint64_t size;
    uint32_t n_hashes;
    vector<bool> bit_set;
    
    uint64_t needed_bits(double f_pos_rate, uint64_t items);
    uint64_t optimal_nhashes(uint64_t bits, uint64_t items);
public:
    Bloomfilter();
    Bloomfilter(uint64_t n, uint8_t hashes);
    Bloomfilter(double rate, uint64_t expected_items);
    void add_value(uint64_t x);
    bool is_present(uint64_t x);
};

#endif