#ifndef BLOOMFILTER_HPP
#define BLOOMFILTER_HPP

#include <vector>
#include "hash.hpp"

using namespace std;

class Bloomfilter {
    uint32_t size;
    uint8_t n_hashes;
    vector<bool> bit_set;
public:
    Bloomfilter(uint32_t n, uint8_t hashes);
    void add_value(uint64_t x);
    bool is_present(uint64_t x);
};

#endif