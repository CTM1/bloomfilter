#include "Bloomfilter.hpp"

using namespace std;

Bloomfilter::Bloomfilter() {
}

Bloomfilter::Bloomfilter(uint32_t n, uint8_t hashes) {
    size = n;
    n_hashes = hashes;
    bit_set = vector<bool> (n);
}

Bloomfilter::Bloomfilter(double rate, uint64_t expected_items) {
    size = needed_bits(rate, expected_items);
    n_hashes = optimal_nhashes(size, expected_items);
    bit_set = vector<bool> (size);
}

void Bloomfilter::add_value(uint64_t x) {
    uint64_t hashes[(int) n_hashes];

    multihash(x, hashes, n_hashes, size);

    
    for (uint64_t k : hashes) {
         bit_set.at(k) = true;
    }
}

bool Bloomfilter::is_present(uint64_t x) {
    uint64_t hashes[(int) n_hashes];
    multihash(x, hashes, n_hashes, size);

    for (uint64_t k : hashes) {
        if (!bit_set.at(k))
            return (false);
    }

    return (true);
}

// With the help of https://hur.st/bloomfilter/

uint32_t Bloomfilter::needed_bits(double f_pos_rate, uint64_t items) {
    return ceil((items * log(f_pos_rate)) / (log(1 / pow(2, log(2)))));
}

uint32_t Bloomfilter::optimal_nhashes(uint32_t bits, uint64_t items) {
    return ((bits / items) * log(2));
}