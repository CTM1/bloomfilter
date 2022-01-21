#include "Bloomfilter.hpp"

using namespace std;

Bloomfilter::Bloomfilter(uint32_t n, uint8_t hashes) {
    size = n;
    n_hashes = hashes;
    bit_set = vector<bool> (n);
}

void Bloomfilter::add_value(uint64_t x) {
    uint64_t hashes[(int) n_hashes];

    multihash(x, hashes, n_hashes, size);

    for (uint64_t k : hashes) {
         bit_set.at(k % size) = true;
    }
}

bool Bloomfilter::is_present(uint64_t x) {
    uint64_t hashes[(int) n_hashes];
    multihash(x, hashes, n_hashes, size);

    for (uint64_t k : hashes) {
        if (!bit_set.at(k % size))
            return (false);
    }

    return (true);
}