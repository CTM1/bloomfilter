#ifndef HASH_HPP
#define HASH_HPP

#include <cstdint>
using namespace std;

/** Generate multiple hashes of the same value using sequences of xorshifts.
 * @param x The value to hash
 * @param hashes An array already allocated to fill with hash values
 * @param nb_hashes The number of hash values needed
 * @param max_val The maximum number that a hashed value can be.
 */
void multihash(uint64_t x, uint64_t * hashes, uint64_t nb_hashes, uint64_t max_val);

#endif