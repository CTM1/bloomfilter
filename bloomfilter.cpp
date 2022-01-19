#include <bitset>
#include <sys/types.h>
#include <vector>
#include "hash.cpp"

using namespace std;

class Bloomfilter {
private:
    u_int64_t size;
    u_int32_t n_hashes;
    vector<bool> bit_set;
public:
    Bloomfilter(uint64_t n, u_int32_t hashes) {
        size = n;
        n_hashes = hashes;
        bit_set = vector<bool> (n);
    }

    void add_value(u_int64_t x) {
        u_int64_t hashes[(int) n_hashes];

        multihash(x, hashes, n_hashes, size);

        for (u_int64_t k : hashes) {
            bit_set.at(k % size) = true;
        }
    }

    bool is_present(u_int64_t x) {
        u_int64_t hashes[(int) n_hashes];

         multihash(x, hashes, n_hashes, size);

         for (u_int64_t k : hashes) {
             if (!bit_set.at(k))
                return (false);
         }

         return (true);
    }
};