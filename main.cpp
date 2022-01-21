#include <cstdint>
#include <cstdio>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <cmath>
#include <string.h>
#include "hash.h"
#include <deque>
#include <iterator>
#include "bloomfilter.cpp"

using namespace std;


struct Parameters {
    char * filename;
    uint16_t k;  // Size of kmers
    uint64_t n;  // Size of bloom filter
    uint64_t nf; // Number of hashing functions
    uint32_t r;  // Number of requests
};

void help() {
    printf("Usage:\n\t./bloomfilter <filepath> <kmerSize> <bitvecSize> <hashes> <requests>\n\t");
    printf("Example : ./bloomfilter data/ecoli.fasta 31 456637 3 10000");
    exit(0);
}

deque<u_int8_t> random_kmer(uint16_t kmersize) {
    deque<u_int8_t> kmer;

    for (int i = 0; i < kmersize; i++) {
        kmer[i] = rand() % 4;
    }

    return (kmer);
}
/** Returns the next nucleotide in the DNA sequence,
skipping over headers and N chars */
char next_nucl(ifstream &f) {
    char c;

    c = f.get();

    if (c == '>') {
        while (f.get() != '\n')
            continue;
        c = f.get();
    }

    if (c == 'N') {
        while (f.get() == 'N') {
            continue;
        }
    }

    return (c);
}

/** 
*   A = 1000 00 1 
*   C = 1000 01 1
*   G = 1000 11 1 
*   T = 1010 10 0
*/

/** Encodes A,C,G,T ASCII characters into 0,1,3,2 respectively.
    works the same for a,c,g,t */
u_int8_t nucltoi(char n) {
    return ((u_int8_t) (n >> 1 & 0b11));
}

deque<u_int8_t> next_kmer(deque<u_int8_t> currkmer, ifstream &f) {
    currkmer.pop_front();
    currkmer.push_back(nucltoi(next_nucl(f)));

    return(currkmer);
}

/** Compares a kmer and it's reverse complement,
returns wether the kmer is smaller lexicographically
than it's complement*/
bool comp_kmer(uint64_t kmer, uint64_t rev, uint16_t kmersize) {
    uint8_t i = (kmersize * 2) - 2;
    uint64_t kmercopy = kmer;
    uint64_t revcopy = rev;

    while (i > 0) {
        kmercopy = kmer;
        revcopy = rev;

        kmercopy = (kmercopy >> i & 0b11);
        revcopy = (revcopy >> i & 0b11);

        if (kmercopy == revcopy)
        {
            i -= 2;
        }
        else 
        {
            //if kmerG and revT return kmer
            if ((kmercopy) == 0b11 && (revcopy) == 0b10) {
                return (true); //kmer
            }
            //if kmerT and revG return rev
            if ((revcopy) == 0b11 && (kmercopy) == 0b10) {
                return (false); //rev
            }
            else {
                return (kmercopy < revcopy);
            }
        }
    }

    return (true);
}

/*
*  000000 00 = A   000000 11 = G 
* -2              -2
*  111111 10 = T   000000 01 = C
*
*  000000 01 = C   000000 10 = T
* -2              -2
*  111111 11 = G   000000 00 = A 
*/ 

uint64_t encode_kmer(deque<u_int8_t> currkmer, uint16_t kmersize) {
    uint64_t encoded = 0;
    uint64_t encoded_rev = 0;

    for (uint8_t c : currkmer) {
        encoded <<= 2;
        encoded |= c;
        c = (c - 2) & 0b11;
        encoded_rev |= (uint64_t) c << 62;
        encoded_rev >>= 2;
    }

    encoded_rev >>= 2 * (31 - kmersize);

    if (comp_kmer(encoded, encoded_rev, kmersize)) 
        return encoded;
    return encoded_rev;
}

/*TODO: Test getting next hash from previous hash
H(“bcd”)=H(“abc”)−H(“a”)+H(“d”)*/

/** Rolling hash function
uint64_t hash_kmer(deque<u_int8_t> currkmer, uint16_t kmersize, uint64_t bloomsize) {
    uint64_t kmerint = 0;
    deque<u_int8_t>::iterator it = currkmer.begin();
    int i = 0;

    while (it != currkmer.end()) {
        kmerint += (u_int64_t) (*it++ * pow(4, i));
        i++;
    }
    return (kmerint);
}*/

int main(int argc, char ** argv) {
    if (argc != 6) help();

    Parameters params = {};
    params.filename = argv[1];
    params.k = atoi(argv[2]);
    params.n = atoi(argv[3]);
    params.nf = atoi(argv[4]);
    params.r = atoi(argv[5]);

    ifstream fasta_stream(params.filename);
    
    if (!fasta_stream.is_open()) { 
        perror("Error opening file"); 
        exit(1);
    }

    deque<u_int8_t> kmer;
    deque<u_int8_t>::iterator kmer_it;
    
    Bloomfilter bf = Bloomfilter(params.n, params.nf);

    for (int i = 0; i < params.k; i++) {
        kmer.push_back(nucltoi(next_nucl(fasta_stream)));
    }

    while (fasta_stream.peek() != EOF) {
        bf.add_value(encode_kmer(kmer, params.k));

        kmer = next_kmer(kmer, fasta_stream);
    }

    for (int i = 0; i < params.r; i++) {
        deque<u_int8_t> randkmer = random_kmer(params.k);

        printf("Testing kmer: ");
        
        for (int i = 0; i < params.k; i++) {
            printf("%d ", randkmer[i]);
        }

        u_int64_t x = encode_kmer(randkmer, params.k);

        printf("\n");
        printf("Is present: %d\n", bf.is_present(x));
    }

    fasta_stream.close();

    return(0);
}