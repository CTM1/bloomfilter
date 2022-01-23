#include "hash.hpp"
#include "Bloomfilter.hpp"
#include "kmer_utils.hpp"
#include "fasta_utils.hpp"
#include <unistd.h> 

using namespace std;

struct Parameters {
    char * filename;
    uint8_t k;      // Size of kmers, 3 to 31
    uint8_t nf;     // Number of hashing functions
    uint64_t n;     // Size of bloom filter, max 2^34
    uint32_t r;     // Number of requests
    
    float fp;       // False-positive rate
    uint64_t elems; // Estimate of items to be stored in bloom filter
};

void help() {
    printf("Takes all parameters and builds according bloom filter.\n\n");
    printf("Usage:\n\t./bloomfilter <filepath> <kmerSize> <bitvecSize> <hashes> <requests>\n");
    printf("Example:\n\t./bloomfilter data/ecoli.fasta 31 456637 3 10000\n\n\n\n");

    printf("Takes a false-positive rate and an estimate of the number\n");
    printf("of elements to be stored in the bloom filter, and builds filter accordingly.\n\n");

    printf("Usage:\n\t./bloomfilter -p <filepath> <kmerSize> <fp_rate> <elems> <requests>\n");
    printf("Example:\n\t./bloomfilter -p data/sarscov.fasta 31 0.001 762342 1000\n");
    exit(0);
}

void verify_params(Parameters p) {
    if (p.k > 31 || p.k < 3) {
        cerr << "Invalid kmer value, please choose a value between 3 - 31.\n";
        exit(1);
    }

    // The bitvec could be much bigger if needed
    if (p.n > 17179869185 || p.k < 1) {
        cerr << "Invalid bitvector size, please choose a value between 1 - 2^34.\n";
        exit(2);
    }

    if (p.nf < 0 || p.nf > 64) {
        cerr << "Invalid number of hashes, please choose a value between 0 - 64.\n";
        exit(3);
    }

    if (p.r < 0) {
        cerr << "Invalid number of requests, please choose a positive number.\n";
        exit(4);
    }

    if (p.fp < 0 || p.fp > 1) {
        cerr << "Invalid false-positive rate, please choose a value between 0 - 1.\n";
        exit(5);
    }

    if (p.elems < 0) {
        cerr << "Can't estimate negative number of elements, please choose a positive number.\n";
        exit(6);
    }
}

int main(int argc, char ** argv) {
    Parameters p = {};
    Bloomfilter bf = Bloomfilter();
    if (argc != 7 && argc != 6) help();

    if (argc == 7) {
        p.filename = argv[2];
        p.k = atoi(argv[3]);
        p.fp = strtold(argv[4], NULL);
        p.elems = strtoul(argv[5], NULL, 10);        
        p.r = atoi(argv[6]);

        verify_params(p);

        bf = Bloomfilter(p.fp, p.elems);
    }
    else {
        p.filename = argv[1];
        p.k = atoi(argv[2]);
        p.n = strtoul(argv[3], NULL, 10);
        p.nf = atoi(argv[4]);;
        p.r = atoi(argv[5]);

        verify_params(p);

        bf = Bloomfilter(p.n, p.nf);
    }

    ifstream fasta_stream(p.filename);
    
    if (!fasta_stream.is_open()) { 
        perror("Error opening file"); 
        exit(7);
    }

    // Skipping header
    skip_line(fasta_stream);
    
    uint64_t kmer = 0;
    
    for (uint32_t i = 0; i < p.k; i++) {
        kmer = next_kmer(kmer, fasta_stream, p.k);
    }

    while (fasta_stream.peek() != EOF) {
        bf.add_value(choose_kmer_or_rev(kmer, p.k));

        kmer = next_kmer(kmer, fasta_stream, p.k);
    }

    for (uint32_t i = 0; i < p.r; i++) {
        bf.is_present(random_kmer(p.k));
    }

    fasta_stream.close();
    return(0);
}

