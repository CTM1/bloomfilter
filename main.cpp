#include "hash.hpp"
#include "Bloomfilter.hpp"
#include "kmer_utils.hpp"
#include "fasta_utils.hpp"

using namespace std;

struct Parameters {
    char * filename;
    uint8_t k;  // Size of kmers
    uint64_t n;  // Size of bloom filter
    uint32_t nf; // Number of hashing functions
    uint64_t r;  // Number of requests
};

void help() {
    printf("Usage:\n\t./bloomfilter <filepath> <kmerSize> <bitvecSize> <hashes> <requests>\n\t");
    printf("Example : ./bloomfilter data/ecoli.fasta 31 456637 3 10000");
    exit(0);
}

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

    uint64_t kmer = 0;
    
    Bloomfilter bf = Bloomfilter(params.n, params.nf);
    
    for (int i = 0; i < params.k; i++) {
        kmer = next_kmer(kmer, fasta_stream);
    }

    while (fasta_stream.peek() != EOF) {
        bf.add_value(choose_kmer_or_rev(kmer, params.k));

        kmer = next_kmer(kmer, fasta_stream);
    }

    for (int i = 0; i < params.r; i++) {
        uint64_t randkmer = random_kmer(params.k);

        printf("Testing presence for: ");
        
        print_kmer(randkmer, params.k);

        printf("\nIs present: %d\n\n", bf.is_present(randkmer));
    }

    fasta_stream.close();

    return(0);
}