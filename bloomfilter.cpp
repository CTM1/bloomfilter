#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/** Returns the next nucleotide in the DNA sequence,
skipping over headers **/
char next_nucl(fstream& f) {
    char c;

    c = f.get();

    if (c == '>') {
        while (f.get() != '\n')
            continue;
        c = f.get();
    }

    while (c == 'N') {
        f.get();
    }

    return (c);
}

int main(int argc, char ** argv) {
    string filename(argv[1]);
    fstream fasta_stream(filename);
    
    if (fasta_stream.is_open()) {
        for (int i = 0; i < 200; i++) {
            printf("%c", next_nucl(fasta_stream));
        }
    }
}