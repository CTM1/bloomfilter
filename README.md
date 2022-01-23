# Bloomfilter
Intern candidate test for the Pasteur Institue. Bloom filter implementation.

This implementation creates a [Bloom filter](https://en.wikipedia.org/wiki/Bloom_filter) from a [FASTA](https://en.wikipedia.org/wiki/FASTA_format) file containing a 
single header, then executes a given number of requests, each one veryfing
if a random k-mer is present in the Bloom filter.

The program may take two sets of arguments to create the Bloom filter, for more information, see [Execution](#execution)

1. [Compiling](#compiling)
2. [Execution](#execution)
3. [Documentation](#documentation)
4. [Extra Notes](#extra-notes)

# Compiling
You may compile the code using `g++` by typing `make`. Other options are
available:

> `make clang` :
    Same compilation process as g++, but for the `clang++` compiler.

> `make lto_off` :
    Disables g++ [Link Time Optimization](https://gcc.gnu.org/wiki/LinkTimeOptimization).
    Small negative impact to performance, slightly smaller binary.

> `make unoptimized` :
    Disables g++ [Optimization Options](https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html).
    Big impact to performance, biggest binary.

> `make clang_lto_off` :
    Disables clang `Link Time Optimization`
    Negligeable impact to performance, slightly smaller binary.

> `make clang_unoptimized` :
    Disables clang `Optimization Options`.
    Big impact to performance, biggest binary.

# Execution

The program will always take a **filepath**, **k**, and **number of requests** (`r`). 

By default, you may specify:

- **Size of your bloom filter bit vector** in bits (max 2^34) (`size`).
- Number of **hashes** (max 64) (`nh`).

As such:
```bash
  #                  file        k   size  nh  r
  ./bloomfilter data/ecoli.fasta 31 456637 3 1000000
```

The `-p`\* option will take:
- **A wanted probability of false positives** as a floating point number (`fp`).
- **An estimate of the numbers of items** you would like to store in your bloom filter (`n`).
```bash
  #                     file        k   fp      n       r 
  ./bloomfilter -p data/ecoli.fasta 31 0.0001 4967535 1000000
```
Notes: 
- **This option does not respect the boundaries set by this [README](https://github.com/yoann-dufresne/bloomtest/)**,

Specifically the size of the bloom filter may be over 2^34 bits and there may be more than 64
hashing functions to accomodate for the false positive rate.

- **Unexpected behaviour may occur if the estimated size of the bit vector goes above $2^{63} - 1$**. 

Since this is the [maximum vector size in C++](https://en.cppreference.com/w/c/types/size_t) and
also of an uint64_t.

I estimated you would need $10^{13}$ expected elements with a false positive rate $10^{-17}$, 
so it is not accounted for.


\* As a matter of fact, the command line parsing isn't very robust, you can
just add any argument here and it will execute this option.

# Documentation

You may generate documentation files by using `doxygen`, to do so, install doxygen and type

```bash
doxygen Doxyfile
```
`html` and `latex` folders will be created.

To view the documentation open the `html/index.html` file in your browser. In GNOME desktops:

```bash
xdg-open html/index.html
```

`make clean` will remove both folders, and you will have to create it again.

# Extra Notes

### Small quirk
The original [README](https://github.com/yoann-dufresne/bloomtest/) specified the xorshift
function would always hash the 0 value to 0. 

For any other number, [it does not return zero](https://stackoverflow.com/questions/44753463/can-xorshift-return-zero).

Since each `uint64_t kmer` is unique in this implementation, changing the kmer value is 
equivalent to making it the same as another kmer, creating an unwanted collision.

In this implementation, every `k-mer` only containing the nucleotide `'A'` will aways hash to zero,
therefore in the bloom filter, the only bit representing it's presence will be the first one.

I could not find an efficient way to circumvent this quirk. However,
there is no chance of collision.

---

### Things I would've liked to add with more time
- A more robust command-line (please, be gentle, the option parsing is basically non-existent).
- Parsing a FASTA file with multiple headers, options to only parse specific ones and separate bloomfilters for each one.
- A way of storing kmers with $k > 31$ in the same [Huffman-encoding](https://en.wikipedia.org/wiki/Huffman_coding) fashion.

---
### Thanks
This was pretty fun.