
## grl-BWT : Efficient construction of the BWT using string compression 

This repository contains the implementation of grl-BWT, a linear-time
semi-external algorithm for building the BCR BWT of a string collection.
The novelty of grl-BWT is that it resorts to compression to maintain the
intermediate data of the BWT construction in compact form. The compact
format used by grl-BWT reduces working memory and avoids redundant computations,
thus decreasing computing time too.

## Third-party libraries

1. [SDSL-lite](https://github.com/simongog/sdsl-lite)
2. [xxHash](https://github.com/Cyan4973/xxHash)

## Prerequisites

1. C++ >= 17 
2. CMake >= 3.7
3. SDSL-lite

The xxHash library is already included in the source files. We include a CMake module that will search for the local
installation of the SDSL-library. No need to indicate the path during the compilation.

## Installation

Clone repository, enter the project folder and execute
the following commands:

```
mkdir build
cd build
cmake ..
make
```

## Computing the BCR BWT 
```
./grlBWT input_file.txt -o output_bcr_bwt
```

For the moment, grlBWT only supports input files in one-string-per-line format.
We plan to expand to FASTA/Q files for genomic collections.

## Output

Our tool outputs the BCR BWT as a run-length compressed array. Concretely, it produces a sequence of equal-symbol runs
encoded as pairs (a,l), where a is the run symbol and l is its length. To reduce the space, we represent the run symbols
and the run lengths using cells of different widths. The width for the symbols is the smallest number of bytes that fits
the alphabet. On the other hand, the width for the lengths is the smallest number of bytes that fits the length of the
longest equal-symbol run in the BCR BWT.

We plan to write a parser in the future to produce the BCR BWT in plain format. 

## Multithreading

Our algorithm has a hashing step that supports multithreading to some extent. If you want to use it, pass the parameter
-t to the execution of grlBWT.

## Printing the original text

If you want to print the original strings, please use the print_seqs