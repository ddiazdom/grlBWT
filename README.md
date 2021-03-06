## grlBWT : Efficient construction of the BWT using string compression

This repository contains the implementation of **grlBWT**, a linear-time semi-external algorithm for building the BCR BWT
of a string collection. The novelty of grlBWT is that it resorts to compression to maintain the intermediate data of
the BWT construction in compact form. The compact format used by grlBWT reduces working memory and avoids redundant
computations, thus decreasing computing time too.

## Third-party libraries

1. [SDSL-lite](https://github.com/simongog/sdsl-lite)
2. [xxHash](https://github.com/Cyan4973/xxHash)
3. [CLI11](https://github.com/CLIUtils/CLI11)

## Prerequisites

1. C++ >= 17
2. CMake >= 3.7
3. SDSL-lite

The xxHash and CLI11 libraries are already included in the source files of this repository. The SDSL-lite has to be installed beforehand.
However, we include a CMake module that will search for its local
installation. No need to indicate the path during the compilation.

## Installation

Clone repository, enter the project folder and execute the following commands:

```
mkdir build
cd build
cmake ..
make
```

## Example execution

```
./grlbwt ../tests/sample_file.txt
```

Our tool currently supports string collections in one-string-per-line format or the classical FASTA and FASTQ formats.
It also handles gzipped inputs. In principle, grlBWT automatically detects the format. However,
we have not tested the automatic detection functionality exhaustively yet.

## Output

Our tool outputs the BCR BWT as a run-length compressed array. Concretely, it produces a sequence of equal-symbol runs
encoded as pairs (*a*,*l*), where *a* is the run symbol and *l* is its length. We represent the run symbols and the run
lengths using cells of different widths to reduce space usage. The width for the symbols is the smallest number of bytes
that fits the alphabet. On the other hand, the width for the lengths is the smallest number of bytes that fits the
length of the longest equal-symbol run in the BCR BWT.

We plan to write a parser in the future to produce the BCR BWT in plain format.

## Printing the original text

If you want to print the original strings, please use the *print_seqs* binary:

```
./print_seqss sample_file.txt.rl_bwt 10
```

This command will print the first 10 strings of the original file used in example execution above.

## Multithreading

Our algorithm has a hashing step that supports multithreading to some extent. If you want to use it, pass the parameter
-t to the execution of grlBWT.

## How to cite

We will put a source to cite soon.

## Bugs

This tool still contains experimental code, so you will probably find bugs. Please report them in this repository.

## Author

This implementation was written by [Diego D??az](https://github.com/ddiazdom) .
