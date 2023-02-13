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

## Reversing the BCR BWT to the original text

If you want to print the original strings, please use the *print_seqs* binary:

```
./print_seqss sample_file.txt.rl_bwt 10
```

This command will print the first 10 strings of the original file used in example execution above. If you don't specify
a number of strings, the script will decompress all the strings.

## Store the output BCR BWT in other formats  

The output of our tool is in a custom run-length-encoded format. However, you can easily transform it to
plain text.  

```
./grl2plain sample_file.txt.rl_bwt sample_file_bwt_plain.txt 
```

It is also possible to have the BWT in standard run-length-encoding:

```
./grl2rle sample_file.txt.rl_bwt sample_file_bwt_plain.rle 
```

The file `sample_file_bwt_plain.rle` is a sequence of pairs (a_1, l_1), (a_2, l_2), .... , (a_r, l_r) that encodes the 
r runs of the BCR BWT. Each pair (a_i, l_i) uses 16 bytes, 8 for a_i and 8 for l_i. 

## Multithreading

Our algorithm has a hashing step that supports multithreading to some extent. If you want to use it, pass the parameter
`-t` to the execution of grlBWT.

## Controlling the space usage of the runs

Our algorithm produces partial BCR BWTs so that each new BWT induces the next one. A relevant issue in this scheme is how to
encode the runs of those BWTs. Our default strategy is as follows: if a BWT has an alphabet of length sigma and the longest
run has l symbols, then we use ceil(log(sigma))+ ceil(log(2)) bits for every run. However, in most of the cases, this
value is an overestimation. In non-repetitive instances,
over ~97% of the runs are shorter than 256 (1 byte). This percentage might decrease to 80-60% in highly repetitive datasets,
but ~95% of the run lengths still fit two bytes. In contrast, the number of bytes to encode l can easily exceed the 4 bytes
in any kind of dataset.

We now have an alternative strategy to control the space usage of the runs. You can pass a fix value x (1-4) by the command
line (`-b` parameter) to tell our method that it can't use more than x bytes to encode run lengths. The algorithm will keep the runs with overflow
in a separate data structure. In most of the cases, choosing 1 or 2 should greatly reduce the space usage
during the induction of the BWTs, with almost not impact in the performance.


The general rule is simple: 

- non-repetitive input = `-b 1`
- highly repetitive input = `-b 2`
- insanely repetitive and massive input = `-b 3` (we suppose, but we haven't found an input for this)
- you don't know if the input is repetitive = don't use `-b`. The algorithm will use the default strategy

Example:

```
./grlbwt ../tests/sample_file.txt -b 1
```

If you chose a bad number, don't worry, it won't crash. It might run a bit slower or use more RAM than
it actually needed.

## How to cite

We will put a source to cite soon.

## Bugs

This tool still contains experimental code, so you will probably find bugs. Please report them in this repository.

## Author

This implementation was written by [Diego Díaz](https://github.com/ddiazdom) .
