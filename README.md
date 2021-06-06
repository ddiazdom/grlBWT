
## A grammar self-index based on LMS substrings

A grammar self-index of a text $T$ (Claude et al. 2012) consists of a grammar $\mathcal{G}$ that only produces $T$ and a geometric data structure that indexes the string cuts of the right-hand sides of $\mathcal{G}$'s rules. This representation uses space proportional to $G$, the size of the grammar, which is small when the text is repetitive. However, the index is slow for pattern matching; it finds the $occ$ occurrences of a pattern $P[1..m]$ in $O((m^{2}+occ)\log G)$ time. The most expensive part is a set of binary searches for the different cuts $P[1..j]P[j+1..m]$ in the geometric data structure. Christiansen et al. (2019) solved this problem by building a locally consistent grammar that only searches for $O(\log m)$ cuts of $P$. Their representation, however, requires significant extra space (tough still in $O(G)$) to store a set of permutations for the nonterminal symbols. In this work, we propose another locally consistent grammar that builds on the idea of LMS substrings (Nong et al. 2009). Our grammar also requires to try $O(\log m)$ cuts when searching for $P$, but it does not need to store permutations. 
As a result, we obtain a self-index that searches in time $O((m\log m+occ) \log G)$ and is of practical size. Our experiments showed that our index is faster than previous grammar-indexes at the price of increasing the space by a 1.8x factor on average. Other experimental results showed that our data structure becomes convenient when the patterns to search for are long.

## Third-party libraries

1. [SDSL-lite](https://github.com/simongog/sdsl-lite)
2. [xxHash](https://github.com/Cyan4973/xxHash)

## Prerequisites

1. C++ >= 14
2. CMake >= 3.7
3. SDSL-lite

The xxHash library is already included in the source files. We include a CMake module that will search for
the local installation of the SDSL-library. No need to indicate the path during the compilation.

**Important** : before installing the SDSL-lite, you need to make the following classes public:


## Installation

Clone repository, enter the project folder and execute
the following commands:

```
mkdir build
cd build
cmake ..
make
```

## Creating the index
```
./lpg index tests/sample_file.txt
```

## Search for a pattern

```
./lpg search sample_file.txt.idx -F pattern_list.txt -p "test pattern"
```

The variable **sample_file.txt.idx** is the index created in the previous step and **pattern_list.txt** is a plain text file with a list of patterns separated by a new line

