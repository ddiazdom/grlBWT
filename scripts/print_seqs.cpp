#include <iostream>
#include "fm_index.h"

int main(int argc, char** argv){

    if(argc!=3){
        std::cout<<"usage: ./test_bcrbwt file.rlbwt print n_seqs\n"
                   "file is the bcr bwt file\n"
                   "n_seqs is the number of sequences to be decompressed"<<std::endl;
        exit(0);
    }

    errno = 0;
    char *endptr;
    long n_seqs = strtol(argv[2], &endptr, 10);
    std::string input_file = std::string(argv[1]);
    std::cout<<"decompressing "<<n_seqs<<" sequences"<<std::endl;
    fm_index fmi(input_file);

    std::cout<<"decompressing the first "<<n_seqs<<" sequences "<<std::endl;
    std::string seq;
    for(size_t i=0;i<n_seqs;i++){
        seq.clear();
        auto res = fmi.lf(i);
        size_t next_pos = res.second;
        uint8_t sym = res.first;
        while(sym!=fmi.get_dummy()){
            seq.push_back((char)sym);
            res = fmi.lf(next_pos);
            next_pos = res.second;
            sym = res.first;
        }
        std::reverse(seq.begin(), seq.end());
        std::cout<<seq<<std::endl;
    }
}