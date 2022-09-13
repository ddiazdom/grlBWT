//
// Created by Diaz, Diego on 5.2.2022.
//
#include <iostream>
#include "bwt_io.h"

int main(int argc, char** argv){

    if(argc!=2){
        std::cout<<"usage: ./stats file.rlbwt\n"
                   "file.rlbwt is the bcr bwt file\n"<<std::endl;
        exit(0);
    }

    errno = 0;
    std::string input_file = std::string(argv[1]);
    size_t min_run=std::numeric_limits<size_t>::max();
    size_t max_run=0;
    size_t median=0;

    bwt_buff_reader bwt_reader(input_file);
    size_t sym=0, len=0, acc=0;
    std::cout<<"Reading the input BWT"<<std::endl;
    for(size_t i=0;i<bwt_reader.size();i++){
        bwt_reader.read_run(i, sym, len);
        if(len<min_run) min_run = len;
        if(len>max_run) max_run = len;
        acc+=len;
    }
    std::cout<<"Number of runs: "<<bwt_reader.size()<<std::endl;
    std::cout<<"Text size: "<<acc<<std::endl;
    std::cout<<"n/r: "<<double(acc)/double(bwt_reader.size())<<std::endl;
    std::cout<<"Longest run: "<<max_run<<std::endl;
    std::cout<<"Shortest run: "<<min_run<<std::endl;

}

