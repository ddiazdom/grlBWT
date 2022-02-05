//
// Created by Diaz, Diego on 5.2.2022.
//
#include <iostream>
#include "bwt_io.h"

int main(int argc, char** argv){

    if(argc!=2){
        std::cout<<"usage: ./print_n_runs file.rlbwt\n"
                   "file.rlbwt is the bcr bwt file\n"<<std::endl;
        exit(0);
    }

    errno = 0;
    std::string input_file = std::string(argv[1]);
    bwt_buff_reader bwt_reader(input_file);
    std::cout<<"Number of runs "<<bwt_reader.size()<<std::endl;
}

