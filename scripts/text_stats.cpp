//
// Created by Diaz, Diego on 24.2.2023.
//

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <cstring>

int main (int argc, char** argv){

    if(argc!=2){
        std::cout<<"usage: ./text_stats text"<<std::endl;
        std::cout<<"text_stats: a plain text in one-string-per-line format"<<std::endl;
        exit(0);
    }

    std::string input_file = std::string(argv[1]);
    std::cout<<"Computing text stats for the file "<<input_file<<std::endl;

    std::ifstream ifs(input_file, std::ios::binary);
    assert(ifs.good());
    uint8_t sep_sym;
    size_t sym_frq[256] = {0};
    ifs.seekg(-1, std::ios::end);
    ifs.read((char *)&sep_sym, 1);//we assume the last symbol is the sep symbol
    ifs.seekg(0, std::ios::beg);
    size_t pos = 0, cont=0;

    std::streampos buff_size=8*1024*1024;
    auto buffer = (uint8_t*)malloc(buff_size);
    memset(buffer,0, buff_size);
    size_t read_bytes, str_len;
    uint8_t sym;
    size_t longest_string=0;
    size_t shortest_string=std::numeric_limits<size_t>::max();

    while(true){
        ifs.read((char *)buffer, buff_size);
        read_bytes = ifs.gcount();
        if(read_bytes>0){
            for(size_t i=0;i<read_bytes;i++){
                sym = buffer[i];
                sym_frq[sym]++;
                cont++;
                if(sym==sep_sym){
                    str_len = cont - pos;
                    if(str_len>longest_string) longest_string=str_len;
                    if(str_len<shortest_string) shortest_string=str_len;
                    pos = cont;
                }
            }
        }else{
            break;
        }
    }

    ifs.close();
    std::vector<uint8_t> alphabet;
    std::vector<size_t> sym_freqs;
    size_t n_symbols=0;
    for(size_t i=0;i<256;i++){
        if(sym_frq[i]!=0){
            alphabet.push_back(i);
            sym_freqs.push_back(sym_frq[i]);
            n_symbols+=sym_frq[i];
        }
    }

    std::cout<<"Total symbols:         "<<n_symbols<<std::endl;
    std::cout<<"Alphabet:              "<<alphabet.size()<<std::endl;
    for(size_t i=0;i<alphabet.size();i++){
        std::cout<<"    symbol: "<<(int)alphabet[i]<<" freq: "<<sym_freqs[i]<<std::endl;
    }
    std::cout<<"Smallest symbol:       "<<(int)alphabet[0]<<std::endl;
    std::cout<<"Greatest symbol:       "<<(int)alphabet.back()<<std::endl;
    std::cout<<"Number of strings:     "<<sym_freqs[0]<<std::endl;
    std::cout<<"Longest string:        "<<longest_string<<std::endl;
    std::cout<<"Shortest string:       "<<shortest_string<<std::endl;
    std::cout<<"Average string length: "<<n_symbols/sym_freqs[0]<<std::endl;
    std::cout<<"Sep. symbol:           "<<(int)sep_sym<<std::endl;

    if(sep_sym!=alphabet[0]){
        std::cerr<<"Warning: the separator symbol is not the smallest symbol in the text"<<std::endl;
    }
    free(buffer);
    return 0;
}