//
// Created by Diaz, Diego on 27.10.2021.
//

#ifndef CDS_UTILS_H
#define CDS_UTILS_H

#include <chrono>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <cassert>
#include <sstream>
#include <random>
#include <iomanip>
#include "memory_handler.hpp"

#ifdef __APPLE__
#include <unistd.h>
#endif

struct str_collection {
    std::vector<long> str_ptrs;
    size_t n_strings=0;
    size_t longest_string=0;
    size_t max_sym=0;
    size_t min_sym=std::numeric_limits<size_t>::max();
    size_t n_syms=0;
    size_t max_sym_freq=0;
};

template<class sym_type>
str_collection collection_stats(std::string& input_file){

    //We assume that the smallest symbol in the collection is the separator symbol.
    // This should be the last character in the file
    str_collection str_coll;

    size_t sym_frq[256] = {};

    std::ifstream ifs(input_file, std::ios::binary);

    sym_type sep_sym;
    ifs.seekg(-1*sizeof(sym_type), std::ios::end);

    ifs.read(reinterpret_cast<char *>(&sep_sym), 1*sizeof(sym_type));
    str_coll.n_syms = ifs.tellg()/sizeof(sym_type);
    assert(str_coll.n_syms>0);
    str_coll.max_sym_freq = str_coll.n_syms;
    ifs.seekg(0, std::ios::beg);
    size_t pos = 0, cont=0;

    std::streampos buff_size = 8388608;
    std::streampos n_elms = buff_size/sizeof(sym_type);
    auto buffer = mem::allocate<sym_type>(n_elms);
    memset(buffer, 0, buff_size);

    size_t read_bytes, str_len, read_syms;
    sym_type sym=0, min_sym=std::numeric_limits<sym_type>::max(), max_sym=0;

    while(true){
        ifs.read(reinterpret_cast<char *>(buffer), buff_size);

        read_bytes = ifs.gcount();

        if(read_bytes>0){

            read_syms = read_bytes/sizeof(sym_type);
            for(size_t i=0;i<read_syms;i++){
                sym = buffer[i];

                if constexpr (std::is_same<sym_type, uint8_t>::value){
                    ++sym_frq[sym];
                }else{
                    if(sym>max_sym) max_sym = sym;
                    if(sym<min_sym) min_sym = sym;
                }

                cont++;
                if(sym==sep_sym){
                    str_coll.str_ptrs.push_back((long)pos);
                    str_len = cont - pos;
                    if(str_len>str_coll.longest_string) str_coll.longest_string=str_len;
                    pos = cont;
                }
            }
        } else {
            break;
        }
    }

    str_coll.str_ptrs.shrink_to_fit();

    if constexpr (std::is_same<sym_type, uint8_t>::value) {
        size_t max_sym_freq = 0;
        for (size_t i = 0; i < 256; ++i) {
            if (sym_frq[i] > 0) {
                if (min_sym>i) min_sym = i;
                max_sym = i;
                if (sym_frq[i] > max_sym_freq) max_sym_freq = sym_frq[i];
            }
        }
        str_coll.max_sym_freq = max_sym_freq;
    }

    if(sep_sym!=min_sym || sym!=sep_sym){
        std::cerr<<"Error: the file is ill formed"<<std::endl;
        exit(1);
    }

    str_coll.min_sym = min_sym;
    str_coll.max_sym = max_sym;
    str_coll.n_strings = str_coll.str_ptrs.size();

    ifs.close();
    mem::deallocate(buffer);
    return str_coll;
}

//from https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
    std::ostringstream out;
    out.precision(n);
    out << a_value;
    return std::move(out).str();
}

inline std::string format_time(double seconds) {
    std::ostringstream os;
    if (seconds < 60) {
        os << std::fixed << std::setprecision(3)<< seconds << "s";
    } else if (seconds < 3600) {
        int m = static_cast<int>(seconds) / 60;
        int s = static_cast<int>(seconds) % 60;
        os << m << "m"<< std::setw(2) << std::setfill('0') << s << "s";
    }
    else {
        int h = static_cast<int>(seconds) / 3600;
        int m = (static_cast<int>(seconds) % 3600) / 60;
        os << h << "h"<< std::setw(2) << std::setfill('0') << m << "m";
    }
    return os.str();
}

inline std::string format_space(off_t bytes){
    if(bytes<1000){
        return std::to_string(bytes)+" bytes";
    }

    if(bytes<1000000){
        float b = static_cast<float>(bytes)/1000;
        return to_string_with_precision(b, 3)+" KBs";
    }

    if(bytes < 1000000000L){
        float b = static_cast<float>(bytes)/1000000;
        return to_string_with_precision(b, 3)+" MBs";
    }

    if(bytes < 1000000000000L){
        double b = static_cast<double>(bytes)/1000000000L;
        return to_string_with_precision(b, 3)+" GBs";
    }

    double b = static_cast<double>(bytes)/1000000000000L;
    return to_string_with_precision(b, 3)+" TBs";
}
#endif //CDS_UTILS_H
