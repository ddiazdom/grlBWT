//
// Created by Diaz, Diego on 22.11.2022.
//

#include "int_array.h"
#include "rank_support.h"
#include "sdsl/bit_vectors.hpp"

int main(int argc, char** argv) {
    /*bit_array bv;
    rank_support<> rs;
    bv.reserve(2000);

    sdsl::bit_vector bv_sdsl;
    bv_sdsl.resize(2000);

    for(size_t i=0;i<2000;i++){
        bool bit = rand() & 1UL;
        bv.write(i, bit);
        bv_sdsl[i] = bit;
    }

    rs.set_bit_vector(&bv);
    sdsl::bit_vector::rank_1_type rs_sdsl(&bv_sdsl);

    for(size_t i=0;i<2000;i++){
        //assert(bv[i]==bv_sdsl[i]);
        //std::cout<<bv_sdsl[i]<<" "<<rs_sdsl(i)<<std::endl;
        std::cout<<bv[i]<<" "<<bv_sdsl[i]<<" "<<rs(i)<<" "<<rs_sdsl(i)<<std::endl;
        assert(rs(i)==rs_sdsl(i));
    }*/
    int_array<size_t> vector(2000, 33);
    for(size_t i=0;i<2000;i++){
        vector[i] = rand() % 20000;
        size_t elm = vector[i]+1;

        std::cout<<vector[i]<<" "<<elm<<std::endl;
        vector[i]++;
        vector[i]--;
        vector[i] = vector[i]>>1UL;
        vector[i] = vector[i] | 1UL;
        if(i>0){
            std::cout<<(vector[i]==vector[i-1])<<std::endl;
            vector[i] = vector[i-1];
        }
    }
}
