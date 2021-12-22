//
// Created by Diaz, Diego on 15.12.2021.
//

#ifndef LPG_COMPRESSOR_COMPUTE_BWT_H
#define LPG_COMPRESSOR_COMPUTE_BWT_H
#include "grammar.hpp"

template<class vector_t>
struct rank_support{
    size_t w, s;
    const vector_t& vector;
    sdsl::int_vector<> big_block;
    sdsl::int_vector<> small_block;

    explicit rank_support(const vector_t& _vector): w(sizeof(size_t)*8),
                                                    s(w*w),
                                                    vector(_vector),
                                                    big_block(INT_CEIL(vector[vector.size()-1],s)+1, 0),
                                                    small_block(INT_CEIL(vector[vector.size()-1],w)+1, 0, sdsl::bits::hi(s)+1){
        for(unsigned long long i : vector){
            big_block[(i/s)]++;
            small_block[(i/w)]++;
        }

        size_t acc=0, tmp;
        for(auto && i : big_block){
            tmp = i;
            i = acc;
            acc+=tmp;
        }

        acc=0; size_t j=0;
        for(auto && i : small_block){
            tmp = i;
            if((j & (w-1))==0) acc = 0;
            i = acc;
            acc+=tmp;
            j++;
        }
    };

    [[nodiscard]] inline size_t binary_search(size_t value, size_t left, size_t right) const {
        size_t mid, val;
        while(true){
            mid = left + ((right-left)>>1UL);
            val = vector[mid];

            if(val<=value && ((mid+1)==vector.size() || value<vector[mid+1])) break;

            if(value<val){
                right = mid-1;
            }else{
                left = mid+1;
            }
        }
        return mid;
    }

    [[nodiscard]] size_t successor(size_t idx) const {
        size_t bb = idx/s;
        size_t sb = idx/w;
        size_t l_bound = big_block[bb] + small_block[sb];
        if(vector[l_bound]==idx){
            return l_bound==(vector.size()-1) ? vector[l_bound] : vector[l_bound+1];
        } if(vector[l_bound]>idx) {
            return vector[l_bound];
        } else{
            size_t r_bound = ((sb+1) & (w-1)) ==0 ? big_block[bb+1] : big_block[bb] + small_block[sb+1];
            size_t res = binary_search(idx, l_bound, r_bound);
            return res==(vector.size()-1) ? vector[res] : vector[res+1];
        }
    }
};

template<class grammar_t>
struct gramm_extra_feat{
    rank_support<sdsl::int_vector<>> rs_rp;
    rank_support<sdsl::int_vector<>> rs_sp;
    size_t                           c_start;
    explicit gramm_extra_feat(grammar_t& gram): rs_rp(gram.nter_ptr),
                                                rs_sp(gram.seq_ptr),
                                                c_start(gram.nter_ptr[gram.symbols() - 1]){};

    [[nodiscard]] inline size_t rb(size_t idx) const {
        if(idx>=c_start){
            return c_start+(rs_sp.successor(idx-c_start)-1);
        }else{
            return rs_rp.successor(idx)-1;
        }
    }
};

template<class grammar_t>
void gram2bwt(grammar_t& gram);

#endif //LPG_COMPRESSOR_COMPUTE_BWT_H
