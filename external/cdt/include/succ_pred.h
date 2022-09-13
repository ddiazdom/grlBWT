//
// Created by Diaz, Diego on 3.3.2022.
//

#ifndef CDT_SUCC_PRED_H
#define CDT_SUCC_PRED_H

#include <iostream>
#include "int_array.h"

//a successor/predecessor data structure that relies on the two-level rank data structure of
template<class sp_bv_type>
struct succ_pred{

    size_t w{};//sampling rate for the small blocks
    size_t s{};//sampling rate for the big blocks
    //pointer to the original bit vector, which we assume is encoded using a representation for sparse bit vectors.
    sp_bv_type * vector=nullptr;
    std::vector<size_t> big_block; //samples every s positions
    int_array<size_t> small_block; //samples every w positions
    typename sp_bv_type::iterator_t it_end; //a static iterator pointing to the end of vector (it has a practical purpose)

    succ_pred()=default;

    explicit succ_pred(sp_bv_type& _vector): w(sizeof(size_t)*8),
                                                 s(w*w),
                                                 vector(&_vector),
                                                 big_block(INT_CEIL(_vector.back(),s)+1, 0),
                                                 small_block(INT_CEIL(_vector.back(),w)+1, sym_width(s)+1),
                                                 it_end((*vector).end()){

        small_block.initialize(0, INT_CEIL(_vector.back(),w)+1);

        for(auto const i : _vector){
            big_block[(i/s)]++;
            small_block.write(i/w, small_block[(i/w)]+1);
        }

        size_t acc=0, tmp;
        for(unsigned long & i : big_block){
            tmp = i;
            i = acc;
            acc+=tmp;
        }

        acc=0; size_t j=0;
        for(size_t i=0;i<small_block.size();i++){
            tmp = small_block[i];
            if((j & (w-1))==0) acc = 0;
            small_block.write(i, acc);
            acc+=tmp;
            j++;
        }
    };

    void swap(succ_pred& other){
        std::swap(w, other.w);
        std::swap(s, other.s);
        std::swap(vector, other.vector);
        big_block.swap(other.big_block);
        small_block = std::move(other.small_block);
    }

    [[nodiscard]] inline size_t binary_search(size_t value, size_t left, size_t right, size_t& val) const {
        size_t mid;
        while(true){
            mid = left + ((right-left)>>1UL);
            val = (*vector)[mid];

            if(val<=value && ((mid+1)==(*vector).size() || value<(*vector)[mid+1])) break;

            if(value<val){
                right = mid-1;
            }else{
                left = mid+1;
            }
        }
        return mid;
    }

    //returns a pair {rank, j} where j is the rightmost index such that bv[j]=1 and j<=idx.
    // The value for rank is one-based, and it is the rank for bv[j].
    [[nodiscard]] std::pair<size_t, size_t> pred(size_t idx) const {

        if(idx>vector->back()){
            return {vector->size(), vector->back()};
        }

        size_t bb = idx/s;
        size_t sb = idx/w;
        size_t l_bound = big_block[bb] + small_block[sb];
        size_t pred=0, rank;

        if(l_bound>0) l_bound--;

        auto it = vector->begin(l_bound);
        rank = l_bound;

        //The expensive part of pred is in this loop (I already profiled the code).
        // The iteration requires us to decode up to w=64 consecutive positions of vector.
        // In every operation ++it, I have to decode one elias-delta code.
        // The binary search does not improve the performance as I don't know how to
        // move backward in a delta-encoded vector. Hence, although I require log(w) = log(64)=6 access to
        // vector during the binary search (in the worst case), every operation vector[i] costs me
        // up to 64 elias-delta decodes (64 is the sampling rate I am using in vector).
        // Using elias-fano does not improve the situation.
        while(it!=it_end && *it<=idx){
            pred = *it;
            rank++;
            ++it;
        }
        return {rank, pred};

        /*size_t l_idx, prev_val;
        if(l_idx==idx){
            return {l_bound, idx};
        }else if(l_idx>idx){
            if(l_bound==0) return {0,0};
            return {l_bound-1, prev_val};
        }else{
            size_t r_bound = ((sb+1) & (w-1)) ==0 ? big_block[bb+1] : big_block[bb] + small_block[sb+1];
            size_t val;
            size_t res = binary_search(idx, l_bound, r_bound, val);
            return {res, val};
        }*/
    }

    /*[[nodiscard]] size_t succ(size_t idx) const {
        size_t bb = idx/s;
        size_t sb = idx/w;
        size_t l_bound = big_block[bb] + small_block[sb];
        if((*vector)[l_bound]==idx){
            return l_bound==((*vector).size()-1) ? (*vector)[l_bound] : (*vector)[l_bound+1];
        } if((*vector)[l_bound]>idx) {
            return (*vector)[l_bound];
        } else{
            size_t r_bound = ((sb+1) & (w-1)) ==0 ? big_block[bb+1] : big_block[bb] + small_block[sb+1];
            size_t val;
            size_t res = binary_search(idx, l_bound, r_bound, val);
            return res==((*vector).size()-1) ? val : (*vector)[res+1];
        }
    }

    [[nodiscard]] size_t rank(size_t idx) const {
        size_t bb = idx/s;
        size_t sb = idx/w;
        size_t l_bound = big_block[bb] + small_block[sb];
        if((*vector)[l_bound]==idx){
            return l_bound;
        } if((*vector)[l_bound]>idx) {
            return l_bound-1;
        } else{
            size_t r_bound = ((sb+1) & (w-1)) ==0 ? big_block[bb+1] : big_block[bb] + small_block[sb+1];
            size_t val;
            size_t res = binary_search(idx, l_bound, r_bound, val);
            return res;
        }
    }

    [[nodiscard]] bool is_set(size_t idx) const {
        size_t bb = idx/s;
        size_t sb = idx/w;
        size_t l_bound = big_block[bb] + small_block[sb];
        if((*vector)[l_bound]==idx){
            return true;
        } if((*vector)[l_bound]>idx) {
            return false;
        } else{
            size_t r_bound = ((sb+1) & (w-1)) ==0 ? big_block[bb+1] : big_block[bb] + small_block[sb+1];
            size_t val;
            binary_search(idx, l_bound, r_bound, val);
            return  val == idx;
        }
    }*/

    size_t serialize(std::ostream &out) const{
        size_t written_bytes= serialize_elm(out, w);
        written_bytes += serialize_elm(out, s);
        written_bytes += serialize_plain_vector(out, big_block);
        written_bytes += small_block.serialize(out);
        return written_bytes;
    }

    void load(std::istream &in){
        load_elm(in, w);
        load_elm(in, s);
        load_plain_vector(in, big_block);
        small_block.load(in);
    }

    void set_vector(sp_bv_type* _vector){
        vector = _vector;
        it_end = vector->end();
    }

};
#endif //CDT_SUCC_PRED_H
