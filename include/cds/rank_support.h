//
// Created by Diaz, Diego on 22.11.2022.
//

#ifndef CDT_RANK_SUPPORT_H
#define CDT_RANK_SUPPORT_H
#include <vector>
#include "int_array.h"

//TODO k>W is supported just if k is a power of two (we can support Rank for any K but I have to modify rank operation)
template<size_t bl=64, size_t s_bl=bl*4>
class rank_support {
    static_assert(bl >= 64 && bl%64==0 && bl<s_bl && s_bl%bl==0);
    constexpr static uint8_t shift = 6;
    constexpr static size_t words_per_block = bl/64;
    constexpr static size_t blocks_per_super = s_bl/bl;

    const bit_array *m_bv = nullptr;
    std::vector<size_t> R;
    int_array<> r;

    void populate_rank_samples(){
        size_t R_acc=0, r_acc=0, words_in_block=0, blocks_in_super=0;

        R.clear();
        r.clear();

        R.reserve(m_bv->size()/s_bl + 1);
        r.set_width(sym_width(s_bl-bl));
        r.reserve(m_bv->size()/bl + 1);

        R.push_back(0);
        r.push_back(0);

        for(size_t i=0;i<m_bv->n_words();i++){

            size_t pc = __builtin_popcountll(m_bv->bits.stream[i]);
            R_acc+=pc;
            r_acc+=pc;

            if (++words_in_block==words_per_block) {
                words_in_block=0;
                if (++blocks_in_super==blocks_per_super) {
                    blocks_in_super=0;
                    R.push_back(R_acc);
                    r_acc = 0;
                }
                r.push_back(r_acc);
            }
        }
    }

    void copy(const rank_support &other){
        m_bv = other.m_bv;
        R = other.R;
        r = other.r;
    }

public:

    explicit rank_support(const bit_array* bv): m_bv(bv) {
        r.set_width(sym_width(s_bl-bl));
        populate_rank_samples();
    };

    rank_support()= default;

    void set_bit_vector(const bit_array *bv){
        m_bv = bv;
        populate_rank_samples();
    };

    size_t operator()(size_t pos) const{
        assert(pos<=m_bv->size());
        if(pos==0) return 0;

        pos--;
        size_t bl_index = pos/bl;//block index
        size_t word_index = (pos%bl)>>shift;//word index inside the block
        size_t start = bl_index*words_per_block;
        size_t end = start + word_index;

        size_t acc = R[pos/s_bl] + r[bl_index];
        for(size_t i=start;i<end;++i){
            acc += __builtin_popcountll(m_bv->bits.stream[i]);
        }

        return  acc + m_bv->bits.pop_count((bl_index*bl) + (word_index<<shift), pos);
    }

    rank_support& swap(rank_support& other) noexcept {
        std::swap(m_bv, other.m_bv);
        std::swap(R, other.R);
        r.swap(other.r);
        return *this;
    }

    rank_support& operator=(const rank_support& other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    };
};
#endif //CDT_RANK_SUPPORT_H
