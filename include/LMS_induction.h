//
// Created by Diaz, Diego on 4.1.2022.
//

#ifndef LPG_COMPRESSOR_LMS_INDUCTION_H
#define LPG_COMPRESSOR_LMS_INDUCTION_H
#include "common.h"
#include "lc_gram_algo.hpp"
#include "bwt_io.h"

#define INSERT_BWT_SYM_FOR_SUFFIX(pos, bwt_sym) \
size_t phrase = d_lim_rs(pos); \
size_t freq = dict.freqs.read(phrase); \
size_t phrase_flag = freq & 3UL; \
assert(phrase_flag!=0 && phrase_flag<=3); \
freq>>=2UL;\
if(bwt_sym==dict.bwt_dummy){\
    if(bwt_sym==dict.bwt_dummy && phrase_flag==2){\
        bwt_sym = dict.end_str_dummy;       \
    }else if(bwt_sym==dict.bwt_dummy && phrase_flag==3) /*[[unlikely]]*/ { \
        /*TODO this is the case when a phrase occurs as a suffix but also as an entire string*/\
        exit(0);\
    }\
    meta_sym_list.write(phrase, (met_sym<<1UL) | (freq>1)); \
    met_sym++;                                     \
}\
if(pre_bwt.size()>0 && pre_bwt.last_sym() == bwt_sym){\
    pre_bwt.inc_freq_last(freq);\
}else{\
    pre_bwt.push_back(bwt_sym, freq);\
}

#define INCREMENT_BWT(start, end,  bwt_sym)\
if(is_gr_phrase) { \
    bwt_sym = dict.bwt_dummy+((end-start+1)>1);\
    for(size_t j=start;j<=end;j++){       \
    }                                      \
    if(n_phrases==1){                     \
        meta_sym_list.write(curr_phrase, (met_sym<<1UL) | (curr_phrase_freq>1));    \
    }                                      \
    met_sym++;                              \
}\
if(pre_bwt.size()>0 && pre_bwt.last_sym() == bwt_sym){\
    pre_bwt.inc_freq_last(acc_freq);\
}else{\
    pre_bwt.push_back(bwt_sym, acc_freq);\
}\

template <class value_type>
value_type * suffix_induction(const dictionary &dict, tmp_workspace& ws);

template <class value_type>
void induce_L_type(value_type *sa, size_t sa_size, const dictionary &dict, value_type *buckets, bv_t& solved_sym);

template <class value_type>
void induce_S_type(value_type *sa, size_t sa_size, const dictionary& dict, value_type *buckets, bv_t& solved_sym, tmp_workspace& ws);

//template <class value_type>
//void increment_bwt(size_t start, size_t end, value_type *sa, const dictionary& dict, bool is_gr_phrase, size_t acc_freq, size_t bwt_sym, bwt_buff_writer& bwt_writer);
//void insert_bwt_sym_for_suffix(size_t pos, size_t bwt_sym, const dictionary& dict,bv_rs_t& d_lim_rs, bwt_buff_writer& bwt_writer);
void invert_data(tmp_workspace& ws, size_t n_phrases, vector_t& phrase_metasymbols);
#endif //LPG_COMPRESSOR_LMS_INDUCTION_H
