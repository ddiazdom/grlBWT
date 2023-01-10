//
// Created by Diaz, Diego on 4.1.2022.
//

#ifndef LPG_COMPRESSOR_LMS_INDUCTION_H
#define LPG_COMPRESSOR_LMS_INDUCTION_H
#include "common.h"
#include "lc_gram_algo.hpp"
#include "bwt_io.h"

template <class value_type>
value_type * suffix_induction(const dictionary &dict, tmp_workspace& ws);

template <class value_type>
void induce_L_type(value_type *sa, size_t sa_size, const dictionary &dict, value_type *buckets, bv_t& solved_sym);

template <class value_type>
void induce_S_type(value_type *sa, size_t sa_size, const dictionary& dict, value_type *buckets, bv_t& solved_sym, tmp_workspace& ws);

template <class value_type>
void increment_bwt(size_t start, size_t end, value_type *sa, const dictionary& dict, bool is_gr_phrase, size_t acc_freq, size_t bwt_sym, bwt_buff_writer& bwt_writer);
void insert_bwt_sym_for_suffix(size_t pos, size_t bwt_sym, const dictionary& dict,bv_rs_t& d_lim_rs, bwt_buff_writer& bwt_writer);
void invert_data(tmp_workspace& ws);
#endif //LPG_COMPRESSOR_LMS_INDUCTION_H
