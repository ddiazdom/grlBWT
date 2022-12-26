//
// Created by Diaz, Diego on 4.1.2022.
//

#ifndef LPG_COMPRESSOR_LMS_INDUCTION_H
#define LPG_COMPRESSOR_LMS_INDUCTION_H
#include "common.h"
#include "lc_gram_algo.hpp"

template <class value_type>
value_type * suffix_induction(const dictionary &dict);

template <class value_type>
void induce_L_type(value_type *sa, size_t sa_size, const dictionary &dict, vector_t &buckets, bv_t& solved_sym);
template <class value_type>
void induce_S_type(value_type *sa, size_t sa_size, const dictionary& dict, vector_t &buckets, bv_t& solved_sym);

#endif //LPG_COMPRESSOR_LMS_INDUCTION_H
