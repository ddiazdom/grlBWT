//
// Created by Diaz, Diego on 4.1.2022.
//

#ifndef LPG_COMPRESSOR_LMS_INDUCTION_H
#define LPG_COMPRESSOR_LMS_INDUCTION_H
#include "common.h"

void suffix_induction(vector_t& sa, vector_t &dict, bv_t &d_lim, size_t alphabet, sdsl::cache_config& config);
void induce_L_type(vector_t& sa, vector_t& dict, bv_t& d_lim, vector_t& buckets);
void induce_S_type(vector_t& sa, vector_t& dict, bv_t& d_lim, vector_t& buckets);

#endif //LPG_COMPRESSOR_LMS_INDUCTION_H
