//
// Created by Diaz, Diego on 4.1.2022.
//

#ifndef LPG_COMPRESSOR_LMS_INDUCTION_H
#define LPG_COMPRESSOR_LMS_INDUCTION_H
#include "common.h"
#include "lc_gram_algo.hpp"

void suffix_induction(vector_t &sa, const dictionary &dict);
void induce_L_type(vector_t &sa, const dictionary &dict, vector_t &buckets);
void induce_S_type(vector_t &sa, const dictionary& dict, vector_t &buckets);

#endif //LPG_COMPRESSOR_LMS_INDUCTION_H
