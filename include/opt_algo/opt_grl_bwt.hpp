//
// Created by Diaz, Diego on 26.1.2023.
//

#ifndef OPT_GRLBWT_HPP
#define OPT_GRLBWT_HPP

#include "bwt_io.h"
#include "utils.h"
#include "common.h"
#include "lc_gram_algo.hpp"

size_t opt_compute_hocc_size(dictionary& dict, bv_rs_t& hocc_rs, vector_t& hocc_buckets, size_t p_round, tmp_workspace & ws);

void opt_infer_lvl_bwt(tmp_workspace& ws, size_t p_round);

template<uint8_t b_f_r>
void opt_infer_lvl_bwt(tmp_workspace& ws, size_t p_round);

template<uint8_t b_f_r>
void opt_ind_phase(tmp_workspace& ws, size_t p_round);
#endif //OPT_GRLBWT_HPP
