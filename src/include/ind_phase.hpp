//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef EXACT_GRLBWT_HPP
#define EXACT_GRLBWT_HPP

#include "par_phase.hpp"
#include "utils.h"

size_t compute_hocc_size(dictionary &dict, bv_rs_t &hocc_rs, vector_t &hocc_buckets,
                         size_t p_round, tmp_workspace &ws);

template<uint8_t b_f_r>
void infer_lvl_bwt(tmp_workspace &ws, size_t p_round);
template<uint8_t b_f_r>
void ind_phase(tmp_workspace &ws, size_t p_round, logger& log);

void infer_lvl_bwt(tmp_workspace &ws, size_t p_round, logger& log);
void parse2bwt(tmp_workspace &ws, size_t &p_round, logger& log);

#endif //EXACT_GRLBWT_HPP
