//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_HPP
#define GRLBWT_HPP

#include "lc_gram_algo.hpp"
#include "utils.h"

size_t compute_hocc_size(dictionary& dict, bv_rs_t& hocc_rs, vector_t& hocc_buckets, size_t p_round, tmp_workspace & ws);
void infer_lvl_bwt(tmp_workspace& ws, size_t p_round);
template<uint8_t b_f_r>
void infer_lvl_bwt_ht(tmp_workspace& ws, size_t p_round);
template<uint8_t b_f_r>
void ind_phase(tmp_workspace& ws, size_t p_round);

void parse2bwt(tmp_workspace& ws, size_t p_round);

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
void grl_bwt_algo(std::string &i_file, std::string& o_file, tmp_workspace & tmp_ws,
                  size_t n_threads, str_collection& str_coll, float hbuff_frac, uint8_t b_f_r);
#endif //GRLBWT_HPP
