//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_BUILD_HPP
#define LPG_COMPRESSOR_GRAMMAR_BUILD_HPP

#include "lc_gram_algo.hpp"

size_t compute_hocc_size(ivb_t& bwt, dictionary& dict, bv_rs_t& hocc_rs,
                         vector_t& hocc_buckets, size_t p_round, sdsl::cache_config & config);
void infer_lvl_bwt(sdsl::cache_config& config, size_t p_round);
void parse2bwt(sdsl::cache_config& config, size_t p_round);
void infer_bwt(sdsl::cache_config& config, size_t p_round);

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
void grl_bwt_algo(std::string &i_file, std::string& o_file, std::string& tmp_folder, size_t n_threads, float hbuff_frac);
alpha_t get_alphabet(std::string &i_file);

#endif //LPG_COMPRESSOR_GRAMMAR_BUILD_HPP
