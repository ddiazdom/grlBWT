//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_BUILD_HPP
#define LPG_COMPRESSOR_GRAMMAR_BUILD_HPP

#include "lc_gram_algo.hpp"
#include "suffpair_algo.hpp"

//this source code is for debugging
void decomp(size_t nt, sdsl::int_vector<> &rules,
            bv_ss_t &rlim_ss, bv_t &rem_nt,
            bv_rs_t &rem_nt_rs, ivb_t &dec);

//mark the nonterminals that can be removed from the grammar
bv_t mark_unique_nonterminals(gram_info_t& p_gram);
void simplify_grammar(gram_info_t &p_gram, bv_t &rem_nts, bv_rs_t &rem_nts_rs);
void run_length_compress(gram_info_t& p_gram, sdsl::cache_config& config);

/***
 *
 * @param i_file : input text file
 * @param p_gram_file : file where the plain grammar will be stored
 * @param n_threads : number of working threads
 * @param config : temporal files handler
 * @param hbuff_size : buffer size for the hashing step
 */
void build_gram(std::string &i_file, std::string &p_gram_file, std::string& tmp_folder, size_t n_threads, float hbuff_frac);

/***
 * check if the grammar is correct
 * @param g_file : file with the plain grammar
 * @param uncomp_file : original input text
 */
void check_plain_grammar(gram_info_t& p_gram, std::string& uncomp_file);
alpha_t get_alphabet(std::string &i_file);

#endif //LPG_COMPRESSOR_GRAMMAR_BUILD_HPP
