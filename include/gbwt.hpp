//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_BUILD_HPP
#define LPG_COMPRESSOR_GRAMMAR_BUILD_HPP

#include "lc_gram_algo.hpp"

//this source code is for debugging
/***
 * check if the grammar is correct
 * @param g_file : file with the plain grammar
 * @param uncomp_file : original input text
 */
void check_plain_grammar(gram_info_t& p_gram, std::string& uncomp_file);

void sp_sa2bwt(std::string& sp_sa, std::string& output);
void induce_sym_order(gram_info_t& p_gram);
void gram2sp_sa(gram_info_t& p_gram);

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
void g_bwt_algo(std::string &i_file, std::string& o_file, std::string& tmp_folder, size_t n_threads, float hbuff_frac);

alpha_t get_alphabet(std::string &i_file);

#endif //LPG_COMPRESSOR_GRAMMAR_BUILD_HPP
