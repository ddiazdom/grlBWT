//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_HPP
#define GRLBWT_HPP

#include "utils.h"
#include "lc_gram_algo.hpp"
#include "exact_grl_bwt.hpp"
#include "opt_grl_bwt.hpp"

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
template<bool opt_bwt>
void grl_bwt_algo(std::string &i_file, std::string& o_file, tmp_workspace & tmp_ws,
                  size_t n_threads, str_collection& str_coll, float hbuff_frac, uint8_t b_p_r){

    auto hbuff_size = std::max<size_t>(64 * n_threads, size_t(std::ceil(float(str_coll.n_char) * hbuff_frac)));

    if constexpr (opt_bwt){
        std::cout<<"Constructing the optimal BCR BWT of "<<i_file<<std::endl;
        size_t p_rounds = build_lc_gram(i_file, n_threads, hbuff_size, str_coll, tmp_ws);
        switch (b_p_r) {
            case 0:
                opt_ind_phase<0>(tmp_ws, p_rounds);
                break;
            case 1:
                opt_ind_phase<1>(tmp_ws, p_rounds);
                break;
            case 2:
                opt_ind_phase<2>(tmp_ws, p_rounds);
                break;
            case 3:
                opt_ind_phase<3>(tmp_ws, p_rounds);
                break;
            case 4:
                opt_ind_phase<4>(tmp_ws, p_rounds);
                break;
            default:
                opt_ind_phase<5>(tmp_ws, p_rounds);
        };
    }else{
        std::cout<<"Constructing the exact BCR BWT of "<<i_file<<std::endl;
        size_t p_rounds = build_lc_gram(i_file, n_threads, hbuff_size, str_coll, tmp_ws);
        switch (b_p_r) {
            case 0:
                exact_ind_phase<0>(tmp_ws, p_rounds);
                break;
            case 1:
                exact_ind_phase<1>(tmp_ws, p_rounds);
                break;
            case 2:
                exact_ind_phase<2>(tmp_ws, p_rounds);
                break;
            case 3:
                exact_ind_phase<3>(tmp_ws, p_rounds);
                break;
            case 4:
                exact_ind_phase<4>(tmp_ws, p_rounds);
                break;
            default:
                exact_ind_phase<5>(tmp_ws, p_rounds);
        };
    }
    std::filesystem::rename(tmp_ws.get_file("bwt_lev_0"), o_file);
    std::cout<<"The resulting BCR BWT was stored in "<<o_file<<std::endl;
}
#endif //GRLBWT_HPP
