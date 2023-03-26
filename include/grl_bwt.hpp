//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_HPP
#define GRLBWT_HPP

#include "utils.h"

#include "exact_ind_phase.hpp"
#include "exact_par_phase.hpp"

#include "opt_par_phase.hpp"
#include "opt_ind_phase.hpp"


/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
template<class sym_type,
         bool opt_bwt>
void grl_bwt_algo(std::string &i_file, std::string& o_file, tmp_workspace & tmp_ws,
                  size_t n_threads, float hbuff_frac, uint8_t b_p_r){


    if constexpr (opt_bwt) {
        std::cout<<"This option is broken"<<std::endl;
        exit(0);
        /*std::cout<<"Constructing the optimal BCR BWT of "<<i_file<<std::endl;
        size_t p_rounds = opt_algo::par_phase(i_file, n_threads, hbuff_size, str_coll, tmp_ws);
        switch (b_p_r) {
            case 0:
                opt_algo::ind_phase<0>(tmp_ws, p_rounds);
                break;
            case 1:
                opt_algo::ind_phase<1>(tmp_ws, p_rounds);
                break;
            case 2:
                opt_algo::ind_phase<2>(tmp_ws, p_rounds);
                break;
            case 3:
                opt_algo::ind_phase<3>(tmp_ws, p_rounds);
                break;
            case 4:
                opt_algo::ind_phase<4>(tmp_ws, p_rounds);
                break;
            default:
                opt_algo::ind_phase<5>(tmp_ws, p_rounds);
        };*/
    }else{
        size_t p_rounds = exact_algo::par_phase<sym_type>(i_file, n_threads, hbuff_frac,  tmp_ws);

        switch (b_p_r) {
            case 0:
                exact_algo::ind_phase<0>(tmp_ws, p_rounds);
                break;
            case 1:
                exact_algo::ind_phase<1>(tmp_ws, p_rounds);
                break;
            case 2:
                exact_algo::ind_phase<2>(tmp_ws, p_rounds);
                break;
            case 3:
                exact_algo::ind_phase<3>(tmp_ws, p_rounds);
                break;
            case 4:
                exact_algo::ind_phase<4>(tmp_ws, p_rounds);
                break;
            default:
                exact_algo::ind_phase<5>(tmp_ws, p_rounds);
        };
    }

    std::filesystem::rename(tmp_ws.get_file("bwt_lev_0"), o_file);
    std::cout<<"The resulting BCR BWT was stored in "<<o_file<<std::endl;
}
#endif //GRLBWT_HPP
