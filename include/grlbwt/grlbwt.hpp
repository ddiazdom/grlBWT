//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_HPP
#define GRLBWT_HPP

#include "ind_phase.hpp"
#include "par_phase.hpp"

/***
 *
 * @param i_file : input text file
 * @param o_file : output file where the BCR will be stored
 * @param n_threads : number of working threads
 * @param hbuff_frac : memory used during the parsing
 * @param b_p_r : bytes per run length (used during the induction phase)
 */
template<class sym_type>
void grl_bwt_algo(std::string &i_file, std::string& o_file, tmp_workspace & tmp_ws,
                  size_t n_threads, float hbuff_frac, uint8_t b_p_r){

    size_t p_rounds = par_phase<sym_type>(i_file, n_threads, hbuff_frac,  tmp_ws);
    switch (b_p_r) {
        case 0:
            ind_phase<0>(tmp_ws, p_rounds);
            break;
        case 1:
            ind_phase<1>(tmp_ws, p_rounds);
            break;
        case 2:
            ind_phase<2>(tmp_ws, p_rounds);
            break;
        case 3:
            ind_phase<3>(tmp_ws, p_rounds);
            break;
        case 4:
            ind_phase<4>(tmp_ws, p_rounds);
            break;
        default:
            ind_phase<5>(tmp_ws, p_rounds);
    };

    std::filesystem::rename(tmp_ws.get_file("bwt_lev_0"), o_file);
    std::cout<<"The resulting BCR BWT was stored in "<<o_file<<std::endl;
}
#endif //GRLBWT_HPP
