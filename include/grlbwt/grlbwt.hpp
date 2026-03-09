//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_HPP
#define GRLBWT_HPP

#include "ind_phase.hpp"
#include "par_phase.hpp"
#include "cdt/workspace.h"

/***
 *
 * @param i_file : input text file
 * @param o_file : output file where the BCR will be stored
 * @param n_threads : number of working threads
 * @param hbuff_frac : memory used during the parsing
 * @param b_p_r : bytes per run length (used during the induction phase)
 */
template<class sym_type>
void grl_bwt_algo(std::string &i_file, std::string& o_file, std::string& tmp_ws,
                  size_t n_threads, float hbuff_frac, uint8_t b_p_r){

    INIT_TMP_FOLDER(tmp_ws, "grl.bwt", true);
    LOG_INFO("Temporary folder: "+TMP_WORKSPACE);
    size_t p_rounds = par_phase<sym_type>(i_file, n_threads, hbuff_frac);
    switch (b_p_r) {
        case 0:
            ind_phase<0>(p_rounds);
            break;
        case 1:
            ind_phase<1>(p_rounds);
            break;
        case 2:
            ind_phase<2>(p_rounds);
            break;
        case 3:
            ind_phase<3>(p_rounds);
            break;
        case 4:
            ind_phase<4>(p_rounds);
            break;
        default:
            ind_phase<5>(p_rounds);
    };
    std::filesystem::rename(TMP_FILE_NAME("bwt_lev_0"), o_file);
}
#endif //GRLBWT_HPP