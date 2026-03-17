//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_HPP
#define GRLBWT_HPP

#include "ind_phase.hpp"
#include "par_phase.hpp"
#include <cdt/workspace.h>

/***
 *
 * @param i_file : input text file
 * @param o_file : output file where the BCR will be stored
 * @param n_threads : number of working threads
 * @param hbuff_frac : memory used during the parsing
 * @param tmp_ws: temporary workspace
 * @param log_lvl: verbosity level
 */
template<class sym_type>
void grl_bwt_algo(std::string &i_file, std::string& o_file, size_t n_threads, float hbuff_frac,
                  const std::string& tmp_ws= SYSTEM_TMP,
                  const log_level log_lvl=log_level::INFO){
    Logger::level = log_lvl;
    INIT_TMP_FOLDER(tmp_ws, "grl.bwt", true);
    LOG_INFO("Temporary folder: "+TMP_WORKSPACE);
    size_t p_rounds = par_phase<sym_type>(i_file, n_threads, hbuff_frac);
    ind_phase(p_rounds);
    std::filesystem::rename(TMP_FILE_NAME("bwt_lev_0"), o_file);
}
#endif //GRLBWT_HPP