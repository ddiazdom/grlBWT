//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_HPP
#define GRLBWT_HPP

#include "utils.h"
#include "ind_phase.hpp"
#include "par_phase.hpp"

/***
 *
 * @param i_file : input text file
 * @param n_threads : number of working threads
 * @param hbuff_size : buffer size for the hashing step
 */
template<class sym_type, uint8_t bytes_per_run>
void grl_bwt_algo(std::string &i_file, std::string& o_file, size_t n_threads, const std::string& tmp_dir=""){

    logger log(LOG_INFO);

    tmp_workspace tmp_ws(tmp_dir, true, "grl.bwt");
    log.info("Temporary folder: "+tmp_ws.folder());

    size_t p_rounds = par_phase<sym_type>(i_file, n_threads, tmp_ws, log);
    ind_phase<bytes_per_run>(tmp_ws, p_rounds);

    std::filesystem::rename(tmp_ws.get_file("bwt_lev_0"), o_file);
    log.info("The resulting BCR BWT was stored in "+o_file);
}
#endif //GRLBWT_HPP
