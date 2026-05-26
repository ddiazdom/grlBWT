//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_IND_PHASE_HPP
#define GRLBWT_IND_PHASE_HPP

#include "par_phase.hpp"
#include <cds/logger.h>
#include <cds/workspace.h>
#include <cds/memory_handler.hpp>
#include <cds/rle_vbyte_streams.h>

#ifdef __linux__
#include <malloc.h>
#endif

//extract freq symbols from bwt[j] onwards and put them in new_bwt
inline void extract_rl_syms(rle_vbyte_buff_reader &bwt_buff, rle_vbyte_buff_writer &new_bwt_buff, uint64_t freq) {

    uint64_t tmp_freq, sym;
    while (freq > 0) {
        bwt_buff.read_run(sym, tmp_freq);
        if (tmp_freq <= freq) {
            freq -= tmp_freq;
            bwt_buff.next_run();
        } else {
            bwt_buff.decrease_len(freq);
            tmp_freq = std::exchange(freq, 0);
        }
        new_bwt_buff.push_back(sym, tmp_freq);
    }
}

inline size_t compute_hocc_size(dictionary &dict, bv_rs_t &hocc_rs, vector_t &hocc_buckets, size_t p_round) {

    LOG_DEBUG("Computing the number of induced symbols");
    std::string prev_bwt_f = TMP_FILE_NAME("bwt_lev_" + std::to_string(p_round + 1));
    rle_vbyte_buff_reader bwt_buff(prev_bwt_f);

    size_t pos, dummy_sym = dict.bwt_dummy + 1, left_sym, freq = 0, al_b, fr_b, bps, n_runs = 0;

    al_b = INT_CEIL(sym_width(dict.alphabet), 8);
    fr_b = INT_CEIL(sym_width(dict.max_sym_freq), 8);

    bps = al_b + fr_b;
    size_t tot_bytes = bps * hocc_rs(dict.phrases_has_hocc.size());
    auto *hocc_counts = mem::allocate<char>(tot_bytes);
    memset(hocc_counts, 0, tot_bytes);

    uint64_t sym, len;
    char *ptr;
    while(bwt_buff.has_next()) {
        bwt_buff.next_run(sym, len);
        if (dict.phrases_has_hocc[sym]) {
            ptr = hocc_counts + hocc_rs(sym) * bps;
            if (memcmp(ptr, &dummy_sym, al_b) != 0) {
                memcpy(ptr, &dummy_sym, al_b);
                memcpy(&freq, ptr + al_b, fr_b);
                freq++;
                n_runs++;
                memcpy(ptr + al_b, &freq, fr_b);
            }
        }

        pos = 2 * sym + 1;
        sym = dict.dict.read(pos);

        while (sym >= dict.alphabet) {
            sym -= dict.alphabet;
            left_sym = dict.dict.read(pos - 1) + 1;
            assert(sym < dict.phrases_has_hocc.size());
            ptr = hocc_counts + hocc_rs(sym) * bps;
            if (memcmp(ptr, &left_sym, al_b) != 0) {
                memcpy(ptr, &left_sym, al_b);
                memcpy(&freq, ptr + al_b, fr_b);
                freq++;
                n_runs++;
                memcpy(ptr + al_b, &freq, fr_b);
            }
            pos = 2 * sym + 1;
            sym = dict.dict.read(pos);
        }
    }

    hocc_buckets.set_width(sym_width(std::max<size_t>(n_runs, 1)));
    hocc_buckets.resize(hocc_rs(dict.phrases_has_hocc.size()) + 1);

    size_t acc = 0;
    ptr = hocc_counts;
    for (size_t i = 0; i < hocc_buckets.size() - 1; i++) {
        hocc_buckets.write(i, acc);
        memcpy(&freq, ptr + al_b, fr_b);
        acc += freq;
        ptr += bps;
    }
    assert(acc == n_runs);
    hocc_buckets.write(hocc_buckets.size() - 1, acc);
    mem::deallocate(hocc_counts);
    return acc;
}

inline void induce_from_prev_bwt(size_t p_round, dictionary &dict, char* hocc, bv_rs_t &hocc_rs, vector_t &hocc_buckets) {

    std::string prev_bwt_f = TMP_FILE_NAME("bwt_lev_" + std::to_string(p_round + 1));
    rle_vbyte_buff_updater bwt_buff(prev_bwt_f);

    size_t rank;
    const size_t dummy_sym = dict.bwt_dummy + 1;
    const size_t al_b = INT_CEIL(sym_width(dict.alphabet), 8);
    const size_t fr_b = INT_CEIL(sym_width(dict.max_sym_freq), 8);
    const size_t bps = al_b + fr_b;
    char *hocc_ptr;
    uint64_t sym, freq;

    LOG_DEBUG("Inducing from the previous BWT");
    while(bwt_buff.has_next()) {

        bwt_buff.next_run(sym, freq);

        if (dict.phrases_has_hocc[sym]) {
            rank = hocc_rs(sym);
            hocc_ptr = hocc + hocc_buckets.read(rank) * bps;

            if (hocc_ptr != hocc && memcmp(hocc_ptr - bps, &dummy_sym, al_b) == 0) {
                size_t new_freq = 0;
                memcpy(&new_freq, hocc_ptr - fr_b, fr_b);
                new_freq += freq;
                memcpy(hocc_ptr - fr_b, &new_freq, fr_b);
            } else {
                memcpy(hocc_ptr, &dummy_sym, al_b);
                memcpy(hocc_ptr + al_b, &freq, fr_b);
                ++hocc_buckets[rank];
            }
        }

        size_t pos = 2 * sym + 1;
        sym = dict.dict.read(pos);
        while (sym >= dict.alphabet) {
            sym -= dict.alphabet;
            size_t left_sym = dict.dict.read(pos - 1);
            assert(left_sym < dict.alphabet && dict.phrases_has_hocc[sym]);

            rank = hocc_rs(sym);
            hocc_ptr = hocc + hocc_buckets[rank] * bps;
            left_sym++;

            if (hocc_ptr != hocc && memcmp(hocc_ptr - bps, &left_sym, al_b) == 0) {
                size_t new_freq = 0;
                memcpy(&new_freq, hocc_ptr - fr_b, fr_b);
                new_freq += freq;
                memcpy(hocc_ptr - fr_b, &new_freq, fr_b);
            } else {
                memcpy(hocc_ptr, &left_sym, al_b);
                memcpy(hocc_ptr + al_b, &freq, fr_b);
                ++hocc_buckets[rank];
            }
            pos = 2 * sym + 1;
            sym = dict.dict.read(pos);
        }
        assert(dict.dict.read(pos - 1) == dict.metasym_dummy);
        bwt_buff.update_sym(sym);
    }

    dict.dict.erase();
    hocc_buckets.erase();
    clear(hocc_rs);
    bwt_buff.close();

#ifdef __linux__
     malloc_trim(0);
#endif
}

inline void infer_lvl_bwt(size_t p_round) {

    TRACE_SCOPE();

    dictionary dict;
    std::string dict_file = "dict_lev_" + std::to_string(p_round);
    load_from_file(TMP_FILE_NAME(dict_file), dict);

    bv_rs_t hocc_rs(&dict.phrases_has_hocc);
    vector_t hocc_buckets;
    size_t n_runs = compute_hocc_size(dict, hocc_rs, hocc_buckets, p_round);

    size_t al_b = INT_CEIL(sym_width(dict.alphabet), 8);//number of bytes for the dictionary alphabet
    size_t fr_b = INT_CEIL(sym_width(dict.max_sym_freq), 8);//number of bytes for the dictionary max freq
    size_t bps = al_b + fr_b;
    auto *hocc = mem::allocate<char>(n_runs * bps);
    memset(hocc, 0, n_runs * bps);

    induce_from_prev_bwt(p_round, dict, hocc, hocc_rs, hocc_buckets);

    LOG_DEBUG("Assembling the new BWT");

    //BWT from the previous round, now updated for the current round
    std::string prev_bwt_f = "bwt_lev_" + std::to_string(p_round + 1);
    rle_vbyte_buff_reader prev_bwt(TMP_FILE_NAME(prev_bwt_f));

    //Preliminary BWT from the parsing phase
    std::string prelim_bwt_f = "pre_bwt_lev_" + std::to_string(p_round);
    rle_vbyte_buff_reader prelim_bwt(TMP_FILE_NAME(prelim_bwt_f));

    //new BWT
    std::string new_bwt_f = "bwt_lev_" + std::to_string(p_round);
    rle_vbyte_buff_writer new_bwt(TMP_FILE_NAME(new_bwt_f));

    size_t freq=0;
    char *hocc_ptr = hocc;
    char *last = hocc + n_runs * bps;
    uint64_t pbwt_freq, sym;

    while (prelim_bwt.has_next()) {

        prelim_bwt.next_run(sym, pbwt_freq);

        if (sym >= dict.bwt_dummy) {// an unsolved segment of the preliminary BWT

            if (sym == dict.hocc_dummy) {// copy from the buffer of induced symbols

                //copy from hocc+bwt
                while (pbwt_freq > 0) {
                    memcpy(&sym, hocc_ptr, al_b);
                    memcpy(&freq, hocc_ptr + al_b, fr_b);

                    while (sym == 0 && hocc_ptr != last) {
                        hocc_ptr += bps;
                        memcpy(&sym, hocc_ptr, al_b);
                        memcpy(&freq, hocc_ptr + al_b, fr_b);
                    }

                    assert(sym > 0);

                    if (freq <= pbwt_freq) {
                        pbwt_freq -= freq;
                        hocc_ptr += bps;
                    } else {
                        size_t tmp = 0;
                        memcpy(&tmp, hocc_ptr + al_b, fr_b);
                        assert(tmp == freq);
                        assert(tmp > pbwt_freq);
                        tmp -= pbwt_freq;
                        memcpy(hocc_ptr + al_b, &tmp, fr_b);
                        freq = std::exchange(pbwt_freq, 0);
                    }

                    sym--;
                    if (sym == dict.bwt_dummy) {//from the bwt i+1
                        extract_rl_syms(prev_bwt, new_bwt, freq);
                    } else {
                        new_bwt.push_back(sym, freq);
                    }
                }
            } else {//copy from bwt i+1
                extract_rl_syms(prev_bwt, new_bwt, pbwt_freq);
            }
        } else {// a segment in the preliminary BWT that was already solved
            new_bwt.push_back(sym, pbwt_freq);
        }
    }

    prelim_bwt.close();
    new_bwt.close();
    prev_bwt.close();

    TMP_REMOVE_FILE(dict_file);
    TMP_REMOVE_FILE(prev_bwt_f);
    TMP_REMOVE_FILE(prelim_bwt_f);

    LOG_DEBUG("Stats:");
    LOG_DEBUG("  BWT size (n):         "+std::to_string(new_bwt.size()));
    LOG_DEBUG("  Number of runs (r):   "+std::to_string(new_bwt.tot_runs()));
    LOG_DEBUG("  n/r:                  "+to_string_with_precision(new_bwt.avg_len(), 3));
    mem::deallocate(hocc);
}

template<class sym_type>
void parse2bwt_int(size_t& p_round) {

    std::string parse_file = TMP_FILE_NAME("tmp_input");
    std::ifstream c_vec(parse_file, std::ifstream::binary);
    c_vec.seekg(0, std::ifstream::end);
    size_t tot_bytes = c_vec.tellg();
    c_vec.seekg(0, std::ifstream::beg);

    auto *buffer = mem::allocate<sym_type>(BUFFER_SIZE/sizeof(sym_type));
    size_t read_bytes = 0;

    std::string bwt_lev_file = TMP_FILE_NAME("bwt_lev_" + std::to_string(p_round+1));
    rle_vbyte_buff_writer bwt_buff(bwt_lev_file);

    while (read_bytes < tot_bytes) {
        c_vec.read(reinterpret_cast<char *>(buffer), BUFFER_SIZE);
        read_bytes += c_vec.gcount();
        assert((c_vec.gcount() % sizeof(sym_type)) == 0);

        for (size_t i = 0; i < c_vec.gcount() / sizeof(sym_type); i++) {
            buffer[i]>>=1UL;//the symbols are shifted by one as I use the first bit to mark repeated elements
            bwt_buff.push_back(buffer[i], 1);
        }
    }
    c_vec.close();
    mem::deallocate(buffer);
    bwt_buff.close();

    TMP_REMOVE_FILE("tmp_input");

    LOG_DEBUG("Stats:");
    LOG_DEBUG("  BWT size (n):         "+std::to_string(bwt_buff.size()));
    LOG_DEBUG("  Number of runs (r):   "+std::to_string(bwt_buff.tot_runs()));
    LOG_DEBUG("  n/r:                  "+to_string_with_precision(bwt_buff.avg_len(), 3));
}

inline void parse2bwt(size_t& p_round) {

    LOG_INFO("Computing the deepest recursive BWT");

    SCOPE_INFO();
    TRACE_SCOPE();
    std::string dict_file = "dict_lev_" + std::to_string(p_round);
    dictionary dict;
    load_from_file(TMP_FILE_NAME(dict_file), dict);
    size_t bps = sym_width(dict.n_phrases)+1;

    if(bps<=8){
        parse2bwt_int<uint8_t>(p_round);
    }else if(bps<=16){
        parse2bwt_int<uint16_t>(p_round);
    } else if(bps<=32){
        parse2bwt_int<uint32_t>(p_round);
    } else{
        parse2bwt_int<uint64_t>(p_round);
    }
    p_round++;
}

inline void ind_phase(size_t p_round) {
    LOG_INFO("Completing the BWT");
    SCOPE_INFO();
    parse2bwt(p_round);
    while (p_round-- > 0) {
        LOG_INFO("Level "+std::to_string(p_round + 1));
        SCOPE_INFO();
        infer_lvl_bwt(p_round);
    }
}
#endif //GRLBWT_IND_PHASE_HPP
