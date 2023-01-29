//
// Created by Diaz, Diego on 23.11.2021.
//
#include "exact_ind_phase.hpp"
#include "bwt_io.h"
#include <filesystem>
#ifdef __linux__
#include <malloc.h>
#endif
//#include "malloc_count.h"

namespace exact_algo {

//extract freq symbols from bwt[j] onwards and put them in new_bwt
    void extract_rl_syms(bwt_buff_writer &bwt_buff, bwt_buff_writer &new_bwt_buff, size_t &j, size_t freq) {

        size_t tmp_freq, sym;
        while (freq > 0) {

            bwt_buff.read_run(j, sym, tmp_freq);

            if (tmp_freq <= freq) {
                freq -= tmp_freq;
                j++;
            } else {
                bwt_buff.dec_freq(j, freq);
                tmp_freq = std::exchange(freq, 0);
            }

            if (new_bwt_buff.size() > 0 && new_bwt_buff.last_sym() == sym) {
                new_bwt_buff.inc_freq_last(tmp_freq);
            } else {
                new_bwt_buff.push_back(sym, tmp_freq);
            }
        }
    }

    size_t compute_hocc_size(dictionary &dict, bv_rs_t &hocc_rs, vector_t &hocc_buckets,
                                   size_t p_round, tmp_workspace &ws) {

        std::string prev_bwt_f = ws.get_file("bwt_lev_" + std::to_string(p_round + 1));
        bwt_buff_reader bwt_buff(prev_bwt_f);

        size_t sym, pos, dummy_sym = dict.bwt_dummy + 1, left_sym, freq = 0, al_b, fr_b, bps, n_runs = 0;
        al_b = INT_CEIL(sym_width(dict.alphabet), 8);
        fr_b = INT_CEIL(sym_width(dict.max_sym_freq), 8);
        bps = al_b + fr_b;
        size_t tot_bytes = bps * hocc_rs(dict.phrases_has_hocc.size());
        auto *hocc_counts = (char *) malloc(tot_bytes);
        memset(hocc_counts, 0, tot_bytes);

        char *ptr;
        for (size_t i = 0; i < bwt_buff.size(); i++) {

            sym = bwt_buff.read_sym(i);

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

        bwt_buff.close();
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
        free(hocc_counts);
        return acc;
    }

    template<uint8_t b_f_r>
    void infer_lvl_bwt(tmp_workspace &ws, size_t p_round) {

        dictionary dict;
        std::string dict_file = ws.get_file("dict_lev_" + std::to_string(p_round));
        sdsl::load_from_file(dict, dict_file);
        bv_rs_t hocc_rs(&dict.phrases_has_hocc);

        size_t sym, left_sym, pos, freq, rank, dummy_sym = dict.bwt_dummy + 1, max_run_len = (1UL << (b_f_r * 8)) - 1;

        std::cout << "    Computing the number of induced symbols" << std::flush;
        auto start = std::chrono::steady_clock::now();
        vector_t hocc_buckets;
        size_t n_runs = compute_hocc_size(dict, hocc_rs, hocc_buckets, p_round, ws);

        size_t al_b = INT_CEIL(sym_width(dict.alphabet), 8);
        size_t fr_b = 1;
        size_t bps = al_b + fr_b;

        bit_hash_table<uintptr_t, sizeof(uintptr_t) * 8, size_t, 8, true> ht;

        auto *hocc = (char *) malloc(n_runs * bps);
        memset(hocc, 0, n_runs * bps);
        char *hocc_ptr;

        std::string prev_bwt_f = ws.get_file("bwt_lev_" + std::to_string(p_round + 1));
        bwt_buff_writer bwt_buff(prev_bwt_f, std::ios::in);
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 9);

        std::cout << "    Performing the induction from the previous BWT" << std::flush;
        start = std::chrono::steady_clock::now();
        for (size_t i = 0; i < bwt_buff.size(); i++) {

            bwt_buff.read_run(i, sym, freq);

            if (dict.phrases_has_hocc[sym]) {
                rank = hocc_rs(sym);
                hocc_ptr = hocc + hocc_buckets.read(rank) * bps;

                if (hocc_ptr != hocc && memcmp(hocc_ptr - bps, &dummy_sym, al_b) == 0) {
                    size_t new_freq = 0;
                    memcpy(&new_freq, hocc_ptr - fr_b, fr_b);

                    if (new_freq == 0) {
                        auto ptr_addr = reinterpret_cast<uintptr_t>(hocc_ptr - bps);
                        auto res = ht.find(&ptr_addr, sizeof(ptr_addr) * 8);
                        assert(res.second);
                        new_freq = 0;
                        ht.get_value_from(res.first, new_freq);
                        new_freq += freq;
                        ht.insert_value_at(res.first, new_freq);
                    } else {
                        new_freq += freq;
                        if (new_freq > max_run_len) {
                            auto ptr_addr = reinterpret_cast<uintptr_t>(hocc_ptr - bps);
                            auto res = ht.insert(&ptr_addr, sizeof(ptr_addr) * 8, new_freq);
                            if (!res.second) {
                                size_t tmp = 0;
                                ht.get_value_from(res.first, tmp);
                                new_freq += tmp;
                                ht.insert_value_at(res.first, new_freq);
                            }
                            new_freq = 0;
                        }
                        memcpy(hocc_ptr - fr_b, &new_freq, fr_b);
                    }
                } else {
                    memcpy(hocc_ptr, &dummy_sym, al_b);
                    size_t new_freq = freq;
                    if (new_freq > max_run_len) {
                        auto ptr_addr = reinterpret_cast<uintptr_t>(hocc_ptr);
                        auto res = ht.insert(&ptr_addr, sizeof(ptr_addr) * 8, new_freq);
                        assert(res.second);
                        new_freq = 0;
                        //std::cout<<ptr_addr<<std::endl;
                    }
                    memcpy(hocc_ptr + al_b, &new_freq, fr_b);
                    hocc_buckets[rank]++;
                }
            }

            pos = 2 * sym + 1;
            sym = dict.dict.read(pos);
            while (sym >= dict.alphabet) {
                sym -= dict.alphabet;
                left_sym = dict.dict.read(pos - 1);
                assert(left_sym < dict.alphabet && dict.phrases_has_hocc[sym]);

                rank = hocc_rs(sym);
                hocc_ptr = hocc + hocc_buckets[rank] * bps;
                left_sym++;

                if (hocc_ptr != hocc && memcmp(hocc_ptr - bps, &left_sym, al_b) == 0) {
                    size_t new_freq = 0;
                    memcpy(&new_freq, hocc_ptr - fr_b, fr_b);

                    if (new_freq == 0) /*[[unlikely]]*/ {
                        auto ptr_addr = reinterpret_cast<uintptr_t>(hocc_ptr - bps);
                        auto res = ht.find(&ptr_addr, sizeof(ptr_addr) * 8);
                        assert(res.second);
                        new_freq = 0;
                        ht.get_value_from(res.first, new_freq);
                        new_freq += freq;
                        ht.insert_value_at(res.first, new_freq);
                    } else {
                        new_freq += freq;
                        if (new_freq > max_run_len) /*[[unlikely]]*/ {
                            auto ptr_addr = reinterpret_cast<uintptr_t>(hocc_ptr - bps);
                            auto res = ht.insert(&ptr_addr, sizeof(ptr_addr) * 8, new_freq);
                            if (!res.second) {
                                size_t tmp = 0;
                                ht.get_value_from(res.first, tmp);
                                new_freq += tmp;
                                ht.insert_value_at(res.first, new_freq);
                            }
                            new_freq = 0;
                        }
                        memcpy(hocc_ptr - fr_b, &new_freq, fr_b);
                    }
                } else {
                    memcpy(hocc_ptr, &left_sym, al_b);
                    size_t new_freq = freq;
                    if (new_freq > max_run_len) /*[[unlikely]]*/ {
                        auto ptr_addr = reinterpret_cast<uintptr_t>(hocc_ptr);
                        auto res = ht.insert(&ptr_addr, sizeof(ptr_addr) * 8, new_freq);
                        assert(res.second);
                        new_freq = 0;
                    }
                    memcpy(hocc_ptr + al_b, &new_freq, fr_b);
                    hocc_buckets[rank]++;
                }
                pos = 2 * sym + 1;
                sym = dict.dict.read(pos);
            }
            assert(dict.dict.read(pos - 1) == dict.metasym_dummy);
            bwt_buff.write_sym(i, sym);
        }
        dict.dict.erase();
        hocc_buckets.erase();
        sdsl::util::clear(hocc_rs);

#ifdef __linux__
        malloc_trim(0);
#endif

        end = std::chrono::steady_clock::now();
        report_time(start, end, 2);

        std::cout << "    Assembling the new BWT" << std::flush;
        start = std::chrono::steady_clock::now();
        std::string new_bwt_f = ws.get_file("bwt_lev_" + std::to_string(p_round));

        uint8_t new_al_b = INT_CEIL(sym_width(std::max(dict.alphabet, dict.prev_alphabet)), 8);
        uint8_t new_fr_b = INT_CEIL(sym_width(dict.max_sym_freq), 8);
        bwt_buff_writer new_bwt_buff(new_bwt_f, std::ios::out, new_al_b, new_fr_b);

        std::string p_bwt_file = ws.get_file("pre_bwt_lev_" + std::to_string(p_round));
        bwt_buff_reader p_bwt(p_bwt_file);

        size_t i = 0, j = 0, pbwt_freq, new_bwt_size = 0, ht_addr = 0;
        hocc_ptr = hocc;
        char *last = hocc + n_runs * bps;
        freq = 0;
        bool freq_in_ht;

        while (i < p_bwt.size()) {

            p_bwt.read_run(i, sym, pbwt_freq);

            if (sym >= dict.bwt_dummy) {// an unsolved segment of the preliminary BWT

                if (sym == dict.hocc_dummy) {// copy from the buffer of induced symbols

                    //copy from hocc+bwt
                    while (pbwt_freq > 0) {
                        freq = 0;//this clean the previous bits
                        memcpy(&sym, hocc_ptr, al_b);
                        memcpy(&freq, hocc_ptr + al_b, fr_b);

                        while (sym == 0 && hocc_ptr != last) {
                            hocc_ptr += bps;
                            memcpy(&sym, hocc_ptr, al_b);
                            memcpy(&freq, hocc_ptr + al_b, fr_b);
                        }

                        assert(sym > 0);

                        freq_in_ht = freq == 0;
                        if (freq_in_ht) /*[[unlikely]]*/ {
                            auto ptr_addr = reinterpret_cast<uintptr_t>(hocc_ptr);
                            auto res = ht.find(&ptr_addr, sizeof(ptr_addr) * 8);
                            assert(res.second);
                            ht_addr = res.first;
                            ht.get_value_from(ht_addr, freq);
                        }

                        if (freq <= pbwt_freq) {
                            pbwt_freq -= freq;
                            hocc_ptr += bps;
                        } else {
                            freq -= pbwt_freq;
                            if (freq_in_ht) /*[[unlikely]]*/ {
                                ht.insert_value_at(ht_addr, freq);
                            } else {
                                assert(freq <= max_run_len);
                                memcpy(hocc_ptr + al_b, &freq, fr_b);
                            }
                            freq = std::exchange(pbwt_freq, 0);
                        }

                        sym--;
                        if (sym == dict.bwt_dummy) {//from the bwt i+1
                            extract_rl_syms(bwt_buff, new_bwt_buff, j, freq);
                            new_bwt_size += freq;
                        } else {
                            if (new_bwt_buff.size() > 0 && new_bwt_buff.last_sym() == sym) {
                                new_bwt_buff.inc_freq_last(freq);
                            } else {
                                assert(sym < dict.alphabet);
                                new_bwt_buff.push_back(sym, freq);
                            }
                            new_bwt_size += freq;
                        }
                    }
                } else {//copy from the bwt i+1
                    assert(sym == dict.bwt_dummy);
                    extract_rl_syms(bwt_buff, new_bwt_buff, j, pbwt_freq);
                    new_bwt_size += pbwt_freq;
                }
            } else {// a segment in the preliminary BWT that was already solved
                if (new_bwt_buff.size() > 0 && new_bwt_buff.last_sym() == sym) {
                    new_bwt_buff.inc_freq_last(pbwt_freq);
                } else {
                    assert(sym < dict.alphabet);
                    new_bwt_buff.push_back(sym, pbwt_freq);
                }
                new_bwt_size += pbwt_freq;
            }
            i++;
        }

        p_bwt.close(true);
        bwt_buff.close(true);
        new_bwt_buff.close();
        if (remove(dict_file.c_str())) {
            std::cout << "Error trying to remove file " << dict_file << std::endl;
        }
        end = std::chrono::steady_clock::now();
        report_time(start, end, 26);

        std::cout << "    Stats:       " << std::endl;
        std::cout << "      BWT size (n):                  " << new_bwt_size << std::endl;
        std::cout << "      Number of runs (r):            " << new_bwt_buff.size() << std::endl;
        std::cout << "      n/r:                           " << double(new_bwt_size) / double(new_bwt_buff.size())<< std::endl;
        std::cout << "      Run stats:                     " << std::endl;
        std::cout << "        Bytes for the symbol:        " << (int) new_al_b << std::endl;
        std::cout << "        Bytes for the length:        " << (int) new_fr_b << std::endl;
        if (ht.size() > 0) {
            std::cout << "        Runs with length overflow:   " << ht.size() << " ("
                      << (double(ht.size()) / double(new_bwt_buff.size())) * 100 << "%)" << std::endl;
        } else {
            std::cout << "        Runs with length overflow:   0" << std::endl;
        }
        free(hocc);
    }

    void infer_lvl_bwt(tmp_workspace &ws, size_t p_round) {

        dictionary dict;
        std::string dict_file = ws.get_file("dict_lev_" + std::to_string(p_round));
        sdsl::load_from_file(dict, dict_file);
        bv_rs_t hocc_rs(&dict.phrases_has_hocc);

        //TODO testing
        /*{
            std::string tmp1 = ws.get_file("bwt_lev_" + std::to_string(p_round + 1));
            std::string tmp2 = ws.get_file("pre_bwt_lev_" + std::to_string(p_round));
            bwt_buff_reader bwt1(tmp1);
            bwt_buff_reader bwt2(tmp2);
            for (size_t i = 0; i < bwt1.size(); i++) {
                size_t a, b;
                bwt1.read_run(i, a, b);
                std::cout << a << " -> " << b << std::endl;
            }
            std::cout <<" "<< std::endl;
            for (size_t i = 0; i < bwt2.size(); i++) {
                size_t a, b;
                bwt2.read_run(i, a, b);
                if(a==dict.sym_dummy){
                    std::cout <<  "*  -> " << b << std::endl;
                } else if(a==dict.sym_end_string)
                    std::cout <<  "$  -> " << b << std::endl;
                else{
                    std::cout << a << " -> " << b << std::endl;
                }
            }
            std::cout << " fin " << std::endl;
            bwt1.close();
            bwt2.close();
        }*/
        //

        size_t sym, left_sym, pos, freq, rank, dummy_sym = dict.bwt_dummy + 1;

        std::cout << "    Computing the number of induced symbols" << std::flush;
        auto start = std::chrono::steady_clock::now();
        vector_t hocc_buckets;
        size_t n_runs = compute_hocc_size(dict, hocc_rs, hocc_buckets, p_round, ws);

        size_t al_b = INT_CEIL(sym_width(dict.alphabet), 8);
        size_t fr_b = INT_CEIL(sym_width(dict.max_sym_freq), 8);
        size_t bps = al_b + fr_b;

        auto *hocc = (char *) malloc(n_runs * bps);
        memset(hocc, 0, n_runs * bps);
        char *hocc_ptr;

        std::string prev_bwt_f = ws.get_file("bwt_lev_" + std::to_string(p_round + 1));
        bwt_buff_writer bwt_buff(prev_bwt_f, std::ios::in);
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 9);

        std::cout << "    Performing the induction from the previous BWT" << std::flush;
        start = std::chrono::steady_clock::now();
        for (size_t i = 0; i < bwt_buff.size(); i++) {

            bwt_buff.read_run(i, sym, freq);

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
                    hocc_buckets[rank]++;
                }
            }

            pos = 2 * sym + 1;
            sym = dict.dict.read(pos);
            while (sym >= dict.alphabet) {
                sym -= dict.alphabet;
                left_sym = dict.dict.read(pos - 1);
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
                    hocc_buckets[rank]++;
                }
                pos = 2 * sym + 1;
                sym = dict.dict.read(pos);
            }
            assert(dict.dict.read(pos - 1) == dict.metasym_dummy);
            bwt_buff.write_sym(i, sym);
        }

        dict.dict.erase();
        hocc_buckets.erase();
        sdsl::util::clear(hocc_rs);

#ifdef __linux__
        malloc_trim(0);
#endif

        end = std::chrono::steady_clock::now();
        report_time(start, end, 2);

        std::cout << "    Assembling the new BWT" << std::flush;
        start = std::chrono::steady_clock::now();
        std::string new_bwt_f = ws.get_file("bwt_lev_" + std::to_string(p_round));

        uint8_t new_al_b = INT_CEIL(sym_width(std::max(dict.alphabet, dict.prev_alphabet)), 8);
        bwt_buff_writer new_bwt_buff(new_bwt_f, std::ios::out, new_al_b, fr_b);

        std::string p_bwt_file = ws.get_file("pre_bwt_lev_" + std::to_string(p_round));
        bwt_buff_reader p_bwt(p_bwt_file);

        size_t i = 0, j = 0, pbwt_freq, new_bwt_size = 0;
        hocc_ptr = hocc;
        char *last = hocc + n_runs * bps;
        freq = 0;

        while (i < p_bwt.size()) {

            p_bwt.read_run(i, sym, pbwt_freq);

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
                            extract_rl_syms(bwt_buff, new_bwt_buff, j, freq);
                        } else {
                            if (new_bwt_buff.size() > 0 && new_bwt_buff.last_sym() == sym) {
                                new_bwt_buff.inc_freq_last(freq);
                            } else {
                                assert(sym < dict.alphabet);
                                new_bwt_buff.push_back(sym, freq);
                            }
                        }
                        new_bwt_size += freq;
                    }
                } else {//copy from the bwt i+1
                    extract_rl_syms(bwt_buff, new_bwt_buff, j, pbwt_freq);
                    new_bwt_size += pbwt_freq;
                }
            } else {// a segment in the preliminary BWT that was already solved
                if (new_bwt_buff.size() > 0 && new_bwt_buff.last_sym() == sym) {
                    new_bwt_buff.inc_freq_last(pbwt_freq);
                } else {
                    assert(sym < dict.alphabet);
                    new_bwt_buff.push_back(sym, pbwt_freq);
                }
                new_bwt_size += pbwt_freq;
            }
            i++;
        }

        p_bwt.close(true);
        bwt_buff.close(true);
        new_bwt_buff.close();
        if (remove(dict_file.c_str())) {
            std::cout << "Error trying to remove file " << dict_file << std::endl;
        }
        end = std::chrono::steady_clock::now();
        report_time(start, end, 26);

        std::cout << "    Stats:       " << std::endl;
        std::cout << "      BWT size (n):       " << new_bwt_size << std::endl;
        std::cout << "      Number of runs (r): " << new_bwt_buff.size() << std::endl;
        std::cout << "      n/r:                " << double(new_bwt_size) / double(new_bwt_buff.size()) << std::endl;
        std::cout << "      Bytes per run:      " << bps << std::endl;
        std::cout << "        Bytes per symbol: " << al_b << std::endl;
        std::cout << "        Bytes per freq.:  " << fr_b << std::endl;
        free(hocc);
    }

    template<class sym_type>
    void parse2bwt_int(tmp_workspace &ws, dictionary& dict, size_t& p_round) {

        std::string parse_file = ws.get_file("tmp_input");
        std::ifstream c_vec(parse_file, std::ifstream::binary);
        c_vec.seekg(0, std::ifstream::end);
        size_t tot_bytes = c_vec.tellg();
        c_vec.seekg(0, std::ifstream::beg);
        auto *buffer = reinterpret_cast<sym_type *>(malloc(BUFFER_SIZE));
        size_t read_bytes = 0;
        size_t len = tot_bytes / sizeof(sym_type);

        size_t sb = INT_CEIL(sym_width(std::max(dict.prev_alphabet, dict.alphabet)) + 1, 8);
        size_t fb = INT_CEIL(sym_width(len) + 1, 8);

        std::string bwt_lev_file = ws.get_file("bwt_lev_" + std::to_string(p_round+1));
        bwt_buff_writer bwt_buff(bwt_lev_file, std::ios::out, sb, fb);

        while (read_bytes < tot_bytes) {
            c_vec.read((char *) buffer, BUFFER_SIZE);
            read_bytes += c_vec.gcount();
            assert((c_vec.gcount() % sizeof(sym_type)) == 0);

            for (size_t i = 0; i < c_vec.gcount() / sizeof(sym_type); i++) {
                buffer[i]>>=1UL;//the symbols are shifted by one as I use the first bit to mark repeated elements
                if (bwt_buff.size() == 0 || buffer[i] != bwt_buff.last_sym()) {
                    bwt_buff.push_back(buffer[i], 1);
                } else {
                    assert(buffer[i] == bwt_buff.last_sym());
                    bwt_buff.inc_freq_last(1);
                }
            }
        }
        c_vec.close();
        free(buffer);
        bwt_buff.close();

        if (remove(parse_file.c_str())) {
            std::cout << "Error trying to delete file " << parse_file << std::endl;
        }

        std::cout << "    Stats:       " << std::endl;
        std::cout << "      BWT size (n):       " << len << std::endl;
        std::cout << "      Number of runs (r): " << bwt_buff.size() << std::endl;
        std::cout << "      n/r:                " << double(len) / double(bwt_buff.size()) << std::endl;
        std::cout << "      Bytes per run:      " << (sb + fb) << std::endl;
        std::cout << "        Bytes per symbol: " << sb << std::endl;
        std::cout << "        Bytes per freq.:  " << fb << std::endl;
    }

    void parse2bwt(tmp_workspace &ws, size_t& p_round) {

        std::string dict_file = ws.get_file("dict_lev_" + std::to_string(p_round));
        dictionary dict;
        sdsl::load_from_file(dict, dict_file);
        size_t bps = sym_width(dict.n_phrases)+1;

        if(bps<=8){
            parse2bwt_int<uint8_t>(ws, dict, p_round);
        }else if(bps<=16){
            parse2bwt_int<uint16_t>(ws, dict, p_round);
        } else if(bps<=32){
            parse2bwt_int<uint32_t>(ws, dict, p_round);
        } else{
            parse2bwt_int<uint64_t>(ws, dict, p_round);
        }

        p_round++;
    }

    template<uint8_t b_f_r> //b_f_r = number of bytes to encode the length of a BWT run
    void ind_phase(tmp_workspace &ws, size_t p_round) {

        static_assert(b_f_r >= 0 && b_f_r <= 5);

        std::cout << "Inferring the BWT" << std::endl;
        parse2bwt(ws, p_round);

        while (p_round-- > 0) {
            std::cout << "  Inducing the BWT for parse " << (p_round + 1) << std::endl;
            auto start = std::chrono::steady_clock::now();
            if constexpr (b_f_r == 0) {
                infer_lvl_bwt(ws, p_round);
            } else {
                infer_lvl_bwt<b_f_r>(ws, p_round);
            }
            auto end = std::chrono::steady_clock::now();
            report_time(start, end, 4);
            report_mem_peak();
            //malloc_count_print_status();
            //malloc_count_reset_peak();
        }
    }

    template void ind_phase<0>(tmp_workspace &ws, size_t p_round);
    template void ind_phase<1>(tmp_workspace &ws, size_t p_round);
    template void ind_phase<2>(tmp_workspace &ws, size_t p_round);
    template void ind_phase<3>(tmp_workspace &ws, size_t p_round);
    template void ind_phase<4>(tmp_workspace &ws, size_t p_round);
    template void ind_phase<5>(tmp_workspace &ws, size_t p_round);
}