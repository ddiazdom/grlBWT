//
// Created by diego on 17-09-20.
//

#ifndef LPG_COMPRESSOR_SUFFPAIR_ALGO_H
#define LPG_COMPRESSOR_SUFFPAIR_ALGO_H
#include <iostream>
//#include <mutex>
//#include <condition_variable>

#include "grammar.hpp"
#include "common.h"

//struct pairing_data{
//    size_t                    next_av_rule;
//    size_t                    lim_id;
//    size_t                    n_threads;
//    size_t                    hbuff_size;
//    size_t                    tot_lms_rules;
//    size_t                    gsyms;
//    uint8_t                   s_width;
//
//    std::string               r_file;
//    std::string               nr_file;
//    ivb_t                     new_rules;
//
//    bv_t                      r_lim;
//    bv_t::rank_1_type         r_lim_rs;
//    const bv_t&               rep_sym;
//    const sdsl::cache_config& config;
//
//    pairing_data(gram_info_t& gram_info, bv_t& rep_syms_, size_t n_treads_, size_t hbuff_size_, sdsl::cache_config& config_):
//            rep_sym(rep_syms_),
//            config(config_){
//
//        hbuff_size = hbuff_size_;
//        n_threads = n_treads_;
//        next_av_rule = gram_info.r - 1; //I subtract one as the last rule is the compressed string
//        tot_lms_rules = gram_info.rules_breaks[gram_info.n_p_rounds]-(gram_info.max_tsym+1); // the same
//        lim_id = 2*gram_info.g - gram_info.r;
//        s_width = sdsl::bits::hi(lim_id)+1;
//        r_file = gram_info.rules_file;
//        nr_file = sdsl::cache_file_name("nr_file", config);
//        new_rules = ivb_t(nr_file, std::ios::out, BUFFER_SIZE);
//        gsyms = gram_info.g - gram_info.c;
//
//        sdsl::load_from_file(r_lim, gram_info.rules_lim_file);
//        sdsl::util::init_support(r_lim_rs, &r_lim);
//    }
//
//    ~pairing_data(){
//        if(new_rules.is_open()) new_rules.close(true);
//    }
//};
//struct thread_data{
//    const pairing_data& p_data;
//    const size_t        start;
//    const size_t        end;
//
//    bool                finished;
//    std::string         pl_chunk_file;
//    std::string         ht_file;
//    phrase_map_t&       ht;
//    void *              hb_addr;
//    size_t              hb_size;
//
//    thread_data(const pairing_data& p_data_, const size_t start_, const size_t end_, phrase_map_t& ht_, size_t hb_size_, void * hb_addr_):
//                p_data(p_data_),
//                start(start_),
//                end(end_),
//                ht(ht_),
//                hb_addr(hb_addr_),
//                hb_size(hb_size_){
//        std::stringstream ss;
//        ss << p_data.config.dir<< "/range_" << start << "_" << end;
//        std::string prefix = ss.str();
//
//        ht_file = prefix + "_pairs";
//        pl_chunk_file = prefix+"_pt_chunk_file";
//        finished = false;
//    }
//};
//void create_new_rules(phrase_map_t& ht, pairing_data& p_data);
//void update_grammar(pairing_data& p_data, gram_info_t& gram);
//void prepare_input(gram_info_t& gram_info, bv_t& rep_syms, sdsl::cache_config& config);
//void collect_pairs(thread_data* d, i_file_stream<size_t>& p_list, o_file_stream<size_t>& r);
//void replace_pairs(const phrase_map_t& ht, const pairing_data& p_data, std::string& pl_file, o_file_stream<size_t>& r);
//void * suffpair_thread(void * data);
//void merge_ht_data(std::vector<thread_data>& t_data);
//void suffixpair_int(pairing_data& p_data)
size_t
create_new_rules(phrase_map_t &ht, gram_info_t &p_gram, bv_t &rep_sym, sdsl::cache_config &config, vector_t &rules,
                 bv_t &r_lim);
void compress_grammar(gram_info_t& p_gram, vector_t& rules, bv_t& r_lim, bv_t& rep_sym,
                      phrase_map_t& ht, sdsl::cache_config &config);
void mark_repeated(bv_t& rep_sym, vector_t& rules, bv_t& r_lim, gram_info_t& p_gram);
void suffpair(gram_info_t& p_gram, sdsl::cache_config &config);

#endif //LPG_COMPRESSOR_SUFFPAIR_ALGO_H
