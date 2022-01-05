//
// Created by diego on 06-03-20.
//

#ifndef LG_COMPRESSOR_LMS_ALGO_H
#define LG_COMPRESSOR_LMS_ALGO_H

#include "grammar.hpp"
#include "common.h"
#include "lc_parsers.hpp"
#include "LMS_induction.h"

template<class istream_t,
         class out_sym_t=size_t,
         class ostream_t=o_file_stream<out_sym_t>>
struct parse_data_t {

    istream_t             ifs;
    ostream_t             ofs;

    phrase_map_t&         m_map;
    size_t                start;
    size_t                end;
    dict_t                thread_dict;

    parse_data_t(std::string &i_file_, std::string &o_file_, phrase_map_t &m_map_,
                 size_t start_, size_t end_,
                 const size_t &hb_size, void *hb_addr) : ifs(i_file_, BUFFER_SIZE),
                                                         ofs(o_file_, BUFFER_SIZE, std::ios::out),
                                                         m_map(m_map_),
                                                         start(start_),
                                                         end(end_),
                                                         thread_dict(hb_size, o_file_ + "_phrases", 0.7, hb_addr) {
        //TODO for the moment the input string has to have a sep_symbol appended at the end
        //TODO assertion : sep_symbols cannot be consecutive
    };
};

template<typename parse_data_t,
         typename parser_t>
struct hash_functor{
    void operator()(parse_data_t& data, parser_t& parser){
        auto task = [&](string_t& phrase, bool is_full_str){
            phrase.mask_tail();
            if(!is_full_str){
                data.thread_dict.insert(phrase.data(), phrase.n_bits(), false);
            }
        };
        parser(data.ifs, data.start, data.end, task);
        data.thread_dict.flush();
        pthread_exit(nullptr);
    };
};

template<typename parse_data_t,
         typename parser_t>
struct parse_functor{
    void operator()(parse_data_t& data, parser_t& parser){
        auto task = [&](string_t& phrase, bool is_full_str){
            phrase.mask_tail();
            if(!is_full_str){
                auto res = data.m_map.find(phrase.data(), phrase.n_bits());
                assert(res.second);
                size_t id = 0;
                data.m_map.get_value_from(res.first, id);
                data.ofs.push_back(id>>1UL);
            }else{
                data.ofs.push_back(phrase[0]);
            }
        };
        parser(data.ifs, data.start, data.end, task);
        data.ofs.close();
        data.ifs.close();
        pthread_exit(nullptr);
    };
};

template<class parser_t> std::vector<std::pair<size_t, size_t>>
compute_thread_ranges(size_t n_threads, std::string& i_file, parser_t& phrase_desct);

template<template<class, class> class lc_parser_t>
void build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size,
                    gram_info_t &p_gram, alpha_t alphabet, sdsl::cache_config &config);
template<class parser_t, class out_sym_t=size_t>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                         gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                         sdsl::int_vector<2> &phrase_desc, sdsl::cache_config &config);
void join_parse_chunks(const std::string &output_file,
                       std::vector<std::string> &chunk_files);
size_t join_thread_phrases(phrase_map_t& map, std::vector<std::string> &files);

void assign_ids(phrase_map_t &mp_map, size_t max_sym, size_t min_sym, key_wrapper &key_w, bvb_t &r_lim,
                sdsl::cache_config &config, ivb_t &r, size_t dict_syms);


#endif //LG_COMPRESSOR_LMS_ALGO_H