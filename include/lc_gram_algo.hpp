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
    phrase_map_t          thread_dict;

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

struct dictionary {
    typedef size_t size_type;
    size_t min_sym{};
    size_t max_sym{};
    size_t alphabet{};
    size_t n_phrases{};
    size_t t_size{};
    vector_t dict;
    vector_t freqs;
    bv_t     d_lim;
    vector_t phrases_ptr;
    bv_t phrases_has_hocc;
    vector_t hocc_buckets;
    bv_t* desc_bv=nullptr;

    dictionary()=default;
    dictionary(phrase_map_t &mp_map, size_t _min_sym, size_t _max_sym,
               key_wrapper &key_w, size_t dict_syms, size_t max_freq,
               bv_t& is_suffix_bv, size_t _t_size): min_sym(_min_sym),
                                                    max_sym(_max_sym),
                                                    alphabet(max_sym-min_sym+1),
                                                    n_phrases(mp_map.size()),
                                                    t_size(_t_size),
                                                    dict(dict_syms, 0, sdsl::bits::hi(alphabet)+1),
                                                    freqs(n_phrases, 0, sdsl::bits::hi(max_freq)+1),
                                                    d_lim(dict_syms, false),
                                                    phrases_has_hocc(dict.size(), false),
                                                    desc_bv(&is_suffix_bv){
        size_t j=0, k=0, freq;
        for (auto const &ptr : mp_map) {
            for(size_t i=key_w.size(ptr);i-->0;){
                dict[j] = key_w.read(ptr, i);
                assert(key_w.read(ptr,i)==dict[j]);
                d_lim[j++] = false;
            }
            d_lim[j-1] = true;

            freq = 0;
            mp_map.get_value_from(ptr, freq);
            freqs[k++] = freq;
        }
        assert(j==dict_syms);
    }

    [[nodiscard]] inline bool is_suffix(size_t sym) const{
        return (*desc_bv)[sym];
    };

    size_type serialize(std::ostream& out, sdsl::structure_tree_node * v=nullptr, std::string name="") const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child( v, name, sdsl::util::class_name(*this));
        size_type written_bytes = sdsl::write_member(alphabet, out, child, "alphabet");
        written_bytes+= sdsl::write_member(n_phrases, out, child, "n_phrases");
        written_bytes+= sdsl::write_member(t_size, out, child, "t_size");
        dict.serialize(out, child);
        phrases_ptr.serialize(out, child);
        phrases_has_hocc.serialize(out, child);
        hocc_buckets.serialize(out, child);
        return written_bytes;
    }

    void load(std::istream& in){
        sdsl::read_member(alphabet, in);
        sdsl::read_member(n_phrases, in);
        sdsl::read_member(t_size, in);
        dict.load(in);
        phrases_ptr.load(in);
        phrases_has_hocc.load(in);
        hocc_buckets.load(in);
    }
};

template<typename parse_data_t,
         typename parser_t>
struct hash_functor{
    void operator()(parse_data_t& data, parser_t& parser){
        auto task = [&](string_t& phrase){
            phrase.mask_tail();
            auto res = data.thread_dict.insert(phrase.data(), phrase.n_bits(), 1);
            if(!res.second){
                size_t val;
                data.thread_dict.get_value_from(res.first, val);
                val++;
                data.thread_dict.insert_value_at(res.first, val);
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
        auto task = [&](string_t& phrase){
            phrase.mask_tail();
            auto res = data.m_map.find(phrase.data(), phrase.n_bits());
            assert(res.second);
            size_t id = 0;
            data.m_map.get_value_from(res.first, id);
            data.ofs.push_back(id);
        };
        parser(data.ifs, data.start, data.end, task);
        data.ofs.close();
        data.ifs.close();
        pthread_exit(nullptr);
    };
};

template<template<class, class> class lc_parser_t>
void build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size,
                    gram_info_t &p_gram, alpha_t alphabet, sdsl::cache_config &config);
template<class parser_t, class out_sym_t=size_t>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                         gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                         bv_t &phrase_desc, sdsl::cache_config &config);
void join_parse_chunks(const std::string &output_file,
                       std::vector<std::string> &chunk_files);
std::pair<size_t, size_t> join_thread_phrases(phrase_map_t& map, std::vector<std::string> &files);

void assign_ids(phrase_map_t &mp_map, ivb_t &r, bvb_t &r_lim, dictionary &dict, gram_info_t &p_gram,
                sdsl::cache_config &config);


#endif //LG_COMPRESSOR_LMS_ALGO_H