//
// Created by diego on 06-03-20.
//

#ifndef LG_COMPRESSOR_LMS_ALGO_H
#define LG_COMPRESSOR_LMS_ALGO_H

#include "common.h"
#include "utils.h"
#include "lc_parsers.hpp"

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

struct parsing_info{
    size_t lms_phrases=0; //number of LMS phrases in the current parsing round
    size_t tot_phrases=0; //total number of phrases generated in the current parsing round
    size_t p_round=0; //parsing round
    size_t prev_alph=0; //size of the alphabet in the previous parsing round
    size_t max_sym_freq=0; // most repeated symbol in the input text of the round
};

struct dictionary {
    typedef size_t size_type;
    size_t alphabet{};       //alphabet of the dictionary
    size_t p_alpha_size{};   //size of the previous alphabet
    size_t n_phrases{};      //number of LMS phrases in the dictionary
    size_t t_size{};         //size of the text from the dictionary was generated
    size_t max_sym_freq{};   //maximum symbol frequency in the parse from which the dictionary was created
    size_t dict_dummy{};     //dummy symbol for the dictionary
    vector_t dict;           //list of phrases in the dictionary
    vector_t freqs;          //frequency of every dictionary phrase in the original text
    bv_t     d_lim;
    bv_t phrases_has_hocc;   //mark the phrases with hidden occurrences
    bv_t* desc_bv=nullptr;

    dictionary()=default;
    dictionary(phrase_map_t &mp_map, size_t dict_syms,
               size_t max_freq, bv_t& is_suffix_bv, size_t _t_size,
               size_t _p_alph_size, size_t _max_sym_freq): alphabet(is_suffix_bv.size()),
                                                           p_alpha_size(_p_alph_size),
                                                           n_phrases(mp_map.size()),
                                                           t_size(_t_size),
                                                           max_sym_freq(_max_sym_freq),
                                                           dict(dict_syms, 0, sdsl::bits::hi(alphabet)+1),
                                                           freqs(n_phrases, 0, sdsl::bits::hi(max_freq)+1),
                                                           d_lim(dict_syms, false),
                                                           phrases_has_hocc(dict.size(), false),
                                                           desc_bv(&is_suffix_bv){

        key_wrapper key_w{dict.width(), mp_map.description_bits(), mp_map.get_data()};
        size_t j=0, k=0, freq;
        for (auto const &ptr : mp_map) {
            for(size_t i=key_w.size(ptr);i-->0;){
                dict[j] = key_w.read(ptr, i);
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
        written_bytes+= sdsl::write_member(p_alpha_size, out, child, "p_alpha_size");
        written_bytes+= sdsl::write_member(n_phrases, out, child, "n_phrases");
        written_bytes+= sdsl::write_member(t_size, out, child, "t_size");
        written_bytes+= sdsl::write_member(max_sym_freq, out, child, "max_sym_freq");
        written_bytes+= sdsl::write_member(dict_dummy, out, child, "dummy_sym");
        dict.serialize(out, child);
        phrases_has_hocc.serialize(out, child);
        return written_bytes;
    }

    void load(std::istream& in){
        sdsl::read_member(alphabet, in);
        sdsl::read_member(p_alpha_size, in);
        sdsl::read_member(n_phrases, in);
        sdsl::read_member(t_size, in);
        sdsl::read_member(max_sym_freq, in);
        sdsl::read_member(dict_dummy, in);
        dict.load(in);
        phrases_has_hocc.load(in);
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

            //auto res = data.m_map.find(phrase.data(), phrase.n_bits());
            auto res = data.m_map.new_find(phrase.data(), phrase.n_bits());
            //assert(res.second==res2.second);

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

void dict2gram(dictionary &dict, phrase_map_t& phrases_ht, vector_t& s_sa, bv_t& phr_marks,
               parsing_info& p_info, tmp_workspace& ws);
void get_pre_bwt(dictionary &dict, vector_t &sa, parsing_info& p_info, bv_t& phr_marks,
                 phrase_map_t& new_phrases_ht, tmp_workspace& ws);
template<template<class, class> class lc_parser_t>
size_t build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size,
                     str_collection& str_coll, tmp_workspace &ws);
template<class parser_t, class out_sym_t=size_t>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                         parsing_info &p_info, bv_t &phrase_desc, tmp_workspace &ws);
void join_parse_chunks(const std::string &output_file, std::vector<std::string> &chunk_files);
std::pair<size_t, size_t> join_thread_phrases(phrase_map_t& map, std::vector<std::string> &files);

size_t process_dictionary(dictionary &dict, parsing_info &p_info, tmp_workspace &config);


#endif //LG_COMPRESSOR_LMS_ALGO_H