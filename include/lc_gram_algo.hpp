//
// Created by diego on 06-03-20.
//

#ifndef LG_COMPRESSOR_LMS_ALGO_H
#define LG_COMPRESSOR_LMS_ALGO_H

#include "grammar.hpp"
#include "common.h"
#include "lc_parsings.hpp"

template<class sym_t>
struct parse_data_t {

    i_file_stream<sym_t>       ifs;
    o_file_stream<size_t>      ofs;
    const sdsl::int_vector<2>& phrase_desc;

    phrase_map_t&              m_map;
    size_t                     start;
    size_t                     end;
    const uint8_t              sym_width;
    dict_t                     thread_dict;

    parse_data_t(std::string &i_file_, std::string &o_file_, phrase_map_t &m_map_,
                 size_t start_, size_t end_,
                 const size_t &alph,
                 const size_t &hb_size, void *hb_addr,
                 const sdsl::int_vector<2> &phrase_desc_) : ifs(i_file_, BUFFER_SIZE),
                                                            ofs(o_file_, BUFFER_SIZE, std::ios::out),
                                                            phrase_desc(phrase_desc_),
                                                            m_map(m_map_),
                                                            start(start_),
                                                            end(end_),
                                                            sym_width(sdsl::bits::hi(alph)+1),
                                                            thread_dict(hb_size, o_file_ + "_phrases", 0.8, hb_addr) {
        //TODO for the moment the input string has to have a sep_symbol appended at the end
        //TODO assertion : sep_symbols cannot be consecutive
    };

    /*inline void hash_phrase(string_t& phrase, bool is_full_str) {
        phrase.mask_tail();
        if(!is_full_str){
            thread_dict.insert(phrase.data(), phrase.n_bits(), false);
        }
    };

    inline void store_phrase(string_t& phrase, bool is_full_str){
        phrase.mask_tail();
        if(!is_full_str){
            auto res = m_map.find(phrase.data(), phrase.n_bits());
            assert(res.second);
            ofs.push_back(res.first.value()>>1UL);
        }else{
            assert(phrase.size()==1 && is_suffix(phrase[0]));
            ofs.push_back(phrase[0]);
        }
    };*/

    inline bool is_suffix(sym_t symbol) const{
        return phrase_desc[symbol] & 2;
    }
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
                data.ofs.push_back(res.first.value()>>1UL);
            }else{
                //assert(phrase.size()==1 && is_suffix(phrase[0]));
                data.ofs.push_back(phrase[0]);
            }
        };
        parser(data.ifs, data.start, data.end, task);
        data.ofs.close();
        data.ifs.close();
        pthread_exit(nullptr);
    };
};

/***
 * Find the previous S* suffix before S[idx]
 * @tparam sym_t : type of symbol in S
 * @param idx : the scan starts from this position
 * @param ifs : input file stream for S
 * @param phrase_desc : input array with the class of each symbol S
 * @return
 */
template<class sym_t>
long long prev_lms_sym(long long idx, i_file_stream<sym_t>& ifs, sdsl::int_vector<2>& phrase_desc) {
    bool type, prev_type;
    size_t sym, prev_sym = ifs.read(idx);

    //get the suffix type of the symbol at idx
    if ((phrase_desc[prev_sym] & 2UL)) {
        return idx;
    } else {
        long long pos = idx;
        sym = ifs.read(++pos);
        while (!(phrase_desc[sym] & 2UL) && prev_sym == sym) {
            sym = ifs.read(++pos);
        }
        prev_type = prev_sym < sym;
    }

    //reach the next LMS-symbol
    while(true) {
        idx--;
        if(idx<0) return idx;

        sym = ifs.read(idx);
        type = sym==prev_sym? prev_type : sym < prev_sym;

        if(phrase_desc[sym] & 2U) return idx;
        if(type==L_TYPE && prev_type==S_TYPE) return idx+1;

        prev_sym = sym;
        prev_type = type;
    }
}

template<class sym_type> std::vector<std::pair<size_t, size_t>>
compute_thread_ranges(size_t n_threads, std::string& i_file, sdsl::int_vector<2>& phrase_desct);

/*template<class sym_t>
void * hash_phrases(void * data);
template<class sym_t>
void * record_phrases(void *data);*/

void build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size,
                      gram_info_t &p_gram, alpha_t alphabet, sdsl::cache_config &config);
template<class sym_type>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                         gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                         sdsl::int_vector<2> &phrase_desc, sdsl::cache_config &config);
void join_parse_chunks(const std::string &output_file,
                       std::vector<std::string> &chunk_files);
void join_thread_phrases(phrase_map_t& mp_map, std::vector<std::string> &chunk_files);

void assign_ids(phrase_map_t &mp_map, size_t max_sym, key_wrapper &key_w, ivb_t &r,
                bvb_t &r_lim, size_t n_threads, sdsl::cache_config &config);


#endif //LG_COMPRESSOR_LMS_ALGO_H