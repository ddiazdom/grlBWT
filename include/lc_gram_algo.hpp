//
// Created by diego on 06-03-20.
//

#ifndef LG_COMPRESSOR_LMS_ALGO_H
#define LG_COMPRESSOR_LMS_ALGO_H

#include "common.h"
#include "utils.h"
#include "lc_parsers.hpp"
#include "parsing_strategies.h"

struct dictionary {

    typedef size_t size_type;
    size_t alphabet{};       //alphabet of the dictionary
    size_t p_alpha_size{};   //size of the previous alphabet
    size_t n_phrases{};      //number of LMS phrases in the dictionary
    size_t t_size{};         //size of the text from which the dictionary was generated
    size_t max_sym_freq{};   //maximum symbol frequency in the parse from which the dictionary was created
    size_t dict_dummy{};     //dummy symbol for the dictionary
    vector_t dict;           //list of phrases in the dictionary
    vector_t freqs;          //frequency of every dictionary phrase in the original text
    bv_t     d_lim;
    bv_t phrases_has_hocc;   //mark the phrases with hidden occurrences
    bv_t* desc_bv=nullptr;

    dictionary()=default;
    dictionary(phrase_map_t &mp_map, size_t dict_syms,
               size_t max_freq, bv_t& is_suffix_bv, size_t _t_size, size_t _p_alph_size,
               size_t _max_sym_freq, size_t _dummy_sym): alphabet(is_suffix_bv.size()),
                                                         p_alpha_size(_p_alph_size),
                                                         n_phrases(mp_map.size()),
                                                         t_size(_t_size),
                                                         max_sym_freq(_max_sym_freq),
                                                         dict_dummy(_dummy_sym),
                                                         dict(dict_syms, 0, sym_width(alphabet+dict_syms)),
                                                         freqs(n_phrases, 0, sdsl::bits::hi(max_freq)+1),
                                                         d_lim(dict_syms, false),
                                                         phrases_has_hocc(dict.size(), false),
                                                         desc_bv(&is_suffix_bv){

        key_wrapper key_w{sym_width(alphabet), mp_map.description_bits(), mp_map.get_data()};
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

        //TODO testing
        std::vector<size_t> g_freqs(alphabet, 0);
        for(size_t i=0;i<dict.size();i++){
            if(!d_lim[i] || (d_lim[i] && is_suffix(dict[i]))){
                g_freqs[dict[i]]++;
            }
        }
        size_t i=0;
        bool phrase_uniq;
        size_t n_uniq=0, phrase=0, n_rep=0;
        while(i<dict.size()){
            assert(i==0 || d_lim[i-1]);
            phrase_uniq=true;
            while(!d_lim[i]){
                if(g_freqs[dict[i]]>1){
                    phrase_uniq=false;
                }
                i++;
            }
            assert(d_lim[i]);
            if(is_suffix(dict[i]) && g_freqs[dict[i]]>1){
                phrase_uniq = false;
            }
            i++;
            if(phrase_uniq){
                n_uniq++;
                if(freqs[phrase]>1) n_rep++;
                if(n_phrases==8036990){
                    size_t tmp=i-2;
                    while(tmp>0 && !d_lim[tmp]) tmp--;
                    if(tmp>0) tmp++;
                    assert(tmp==0 || d_lim[tmp-1]);
                    do{
                        std::cout<<dict[tmp]<<" ";
                    }while(!d_lim[tmp++]);
                    std::cout<<" - > "<<freqs[phrase]<<std::endl;
                }
            }
            phrase++;
        }
        std::cout<<"\nThere are "<<n_uniq<<" phrases that have a unique context. From them, "<<n_rep<<" are repeated in the text"<<std::endl;
        //
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
        dict.serialize(out);
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


void dict2gram(dictionary &dict, phrase_map_t& phrases_ht, vector_t& s_sa, bv_t& phr_marks,
               parsing_info& p_info, tmp_workspace& ws);

void get_pre_bwt(dictionary &dict, vector_t &sa, parsing_info& p_info, bv_t& phr_marks,
                 phrase_map_t& new_phrases_ht, tmp_workspace& ws);

template<class value_type>
void dict2gram2(dictionary &dict, value_type *s_sa, size_t sa_size, vector_t& first_symbol, parsing_info& p_info, tmp_workspace& ws);

template<class value_type>
size_t get_pre_bwt2(dictionary &dict, value_type * sa, size_t sa_size, parsing_info& p_info, tmp_workspace& ws);

size_t build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size,
                     str_collection& str_coll, tmp_workspace &ws);

template<class parse_strategy>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file, parse_strategy& p_strat,
                         parsing_info &p_info, bv_t &phrase_desc, tmp_workspace &ws);

size_t process_dictionary(dictionary &dict, parsing_info &p_info, tmp_workspace &config);

template<class vector_type>
size_t process_dictionary_int(dictionary &dict, parsing_info &p_info, tmp_workspace &config);

#endif //LG_COMPRESSOR_LMS_ALGO_H