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
    size_t alphabet{};      //alphabet of the dictionary
    size_t prev_alphabet{};  //size of the previous alphabet
    size_t n_phrases{};     //number of LMS phrases in the dictionary
    size_t t_size{};        //size of the text from which the dictionary was generated
    size_t max_sym_freq{};  //maximum symbol frequency in the parse from which the dictionary was created
    size_t end_str_dummy{}; //symbol to mark the end of a string
    size_t bwt_dummy{};     //dummy symbol indicating that the induction is from the BWT i+1
    size_t hocc_dummy{};    //dummy symbol indicating that the induction is from the hocc array
    size_t metasym_dummy{}; //dummy metasymbol for the compressed suffixes
    vector_t dict;          //list of phrases in the dictionary
    vector_t freqs;         //frequency of every dictionary phrase in the original text
    bv_t     d_lim;
    bv_t phrases_has_hocc; //mark the phrases with hidden occurrences
    bv_t* desc_bv=nullptr;

    dictionary()=default;
    dictionary(phrase_map_t &mp_map, size_t dict_syms,
               size_t max_freq, bv_t& is_suffix_bv, size_t _t_size, size_t _p_alph_size,
               size_t _max_sym_freq): alphabet(is_suffix_bv.size()),
                                      prev_alphabet(_p_alph_size),
                                      n_phrases(mp_map.size()),
                                      t_size(_t_size),
                                      max_sym_freq(_max_sym_freq),
                                      end_str_dummy(alphabet),
                                      bwt_dummy(alphabet+1),
                                      hocc_dummy(alphabet+2),
                                      dict(dict_syms, 0, sym_width(alphabet+3+dict_syms)),//the +3 corresponds to the extra symbols I use for the BWT induction
                                      freqs(n_phrases, 0, sym_width(max_freq)),
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

            //For the moment, I won't deal with this problem
            assert((freq & 3UL)!=3);

            if(!(*desc_bv)[dict[j-1]]) freq>>=2UL;
            assert(freq<=max_freq);
            freqs[k++] = freq;
        }
        //std::cout<<"There are "<<cont<<" phrases that cover an entire string"<<std::endl;
        assert(j==dict_syms);

        //TODO testing
        /*std::vector<size_t> g_freqs(alphabet, 0);
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
            }
            phrase++;
        }
        std::cout<<"\nThere are "<<n_uniq<<" phrases that have a unique context. From them, "<<n_rep<<" are repeated in the text"<<std::endl;*/
        //
    }

    [[nodiscard]] inline bool is_suffix(size_t sym) const{
        return (*desc_bv)[sym];
    };

    size_type serialize(std::ostream& out, sdsl::structure_tree_node * v=nullptr, std::string name="") const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child( v, name, sdsl::util::class_name(*this));
        size_type written_bytes = sdsl::write_member(alphabet, out, child, "alphabet");
        written_bytes+= sdsl::write_member(prev_alphabet, out, child, "p_alpha_size");
        written_bytes+= sdsl::write_member(n_phrases, out, child, "n_phrases");
        written_bytes+= sdsl::write_member(t_size, out, child, "t_size");
        written_bytes+= sdsl::write_member(max_sym_freq, out, child, "max_sym_freq");
        written_bytes+= sdsl::write_member(end_str_dummy, out, child, "end_str_dummy");
        written_bytes+= sdsl::write_member(bwt_dummy, out, child, "bwt_dummy");
        written_bytes+= sdsl::write_member(hocc_dummy, out, child, "hocc_dummy");
        written_bytes+= sdsl::write_member(metasym_dummy, out, child, "metasym_dummy");
        dict.serialize(out);
        phrases_has_hocc.serialize(out, child);
        return written_bytes;
    }

    void load(std::istream& in){
        sdsl::read_member(alphabet, in);
        sdsl::read_member(prev_alphabet, in);
        sdsl::read_member(n_phrases, in);
        sdsl::read_member(t_size, in);
        sdsl::read_member(max_sym_freq, in);
        sdsl::read_member(end_str_dummy, in);
        sdsl::read_member(bwt_dummy, in);
        sdsl::read_member(hocc_dummy, in);
        sdsl::read_member(metasym_dummy, in);
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