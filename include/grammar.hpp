//
// Created by Diaz, Diego on 10.11.2021.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_HPP
#define LPG_COMPRESSOR_GRAMMAR_HPP
#include <iostream>
#include <sdsl/int_vector.hpp>

struct gram_info_t{
    uint8_t                            sigma{}; // terminal's alphabet
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the rules' right-hand sides
    size_t                             max_tsym{}; //highest terminal symbol
    std::unordered_map<size_t,uint8_t> sym_map; //map terminal symbols to their original byte symbols
    std::string                        rules_file; // rules are concatenated in this array
    std::string                        rules_lim_file; // bit vector marking the last symbol of every right-hand
    std::vector<size_t>                rules_per_level; // number of rules generated in every LMS parsing round
    //std::string                        is_rl_file; //bit vector [0..r-1] indicating which rules are run-length compressed
    std::string                        lvl_breaks_file; //file with the LMS breaks per level
    size_t                             n_lcp_rules; // number of locally consistent parsing rules
    size_t                             n_rl_rules; // number of run-length rules
    size_t                             n_sp_rules; // number of suffix pair rules

    gram_info_t() = default;
    gram_info_t(std::string& r_file_,
                std::string& r_lim_file_,
                //std::string& rl_rules_file_,
                std::string& lvl_breaks_file_): rules_file(r_file_),
                                                rules_lim_file(r_lim_file_),
                                                //is_rl_file(rl_rules_file_),
                                                lvl_breaks_file(lvl_breaks_file_){}

    void save_to_file(std::string& output_file);
    void load_from_file(std::string &g_file);

    [[nodiscard]] bool is_terminal(const size_t& id) const {
        return sym_map.find(id) != sym_map.end();
    }
};

class grammar {
private:
    size_t orig_size{}; //original size of the text
    size_t n_seqs{}; //original size of the text
    size_t grammar_size;//size of the grammar (i.e., sum of the lengths of the right-hand sides of the rules)
    size_t sigma; //number of terminals
    size_t n_nter; //number of nonterminals
    size_t comp_string_size; //size of the compressed string
    sdsl::int_vector<8> symbols_map; //map a terminal to its original byte symbol
    sdsl::int_vector<> rules; //list of rules
    sdsl::int_vector<> nter_ptr; //pointers of the rules in the rules' array
    //sdsl::bit_vector is_rl; //bit vector indicating if the rule is run-length compressed
    sdsl::int_vector<> seq_pointers; //pointers to the boundaries of the strings in the text

    void mark_str_boundaries();
public:
    typedef size_t size_type;
    explicit grammar(gram_info_t& gram_info, size_t _orig_size, size_t _n_seqs): orig_size(_orig_size),
                                                                                 n_seqs(_n_seqs),
                                                                                 grammar_size(gram_info.g),
                                                                                 sigma(gram_info.sigma),
                                                                                 n_nter(gram_info.r),
                                                                                 comp_string_size(gram_info.c){

        //TODO check if the alphabet is compressed
        symbols_map.resize(gram_info.sigma);
        for(auto const pair: gram_info.sym_map){
            symbols_map[pair.first] = pair.second;
        }

        sdsl::int_vector_buffer<> rules_buff(gram_info.rules_file, std::ios::in);
        rules.width(sdsl::bits::hi(sigma + n_nter)+1);
        rules.resize(rules_buff.size());
        size_t i=0;
        for(auto const& sym : rules_buff) rules[i++] = sym;
        rules_buff.close();

        sdsl::int_vector_buffer<1> rules_lim_buff(gram_info.rules_lim_file, std::ios::in);
        nter_ptr.width(sdsl::bits::hi(grammar_size)+1);
        nter_ptr.resize(n_nter);
        nter_ptr[0] = 0;

        i=1;
        for(size_t j=1;j<rules_lim_buff.size();j++){
            if(rules_lim_buff[j-1]){
                nter_ptr[i++] = j;
            }
        }
        rules_lim_buff.close();
        //sdsl::load_from_file(is_rl, gram_info.is_rl_file);
        mark_str_boundaries();
    }

    std::string decompress_seq(size_t idx);
    size_type serialize(std::ostream& out, sdsl::structure_tree_node * v=nullptr, std::string name="") const;
    void load(std::ifstream& in);
};

#endif //LPG_COMPRESSOR_GRAMMAR_HPP
