//
// Created by Diaz, Diego on 10.11.2021.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_HPP
#define LPG_COMPRESSOR_GRAMMAR_HPP

#include <iostream>
#include <unordered_map>
#include <vector>

struct gram_info_t{
    uint8_t                            sigma{}; // terminal's alphabet
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the m_rules' right-hand sides
    size_t                             max_tsym{}; //highest terminal symbol
    size_t                             last_dict_size{}; //size of the last dictionary
    std::unordered_map<size_t,uint8_t> sym_map; //map terminal symbols to their original byte symbols
    std::string                        rules_file; // m_rules are concatenated in this array
    std::string                        rules_lim_file; // bit vector marking the last symbol of every right-hand
    std::vector<size_t>                rules_breaks; // number of m_rules generated in every LMS parsing round
    size_t                             n_p_rounds{}; // number of parsing rounds

    gram_info_t() = default;
    gram_info_t(std::string& r_file_,
                std::string& r_lim_file_): rules_file(r_file_),
                                           rules_lim_file(r_lim_file_) {}

    void save_to_file(std::string& output_file);
    void load_from_file(std::string &g_file);

    [[nodiscard]] inline long long int parsing_level(size_t symbol) const{
        if(symbol < sigma) return 0;
        if(symbol>=rules_breaks[n_p_rounds]){
            return -1;
        }

        for(long long int i=0;i<int(n_p_rounds);i++){
            if(rules_breaks[i]<=symbol && symbol<rules_breaks[i+1]) return i+1;
        }
        return -1;
    }
};
#endif //LPG_COMPRESSOR_GRAMMAR_HPP
