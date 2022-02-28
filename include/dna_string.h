//
// Created by Diaz, Diego on 13.10.2021.
//

#ifndef LHTIGS_DNA_STRING_H
#define LHTIGS_DNA_STRING_H
#include <string>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>
#include "utils.h"

class dna_string : public std::basic_string<uint8_t> {
private:
    std::vector<uint16_t> run_lengths;
    size_t                e_size=0;

public:
    struct iterator{
        using value_type        = char;
    private:
        const dna_string& seq;
        size_t s_pos{};
        size_t r_pos{};
        size_t run_len{};
        size_t v_syms{};
        value_type sym{};

    public:
        explicit iterator(const dna_string& _seq, size_t v_syms_): seq(_seq),
                                                             s_pos(0),
                                                             r_pos(0),
                                                             run_len(0),
                                                             v_syms(v_syms_),
                                                             sym(0){
            if(v_syms<seq.eff_size()){
                uint8_t tmp = seq[s_pos++];
                run_len = (tmp & 1U) ? seq.get_run(r_pos++)-1 : 0 ;
                sym = char(tmp>>1U);
            }
        };

        iterator& operator++(){
            if(v_syms<seq.eff_size()){
                if(run_len==0){
                    uint8_t tmp = seq[s_pos++];
                    run_len = (tmp & 1U) ? seq.get_run(r_pos++) : 1;
                    sym = char(tmp>>1U);
                }
                run_len--;
            }
            v_syms = std::min(v_syms+1, seq.eff_size());
            return *this;
        }

        inline char operator*() const {
            return sym;
        }

        inline bool operator==(const iterator& a) const {
            return a.v_syms == v_syms;
        };

        inline bool operator!=(const iterator& a) const {
            return a.v_syms != v_syms;
        };
    };

    static uint8_t rev_comp[170];

    using std::basic_string<uint8_t>::basic_string;

    dna_string& reverse_complement(){
        for(auto i=begin(),j=end()-1;i<j;++i,--j) {
            std::swap(*i, *j);
            *i = rev_comp[(uint8_t)*i];
            *j = rev_comp[(uint8_t)*j];
            assert(*j!=0 && *i!=0);
        }

        if(length() & 1){
            (*this)[length()/2] = rev_comp[(uint8_t)(*this)[length()/2]];
        }
        std::reverse(run_lengths.begin(), run_lengths.end());
        return *this;
    }

    dna_string reverse_complement(size_t n){
        dna_string res;
        assert(n<=length());
        size_t tot_symbols = length();
        for(size_t i=1;i<=n;i++){
            res.push_back(rev_comp[(uint8_t)(*this)[tot_symbols-i]]);
        }
        return res;
    }

    inline void push_sym(uint8_t sym){
        assert((sym>0 && sym<=84) || !(sym & 1U));
        push_back(sym);
        e_size++;
    }

    inline void push_sym_run(uint8_t sym, size_t length){
        assert(sym>84 || (sym & 1U));
        assert(length>1);
        push_back(sym);
        run_lengths.push_back(length);
        e_size+=length;
    }

    [[nodiscard]] inline size_t get_run(size_t idx) const{
        assert(idx<run_lengths.size());
        return run_lengths[idx];
    }

    [[nodiscard]] inline size_t eff_size() const{
        return e_size;
    }

    [[nodiscard]] inline iterator ibegin() const {
        return iterator(*this, 0);
    }

    [[nodiscard]] inline iterator iend() const {
        return iterator(*this, eff_size());
    }
};

std::ostream& operator<<(std::ostream& os, const dna_string& seq);

//get the DNA reverse complements
void get_rev_comp(const std::string& plain_file, str_collection& str_coll);

#endif //LHTIGS_DNA_STRING_H
