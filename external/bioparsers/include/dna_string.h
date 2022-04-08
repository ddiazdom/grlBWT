//
// Created by Diaz, Diego on 13.10.2021.
//

#ifndef BIOPARSER_DNA_STRING_H
#define BIOPARSER_DNA_STRING_H

#include <string>
#include <cassert>
#include <vector>
#include <iostream>
#include <algorithm>

class dna_string : public std::basic_string<uint8_t> {
private:
    std::vector<uint16_t> run_lengths;
    size_t                e_size=0;

public:

    const static uint8_t dna_offset = 84;
    static uint8_t rev_comp[170];

    struct iterator{
        using value_type        = uint8_t;
    private:
        const dna_string& seq;
        size_t s_pos{};
        size_t r_pos{};
        size_t run_len{};
        size_t idx{};
        value_type sym{};

    public:
        explicit iterator(const dna_string& _seq, size_t idx_): seq(_seq),
                                                                s_pos(0),
                                                                r_pos(0),
                                                                run_len(0),
                                                                idx(idx_){
            if(idx<seq.eff_size()){
                sym = seq[s_pos++];
                if(sym>dna_offset){
                    run_len = seq.homopolymer_length(r_pos++)-1;
                    sym-=dna_offset;
                }
            }
        };

        iterator& operator++(){
            if(idx<seq.eff_size()){
                if(run_len==0){
                    sym = seq[s_pos++];
                    if(sym>dna_offset){
                        run_len =  seq.homopolymer_length(r_pos++);
                        sym -= dna_offset;
                    }else{
                        run_len=1;
                    }
                }
                run_len--;
            }
            idx = std::min(idx+1, seq.eff_size());
            return *this;
        }

        inline char operator*() const {
            return (char)sym;
        }

        inline bool operator==(const iterator& other) const {
            return other.idx == idx;
        };

        inline bool operator!=(const iterator& other) const {
            return other.idx != idx;
        };
    };

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

    dna_string& reverse(){
        for(auto i=begin(),j=end()-1;i<j;++i,--j) {
            std::swap(*i, *j);
        }
        std::reverse(run_lengths.begin(), run_lengths.end());
        return *this;
    }

    dna_string& complement(){
        for(unsigned char & sym : *this) {
            sym = rev_comp[(uint8_t)sym];
        }
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
        assert(sym>0 && sym<=dna_offset);
        push_back(sym);
        e_size++;
    }

    inline void push_homopolymer(uint8_t sym, size_t length){
        assert(sym>dna_offset && length>1);
        push_back(sym);
        run_lengths.push_back(length);
        e_size+=length;
    }

    inline void expand_homopolymers(){
        dna_string other;
        size_t i=0;
        for(unsigned char & sym : *this) {
            if(sym>dna_offset){
                size_t len = run_lengths[i++];
                uint8_t real_sym = sym-dna_offset;
                for(size_t j=0;j<len;j++){
                    other.push_sym(real_sym);
                }
            }else{
                other.push_sym(sym);
            }
        }
        this->swap(other);
    }

    inline void contract_homopolymers(){

    }

    inline dna_string& operator+(dna_string& other){
        insert(end(), other.begin(), other.end());
        run_lengths.insert(run_lengths.end(), run_lengths.begin(), run_lengths.end());
        e_size +=other.e_size;
        return *this;
    }

    [[nodiscard]] inline size_t homopolymer_length(size_t run_idx) const{
        assert(run_idx<run_lengths.size());
        return run_lengths[run_idx];
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

#endif //BIOPARSER_DNA_STRING_H
