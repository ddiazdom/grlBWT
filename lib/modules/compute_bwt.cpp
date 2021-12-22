//
// Created by Diaz, Diego on 15.12.2021.
//

#include "modules/compute_bwt.h"
#include "sdsl/int_vector.hpp"

template<class grammar_t>
void print_suffix(size_t start, size_t end, grammar_t& gram){
    std::string decomp;
    for(size_t i=start;i<=end;i++){
        gram.decomp_nt(gram.rules[i], decomp);
        if(i<end){
            decomp.push_back(' ');
            decomp.push_back('|');
            decomp.push_back(' ');
        }
    }
    std::cout<<decomp<<std::endl;
}

template<class grammar_t>
void gram2bwt(grammar_t& gram){

    gramm_extra_feat gram_ef(gram);

    sdsl::int_vector<> gram_sa(gram.gram_size(), 0, sdsl::bits::hi(gram.gram_size())+1);
    sdsl::int_vector<> buckets(gram.symbols()+1, 0, sdsl::bits::hi(gram.gram_size())+1);

    //count the frequency of every symbol
    for(size_t i=0;i<gram.rules.size(); i++){
        buckets[gram.rules[i]]++;
    }

    //create the buckets
    size_t acc=0, tmp;
    for(auto && bucket : buckets){
        tmp = bucket;
        bucket = acc;
        acc +=tmp;
    }

    //put the symbols symbol occurrences in their positions
    for(size_t i=0;i<gram.gram_size();i++){
        gram_sa[buckets[gram.rules[i]]++] = i;
    }

    //TODO testing
    for(size_t i=0;i<gram.symbols()-1;i++){
        size_t start = gram.nter_ptr[i];
        size_t end = gram.nter_ptr[i+1]-1;
        for(size_t j=start;j<=end;j++){
            assert(gram_ef.rb(j)==end);
        }
    }
    size_t c_start = gram.nter_ptr[gram.symbols()-1];
    for(size_t i=0;i<gram.strings();i++){
        size_t start = c_start + (i==0? 0 : gram.seq_ptr[i-1]+1);
        size_t end = c_start+gram.seq_ptr[i];
        for(size_t j=start;j<=end;j++){
            std::cout<<start<<":"<<end<<" -> "<<j<<" "<<gram_ef.rb(j)<<std::endl;
            assert(gram_ef.rb(j)==end);
        }
    }
    //

    for(auto const& pos : gram_sa){
        if(gram.in_rl_zone(pos)){
            std::cout<<pos<<" "<<gram_ef.rb(pos)<<std::endl;
            std::cout<<"blablabla"<<std::endl;
        }else{
            //std::cout<<pos<<" "<<gram_ef.rb(pos)<<std::endl;
            //print_suffix(pos, gram_ef.rb(pos), gram);
        }
    }

    auto descend = [&](size_t& sym, size_t& freq, std::stack<size_t>& stack){
        if(sym<gram.ter()){
            sym = stack.top();
            stack.pop();
            freq = 1;
        } else{
            size_t start = gram.nter_ptr[sym];
            if(gram.is_rl(sym)){
                freq = gram.rules[start + 1];
            } else{
                freq = 1;
                size_t end = gram.nter_ptr[sym+1]-1;
                for(size_t i=end; i>start;i--){
                    stack.push(gram.rules[i]);
                }
            }
            sym = gram.rules[start];
        }
    };

    auto comp_positions = [&](size_t idx_a, size_t idx_b) -> bool {

        size_t sym_a, sym_b, bound_a=gram_ef.rb(idx_a), bound_b=gram_ef.rb(idx_b);

        //move to the right as long as they have the same sequence
        while(idx_a<=bound_a && idx_b<=bound_b &&
              gram.rules[idx_a] == gram.rules[idx_b]){
            idx_a++;
            idx_b++;
        }

        if(idx_a==bound_a || idx_b==bound_b){//one is a proper prefix of the other
            return (bound_a-idx_a+1)<(bound_b-idx_b+1);
        }else{
            size_t freq_a, freq_b, lvl_a, lvl_b;

            std::stack<size_t> stack_a;
            std::stack<size_t> stack_b;

            sym_a = gram.rules[idx_a];
            sym_b = gram.rules[idx_b];

            do{
                while(!gram.is_lc(sym_a)){
                    descend(sym_a, freq_a, stack_a);
                }

                while(!gram.is_lc(sym_b)){
                    descend(sym_b, freq_b, stack_b);
                }

                lvl_a = gram.parsing_level(sym_a);
                lvl_b = gram.parsing_level(sym_b);

                while(lvl_b!=lvl_a){
                    if(lvl_b<lvl_a){
                        descend(sym_a, freq_a, stack_a);
                        lvl_a = gram.parsing_level(sym_a);
                    }else {
                        descend(sym_b, freq_b, stack_b);
                        lvl_b = gram.parsing_level(sym_b);
                    }
                }

                if(sym_a!=sym_b){
                    return sym_a<sym_b;
                }else{
                    sym_a = stack_a.top();
                    sym_b = stack_b.top();
                }
            }while(stack_a.empty() && stack_b.empty());
            return stack_a.empty() && !stack_b.empty();
        }
    };

    //comp_positions(10, 2);
}

template void gram2bwt<grammar<sdsl::int_vector<>>>(grammar<sdsl::int_vector<>>& gram);
template void gram2bwt<grammar<huff_vector<>>>(grammar<huff_vector<>>& gram);
