//
// Created by Diaz, Diego on 15.12.2021.
//

#include "modules/compute_bwt.h"
#include "sdsl/int_vector.hpp"

template<class grammar_t>
void gram2bwt(grammar_t& gram){

    rank_support rs(gram.nter_ptr);

    size_t j=0, succ=gram.nter_ptr[j];
    for(size_t i=0;i<gram.nter_ptr[gram.nter_ptr.size()-1];i++){

        if(i==4097){
            std::cout<<"holaa"<<std::endl;
        }
        std::cout<<i<<" -> "<<succ<<" "<<rs.successor(i)<<std::endl;
        assert(succ==rs.successor(i));
        if(i==succ){
            succ = gram.nter_ptr[++j];
        }
    }


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

        size_t sym_a, sym_b, bound_a=3, bound_b=3;

        //move to the right as long as they have the same sequence
        while(idx_a<=bound_a && idx_b<=bound_b &&
              gram.rules[idx_a] == gram.rules[idx_b]){
            idx_a++;
            idx_b++;
        }

        if(idx_a==bound_a || idx_b==bound_b){//one is a proper prefix of the other
            return (bound_a-idx_a)<(bound_b-idx_b);
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

    comp_positions(10, 2);
}

template void gram2bwt<grammar<sdsl::int_vector<>>>(grammar<sdsl::int_vector<>>& gram);
template void gram2bwt<grammar<huff_vector<>>>(grammar<huff_vector<>>& gram);
