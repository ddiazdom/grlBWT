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
        if(decomp.back()=='\n') decomp.back() = '$';
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
    sdsl::int_vector<> buckets(gram.symbols()+1, 0, sdsl::bits::hi(gram.gram_size())+1);
    using g_iterator = typename grammar_t::iterator;

    //count the frequency of every symbol
    {
        size_t curr_pos = gram.ter(), last_pos = gram_ef.rb(curr_pos), sym;
        std::pair<size_t, size_t> rl_zone = gram.rl_zone();
        size_t limit = rl_zone.first, c_start = gram.nter_ptr[gram.symbols()-1];
        while(true){
            while (curr_pos < limit) {
                if (curr_pos == last_pos) {
                    last_pos = gram_ef.rb(last_pos + 1);
                } else {
                    sym = gram.rules[curr_pos];
                    if(!gram.is_sp(sym)){
                        if(gram.is_rl(sym)) sym =  gram.rules[gram.nter_ptr[sym]];
                        buckets[sym]++;
                    }
                }
                curr_pos++;
            }

            if(curr_pos==c_start){
                while(curr_pos<gram.gram_size()){
                    sym = gram.rules[curr_pos];
                    if(gram.is_rl(sym)) sym =  gram.rules[gram.nter_ptr[sym]];
                    buckets[sym]++;
                    curr_pos++;
                }
                break;
            }
            curr_pos = rl_zone.second + 1;
            last_pos = gram_ef.rb(curr_pos);
            limit = c_start;
        }
    }

    //create the buckets
    size_t acc=0, tmp;
    for(auto && bucket : buckets){
        tmp = bucket;
        bucket = acc;
        acc +=tmp;
    }

    //create the empty SA
    sdsl::int_vector<> gram_sa(acc, 0, sdsl::bits::hi(gram.gram_size())+1);

    //put the symbols symbol occurrences in their positions
    {
        size_t curr_pos = gram.ter(), last_pos = gram_ef.rb(curr_pos), sym;
        std::pair<size_t, size_t> rl_zone = gram.rl_zone();
        size_t limit = rl_zone.first, c_start = gram.nter_ptr[gram.symbols()-1];
        while(true){

            while (curr_pos < limit) {
                if (curr_pos == last_pos) {
                    last_pos = gram_ef.rb(last_pos + 1);
                } else {
                    sym = gram.rules[curr_pos];
                    if(!gram.is_sp(sym)){
                        if(gram.is_rl(sym)) sym =  gram.rules[gram.nter_ptr[sym]];
                        gram_sa[buckets[sym]++] = curr_pos;
                    }
                }
                curr_pos++;
            }

            if(curr_pos==c_start){
                while(curr_pos<gram.gram_size()){
                    sym = gram.rules[curr_pos];
                    if(gram.is_rl(sym)) sym =  gram.rules[gram.nter_ptr[sym]];
                    gram_sa[buckets[sym]++] = curr_pos;
                    curr_pos++;
                }
                break;
            }
            curr_pos = rl_zone.second + 1;
            last_pos = gram_ef.rb(curr_pos);
            limit = c_start;
        }
    }

    //TODO testing
    /*for(size_t i=0;i<gram.symbols()-1;i++){
        size_t start = gram.nter_ptr[i];
        size_t end = gram.nter_ptr[i+1]-1;
        for(size_t j=start;j<=end;j++){
            assert(gram_ef.rb(j)==end);
        }
    }
    size_t c_start = gram.nter_ptr[gram.symbols()-1];
    for(size_t i=0;i<gram.strings();i++){
        size_t start = c_start + gram.seq_ptr[i];
        size_t end = c_start + gram.seq_ptr[i+1]-1;
        for(size_t j=start;j<=end;j++){
            assert(gram_ef.rb(j)==end);
        }
    }*/
    //

    size_t k = 0;
    for(auto const& pos : gram_sa){
        assert(!gram.in_rl_zone(pos));
        if(pos>gram.nter_ptr[gram.symbols()-1]){
            std::cout<<k<<"   "<<pos<<":"<<gram_ef.rb(pos)<<" -> ";
            for(size_t j=pos;j<=gram_ef.rb(pos);j++){
                std::cout<<gram.rules[j]<<","<<gram.parsing_level(gram.rules[j])<<","<<gram.is_sp(gram.rules[j])<<" ";
            }
            std::cout<<" "<<std::endl;
            print_suffix(pos, gram_ef.rb(pos), gram);
        }
        k++;
    }

    auto end_it = gram.end();
    auto comp_positions = [&](g_iterator a, g_iterator b) -> bool {

        //move to the right as long as they have the same sequence
        while(a!=end_it && b!= end_it && a==b){
            a.skip_subtree();
            b.skip_subtree();
            ++a;
            ++b;
        }

        if(a==end_it || b==end_it){//one is a proper prefix of the other
            return a==end_it && b!=end_it;
        }else{
            while(a!=end_it && b!=end_it){
                while(!gram.is_lc(*a)) ++a;
                while(!gram.is_lc(*b)) ++b;
                size_t lvl_a = gram.parsing_level(*a);
                size_t lvl_b = gram.parsing_level(*b);
                while(lvl_b!=lvl_a){
                    if(lvl_b<lvl_a){
                        ++a;
                        lvl_a = gram.parsing_level(*a);
                    }else {
                        ++b;
                        lvl_b = gram.parsing_level(*b);
                    }
                }
                if(*a!=*b) return *a<*b;
                ++a;
                ++b;
            }
            return a==end_it && b!=end_it;
        }
    };

    auto a = gram.range(80412, 80419);
    auto b = gram.range(80540, 80543);

    comp_positions(a, b);
}

template void gram2bwt<grammar<sdsl::int_vector<>>>(grammar<sdsl::int_vector<>>& gram);
template void gram2bwt<grammar<huff_vector<>>>(grammar<huff_vector<>>& gram);
