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
void resort_sa(sdsl::int_vector<>& g_sa, size_t start, size_t end,
               sdsl::int_vector_buffer<>& new_g_sa,
               sdsl::int_vector<>& buckets,
               gramm_extra_feat<grammar_t>& gram_ef){
    for(size_t i=start;i<end;i++){
        if(!gram_ef.is_set(g_sa[i])){
            new_g_sa.push_back(g_sa[i]);
        }else{
            size_t sym= gram_ef.parent_nt(g_sa[i]);
            resort_sa(g_sa, buckets[sym], buckets[sym+1], new_g_sa , buckets, gram_ef);
        }
    }
}

template<class grammar_t>
void gram2bwt(grammar_t& gram){

    gramm_extra_feat gram_ef(gram);
    sdsl::int_vector<> buckets(gram.symbols()+1, 0, sdsl::bits::hi(gram.gram_size())+1);

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
                /*while(curr_pos<gram.gram_size()){
                    sym = gram.rules[curr_pos];
                    if(gram.is_rl(sym)) sym =  gram.rules[gram.nter_ptr[sym]];
                    buckets[sym]++;
                    curr_pos++;
                }*/
                break;
            }else{
                curr_pos = rl_zone.second + 1;
                last_pos = gram_ef.rb(curr_pos);
                limit = c_start;
            }
        }
    }

    //create the buckets
    size_t acc=0, tmp, max_freq=0;
    for(auto && bucket : buckets){
        tmp = bucket;
        bucket = acc;
        acc +=tmp;
        if(tmp>max_freq) max_freq = tmp;
    }

    //create the empty SA
    sdsl::int_vector<> gram_sa(acc, 0, sdsl::bits::hi(gram.gram_size())+1);
    sdsl::int_vector<> freqs(gram.symbols()+1, 0, sdsl::bits::hi(max_freq)+1);

    //TODO testing
    /*size_t curr=0;
    while(curr<gram.rules.size()){
        size_t next = gram_ef.rb(curr);
        std::cout<<curr<<" "<<next<<" "<<gram.rules.size()<<std::endl;
        assert(gram_ef.is_set(curr++));
        while(curr<=next){
            assert(!gram_ef.is_set(curr++));
        }
    }*/
    //

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
                        gram_sa[buckets[sym] + freqs[sym]++] = curr_pos;
                    }
                }
                curr_pos++;
            }

            if(curr_pos==c_start){
                /*while(curr_pos<gram.gram_size()){
                    sym = gram.rules[curr_pos];
                    if(gram.is_rl(sym)) sym =  gram.rules[gram.nter_ptr[sym]];
                    gram_sa[buckets[sym]++] = curr_pos;
                    curr_pos++;
                }*/
                break;
            }else{
                curr_pos = rl_zone.second + 1;
                last_pos = gram_ef.rb(curr_pos);
                limit = c_start;
            }
        }
    }

    auto it_a = gram.range(0,0);
    auto it_b = gram.range(0,0);
    auto it_sa_start = gram_sa.begin();
    sdsl::int_vector<>::iterator it_sa_end;
    size_t head;
    do{
        head = gram.rules[*it_sa_start];
        it_sa_end = it_sa_start+1;
        while(it_sa_end!=gram_sa.end() && gram.rules[*it_sa_end]==head){
            ++it_sa_end;
        }
        std::sort(it_sa_start, it_sa_end, [&](size_t a, size_t b){
            it_a.update_it(a, gram_ef.rb(a));
            it_b.update_it(b, gram_ef.rb(b));
            return it_a < it_b;
        });
        it_sa_start = it_sa_end;
    }while(it_sa_start!=gram_sa.end());

    /*size_t k = 0;
    size_t c_start = gram.nter_ptr[gram.symbols()-1];
    for(auto const& pos : gram_sa){
        assert(!gram.in_rl_zone(pos));
        if(pos<c_start){
                //std::cout<<pos<<":"<<c_start<<" "<<(pos>=c_start)<<" -> ";
                //std::cout<<k<<"   "<<pos<<":"<<gram_ef.rb(pos)<<" -> ";
                //for(size_t j=pos;j<=gram_ef.rb(pos);j++){
                //    std::cout<<gram.rules[j]<<","<<gram.parsing_level(gram.rules[j])<<","<<gram.is_sp(gram.rules[j])<<" ";
                //}
                //std::cout<<" "<<std::endl;
                std::cout<<gram.parsing_level(gram.rules[pos])<<" "<<pos<<"-"<<gram_ef.rb(pos)<<" ";
                print_suffix(pos, gram_ef.rb(pos), gram);
        }
        k++;
    }*/
}

template void gram2bwt<grammar<sdsl::int_vector<>>>(grammar<sdsl::int_vector<>>& gram);
template void gram2bwt<grammar<huff_vector<>>>(grammar<huff_vector<>>& gram);
