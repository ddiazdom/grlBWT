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
std::string decomp(size_t start, size_t end, grammar_t& gram){
    std::string decomp;
    for(size_t i=start;i<=end;i++){
        gram.decomp_nt(gram.rules[i], decomp);
        if(decomp.back()=='\n') decomp.back() = '$';
    }
    return decomp;
}

template<class grammar_t>
void resort_sa(sdsl::int_vector<>& g_sa, size_t start, size_t end,
               sdsl::int_vector_buffer<>& new_g_sa,
               sdsl::int_vector<>& buckets,
               gramm_extra_feat<grammar_t>& gram_ef,
               size_t c_start){
    for(size_t i=start;i<=end;i++){
        size_t pos = g_sa[i];
        if(!gram_ef.is_set(pos) || pos>=c_start){
            new_g_sa.push_back(pos);
        }else{
            size_t sym = gram_ef.parent_nt(pos);
            resort_sa(g_sa, buckets[sym], buckets[sym+1]-1, new_g_sa , buckets, gram_ef, c_start);
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

    //TODO testing
    /*size_t curr=0;
    size_t curr_nt=0;
    while(curr<gram.rules.size()){
        size_t next = gram_ef.rb(curr);
        //std::cout<<curr<<" "<<next<<" "<<gram.rules.size()<<" "<<curr_nt<<" "<<gram_ef.parent_nt(curr)<<std::endl;
        assert(gram_ef.is_set(curr));
        assert(curr_nt==gram_ef.parent_nt(curr));
        curr++;
        while(curr<=next){
            //std::cout<<curr<<" "<<next<<" "<<gram.rules.size()<<" "<<curr_nt<<" "<<gram_ef.parent_nt(curr)<<" "<<gram.symbols()-1<<std::endl;
            assert(!gram_ef.is_set(curr));
            assert(curr_nt==gram_ef.parent_nt(curr));
            curr++;
        }
        if(curr_nt<(gram.symbols()-1)){
            curr_nt++;
        }
    }*/
    //

    //put the symbols symbol occurrences in their positions
    {
        sdsl::int_vector<> freqs(gram.symbols()+1, 0, sdsl::bits::hi(max_freq)+1);
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
                    gram_sa[buckets[sym] + freqs[sym]++] = curr_pos;
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

    /*std::sort(gram_sa.begin()+2584, gram_sa.begin()+2615, [&](size_t a, size_t b){

        if(a==b) return false;
        it_a.update_it(a, gram_ef.rb(a));
        it_b.update_it(b, gram_ef.rb(b));

        std::string tmp1 = decomp(a, gram_ef.rb(a), gram);
        std::string tmp2 = decomp(b, gram_ef.rb(b), gram);


        print_suffix(a, gram_ef.rb(a), gram);
        print_suffix(b, gram_ef.rb(b), gram);
        bool found=false;
        size_t i=0;
        while(i<std::min<size_t>(tmp1.size(), tmp2.size())){
            if(tmp1[i]!=tmp2[i]){
                found = true;
                break;
                assert(tmp1[i]<tmp2[i]==it_a<it_b);
            }
            i++;
        }
        if(!found){
            if(i==tmp1.size()){
                assert(it_a<it_b == false);
            }else{
                assert(it_a<it_b == true);
            }
        }

        return it_a < it_b;
    });*/

    /*for(size_t pos=2584; pos<2614;pos++){
        it_a.update_it(gram_sa[pos], gram_ef.rb(gram_sa[pos]));
        it_b.update_it(gram_sa[2614], gram_ef.rb(gram_sa[2614]));
        print_suffix(gram_sa[pos], gram_ef.rb(gram_sa[pos]), gram);
        print_suffix(gram_sa[2614], gram_ef.rb(gram_sa[2614]), gram);
        std::cout<<(it_a<it_b)<<std::endl;
    }*/

    size_t bck = 0;
    while(bck<gram.symbols()){
        auto start = (long)buckets[bck];
        auto end = (long)buckets[bck+1];
        if(start<end){
            std::sort(it_sa_start+start, it_sa_start+end, [&](size_t a, size_t b){
                if(a==b) return false;
                it_a.update_it(a, gram_ef.rb(a));
                it_b.update_it(b, gram_ef.rb(b));
                return it_a < it_b;
            });
        }
        bck++;
    }

    size_t c_start = gram.nter_ptr[gram.symbols()-1];
    sdsl::int_vector_buffer<> new_g_sa("resorted_g_sa",std::ios::out, 1024*1024, sdsl::bits::hi(gram.gram_size())+1);
    resort_sa(gram_sa, 0, buckets[gram.ter()]-1, new_g_sa, buckets, gram_ef, c_start);
    new_g_sa.close();
    sdsl::load_from_file(gram_sa, "resorted_g_sa");

    //size_t k = 0;
    size_t prev_pos;
    std::string dc, prev_dc;
    for(auto const& pos : gram_sa){
        assert(!gram.in_rl_zone(pos));

        dc = decomp(pos, gram_ef.rb(pos), gram);
        //std::cout<<gram_ef.parent_nt(pos)<<"\t"<<dc<<std::endl;
        //std::cout<<pos<<"-"<<gram_ef.rb(pos)<<" "<<gram.rules[pos]<<"\t";
        //print_suffix(pos, gram_ef.rb(pos), gram);
        if(!prev_dc.empty()){
            size_t y=0,x=0;
            while(y<dc.size() && x<prev_dc.size()){
                if(dc[x]!=prev_dc[y]){
                    //assert(dc[x]>prev_dc[y]);
                    if(!(dc[x]>prev_dc[y])){
                        std::string tmp_dc, tmp_dc2;
                        if(!gram_ef.is_set(prev_pos)){
                            gram.decomp_nt(gram.rules[prev_pos-1], tmp_dc);
                        }
                        std::cout<<tmp_dc<<"\t\t\t";
                        print_suffix(prev_pos, gram_ef.rb(prev_pos), gram);

                        if(!gram_ef.is_set(pos)){
                            gram.decomp_nt(gram.rules[pos-1], tmp_dc2);
                        }
                        std::cout<<tmp_dc2<<"\t\t\t";
                        print_suffix(pos, gram_ef.rb(pos), gram);
                        std::cout<<" "<<std::endl;
                    }
                    goto next;
                }
                x++;
                y++;
            }
            if(x==dc.size() || y==prev_dc.size()){
                if(!(y>=x)){

                    std::string tmp_dc, tmp_dc2;
                    if(!gram_ef.is_set(prev_pos)){
                        gram.decomp_nt(gram.rules[prev_pos-1], tmp_dc);
                    }
                    std::cout<<tmp_dc<<"\t\t\t";
                    print_suffix(prev_pos, gram_ef.rb(prev_pos), gram);

                    if(!gram_ef.is_set(pos)){
                        gram.decomp_nt(gram.rules[pos-1], tmp_dc2);
                    }
                    std::cout<<tmp_dc2<<"\t\t\t";
                    print_suffix(pos, gram_ef.rb(pos), gram);
                    std::cout<<" "<<std::endl;
                }
                //assert(y>=x);
            }
        }
        next:
        prev_dc = dc;
        prev_pos = pos;

        //std::cout<<pos<<":"<<c_start<<" "<<(pos>=c_start)<<" -> ";
        //std::cout<<k<<"   "<<pos<<":"<<gram_ef.rb(pos)<<" -> ";
        //for(size_t j=pos;j<=gram_ef.rb(pos);j++){
        //    std::cout<<gram.rules[j]<<","<<gram.parsing_level(gram.rules[j])<<","<<gram.is_sp(gram.rules[j])<<" ";
        //}
        //std::cout<<" "<<std::endl;
        //std::cout<<k<<" "<<gram_ef.is_set(pos)<<" "<<gram.parsing_level(gram.rules[pos])<<" "<<pos<<"-"<<gram_ef.rb(pos)<<" ";
        //std::cout<<k<<" "<<gram.rules[pos]<<" "<<pos<<"-"<<gram_ef.rb(pos)<<" "<<gram_ef.is_set(pos)<<" "<<gram.parsing_level(gram.rules[pos])<<"\t";
        //std::string tmp_dc;
        //if(!gram_ef.is_set(pos)){
        //    gram.decomp_nt(gram.rules[pos-1], tmp_dc);
        //}
        //std::cout<<tmp_dc<<"\t\t\t";
        //print_suffix(pos, gram_ef.rb(pos), gram);
        //k++;
    }
}

template void gram2bwt<grammar<sdsl::int_vector<>>>(grammar<sdsl::int_vector<>>& gram);
template void gram2bwt<grammar<huff_vector<>>>(grammar<huff_vector<>>& gram);
