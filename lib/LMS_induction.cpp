//
// Created by Diaz, Diego on 4.1.2022.
//
#include "LMS_induction.h"

void suffix_induction(std::string &sa_file, vector_t &dict, bv_t &d_lim, size_t alphabet) {

    vector_t buckets(alphabet+1, 0, sdsl::bits::hi(dict.size())+1);

    for(auto && sym : dict){
        buckets[sym]++;
    }

    size_t acc=0, tmp, max_freq=0;
    for(auto && bucket : buckets){
        tmp = bucket;
        bucket = acc;
        acc +=tmp;
        if(tmp>max_freq) max_freq = tmp;
    }
    buckets[buckets.size()-1] = acc;

    vector_t sa(dict.size(), 0, sdsl::bits::hi(dict.size())+1);
    vector_t freq(alphabet, 0, sdsl::bits::hi(max_freq)+1);

    size_t sym;
    for(size_t i=0;i<dict.size();i++){
        if(d_lim[i]){
            sym = dict[i];
            sa[buckets[sym+1]-(freq[sym]++)-1] = i+1;
        }
    }

    sdsl::util::set_to_value(freq, 0);
    induce_L_type(sa, dict, d_lim, buckets, freq);

    for(size_t bck=0;bck<buckets.size()-1;bck++){
        tmp = buckets[bck+1];
        buckets[bck] = buckets[bck]+freq[bck];
        freq[bck] = tmp - buckets[bck]-1;
    }

    induce_S_type(sa, dict, d_lim, buckets, freq);

    /*size_t next_av=0;
    for(size_t i=0;i<sa.size();i++){
        if(sa[i]!=0){
            sa[next_av++]=sa[i];
        }
    }
    sa.resize(next_av);*/

    sdsl::store_to_file(sa, sa_file);

    /*for(auto && i : sa){
        if(i==0)continue;

        size_t pos = i-1;
        std::cout<<pos<<"\t";
        while(!d_lim[pos]){
            std::cout<<dict[pos++]<<" ";
        }
        assert(d_lim[pos]);
        std::cout<<dict[pos]<<std::endl;
    }*/
}

void induce_L_type(vector_t& sa, vector_t& dict, bv_t& d_lim, vector_t& buckets, vector_t& freq){
    size_t sym, prev_sym;
    for(auto && pos : sa){
        if(pos<=1  || d_lim[pos-2]) continue;

        sym = dict[pos-1];
        prev_sym = dict[pos-2];

        if(prev_sym >= sym){
            sa[buckets[prev_sym]+freq[prev_sym]++] = pos-1;
        }
    }
}

void induce_S_type(vector_t& sa, vector_t& dict, bv_t& d_lim, vector_t& buckets, vector_t& freq){
    size_t pos, sym, prev_sym;
    for(size_t i=sa.size();i-->0;){
        pos = sa[i];
        if(pos<=1  || d_lim[pos-2]) continue;

        sym = dict[pos-1];
        prev_sym = dict[pos-2];

        if(prev_sym < sym ||
           (prev_sym==sym && i >= buckets[sym])){
           sa[buckets[prev_sym]+freq[prev_sym]--] = pos-1;
        }
    }
}
