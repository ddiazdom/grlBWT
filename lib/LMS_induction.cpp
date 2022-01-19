//
// Created by Diaz, Diego on 4.1.2022.
//
#include "LMS_induction.h"

void suffix_induction(vector_t& sa, vector_t &dict, bv_t &d_lim, size_t alphabet, sdsl::cache_config& config) {

    vector_t buckets(alphabet+1, 0, sdsl::bits::hi(dict.size())+1);

    for(auto && sym : dict){
        buckets[sym]++;
    }

    size_t acc=0, freq;
    for(auto && bucket : buckets){
        freq = bucket;
        bucket = acc;
        acc +=freq;
    }
    buckets[buckets.size()-1] = acc;
    sdsl::store_to_file(buckets, sdsl::cache_file_name("iss_buckets", config));

    for(size_t i=0;i<dict.size();i++){
        if(d_lim[i]){
            sa[--buckets[dict[i]+1]] = i+1;
        }
    }
    sdsl::util::clear(buckets);
    sdsl::load_from_file(buckets, sdsl::cache_file_name("iss_buckets", config));

    induce_L_type(sa, dict, d_lim, buckets);

    sdsl::util::clear(buckets);
    sdsl::load_from_file(buckets, sdsl::cache_file_name("iss_buckets", config));
    for(size_t bck=0;bck<buckets.size()-1;bck++){
        buckets[bck] = buckets[bck+1]-1;
    }

    induce_S_type(sa, dict, d_lim, buckets);

    /*size_t next_av=0;
    for(size_t i=0;i<sa.size();i++){
        if(sa[i]!=0){
            sa[next_av++]=sa[i];
        }
    }
    sa.resize(next_av);*/

    //sdsl::store_to_file(sa, sa_file);

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

void induce_L_type(vector_t& sa, vector_t& dict, bv_t& d_lim, vector_t& buckets){
    std::cout<<"Inducing L-type"<<std::endl;
    //size_t sym, l_sym, p_sym, p_l_sym, sa_val;
    //vector_t tmp_sa = sa;
    //vector_t tmp_buckets = buckets;
    /*size_t buff_pos=0;
    size_t buffer[2048]={0};

    size_t k=0;
    while(sa[k]==0) k++;
    buffer[buff_pos++] = sa[k]-1;
    p_sym = dict[sa[k]-1];
    p_l_sym = dict[sa[k]-2];

    for(size_t i=k+1;i<sa.size();i++) {
        sa_val = sa[i];
        if(sa_val<=1){//sa[i]=0
            if(buff_pos!=0 && i==buckets[p_l_sym]) {//we reach a position that has not been induced yet
                for(size_t ind_pos = buckets[p_l_sym], j=0; j<buff_pos;j++,ind_pos++){
                    sa[ind_pos] = buffer[j];
                }
                buckets[p_l_sym]+=buff_pos;
                buff_pos=0;
                sa_val = sa[i];
            }else{
                continue;
            }
        }

        if(sa_val>1 && d_lim[sa_val-2]){
            continue;
        }

        sym = dict[sa_val-1];
        l_sym = dict[sa_val-2];
        if(l_sym!=p_l_sym || sym!=p_sym || buff_pos==2048){
            if(p_l_sym>=p_sym && buff_pos!=0){
                for(size_t ind_pos = buckets[p_l_sym], j=0; j<buff_pos;j++,ind_pos++){
                    sa[ind_pos] = buffer[j];
                }
                buckets[p_l_sym]+=buff_pos;
            }
            buff_pos=0;
        }
        buffer[buff_pos++] = sa_val-1;
        p_sym = sym;
        p_l_sym = l_sym;
    }

    if(p_l_sym>=p_sym){
        for(size_t ind_pos = buckets[p_l_sym], j=0; j<buff_pos;j++,ind_pos++){
            sa[ind_pos] = buffer[j];
        }
        buckets[p_l_sym]+=buff_pos;
    }*/

    size_t pos, l_sym, sym, ind_pos;
    for(size_t i=0;i<sa.size();i++){
        pos = sa[i];
        if(pos<=1  || d_lim[pos-2]) continue;
        sym = dict[pos-1];
        l_sym = dict[pos-2];
        if(l_sym>=sym){
            sa[buckets[l_sym]++] = pos-1;
        }
    }
}

void induce_S_type(vector_t& sa, vector_t& dict, bv_t& d_lim, vector_t& buckets){
    std::cout<<"Inducing S-type"<<std::endl;
    size_t pos, sym, l_sym;
    for(size_t i=sa.size();i-->0;){
        pos = sa[i];
        if(pos<=1  || d_lim[pos-2]) continue;

        sym = dict[pos-1];
        l_sym = dict[pos-2];

        if(l_sym < sym ||
           (l_sym==sym && i >= buckets[sym])){
           sa[buckets[l_sym]--] = pos-1;
        }
    }
}
