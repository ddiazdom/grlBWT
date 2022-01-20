//
// Created by Diaz, Diego on 4.1.2022.
//
#include "LMS_induction.h"

void suffix_induction(vector_t& sa, const dictionary& dict, sdsl::cache_config& config) {

    vector_t buckets(dict.alphabet+1, 0, sdsl::bits::hi(dict.dict.size())+1);

    for(auto && sym : dict.dict){
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

    vector_t freqs(buckets.size(), 0, sdsl::bits::hi(dict.n_phrases)+1);
    size_t sym;
    for(size_t i=0;i<dict.dict.size();i++) {
        if(dict.d_lim[i]){
            sym = dict.dict[i];
            sa[buckets[sym+1]-1-freqs[sym]++] = ((i+1)<<1UL) | 1UL;
        }
    }

    size_t pos;
    for(size_t i=0;i<dict.alphabet;i++){
        if(freqs[i]>0){
            pos = buckets[i+1]-freqs[i];
            sa[pos] = sa[pos] & ~1UL;
        }
    }

    sdsl::util::clear(freqs);
    induce_L_type(sa, dict, buckets);

    sdsl::util::clear(buckets);
    sdsl::load_from_file(buckets, sdsl::cache_file_name("iss_buckets", config));
    for(size_t bck=0;bck<buckets.size()-1;bck++){
        buckets[bck] = buckets[bck+1]-1;
    }

    induce_S_type(sa, dict, buckets);

    /*if(sa.size()>120991){
        size_t j = 0;
        for (size_t i=0;i<sa.size();i++) {
            if (sa[i] == 0) {
                //std::cout << "*" << std::endl;
            } else {
                bool lcs = sa[i] & 1UL;
                pos = (sa[i] >> 1U) - 1;
                if(lcs ||
                   (i<sa.size()-1 && sa[i+1] & 1UL)){
                    std::cout << j << " " << pos << " " << lcs << "\t";

                    if (pos == 0 || dict.d_lim[pos - 1]) {
                        std::cout << "$ | ";
                    } else {
                        std::cout << dict.dict[pos - 1] << " | ";
                    }

                    while (!dict.d_lim[pos]) {
                        std::cout << dict.dict[pos++] << " ";
                    }
                    std::cout << dict.dict[pos] << std::endl;
                }
            }
            j++;
        }
    }*/
}

void induce_L_type(vector_t &sa, const dictionary& dict, vector_t &buckets) {
    std::cout<<"Inducing L-type"<<std::endl;

    size_t pos, l_sym, bck, ind_pos, p_bck=dict.alphabet+1;
    vector_t ind_bck(buckets.size(), 0, sdsl::bits::hi(sa.size())+1);
    bool lcs, first_eq=false;
    size_t lb=1, next_lb = sa.size();

    for(size_t i=0;i<sa.size();i++) {

        pos = sa[i];
        lcs = pos & 1UL;
        pos>>=1UL;
        if(pos==0) continue;

        bck = dict.dict[pos-1];

        if(bck!=p_bck || i==next_lb || !lcs) {
            first_eq = false;
            lb = i+1;
            assert(!lcs);
        }

        if(pos==1 || dict.d_lim[pos-2]){
            p_bck = bck;
            continue;
        }

        l_sym = dict.dict[pos-2];
        if(l_sym>=bck){
            ind_pos = buckets[l_sym]++;
            if(!first_eq && l_sym==bck){
                first_eq = true;
                next_lb = ind_pos;
                lcs = false;
            }
            sa[ind_pos] = (pos-1)<<1UL | (ind_bck[l_sym]>=lb && lcs);
            ind_bck[l_sym] = i+1;
        }
        p_bck=bck;
    }
}

void induce_S_type(vector_t &sa, const dictionary &dict, vector_t &buckets) {
    std::cout<<"Inducing S-type"<<std::endl;
    size_t pos, bck, l_sym, ind_pos, p_bck;
    vector_t ind_bck(buckets.size(), 0, sdsl::bits::hi(sa.size())+1);
    bool lcs, first_eq, new_break, p_lcs=false;
    size_t rb=sa.size(), prev_rb=0;

    for(size_t i=sa.size();i-->0;){

        pos = sa[i];
        lcs = pos & 1UL;
        pos>>=1UL;
        if(pos==0){
            sa[i+1] = sa[i+1] & ~1UL;
            p_lcs = false;
            continue;
        }
        bck = dict.dict[pos-1];

        if(bck!=p_bck || i==prev_rb || !p_lcs) {
            first_eq = false;
            rb = i+1;
            assert(!p_lcs);
        }

        if(pos==1  || dict.d_lim[pos-2]){
            p_bck = bck;
            p_lcs = lcs;
            continue;
        }

        l_sym = dict.dict[pos-2];

        if(l_sym < bck ||
           (l_sym==bck && i >= buckets[bck])){

            ind_pos = buckets[l_sym]--;
            sa[ind_pos] = (pos-1)<<1UL | (lcs && ind_pos>0 && sa[ind_pos-1]==0);

            new_break = false;
            if(!first_eq && l_sym==bck){
                first_eq = true;
                prev_rb = ind_pos;
                new_break = true;
            }

            if(ind_bck[l_sym]==0 || ind_bck[l_sym]>rb || new_break){
                sa[ind_pos+1] = sa[ind_pos+1] & ~1UL;
            }
            ind_bck[l_sym] = i+1;
        }
        p_lcs = lcs;
        p_bck = bck;
    }
}
