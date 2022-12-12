//
// Created by Diaz, Diego on 4.1.2022.
//
#include "LMS_induction.h"

template <class value_type>
void suffix_induction(vector_t &sa, const dictionary &dict) {

    auto * buckets = (value_type *) calloc((dict.alphabet+1), sizeof(value_type));

    for(size_t i=0;i<dict.dict.size();i++){
        //TODO testing
        /*if(dict.dict[i]==196170){
            std::cout<<i<<" | "<<dict.dict[i-1]<<" "<<dict.dict[i]<<" "<<dict.dict[i+1]<<" "<<dict.dict[i+2]<<" "<<dict.dict[i+3]<<" "<<dict.dict[i+4]<<" | "<<dict.dict_dummy<<std::endl;
            std::cout<<dict.d_lim[i-1]<<" "<<dict.d_lim[i]<<" "<<dict.d_lim[i+1]<<" "<<dict.d_lim[i+2]<<" "<<dict.d_lim[i+3]<<std::endl;
        }*/
        //
        if(dict.dict[i]==dict.dict_dummy) continue;
        buckets[dict.dict[i]]++;
    }

    size_t acc=0, freq, max_freq=0;
    for(size_t i=0;i<dict.alphabet;i++){
        freq = buckets[i];
        buckets[i] = acc;
        acc +=freq;
        if(freq>max_freq) max_freq = freq;
    }
    buckets[dict.alphabet] = acc;

    vector_t freqs(dict.alphabet+1, 0, sdsl::bits::hi(max_freq)+1);
    bv_t solv_syms(dict.alphabet+1, false);
    size_t sym, l_sym, r_sym;
    for(size_t i=0;i<dict.dict.size();i++) {
        sym = dict.dict[i];

        if(sym==dict.dict_dummy){
            assert(i>0 && i<(dict.dict.size()-1));
            assert(!dict.d_lim[i] && !dict.d_lim[i-1]);//the dummy is nor the first nor the last
            l_sym = dict.dict[i-1];
            r_sym = dict.dict[i+1];

            //todo testing
            //std::cout<<l_sym<<" <- "<<" # "<<" -> "<<r_sym<<std::endl;
            //std::cout<<buckets[l_sym+1]-buckets[l_sym]<<" "<<buckets[r_sym+1]-buckets[r_sym]<<std::endl;
            assert(buckets[l_sym+1]-buckets[l_sym]==1);
            assert(buckets[r_sym+1]-buckets[r_sym]==1);
            //

            sa[buckets[l_sym]] = i<<1UL;
            sa[buckets[r_sym]] = (i+2)<<1UL;
            solv_syms[l_sym]= true;
            solv_syms[r_sym]= true;
        }else if(dict.d_lim[i]){
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
    freqs.erase();

    for(size_t i=0;i<solv_syms.size();i++){
        if(solv_syms[i]) buckets[i]++;
    }
    sdsl::util::clear(solv_syms);

    induce_L_type<value_type>(sa, dict, buckets);

    //move the pointer of every bucket to the last S symbol
    for(size_t bck=0;bck<dict.alphabet;bck++){
        size_t ptr = buckets[bck];
        size_t limit = buckets[bck+1];
        while(ptr<limit && sa[ptr]==0) ptr++;
        ptr--;
        buckets[bck] = ptr;
    }

    induce_S_type<value_type>(sa, dict, buckets);

    //TODO testing
    /*for(size_t i=273024;i<sa.size();i++){
        pos = (sa[i]>>1UL)-1;
        if((sa[i])>>1==0){
            std::cout<<"whut? "<<i<<std::endl;
        }
        std::cout<<(sa[i] & 1UL)<<" | ";
        do{
            if(dict.dict[pos]==dict.dict_dummy){
                std::cout<<"* ";
            }else{
                std::cout<<dict.dict[pos]<<" ";
            }
        }while(!dict.d_lim[pos++]);
        std::cout<<""<<std::endl;
    }*/
    //
    free(buckets);
}

template <class value_type>
void induce_L_type(vector_t &sa, const dictionary &dict, value_type* buckets) {

    size_t pos, l_sym, bck;
    vector_t ind_bck(dict.alphabet+1, 0, sdsl::bits::hi(sa.size())+1);
    bool lcs, first_eq=false;
    size_t lb=1;

    for(size_t i=0;i<sa.size();i++) {

        pos = sa[i];
        lcs = pos & 1UL;
        pos>>=1UL;
        if(pos==0) continue;

        if(!lcs) {
            first_eq = false;
            lb = i+1;
        }
        if(pos==1 || dict.d_lim[pos-2]) continue;

        if(dict.dict[pos-2]==dict.dict_dummy ||
           (pos>2 && dict.dict[pos-3]==dict.dict_dummy)){
            continue;
        }

        bck = dict.dict[pos-1];
        l_sym = dict.dict[pos-2];

        if(l_sym>=bck){
            if(!first_eq && l_sym==bck){
                first_eq = true;
                lcs = false;
            }
            sa[buckets[l_sym]++] = (pos-1)<<1UL | (ind_bck[l_sym]>=lb && lcs);
            ind_bck[l_sym] = i+1;
        }
    }
}

template <class value_type>
void induce_S_type(vector_t &sa, const dictionary &dict, value_type *buckets) {

    size_t pos, bck, l_sym, ind_pos;
    vector_t ind_bck(dict.alphabet+1, 0, sdsl::bits::hi(sa.size())+1);
    bool lcs, first_eq, new_break, p_lcs=false;
    size_t rb=sa.size();

    for(size_t i=sa.size();i-->0;) {

        pos = sa[i];
        lcs = pos & 1UL;
        pos>>=1UL;
        if(pos==0){
            sa[i+1] = sa[i+1] & ~1UL;
            p_lcs = false;
            continue;
        }

        if(!p_lcs) {
            first_eq = false;
            rb = i+1;
        }

        if(pos==1 || dict.d_lim[pos-2]){
            p_lcs = lcs;
            continue;
        }

        if(dict.dict[pos-2]==dict.dict_dummy ||
           (pos>2 && dict.dict[pos-3]==dict.dict_dummy)){
            p_lcs = lcs;
            continue;
        }

        bck = dict.dict[pos-1];
        l_sym = dict.dict[pos-2];

        if(l_sym < bck || (l_sym==bck && i >= buckets[bck])){

            ind_pos = buckets[l_sym]--;

            sa[ind_pos] = (pos-1)<<1UL | (lcs && ind_pos>0 && sa[ind_pos-1]==0);

            new_break = !first_eq && l_sym==bck;
            if(new_break) first_eq = true;

            if(ind_bck[l_sym]==0 || ind_bck[l_sym]>rb || new_break){
                sa[ind_pos+1] = sa[ind_pos+1] & ~1UL;
                if(ind_pos+1==i) lcs=false;
            }
            ind_bck[l_sym] = i+1;
        }

        p_lcs = lcs;
    }
}

template void suffix_induction<uint8_t>(vector_t &sa, const dictionary &dict);
template void suffix_induction<uint16_t>(vector_t &sa, const dictionary &dict);
template void suffix_induction<uint32_t>(vector_t &sa, const dictionary &dict);
template void suffix_induction<uint64_t>(vector_t &sa, const dictionary &dict);
