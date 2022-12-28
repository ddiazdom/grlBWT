//
// Created by Diaz, Diego on 4.1.2022.
//
#include "LMS_induction.h"

template <class value_type>
value_type * suffix_induction(const dictionary &dict) {

    size_t sa_size = dict.dict.size();
    auto * sa = (value_type *)calloc(sa_size, sizeof(value_type));
    auto * buckets = (value_type *) calloc((dict.alphabet+1), sizeof(value_type));

    for(size_t i=0;i<dict.dict.size();i++){
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
    size_t sym, is_suffix=0;

    for(size_t i=dict.dict.size();i-->0;) {
        sym = dict.dict[i];
        freq = buckets[sym+1]-buckets[sym];
        if(dict.d_lim[i]) is_suffix = dict.is_suffix(sym) ? 2 : 0;

        if(freq==1){
            //if it is a suffix of a phrase, we don't need to mark it
            sa[buckets[sym]] = ((i+1)<<2UL) | is_suffix;
            solv_syms[sym]= true;
        }else if(dict.d_lim[i]){
            sa[buckets[sym+1]-1-freqs[sym]++] = ((i+1)<<2UL) | (is_suffix | 1UL);
        }
    }

    size_t pos;
    for(size_t i=0;i<dict.alphabet;i++){
        if(freqs[i]>0){
            pos = buckets[i+1]-freqs[i];
            sa[pos] &= ~1UL;
        }
    }
    freqs.erase();

    for(size_t i=0;i<solv_syms.size();i++){
        if(solv_syms[i]) buckets[i]++;
    }

    induce_L_type<value_type>(sa, sa_size, dict, buckets, solv_syms);

    //move the pointer of every bucket to the last S symbol
    for(size_t bck=0;bck<dict.alphabet;bck++){
        size_t ptr = buckets[bck];
        size_t limit = buckets[bck+1];
        while(ptr<limit && sa[ptr]==0) ptr++;
        ptr--;
        buckets[bck] = ptr;
    }

    induce_S_type<value_type>(sa, sa_size, dict, buckets, solv_syms);

    //TODO testing
    /*for(size_t i=0;i<sa_size;i++){
        pos = (sa[i]>>2UL)-1;
        std::cout<<pos<<" : "<<(sa[i] & 1UL)<<" | "<<((sa[i] & 2UL)!=0)<<" -> ";
        do{
            std::cout<<dict.dict[pos]<<" ";
        }while(!dict.d_lim[pos++]);
        std::cout<<""<<std::endl;
    }*/
    //
    free(buckets);
    return sa;
}

template <class value_type>
void induce_L_type(value_type* sa, size_t sa_size, const dictionary &dict, value_type* buckets, bv_t& solved_sym) {

    size_t pos, l_sym, bck;
    auto * ind_bck = (value_type*) calloc(dict.alphabet+1, sizeof(value_type));
    bool lcs, first_eq=false;
    size_t lb=1, is_suffix, flag;

    for(size_t i=0;i<sa_size;i++) {

        pos = sa[i];
        lcs = pos & 1UL;
        is_suffix = pos & 2UL;
        assert(is_suffix==2 || is_suffix==0);
        pos>>=2UL;
        if(pos==0) continue;

        if(!lcs) {
            first_eq = false;
            lb = i+1;
        }
        if(pos==1 || dict.d_lim[pos-2]) continue;

        bck = dict.dict.read(pos-1);
        l_sym = dict.dict.read(pos-2);

        if(!solved_sym[l_sym] && l_sym>=bck){
            if(!first_eq && l_sym==bck){
                first_eq = true;
                lcs = false;
            }

            flag = is_suffix | (ind_bck[l_sym]>=lb && lcs);
            sa[buckets[l_sym]++] = ((pos-1)<<2UL) | flag ;
            ind_bck[l_sym] = i+1;
        }
    }
    free(ind_bck);
}

template <class value_type>
void induce_S_type(value_type * sa, size_t sa_size, const dictionary &dict, value_type *buckets, bv_t& solved_sym) {

    size_t pos, bck, l_sym, ind_pos;
    auto * ind_bck = (value_type *) calloc(dict.alphabet+1, sizeof(value_type));
    bool lcs, first_eq, new_break, p_lcs=false;
    size_t rb=sa_size, is_suffix, flag;

    for(size_t i=sa_size;i-->0;) {

        pos = sa[i];
        lcs = pos & 1UL;
        is_suffix = pos & 2UL;
        assert(is_suffix==2 || is_suffix==0);

        pos>>=2UL;
        if(pos==0){
            sa[i+1] &=  ~1UL;
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

        bck = dict.dict.read(pos-1);
        l_sym = dict.dict.read(pos-2);

        if(!solved_sym[l_sym] && (l_sym < bck || (l_sym==bck && i >= buckets[bck]))){

            ind_pos = buckets[l_sym]--;
            flag = is_suffix | (lcs && ind_pos>0 && sa[ind_pos-1]==0);
            sa[ind_pos] = ((pos-1)<<2UL) | flag;

            new_break = !first_eq && l_sym==bck;
            if(new_break) first_eq = true;

            if(ind_bck[l_sym]==0 || ind_bck[l_sym]>rb || new_break){
                sa[ind_pos+1] &= ~1UL;
                if(ind_pos+1==i) lcs=false;
            }
            ind_bck[l_sym] = i+1;
        }
        p_lcs = lcs;
    }
    free(ind_bck);
}

template uint8_t* suffix_induction<uint8_t>(const dictionary &dict);
template uint16_t* suffix_induction<uint16_t>(const dictionary &dict);
template uint32_t* suffix_induction<uint32_t>(const dictionary &dict);
template uint64_t* suffix_induction<uint64_t>(const dictionary &dict);
