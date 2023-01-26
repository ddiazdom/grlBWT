//
// Created by Diaz, Diego on 4.1.2022.
//

#ifndef LPG_COMPRESSOR_LMS_INDUCTION_H
#define LPG_COMPRESSOR_LMS_INDUCTION_H
#include "common.h"
#include "bwt_io.h"

namespace exact_algo {

#define INSERT_BWT_SYM_FOR_SUFFIX(pos, bwt_sym) \
size_t phrase = d_lim_rs(pos); \
size_t freq = dict.freqs.read(phrase); \
size_t phrase_flag = freq & 3UL; \
assert(phrase_flag!=0 && phrase_flag<=3); \
freq>>=2UL;\
if(bwt_sym==dict.bwt_dummy){\
    if(bwt_sym==dict.bwt_dummy && phrase_flag==2){\
        bwt_sym = dict.end_str_dummy;       \
    }else if(bwt_sym==dict.bwt_dummy && phrase_flag==3) /*[[unlikely]]*/ { \
        /*TODO this is the case when a phrase occurs as a suffix but also as an entire string*/\
        exit(0);\
    }                                                \
    reduced_sa.push_back(pos);                \
    meta_sym_list.write(phrase, (met_sym<<1UL) | (freq>1)); \
    met_sym++;                                     \
}\
if(pre_bwt.size()>0 && pre_bwt.last_sym() == bwt_sym){\
    pre_bwt.inc_freq_last(freq);\
}else{\
    pre_bwt.push_back(bwt_sym, freq);\
}\
acc_bwt_syms+=freq                   \

#define INCREMENT_BWT(start, end,  bwt_sym)\
if(is_gr_phrase) { \
    bwt_sym = dict.bwt_dummy;\
    size_t samp_pos = (sa[start]>>2UL)-1;  \
    reduced_sa.push_back(samp_pos);           \
    if((end-start+1)>1){                   \
        bwt_sym++;                                   \
        nested_phrases.push_back(met_sym);               \
        while(!dict.d_lim[samp_pos]){                      \
            nested_phrases.push_back(dict.dict[samp_pos++]<<1UL); \
        }\
        nested_phrases.push_back((dict.dict[samp_pos++]<<1UL) | 1UL); \
        for(size_t j=start;j<=end;j++){           \
            meta_sym_marks[(sa[j]>>2UL)-1] = true; \
        }                                         \
    }                                             \
    if(n_phrases==1){                 \
        meta_sym_list.write(curr_phrase, (met_sym<<1UL) | (curr_phrase_freq>1));    \
    }                                      \
    met_sym++;                              \
}\
if(pre_bwt.size()>0 && pre_bwt.last_sym() == bwt_sym){\
    pre_bwt.inc_freq_last(acc_freq);\
}else{\
    pre_bwt.push_back(bwt_sym, acc_freq);\
}                                        \
acc_bwt_syms+=acc_freq                   \

template <class vector_type, class sa_type>
void induce_L_type(sa_type& sa, const dictionary &dict, vector_type& buckets, bv_t& solved_sym){

    size_t pos, l_sym, bck;
    vector_type ind_bck(dict.alphabet+1, 0);
    bool lcs, first_eq=false;
    size_t lb=1, is_suffix, flag;

    for(size_t i=0;i<sa.size();i++) {

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
}

template <class vector_type, class sa_type>
size_t induce_S_type(sa_type& sa, dictionary& dict, vector_type& buckets, bv_t& solved_sym, tmp_workspace& ws){

    size_t pos, bck, l_sym, ind_pos, sa_size = sa.size();
    vector_type ind_bck(dict.alphabet+1, 0);
    bool lcs, first_eq=true, new_break, p_lcs=false;
    size_t rb=sa.size(), is_suffix, flag;

    //TODO testing
    size_t n_breaks=0, n_phrases=0, null_sym = std::numeric_limits<size_t>::max(), pl_sym = null_sym, fe=sa.size()-1, met_sym=0;
    bool is_full_phrase, is_gr_phrase;
    size_t sb = INT_CEIL(sym_width(std::max(dict.alphabet+3, dict.prev_alphabet)), 8);
    size_t fb = INT_CEIL(sym_width(dict.t_size),8);
    bwt_buff_writer pre_bwt(ws.get_file("inv_pre_bwt"), std::ios::out, sb, fb);
    vector_t meta_sym_list(dict.n_phrases, 0, sym_width(dict.dict.size())+1);
    bv_rs_t d_lim_rs(&dict.d_lim);
    size_t prev_pos=dict.dict.size(), prev_is_suffix=false, acc_freq=0, curr_phrase=0, curr_phrase_freq=0;
    o_file_stream<size_t> reduced_sa(ws.get_file("s_sa"), BUFFER_SIZE, std::ios::out);
    o_file_stream<size_t> nested_phrases(ws.get_file("nested_phrases"), BUFFER_SIZE, std::ios::out);
    std::vector<std::pair<size_t, size_t>> same_suffix_range;
    size_t end_range_s_suffix=0, acc_bwt_syms=0;
    bv_t meta_sym_marks(dict.dict.size(), false);
    //

    for(size_t i=sa.size();i-->0;) {

        pos = sa[i];
        lcs = pos & 1UL;
        is_suffix = pos & 2UL;
        assert(is_suffix==2 || is_suffix==0);

        pos>>=2UL;
        assert(pos>0);

        if(!p_lcs) {
            first_eq = false;
            rb = i+1;
        }

        is_full_phrase = pos==1 || dict.d_lim[pos-2];
        l_sym = dict.bwt_dummy;

        if(!is_full_phrase){
            bck = dict.dict.read(pos-1);
            l_sym = dict.dict.read(pos-2);
            assert(l_sym<dict.end_str_dummy);

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
        }

        if(prev_is_suffix){
            assert(pl_sym!=null_sym);
            INSERT_BWT_SYM_FOR_SUFFIX(prev_pos-1, pl_sym);

            if(!(sa[i+1] & 1UL) && !dict.d_lim[prev_pos-1]){
                //std::cout<<"whut? "<<prev_pos-1<<" "<<!(sa[i+1] & 1UL)<<" "<<prev_is_suffix<<" "<<acc_bwt_syms<<" "<<end_range_s_suffix<<std::endl;
                same_suffix_range.emplace_back(acc_bwt_syms, end_range_s_suffix);
            }
        }

        if(i<(sa_size-1) && !(sa[i+1] & 1UL)){

            if(!prev_is_suffix && !dict.d_lim[prev_pos-1]){
                /*[[likely]]*/
                assert(n_phrases<=1);
                is_gr_phrase = n_phrases==1 || (n_breaks>1);
                INCREMENT_BWT((i+1), fe, pl_sym);
            }

            fe = i;
            n_phrases = 0;
            n_breaks = 0;
            acc_freq = 0;
            pl_sym = null_sym;
            end_range_s_suffix = acc_bwt_syms;
        }

        p_lcs = lcs;

        if(!is_suffix && !dict.d_lim[pos-1]){
            size_t d_rank = d_lim_rs(pos);
            size_t freq = dict.freqs.read(d_rank);
            acc_freq+= freq;
            if(is_full_phrase){
                curr_phrase = d_rank;
                curr_phrase_freq = freq;
            }
        }
        n_breaks +=pl_sym!=l_sym;
        n_phrases += is_full_phrase;

        pl_sym = l_sym;
        prev_pos = pos;
        prev_is_suffix = is_suffix;
    }

    if(prev_is_suffix){
        assert(pl_sym!=null_sym);
        INSERT_BWT_SYM_FOR_SUFFIX(prev_pos-1, pl_sym);

        if(!(sa[0] & 1UL) && !dict.d_lim[prev_pos-1]){
            //std::cout<<"whut? "<<prev_pos-1<<" "<<!(sa[0] & 1UL)<<" "<<prev_is_suffix<<" "<<acc_bwt_syms<<" "<<end_range_s_suffix<<std::endl;
            same_suffix_range.emplace_back(acc_bwt_syms, end_range_s_suffix);
        }
    }

    if(!(sa[0] & 1UL)){
        if(!prev_is_suffix && !dict.d_lim[prev_pos-1]) {
            assert(n_phrases<=1);
            is_gr_phrase = n_phrases==1 || (n_breaks>1);
            INCREMENT_BWT(0, fe, pl_sym);
        }
    }

    //TODO testing
    size_t cont=0, acc_fr=0;
    for(size_t i=0;i<sa.size();i++){
        pos = (sa[i]>>2UL)-1;
        size_t fr = dict.freqs[d_lim_rs(pos)];
        if(sa[i] & 2UL) fr>>=2UL;

        if(!dict.d_lim[pos] || (sa[i] & 2UL)){
            if(pos==0 || dict.d_lim[pos-1]){
                std::cout<<cont++<<" ::: "<<acc_fr<<" "<<fr<<" -> "<<pos<<" : "<<(sa[i] & 1UL)<<" | "<<((sa[i] & 2UL)!=0)<<" -> * | ";
            }else{
                std::cout<<cont++<<" ::: "<<acc_fr<<" "<<fr<<" -> "<<pos<<" : "<<(sa[i] & 1UL)<<" | "<<((sa[i] & 2UL)!=0)<<" -> "<<dict.dict[pos-1]<<" | ";
            }
            do{
                std::cout<<dict.dict[pos]<<" ";
            }while(!dict.d_lim[pos++]);
            std::cout<<" "<<std::endl;
            acc_fr+=fr;
        }
    }
    //

    pre_bwt.close();
    reduced_sa.close();
    nested_phrases.close();
    dict.freqs.erase();

    store_to_file(ws.get_file("phr_ranks"), meta_sym_list);
    store_to_file(ws.get_file("meta_sym_marks"), meta_sym_marks);
    store_pl_vector(ws.get_file("same_suffix_ranges"), same_suffix_range);

    return met_sym;
}

template <class vector_type, class sa_type>
size_t suffix_induction(dictionary &dict, tmp_workspace& ws){

    sa_type sa;
    vector_type buckets(dict.alphabet+1, 0);

    if constexpr (std::is_same<sa_type, vector_t>::value) {
        size_t width = sym_width(dict.dict.size())+2;
        sa.set_width(width);
        sa.resize(dict.dict.size());
        sa.initialize(0, sa.size());
    }else{
        sa.resize(dict.dict.size());
        memset(sa.data(), 0, sa.size()*sizeof(typename sa_type::value_type));
    }

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

    induce_L_type<vector_type>(sa,  dict, buckets, solv_syms);

    //move the pointer of every bucket to the last S symbol
    for(size_t bck=0;bck<dict.alphabet;bck++){
        size_t ptr = buckets[bck];
        size_t limit = buckets[bck+1];
        while(ptr<limit && sa[ptr]==0) ptr++;
        ptr--;
        buckets[bck] = ptr;
    }

    size_t n_phrases = induce_S_type<vector_type>(sa, dict, buckets, solv_syms, ws);

#ifdef __linux__
    malloc_trim(0);
#endif
    return n_phrases;
}

//template <class value_type>
//void increment_bwt(size_t start, size_t end, value_type *sa, const dictionary& dict, bool is_gr_phrase, size_t acc_freq, size_t bwt_sym, bwt_buff_writer& bwt_writer);
//void insert_bwt_sym_for_suffix(size_t pos, size_t bwt_sym, const dictionary& dict,bv_rs_t& d_lim_rs, bwt_buff_writer& bwt_writer);
//void invert_data(tmp_workspace& ws, size_t n_phrases, vector_t& phrase_metasymbols);

}

#endif //LPG_COMPRESSOR_LMS_INDUCTION_H
