//
// Created by Diego Diaz on 4/7/20.
//

#include "lc_gram_algo.hpp"
#include <thread>
#include <chrono>
#include <cstdlib>
#include "utils.hpp"
#ifdef __linux__
#include <malloc.h>
#endif

void comp_dict_int(dictionary &dict, phrase_map_t& phrases_ht, vector_t& s_sa, bv_t& new_is_suffix) {

    //collapse the full dictionary
    bool found;
    dict.dict_dummy = dict.alphabet+s_sa.size()+1;
    size_t pos, len, p_len, em_nt, new_size=0;
    string_t phrase(2, sdsl::bits::hi(dict.alphabet)+1);
    vector_t new_dict(s_sa.size()*2, 0, sdsl::bits::hi(dict.dict_dummy)+1);
    //vector_t phrases_ptrs(s_sa.size(), 0, sdsl::bits::hi(dict.dict.size()+1));

    size_t rank=0;
    for(auto && u : s_sa) {
        pos = u;
        //get the greatest phrase's proper suffix in reverse order
        size_t tmp_pos = pos;
        while(!dict.d_lim[pos]) pos++;
        len = pos-tmp_pos+1;

        if(dict.is_suffix(dict.dict[tmp_pos+len-1])) new_is_suffix[rank] = true;

        phrase.clear();
        p_len = len;
        while(pos>tmp_pos){
            phrase.push_back(dict.dict[pos--]);
            p_len--;
        }


        //check if the proper suffixes exist as phrases
        found = false;
        while(!phrase.empty()){
            phrase.mask_tail();
            auto res = phrases_ht.find(phrase.data(), phrase.n_bits());
            if(res.second){
                em_nt = 0;
                phrases_ht.get_value_from(res.first, em_nt);
                found = true;
                break;
            }
            phrase.pop_back();
            p_len++;
        }

        if(!found) {
            new_dict[new_size++] = dict.dict_dummy;
            new_dict[new_size++] = dict.is_suffix(dict.dict[tmp_pos+len-1]) ?  dict.dict[tmp_pos+len-1] : dict.dict[tmp_pos+len-2];
            //new_dict[new_size++] = len==1?  dict.dict_dummy : dict.dict[tmp_pos+len-2];
            //new_dict[new_size++] = dict.dict[tmp_pos+len-1];
            /*while(!dict.d_lim[tmp_pos]){
                new_dict[new_size++] = dict.dict[tmp_pos++];
            }
            new_dict[new_size++] = dict.dict[tmp_pos];*/
        }else{
            new_dict[new_size++] = dict.dict[tmp_pos+p_len-1];
            new_dict[new_size++] = dict.alphabet+em_nt;

            /*j = 0;
            while(j<p_len){
                new_dict[new_size++] = dict.dict[tmp_pos+j];
                j++;
            }
            new_dict[new_size++] = dict.alphabet+em_nt;*/
        }
        //phrases_ptrs[u] = new_size-1;
        rank++;
    }

    //new_dict.resize(new_size);
    dict.dict.swap(new_dict);
    //dict.phrases_ptr.swap(phrases_ptrs);
    dict.n_phrases = s_sa.size();
    sdsl::util::clear(dict.d_lim);
}

void compress_dictionary(dictionary &dict, vector_t &sa, gram_info_t &p_gram,
                         bvb_t &r_lim, ivb_t &rules, phrase_map_t &mp_map, sdsl::cache_config& config) {

    size_t pos, prev_pos, lcs, l_sym, prev_l_sym, dummy_sym = dict.alphabet+1, rank=0, freq;

    std::string lvl_bwt_file = sdsl::cache_file_name("pbwt_lvl_"+std::to_string(p_gram.rules_breaks.size()), config);
    ivb_t lvl_bwt(lvl_bwt_file, std::ios::out);

    //remove unnecessary entries from the SA to
    // avoid dealing with corner cases
    bv_rs_t d_lim_rs(&dict.d_lim);
    size_t k=0;
    for(size_t j=0;j<sa.size();j++){
        if(sa[j]!=0) sa[k++] = sa[j];
    }
    sa.resize(k);

    size_t width = sdsl::bits::hi(dict.alphabet+dict.dict.size()-dict.n_phrases)+1;
    vector_t ranks(dict.n_phrases, 0, width);
    phrase_map_t new_phrases_ht;
    string_t phrase(2, sdsl::bits::hi(dict.alphabet)+1);

    bool is_maximal, exists_as_rule;
    size_t i=1, run_bg, n_suffixes=0, n_hocc=0;/*, rem_elms=0, new_phrases=0, existing_phrases=0*/;

    //dict.hocc_buckets.width(sdsl::bits::hi(dict.t_size)+1);
    //dict.hocc_buckets.resize(INT_CEIL(dict.n_phrases,4));

    while(i<sa.size()) {

        assert(sa[i]!=0);
        prev_pos = sa[i - 1] - 1;
        if(prev_pos==0 || dict.d_lim[prev_pos - 1]){
            exists_as_rule = true;
            ranks[d_lim_rs(prev_pos)] = rank;
            prev_l_sym = dummy_sym;
        }else{
            exists_as_rule = false;
            prev_l_sym = dict.dict[prev_pos - 1];
        }

        is_maximal = false;
        run_bg = i-1;

        //TODO
        /*size_t k=prev_pos;
        std::cout<<i-1<<" "<<prev_pos<<" "<<exists_as_rule<<" "<<is_maximal<<" "<<prev_l_sym<<"\t";
        while(!dict.d_lim[k]){
            std::cout<<dict.dict[k++]<<" ";
        }
        std::cout<<dict.dict[k]<<std::endl;*/
        //

        do {
            pos = sa[i]-1;
            lcs = 0;
            bool limit_a, limit_b, equal;

            do {
                equal = dict.dict[pos+lcs] == dict.dict[prev_pos+lcs];
                limit_a = dict.d_lim[pos+lcs] || dict.dict[pos+lcs]>dict.max_sym;
                limit_b = dict.d_lim[prev_pos+lcs] || dict.dict[prev_pos+lcs]>dict.max_sym;
                lcs++;
            }while(equal && !limit_a && !limit_b);

            if((equal+limit_a+limit_b)==3) {

                //TODO I think this ifelse statement can be simplified a bit
                if (pos==0 || dict.d_lim[pos - 1]) {
                    ranks[d_lim_rs(pos)] = rank;
                    exists_as_rule = true;
                    l_sym = dummy_sym;
                } else {
                    l_sym = dict.dict[pos - 1];
                    if (!is_maximal && prev_l_sym != dummy_sym) {
                        is_maximal = (dict.is_suffix(dict.dict[pos]) || lcs>1) && (prev_l_sym != l_sym);
                    }
                }
                if (l_sym != dummy_sym) prev_l_sym = l_sym;
                //

                prev_pos = pos;

                //TODO
                /*k=pos;
                std::cout<<i<<" "<<pos<<" "<<exists_as_rule<<" "<<is_maximal<<" "<<l_sym<<"\t";
                while(!dict.d_lim[k]){
                    std::cout<<dict.dict[k++]<<" ";
                }
                std::cout<<dict.dict[k]<<std::endl;*/
                //
            }else{
                lcs=0;
            }
            i++;
        } while (i<sa.size() && lcs!=0);
        //std::cout<<" "<<std::endl;

        assert(run_bg<=i-2);
        if(lcs!=0) assert(i==sa.size());
        size_t run_end = lcs==0 ? i-2 : i-1;
        freq=0;
        for(size_t j=run_bg;j<=run_end;j++){
            freq+=dict.freqs[d_lim_rs(sa[j]-1)];
        }

        if(exists_as_rule || is_maximal) {

            if((run_end-run_bg+1)>1) {

                //put the phrase in reverse order
                size_t tmp_pos = sa[run_bg]-1;
                size_t tmp_pos2 = tmp_pos;
                phrase.clear();
                while(!dict.d_lim[tmp_pos]) tmp_pos++;

                do{
                    phrase.push_back(dict.dict[tmp_pos--]);
                }while(tmp_pos>=tmp_pos2);
                //

                //hash the phrase to check then if it is
                // embedded in another phrase
                /*std::cout<<run_bg<<" "<<run_end<<" "<<rank+1+dict.max_sym<<" -> ";
                for(size_t u=phrase.size();u-->0;){
                    std::cout<<phrase[u]<<" ";
                }
                std::cout<<""<<std::endl;*/

                n_suffixes += phrase.size()==1;
                //hash only those phrases that can be proper
                // suffixes of other phrases
                //assert(phrase.size()>1);
                //if(phrase.size()>1){
                phrase.mask_tail();
                auto res = new_phrases_ht.insert(phrase.data(), phrase.n_bits(), rank);
                assert(res.second);
                //}

                /*std::cout<<sa[run_bg]<<" ";
                for(size_t u=phrase.size();u-->0;){
                    std::cout<<phrase[u]<<" ";
                }
                std::cout<<""<<std::endl;*/

                //invalidate the extra occurrences of the phrase
                /*for(size_t j=run_bg+1;j<=run_end;j++){
                    sa[j] = 0;
                }
                rem_elms+=run_end-run_bg;*/
                dict.phrases_has_hocc[rank] = true;

                //if(n_hocc==dict.hocc_buckets.size()){
                //    auto new_size = size_t(double(n_hocc)*1.15);
                //    dict.hocc_buckets.resize(new_size);
                //}
                //dict.hocc_buckets[n_hocc++] = freq;
            }

            /*k=sa[run_bg]-1;
            while(!dict.d_lim[k]){
                std::cout<<dict.dict[k++]<<" ";
            }
            std::cout<<dict.dict[k]<<" "<<run_bg<<" "<<run_end<<std::endl;*/
            lvl_bwt.push_back(dummy_sym);
            lvl_bwt.push_back(freq);
            sa[rank] = sa[run_bg]-1;
            rank++;
        }else if(!dict.d_lim[sa[run_bg]-1] ||
                  dict.is_suffix(dict.dict[sa[run_bg]-1])){
            assert(!is_maximal);
            size_t sym = dict.dict[sa[run_bg]-2];
            if(lvl_bwt.size()>2 && lvl_bwt[lvl_bwt.size()-2]==sym){
                lvl_bwt[lvl_bwt.size()-1] +=freq;
            }else{
                lvl_bwt.push_back(sym);
                lvl_bwt.push_back(freq);
            }
        }
    }

    //corner case
    if(lcs==0){
        freq = dict.freqs[d_lim_rs(sa[i-1]-1)];
        if(sa[i-1]==0 || dict.d_lim[sa[i-1]-2]){
            ranks[d_lim_rs(sa[i-1]-1)] = rank;
            sa[rank] = sa[i-1]-1;
            lvl_bwt.push_back(dummy_sym);
            lvl_bwt.push_back(freq);
            dict.phrases_has_hocc[rank++] = false;
        }else{
            size_t sym = dict.dict[sa[i-1]-2];
            if(lvl_bwt.size()>2 && lvl_bwt[lvl_bwt.size()-2]==sym){
                lvl_bwt[lvl_bwt.size()-1] +=freq;
            }else{
                lvl_bwt.push_back(sym);
                lvl_bwt.push_back(freq);
            }
        }
        //existing_phrases++;
    }

    sa.resize(rank);
    dict.phrases_has_hocc.resize(rank);

    //dict.hocc_buckets.resize(n_hocc+1);
    //size_t acc = 0, tmp;
    //for(size_t u=0;u<n_hocc;u++){
    //    tmp = dict.hocc_buckets[u];
    //    dict.hocc_buckets[u] = acc;
    //    acc +=tmp;
    //}
    //dict.hocc_buckets[n_hocc] = acc;
    //comp_dict_int(dict, new_phrases_ht, ranks, rank);
    //TODO just checking
    /*freq=0;
    for(size_t u=1;u<lvl_bwt.size();u+=2){
        freq+=lvl_bwt[u];
    }*/
    //

    //assign the ranks to the original phrases
    size_t j=0;
    for(auto const& ptr : mp_map){
        phrase_map_t::val_type val=0;
        mp_map.get_value_from(ptr, val);
        val = ranks[j++];
        mp_map.insert_value_at(ptr, val);
    }
    sdsl::util::clear(ranks);
    sdsl::util::clear(d_lim_rs);
    sdsl::util::clear(dict.freqs);
    new_phrases_ht.shrink_databuff();

#ifdef __linux__
    malloc_trim(0);
#endif

    //TODO this is for testing purposes
    /*bool found;
    size_t p_len, em_nt;
    tmp = rules.size();
    for(auto && u : sa) {

        pos = u;

        //TODO
        if(dict.d_lim[pos]){
            assert(dict.is_suffix(dict.dict[pos]));
        }
        //

        //extract the greatest proper suffix of the phrase
        // in reverse order
        size_t tmp_pos = pos;
        while(!dict.d_lim[pos]){
            pos++;
        }
        p_len = pos-tmp_pos+1;

        phrase.clear();
        while(pos>tmp_pos){
            phrase.push_back(dict.dict[pos--]);
            p_len--;
        }

        //check if the proper suffixes exists as phrases
        found = false;
        while(!phrase.empty()){
            //for(size_t u=phrase.size();u-->0;){
            //    std::cout<<phrase[u]<<" ";
            //}
            //std::cout<<""<<std::endl;
            phrase.mask_tail();
            auto res = new_phrases_ht.find(phrase.data(), phrase.n_bits());
            if(res.second){
                em_nt = 0;
                new_phrases_ht.get_value_from(res.first, em_nt);
                found = true;
                //std::cout<<" found: "<<dict.max_sym+em_nt<<std::endl;
                break;
            }
            phrase.pop_back();
            p_len++;
        }

        //std::cout<<dict.max_sym+rank<<" -> ";
        if(!found){
            while(!dict.d_lim[tmp_pos]){
                //assert(inv_sa_bv[tmp_pos]);
                //sa[inv_sa[inv_sa_bv_rs(tmp_pos)]] = k++;
                //std::cout<<dict.min_sym+dict.dict[tmp_pos]<<" ";
                rules.push_back(dict.min_sym+dict.dict[tmp_pos++]);
                r_lim.push_back(false);
            }

            rules.push_back(dict.min_sym+dict.dict[tmp_pos]);
            //std::cout<<dict.min_sym+dict.dict[tmp_pos]<<std::endl;
            r_lim.push_back(true);
        }else{
            j = 0;
            while(j<p_len){
                //assert(inv_sa_bv[tmp_pos+j]);
                //sa[inv_sa[inv_sa_bv_rs(tmp_pos+j)]] = k++;
                rules.push_back(dict.min_sym+dict.dict[tmp_pos+j]);
                //std::cout<<dict.min_sym+dict.dict[tmp_pos+j]<<" ";
                r_lim.push_back(false);
                j++;
            }
            //std::cout<<" "<<std::endl;

            //assert(em_nt>0);
            //std::cout<<dict.max_sym+em_nt<<std::endl;
            rules.push_back(dict.max_sym+em_nt+1);
            r_lim.push_back(true);
        }
    }*/
    //

    bv_t new_is_suffix(sa.size(), false);
    comp_dict_int(dict, new_phrases_ht, sa, new_is_suffix);
    std::string dict_file = sdsl::cache_file_name("dict_lvl_"+std::to_string(p_gram.rules_breaks.size()), config);
    std::string new_is_suffix_file = sdsl::cache_file_name("new_is_suffix", config);
    sdsl::store_to_file(dict, dict_file);
    sdsl::store_to_file(new_is_suffix, new_is_suffix_file);

    p_gram.rules_breaks.push_back(sa.size());
    p_gram.r += sa.size();

    //TODO this is for testing purposes
    /*for(size_t u=tmp, x=0;u<rules.size();u++, x++){
        if(dict.dict[x]<dict.alphabet){
            assert(dict.min_sym+dict.dict[x]==rules[u]);
        }else{
            assert(dict.max_sym+(dict.dict[x]-dict.alphabet)+1==rules[u]);
        }
    }*/
    //
}

void assign_ids(phrase_map_t &mp_map, ivb_t &r, bvb_t &r_lim,
                dictionary &dict, gram_info_t &p_gram, sdsl::cache_config &config) {

    std::cout<<"Suffix induction"<<std::endl;
    auto start = std::chrono::steady_clock::now();
    vector_t sa(dict.dict.size(), 0, sdsl::bits::hi(dict.dict.size())+1);
    suffix_induction(sa, dict.dict, dict.d_lim, dict.alphabet, config);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 4);

    std::cout<<"Number of phrases : "<<dict.n_phrases<<std::endl;
    std::cout<<"Alphabet : "<<dict.alphabet<<std::endl;

    std::cout<<"dict  "<<double(sdsl::size_in_bytes(dict.dict))/1000000<<std::endl;
    std::cout<<"freqs "<<double(sdsl::size_in_bytes(dict.freqs))/1000000<<std::endl;
    std::cout<<"sa    "<<double(sdsl::size_in_bytes(sa))/1000000<<std::endl;
    std::cout<<"has_hocc    "<<double(sdsl::size_in_bytes(dict.phrases_has_hocc))/1000000<<std::endl;

    /*bv_rs_t d_lim_rs(&dict.d_lim);
    vector_t ranks(dict.n_phrases, 0, sdsl::bits::hi(dict.n_phrases)+1);
    size_t rank = 0, pos;
    for(auto && i : sa){
        if(i==0) continue;
        pos = i-1;
        //write the dictionary in the grammar
        if(pos==0 || dict.d_lim[pos-1]){
            ranks[d_lim_rs(pos)] = rank++;
            size_t j=pos;
            while(!dict.d_lim[j]){
                r.push_back(dict.min_sym+dict.dict[j++]);
                r_lim.push_back(false);
            }
            assert(dict.d_lim[j]);
            r.push_back(dict.min_sym+dict.dict[j]);
            r_lim.push_back(true);
        }
    }*/

    compress_dictionary(dict, sa, p_gram, r_lim, r, mp_map, config);

    //assign the ranks
    /*size_t j=0;
    for(auto const& ptr : mp_map){
        //modify the key value
        phrase_map_t::val_type val=0;
        mp_map.get_value_from(ptr, val);
        val |= (dict.max_sym+ranks[j++]+1)<<1UL;
        mp_map.insert_value_at(ptr, val);
        //
    }
    p_gram.rules_breaks.push_back(rank);
    p_gram.r +=rank;*/
}

void join_parse_chunks(const std::string &output_file, std::vector<std::string> &chunk_files) {

    //concatenate the files
    std::ofstream of(output_file, std::ofstream::binary);
    size_t buff_size = BUFFER_SIZE/sizeof(size_t);
    size_t len, rem, to_read, start, end;
    auto *buffer = new size_t[buff_size];

    for(auto const& file: chunk_files){

        std::ifstream i_file(file, std::ifstream::binary);

        i_file.seekg (0, std::ifstream::end);
        len = i_file.tellg()/sizeof(size_t);
        i_file.seekg (0, std::ifstream::beg);

        rem=len;
        to_read = std::min<size_t>(buff_size, len);

        while(true){

            i_file.seekg( (rem - to_read) * sizeof(size_t));
            i_file.read((char *)buffer, sizeof(size_t)*to_read);
            assert(i_file.good());

            //invert data
            start =0;end=to_read-1;
            while(start<end){
                std::swap(buffer[start++], buffer[end--]);
            }

            of.write((char *)buffer, sizeof(size_t)*to_read);
            assert(of.good());

            rem -= i_file.gcount()/sizeof(size_t);
            to_read = std::min<size_t>(buff_size, rem);
            if(to_read == 0) break;
        }
        i_file.close();

        if(remove(file.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }
    }
    delete[] buffer;
    of.close();
}

std::pair<size_t, size_t> join_thread_phrases(phrase_map_t& map, std::vector<std::string> &files) {

    size_t dic_bits=0, freq, max_freq=0;

    for(auto const& file : files){

        std::ifstream text_i(file, std::ios_base::binary);

        text_i.seekg (0, std::ifstream::end);
        size_t tot_bytes = text_i.tellg();
        text_i.seekg (0, std::ifstream::beg);

        auto buffer = reinterpret_cast<char *>(malloc(tot_bytes));

        text_i.read(buffer, std::streamsize(tot_bytes));

        bitstream<ht_buff_t> bits;
        bits.stream = reinterpret_cast<ht_buff_t*>(buffer);

        size_t next_bit = 32;
        size_t tot_bits = tot_bytes*8;
        size_t key_bits;
        size_t value_bits=sizeof(size_t)*8;
        void* key=nullptr;
        size_t max_key_bits=0;

        while(next_bit<tot_bits){

            key_bits = bits.read(next_bit-32, next_bit-1);

            size_t n_bytes = INT_CEIL(key_bits, bitstream<ht_buff_t>::word_bits)*sizeof(ht_buff_t);
            if(key_bits>max_key_bits){
                if(key==nullptr){
                    key = malloc(n_bytes);
                }else {
                    key = realloc(key, n_bytes);
                }
                max_key_bits = key_bits;
            }

            char *tmp = reinterpret_cast<char*>(key);
            tmp[INT_CEIL(key_bits, 8)-1] = 0;

            bits.read_chunk(key, next_bit, next_bit+key_bits-1);
            next_bit+=key_bits;
            freq = bits.read(next_bit, next_bit+value_bits-1);
            next_bit+=value_bits+32;

            auto res = map.insert(key, key_bits, freq);
            if(!res.second){
                size_t val;
                map.get_value_from(res.first, val);
                val+=freq;
                if(val>max_freq) max_freq = val;
                map.insert_value_at(res.first, val);
            }else{
                if(freq>max_freq) max_freq = freq;
                dic_bits+=key_bits;
            }
        }
        text_i.close();

        if(remove(file.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }
        free(key);
        free(buffer);
    }
    map.shrink_databuff();
    return {dic_bits, max_freq};
}


template<template<class, class> class lc_parser_t>
size_t build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size,
                      gram_info_t &p_gram, alpha_t alphabet, sdsl::cache_config &config) {

    std::cout<<"Generating the grammar from the text:    "<<std::endl;
    std::string output_file = sdsl::cache_file_name("tmp_output", config);
    std::string tmp_i_file = sdsl::cache_file_name("tmp_input", config);

    // mark which symbols represent string boundaries
    bv_t symbol_desc(alphabet.back().first+1,false);
    symbol_desc[alphabet[0].first] = true;

    ivb_t rules(p_gram.rules_file, std::ios::out, BUFFER_SIZE);
    bvb_t rules_lim(p_gram.rules_lim_file, std::ios::out);
    for(size_t i=0;i<p_gram.r; i++){
        rules.push_back(i);
        rules_lim.push_back(true);
    }
    for(auto const& pair : p_gram.sym_map){
        rules[pair.first] = pair.first;
    }

    size_t iter=1;
    size_t rem_phrases;

    typedef lc_parser_t<i_file_stream<uint8_t>, string_t> byte_parser_t;
    typedef lc_parser_t<i_file_stream<size_t>, string_t> int_parser_t;

    std::cout<<"  Parsing round "<<iter++<<std::endl;
    auto start = std::chrono::steady_clock::now();
    rem_phrases = build_lc_gram_int<byte_parser_t>(i_file, tmp_i_file,
                                             n_threads, hbuff_size,
                                             p_gram, rules, rules_lim,
                                             symbol_desc, config);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 4);

    while (rem_phrases > 0) {
        start = std::chrono::steady_clock::now();
        std::cout<<"  Parsing round "<<iter++<<std::endl;
        rem_phrases = build_lc_gram_int<int_parser_t>(tmp_i_file, output_file,
                                                      n_threads, hbuff_size,
                                                      p_gram, rules, rules_lim,
                                                      symbol_desc, config);
        end = std::chrono::steady_clock::now();

        report_time(start, end,4);
        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
    }

    sdsl::util::clear(symbol_desc);

    /*
    {//put the compressed string at end
        std::ifstream c_vec(tmp_i_file, std::ifstream::binary);
        c_vec.seekg(0, std::ifstream::end);
        size_t tot_bytes = c_vec.tellg();
        c_vec.seekg(0, std::ifstream::beg);
        auto *buffer = reinterpret_cast<size_t*>(malloc(BUFFER_SIZE));
        size_t read_bytes =0;
        //size_t min_sym = p_gram.rules_breaks.empty() ? 0 : p_gram.r - p_gram.rules_breaks.back();
        //p_gram.c=0;

        size_t len = tot_bytes/sizeof(size_t);
        vector_t bwt(len*2, 0, sdsl::bits::hi(len)+1);
        size_t pos = 0;
        while(read_bytes<tot_bytes){
            c_vec.read((char *) buffer, BUFFER_SIZE);
            read_bytes+=c_vec.gcount();
            assert((c_vec.gcount() % sizeof(size_t))==0);
            for(size_t i=0;i<c_vec.gcount()/sizeof(size_t);i++){

                if(pos==0 || buffer[i]!=bwt[pos-2]){
                    bwt[pos++] = buffer[i];
                    bwt[pos++] = 1;
                }else{
                    assert(bwt[pos-2]==buffer[i]);
                    bwt[pos-1]++;
                }
                //rules.push_back(min_sym+buffer[i]);
                //rules_lim.push_back(false);
                //p_gram.c++;
            }
        }
        bwt.resize(pos);
        sdsl::store_to_file(bwt, sdsl::cache_file_name("bwt_lvl_"+std::to_string(iter-2), config));
        //rules_lim[rules_lim.size() - 1] = true;
        //p_gram.r++;
        c_vec.close();
        free(buffer);
    }
    p_gram.g = rules.size();
    p_gram.n_p_rounds = p_gram.rules_breaks.size();

    std::vector<size_t> rule_breaks;
    rule_breaks.push_back(p_gram.max_tsym+1);
    for(unsigned long i : p_gram.rules_breaks){
        rule_breaks.push_back(rule_breaks.back()+i);
    }
    std::swap(p_gram.rules_breaks, rule_breaks);*/
    //rules.close();
    //rules_lim.close();
    //if(remove(tmp_i_file.c_str())){
    //    std::cout<<"Error trying to delete file "<<tmp_i_file<<std::endl;
    //}
    return iter-2;
}

template<class parser_t, class out_sym_t>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file,
                         size_t n_threads, size_t hbuff_size,
                         gram_info_t &p_gram, ivb_t &rules,
                         bvb_t &rules_lim, bv_t &phrase_desc,
                         sdsl::cache_config &config) {

#ifdef __linux__
    malloc_trim(0);
#endif
    //std::cout<<"phrases desc: "<<double(sdsl::size_in_bytes(phrase_desc))/1000000<<" MB "<<std::endl;
    //std::cout<<"0. ";report_mem_peak();

    typedef typename parser_t::sym_type          sym_type;
    typedef typename parser_t::stream_type       stream_type;
    typedef parse_data_t<stream_type, out_sym_t> parse_data_type;

    parser_t parser(phrase_desc);
    phrase_map_t mp_table(0, "", 0.8);

    auto thread_ranges = parser.partition_text(n_threads, i_file);

    std::vector<parse_data_type> threads_data;
    threads_data.reserve(thread_ranges.size());

    //how many size_t cells we can fit in the buffer
    size_t buff_cells = hbuff_size/sizeof(size_t);

    //number of bytes per thread
    size_t hb_bytes = (buff_cells / thread_ranges.size()) * sizeof(size_t);

    {
        //std::cout << "1. ";report_mem_peak();
        std::cout << "    Computing the phrases in the text" << std::endl;

        void *buff_addr = malloc(hbuff_size);
        auto tmp_addr = reinterpret_cast<char *>(buff_addr);
        size_t k = 0;
        for (auto &range: thread_ranges) {
            std::string tmp_o_file = o_file.substr(0, o_file.size() - 5);
            tmp_o_file.append("_range_" + std::to_string(range.first) + "_" + std::to_string(range.second));
            threads_data.emplace_back(i_file, tmp_o_file, mp_table, range.first, range.second, hb_bytes,
                                      tmp_addr + (k * hb_bytes));
            k++;
        }
        //std::cout << "1.1 ";report_mem_peak();

        std::vector<std::thread> threads(threads_data.size());
        hash_functor<parse_data_type, parser_t> hf;
        for (size_t i = 0; i < threads_data.size(); i++) {
            threads[i] = std::thread(hf, std::ref(threads_data[i]), std::ref(parser));
        }

        for (size_t i = 0; i < threads_data.size(); i++) {
            threads[i].join();
        }
        free(buff_addr);

#ifdef __linux__
        malloc_trim(0);
#endif
    }
    //std::cout<<"2. ";report_mem_peak();

    //join the different phrase files
    std::vector<std::string> phrases_files;
    for(size_t i=0;i<threads_data.size();i++){
        phrases_files.push_back(threads_data[i].thread_dict.dump_file());
    }
    auto join_res = join_thread_phrases(mp_table, phrases_files);

    size_t psize=0;//<- for the iter stats
    if(mp_table.size()!=p_gram.last_dict_size) {
        size_t width = sdsl::bits::hi(phrase_desc.size())+1;
        size_t min_sym = p_gram.rules_breaks.empty() ? 0 : p_gram.r - p_gram.rules_breaks.back();
        size_t max_sym = p_gram.r-1;
        size_t dict_syms = join_res.first/width;
        size_t max_freq = join_res.second;

        //p_gram.rules_breaks.push_back(mp_table.size());
        const bitstream<ht_buff_t>& stream = mp_table.get_data();
        key_wrapper key_w{width, mp_table.description_bits(), stream};

        //temporal unload of the hash table (not the data)
        std::string st_table = sdsl::cache_file_name("ht_data", config);
        mp_table.unload_table(st_table);
        //std::cout<<"3. ";report_mem_peak();

        {
            //create a dictionary from where the ids will be computed
            std::cout<<" Creating the dictionary from the hash table"<<std::endl;
            dictionary dict(mp_table, min_sym, max_sym, key_w, dict_syms, max_freq, phrase_desc, threads_data[0].ifs.size());
            //std::cout<<"4. ";report_mem_peak();
            //rename phrases according to their lexicographical ranks
            std::cout << "    Assigning identifiers to the phrases" << std::endl;
            assign_ids(mp_table, rules, rules_lim, dict, p_gram, config);
            //std::cout<<"5. ";report_mem_peak();
        }

        //reload the hash table
        mp_table.load_table(st_table);
        if(remove(st_table.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }

        std::cout<<"    Creating the parse of the text"<<std::endl;
        {//store the phrases into a new file
            std::vector<std::thread> threads(threads_data.size());
            parse_functor<parse_data_type, parser_t> pf;
            for(size_t i=0;i<threads_data.size();i++){
                threads[i] = std::thread(pf, std::ref(threads_data[i]), std::ref(parser));
            }

            for(size_t i=0;i<threads_data.size();i++) {
                threads[i].join();
            }
        }
        //std::cout<<"6. ";report_mem_peak();

        std::vector<std::string> chunk_files;
        for(size_t i=0;i<threads_data.size();i++){
            chunk_files.push_back(threads_data[i].ofs.file);
        }
        join_parse_chunks(o_file, chunk_files);

        {// this is just to get the size of the resulting parse
            i_file_stream<size_t> ifs(o_file, BUFFER_SIZE);
            psize = ifs.tot_cells;
        }

        {
            //keep track of the phrases that have to be rephrased
            //std::cout <<"    Updating symbols status" << std::endl;
            //bv_t new_phrase_desc(p_gram.rules_breaks.back(), false);
            bv_t new_is_suffix;
            sdsl::load_from_file(new_is_suffix, sdsl::cache_file_name("new_is_suffix", config));

            //phrase_desc.resize(p_gram.r);
            /*auto it = mp_table.begin();
            auto it_end = mp_table.end();
            size_t sym;
            while (it != it_end) {
                auto val = it.value();
                //read the (reversed) last symbol
                sym = key_w.read(*it, 0);
                new_phrase_desc[val] = phrase_desc[sym];
                assert(new_phrase_desc[val]==new_is_suffix[val]);
                ++it;
            }*/
            phrase_desc.swap(new_is_suffix);
        }

        p_gram.last_dict_size = mp_table.size();
        std::cout<<"    Stats:"<<std::endl;
        std::cout<<"      Parse size:          "<<psize<<std::endl;
        std::cout<<"      New nonterminals:    "<<p_gram.rules_breaks.back()<<std::endl;
        //std::cout<<"7. ";report_mem_peak();
        return mp_table.size();

    }else{ //just copy the input

        std::ifstream in(i_file, std::ios_base::binary);
        std::ofstream out(o_file, std::ios_base::binary);
        auto buffer = reinterpret_cast<char*>(malloc(BUFFER_SIZE));
        do {
            in.read(&buffer[0], BUFFER_SIZE);
            out.write(&buffer[0], in.gcount());
            psize+=in.gcount();
        } while (in.gcount() > 0);
        free(buffer);
        psize/=sizeof(sym_type);

        //remove remaining files
        for(size_t i=0;i<threads_data.size();i++){
            std::string tmp_file =  threads_data[i].ofs.file;
            if(remove(tmp_file.c_str())){
                std::cout<<"Error trying to delete file "<<tmp_file<<std::endl;
            }
        }
        std::cout<<"    No new phrases found"<<std::endl;
        return 0;
    }
}

template size_t build_lc_gram<lms_parsing>(std::string &i_file, size_t n_threads, size_t hbuff_size, gram_info_t &p_gram, alpha_t alphabet, sdsl::cache_config &config);

template size_t build_lc_gram_int<lms_parsing<i_file_stream<uint8_t>, string_t>>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                           gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                                           bv_t &phrase_desc, sdsl::cache_config &config);

template size_t build_lc_gram_int<lms_parsing<i_file_stream<size_t>, string_t>>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                          gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                                          bv_t &phrase_desc, sdsl::cache_config &config);
