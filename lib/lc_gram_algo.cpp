//
// Created by Diego Diaz on 4/7/20.
//
#include "lc_gram_algo.hpp"
#include <thread>
#include <chrono>
#include <cstdlib>
#include "utils.h"
#include "LMS_induction.h"
#include "bwt_io.h"
#ifdef __linux__
#include <malloc.h>
#endif

void dict2gram(dictionary &dict, phrase_map_t& phrases_ht, vector_t& s_sa, bv_t& phr_marks,
               parsing_info& p_info, tmp_workspace& ws) {

    dict.dict_dummy = dict.alphabet+s_sa.size()+1;
    size_t pos, em_nt, new_size=0, sym;
    string_t phrase(2, sdsl::bits::hi(dict.alphabet)+1);
    vector_t new_dict(s_sa.size()*2, 0, sdsl::bits::hi(dict.dict_dummy)+1);

    size_t rank=0;
    for(size_t u=0;u<s_sa.size();u++) {
        pos = s_sa[u];

        if(p_info.p_round==1 && u==1633){
            std::cout<<"holaa "<<dict.dict[pos]<<std::endl;
        }

        if(dict.d_lim[pos]){
            new_dict[new_size++] = dict.dict_dummy;
            new_dict[new_size++] = dict.dict[pos];
        }else{
            pos++;
            while(!phr_marks[pos] && !dict.d_lim[pos]) pos++;
            if(phr_marks[pos]){
                size_t tmp_pos = pos;
                while(!dict.d_lim[pos]) pos++;
                phrase.clear();
                while(pos>=tmp_pos){
                    phrase.push_back(dict.dict[pos--]);
                }
                phrase.mask_tail();

                auto res = phrases_ht.find(phrase.data(), phrase.n_bits());

                assert(res.second);

                em_nt = 0;
                phrases_ht.get_value_from(res.first, em_nt);
                new_dict[new_size++] = dict.dict[tmp_pos-1];
                new_dict[new_size++] = dict.alphabet+em_nt;
            }else{
                sym = dict.dict[pos];
                new_dict[new_size++] = dict.dict_dummy;
                new_dict[new_size++] = dict.is_suffix(sym) ? sym : dict.dict[pos-1];
            }
        }
        rank++;
    }

    dict.dict.swap(new_dict);
    dict.n_phrases = s_sa.size();
    sdsl::util::clear(dict.d_lim);

    std::string dict_file = ws.get_file("dict_lev_"+std::to_string(p_info.p_round));
    sdsl::store_to_file(dict, dict_file);
}

void dict2gram2(dictionary &dict, vector_t& s_sa, vector_t& first_symbol, parsing_info& p_info, tmp_workspace& ws) {

    dict.dict_dummy = dict.alphabet+s_sa.size()+1;
    size_t pos, new_size=0, l_sym, r_sym;
    vector_t new_dict(s_sa.size()*2, 0, sdsl::bits::hi(dict.dict_dummy)+1);

    std::cout<<"bytes compressed dict "<<INT_CEIL(new_dict.size()*new_dict.width(), 8)<<std::endl;

    for(size_t u=0;u<s_sa.size();u++) {
        pos = s_sa.read(u);

        if(dict.d_lim[pos]){
            l_sym = dict.dict.read(pos);
            assert(l_sym>=dict.alphabet);
            l_sym = first_symbol.read(l_sym-dict.alphabet);
            assert(dict.is_suffix(l_sym));

            new_dict.write(new_size++, dict.dict_dummy);
            new_dict.write(new_size++, l_sym);
        }else{
            pos++;
            while(dict.dict.read(pos)<dict.alphabet && !dict.d_lim[pos]) pos++;

            l_sym = dict.dict.read(pos-1);
            if(l_sym >= dict.alphabet) l_sym = first_symbol.read(l_sym-dict.alphabet);
            r_sym = dict.dict.read(pos);

            if(r_sym>=dict.alphabet){
                new_dict.write(new_size++, l_sym);
                new_dict.write(new_size++, r_sym);
            }else{
                if(r_sym >= dict.alphabet) r_sym = first_symbol.read(r_sym-dict.alphabet);
                new_dict.write(new_size++, dict.dict_dummy);
                new_dict.write(new_size++, dict.is_suffix(r_sym) ? r_sym : l_sym);
            }
        }
    }

    dict.dict.swap(new_dict);
    dict.n_phrases = s_sa.size();
    sdsl::util::clear(dict.d_lim);

    std::string dict_file = ws.get_file("dict_lev_"+std::to_string(p_info.p_round));
    sdsl::store_to_file(dict, dict_file);

    //TODO testing
    /*std::cout<<"\nChecking if they are correct .."<<std::endl;
    dictionary orig_dict;
    sdsl::load_from_file(orig_dict, "/Users/ddiaz/version_original_"+std::to_string(p_info.p_round));
    for (size_t i = 0; i < dict.dict.size(); i += 2) {
        if (dict.dict[i] != orig_dict.dict[i] ||
            dict.dict[i + 1] != orig_dict.dict[i + 1]) {
            std::cout << i << " -> (" << dict.dict[i] << ", " << dict.dict[i + 1] << " != " << orig_dict.dict[i]
                      << ", " << orig_dict.dict[i + 1] << ")" << std::endl;
            exit(1);
        }
    }
    assert(dict.alphabet==orig_dict.alphabet);
    assert(dict.p_alpha_size==orig_dict.p_alpha_size);
    assert(dict.n_phrases==orig_dict.n_phrases);
    assert(dict.t_size==orig_dict.t_size);
    assert(dict.max_sym_freq==orig_dict.max_sym_freq);
    for(size_t i=0;i<dict.phrases_has_hocc.size();i++){
        assert(dict.phrases_has_hocc[i]==orig_dict.phrases_has_hocc[i]);
    }
    //*/
#ifdef __linux__
    malloc_trim(0);
#endif

    //TODO testing
    //sdsl::store_to_file(dict, "/Users/ddiaz/version_nueva_"+std::to_string(p_info.p_round));
    //
}

void get_pre_bwt2(dictionary &dict, vector_t &sa, parsing_info& p_info, tmp_workspace& ws) {

    bool is_maximal, exist_as_phrase;
    size_t u=0, d_pos, pl_sym, bg_pos, freq, rank=0, l_sym, dummy_sym = dict.alphabet+1, f_sa_pos, d_rank, phrase_frq;

    bv_rs_t d_lim_rs(&dict.d_lim);

    std::string pre_bwt_file = ws.get_file("pre_bwt_lev_"+std::to_string(p_info.p_round));
    size_t sb = INT_CEIL(sdsl::bits::hi(dummy_sym)+1, 8);
    size_t fb = INT_CEIL(sdsl::bits::hi(dict.t_size)+1,8);
    bwt_buff_writer pre_bwt(pre_bwt_file, std::ios::out, sb, fb);

    size_t width = sdsl::bits::hi(dict.dict.size())+2;
    vector_t ranks(dict.n_phrases, 0, width);
    vector_t first_symbol(dict.n_phrases, sym_width(dict.alphabet));

    while(u<sa.size()) {
        d_pos = (sa.read(u)>>1UL) - 1;

        if(dict.d_lim[d_pos] && !dict.is_suffix(dict.dict.read(d_pos))){
            u++;
            while(u<sa.size() && sa.read(u) & 1UL) u++;
        }else{

            f_sa_pos = d_pos;
            freq = dict.freqs.read(d_lim_rs(d_pos));
            bg_pos = u;
            is_maximal = false;

            if(d_pos==0 || dict.d_lim[d_pos-1]){
                pl_sym = dummy_sym;
            }else{
                pl_sym = dict.dict[d_pos-1];
                if(pl_sym>=dict.alphabet) pl_sym = first_symbol[pl_sym-dict.alphabet];
            }

            if(pl_sym==dummy_sym){
                exist_as_phrase = true;
                d_rank = d_lim_rs(d_pos);
                ranks.write(d_rank, (rank<<1UL) | (dict.freqs[d_rank]>1));
            }else{
                exist_as_phrase = false;
            }

            u++;
            while(u<sa.size() && sa.read(u) & 1UL){
                d_pos = (sa.read(u)>>1UL) - 1;
                phrase_frq = dict.freqs.read(d_lim_rs(d_pos));
                freq += phrase_frq;

                if(d_pos==0 || dict.d_lim[d_pos-1]){
                    l_sym = dummy_sym;
                }else{
                    l_sym = dict.dict.read(d_pos-1);
                    if(l_sym>=dict.alphabet) l_sym = first_symbol.read(l_sym-dict.alphabet);
                }

                if(!is_maximal && l_sym!=pl_sym) is_maximal = true;
                if(!exist_as_phrase && l_sym==dummy_sym){
                    exist_as_phrase = true;
                    d_rank = d_lim_rs(d_pos);
                    ranks.write(d_rank, (rank << 1UL) | (dict.freqs.read(d_rank)>1));
                }
                pl_sym = l_sym;
                u++;
            }

            if(is_maximal || exist_as_phrase) {

                pre_bwt.push_back(dummy_sym, freq);

                size_t f_sym = dict.dict.read(f_sa_pos);
                first_symbol.push_back(f_sym);

                if(u-bg_pos>1){
                    dict.phrases_has_hocc[rank] = true;
                    for(size_t j=bg_pos;j<u;j++){
                        size_t pos = (sa.read(j)>>1UL)-1;
                        assert(dict.dict.read(pos)==f_sym);
                        dict.dict.write(pos, dict.alphabet + rank);
                    }
                }
                sa.write(rank, f_sa_pos);
                rank++;
            }else{
                l_sym = dict.dict.read(f_sa_pos-1);
                if(l_sym>=dict.alphabet) l_sym = first_symbol[l_sym-dict.alphabet];

                assert(l_sym<dict.alphabet);
                if(pre_bwt.size()>1 && pre_bwt.last_sym() == l_sym){
                    pre_bwt.inc_freq_last(freq);
                }else{
                    pre_bwt.push_back(l_sym, freq);
                }
            }
        }
    }

    assert(rank<dict.dict.size());
    pre_bwt.close();
    sa.resize(rank);
    first_symbol.resize(rank);

    dict.phrases_has_hocc.resize(rank);
    sdsl::util::clear(d_lim_rs);
    dict.freqs.erase();
    store_to_file(ws.get_file("phr_ranks"), ranks);
    ranks.erase();

#ifdef __linux__
    malloc_trim(0);
#endif

    std::cout<<"\nbytes SA "<<INT_CEIL(sa.size()*sa.width(), 8)<<std::endl;
    std::cout<<"bytes dict "<<INT_CEIL(dict.dict.size()*dict.dict.width(), 8)<<std::endl;
    std::cout<<"bytes first_symbol "<<INT_CEIL(first_symbol.size()*first_symbol.width(), 8)<<std::endl;
    dict2gram2(dict, sa, first_symbol, p_info, ws);
}

void get_pre_bwt(dictionary &dict, vector_t &sa, parsing_info& p_info, bv_t& phr_marks,
                 phrase_map_t& new_phrases_ht, tmp_workspace& ws) {

    string_t phrase(2, sdsl::bits::hi(dict.alphabet)+1);
    bool is_maximal, exist_as_phrase;
    size_t u=0, d_pos, pl_sym, bg_pos, freq, rank=0, l_sym, dummy_sym = dict.alphabet+1, f_sa_pos, d_rank, phrase_frq;

    bv_rs_t d_lim_rs(&dict.d_lim);

    std::string pre_bwt_file = ws.get_file("pre_bwt_lev_"+std::to_string(p_info.p_round));
    size_t sb = INT_CEIL(sdsl::bits::hi(dummy_sym)+1, 8);
    size_t fb = INT_CEIL(sdsl::bits::hi(dict.t_size)+1,8);
    bwt_buff_writer pre_bwt(pre_bwt_file, std::ios::out, sb, fb);

    size_t width = sdsl::bits::hi(dict.dict.size())+2;
    vector_t ranks(dict.n_phrases, 0, width);

    while(u<sa.size()) {

        d_pos = (sa[u]>>1UL) - 1;
        if(dict.d_lim[d_pos] && !dict.is_suffix(dict.dict[d_pos])){
            u++;
            while(u<sa.size() && sa[u] & 1UL) u++;
        }else{
            f_sa_pos = d_pos;
            freq = dict.freqs[d_lim_rs(d_pos)];
            bg_pos = u;
            is_maximal = false;
            pl_sym = (d_pos==0 || dict.d_lim[d_pos-1]) ? dummy_sym : dict.dict[d_pos-1];
            if(pl_sym==dummy_sym){
                exist_as_phrase = true;
                d_rank = d_lim_rs(d_pos);
                ranks[d_rank] = (rank<<1UL) | (dict.freqs[d_rank]>1);
            }else{
                exist_as_phrase = false;
            }

            u++;
            while(u<sa.size() && sa[u] & 1UL){
                d_pos = (sa[u]>>1UL) - 1;
                phrase_frq = dict.freqs[d_lim_rs(d_pos)];
                freq += phrase_frq;
                l_sym = d_pos==0 || dict.d_lim[d_pos-1] ? dummy_sym : dict.dict[d_pos-1];
                if(!is_maximal && l_sym!=pl_sym) is_maximal = true;
                if(!exist_as_phrase && l_sym==dummy_sym){
                    exist_as_phrase = true;
                    d_rank = d_lim_rs(d_pos);
                    ranks[d_rank] = (rank << 1UL) | (dict.freqs[d_rank]>1);
                }
                pl_sym = l_sym;
                u++;
            }

            if(is_maximal || exist_as_phrase){

                if(u-bg_pos>1) {

                    //put the phrase in reverse order
                    size_t tmp_pos = f_sa_pos;
                    phrase.clear();
                    while(!dict.d_lim[tmp_pos]) tmp_pos++;
                    do{
                        phrase.push_back(dict.dict[tmp_pos--]);
                    }while(tmp_pos>= f_sa_pos);
                    phrase.mask_tail();
                    auto res = new_phrases_ht.insert(phrase.data(), phrase.n_bits(), rank);
                    assert(res.second);
                    //

                    dict.phrases_has_hocc[rank] = true;

                    for(size_t j=bg_pos;j<u;j++){
                        phr_marks[(sa[j]>>1UL)-1] = true;
                    }
                }

                //TODO testing
                if(p_info.p_round==4){
                    std::cout<<rank<<" -> "<<f_sa_pos<<std::endl;
                }
                //

                pre_bwt.push_back(dummy_sym, freq);
                sa[rank] = f_sa_pos;
                rank++;
            }else{
                l_sym = dict.dict[f_sa_pos-1];
                if(pre_bwt.size()>1 && pre_bwt.last_sym() == l_sym){
                    pre_bwt.inc_freq_last(freq);
                }else{
                    pre_bwt.push_back(l_sym, freq);
                }
            }
        }
    }
    assert(rank<dict.dict.size());
    pre_bwt.close();

    //TODO testing
    std::cout<<"\n "<<rank<<" "<<dict.alphabet<<std::endl;
    //

    sa.resize(rank);
    dict.phrases_has_hocc.resize(rank);
    sdsl::util::clear(d_lim_rs);
    dict.freqs.erase();
    store_to_file(ws.get_file("phr_ranks"), ranks);
    ranks.erase();
    new_phrases_ht.shrink_databuff();

    std::cout<<"\n esto usa la tabla de hash "<<new_phrases_ht.hash_table_tot_bytes()<<std::endl;

#ifdef __linux__
    malloc_trim(0);
#endif
}

size_t process_dictionary(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

    vector_t sa(dict.dict.size(), 0, sdsl::bits::hi(dict.dict.size())+2);
    uint8_t width = sdsl::bits::hi(dict.dict.size())+1;

    std::cout<<"    Sorting the dictionary using suffix induction"<<std::flush;
    auto start = std::chrono::steady_clock::now();
    if(width<=8){
        suffix_induction<uint8_t>(sa, dict);
    }else if(width<=16){
        suffix_induction<uint16_t>(sa, dict);
    }else if(width<=32){
        suffix_induction<uint32_t>(sa, dict);
    }else{
        suffix_induction<uint64_t>(sa, dict);
    }
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 4);

    std::cout<<"    Alternative solution"<<std::flush;
    start = std::chrono::steady_clock::now();
    get_pre_bwt2(dict, sa, p_info,ws);
    end = std::chrono::steady_clock::now();
    report_time(start, end, 17);

    /*std::cout<<"    Constructing the preliminary BWT"<<std::flush;
    phrase_map_t new_phrases_ht;
    bv_t phr_marks(dict.dict.size(), false);
    start = std::chrono::steady_clock::now();
    get_pre_bwt(dict, sa, p_info, phr_marks, new_phrases_ht, ws);
    end = std::chrono::steady_clock::now();
    report_time(start, end, 17);

    std::cout<<"    Storing the dictionary as a grammar"<<std::flush;
    start = std::chrono::steady_clock::now();
    dict2gram(dict, new_phrases_ht, sa, phr_marks, p_info, ws);
    end = std::chrono::steady_clock::now();
    report_time(start, end, 14);*/

    return sa.size();
}

size_t build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size, str_collection& str_coll, tmp_workspace &ws) {

    std::cout<<"Parsing the text:    "<<std::endl;
    std::string output_file = ws.get_file("tmp_output");
    std::string tmp_i_file = ws.get_file("tmp_input");

    // mark which symbols represent string boundaries
    bv_t symbol_desc(str_coll.alphabet.back()+1,false);
    symbol_desc[str_coll.alphabet[0]] = true;

    parsing_info p_info;
    for(unsigned long sym_freq : str_coll.sym_freqs){
        if(p_info.max_sym_freq<sym_freq){
            p_info.max_sym_freq = sym_freq;
        }
    }
    p_info.tot_phrases = str_coll.alphabet.back()+1;
    p_info.str_ptrs.swap(str_coll.str_ptrs);
    p_info.str_ptrs.push_back((long)str_coll.n_char);
    p_info.str_ptrs.shrink_to_fit();
    p_info.longest_str = str_coll.longest_string;

    size_t iter=1;
    size_t rem_phrases=0;

    std::cout<<"  Parsing round "<<iter++<<std::endl;
    auto start = std::chrono::steady_clock::now();
    if(n_threads>1) {
        {
            mt_byte_parse_strategy p_strat(i_file, tmp_i_file, p_info, hbuff_size, n_threads);
            rem_phrases = build_lc_gram_int<mt_byte_parse_strategy>(i_file, tmp_i_file, p_strat, p_info, symbol_desc, ws);
        }
    }else{
        {
            st_byte_parse_strategy p_strat(i_file, tmp_i_file, p_info);
            rem_phrases = build_lc_gram_int<st_byte_parse_strategy>(i_file, tmp_i_file, p_strat, p_info, symbol_desc, ws);
        }
    }

    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 4);

    while (rem_phrases > 0) {
        start = std::chrono::steady_clock::now();
        std::cout<<"  Parsing round "<<iter++<<std::endl;
        if(n_threads>1) {
            mt_int_parse_strategy p_strat( tmp_i_file, output_file, p_info, hbuff_size, n_threads);
            rem_phrases = build_lc_gram_int<mt_int_parse_strategy>(tmp_i_file, output_file, p_strat, p_info, symbol_desc, ws);
        }else{
            st_int_parse_strategy p_strat(tmp_i_file, output_file, p_info);
            rem_phrases = build_lc_gram_int<st_int_parse_strategy>(tmp_i_file, output_file, p_strat, p_info, symbol_desc, ws);
        }
        end = std::chrono::steady_clock::now();
        report_time(start, end,4);
        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
    }

    sdsl::util::clear(symbol_desc);
    ws.remove_file("suffix_file");

    return iter-2;
}

template<class parse_strategy_t>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file,
                         parse_strategy_t& p_strategy, parsing_info &p_info,
                         bv_t &phrase_desc, tmp_workspace &ws) {

#ifdef __linux__
    malloc_trim(0);
#endif

    std::cout<< "    Computing the phrases in the text" << std::flush;
    auto start = std::chrono::steady_clock::now();
    auto res = p_strategy.get_phrases();
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 16);

    store_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
    std::vector<long>().swap(p_info.str_ptrs);

#ifdef __linux__
    malloc_trim(0);
#endif

    phrase_map_t & map = p_strategy.map;

    size_t psize=0;//<- for the iter stats
    if(map.size()!=p_info.lms_phrases) {

        //size_t width = sdsl::bits::hi(phrase_desc.size())+1;
        size_t dict_sym = res.first;
        size_t max_freq = res.second;

        //save a copy of the hash table into a file
        std::string ht_file = ws.get_file("ht_data");
        map.store_data_to_file(ht_file);

        //temporal unload of the hash table (not the data)
        map.destroy_table();
        size_t tot_phrases=0;
        {
            //create a dictionary from where the ids will be computed
            std::cout<<"    Creating the dictionary from the hash table"<<std::flush;
            start = std::chrono::steady_clock::now();
            dictionary dict(map, dict_sym, max_freq, phrase_desc,
                            p_strategy.text_size, p_info.prev_alph,
                            p_info.max_sym_freq, p_info.tot_phrases);
            end = std::chrono::steady_clock::now();
            map.destroy_data();
            report_time(start, end, 6);

            //process the dictionary
            tot_phrases = process_dictionary(dict, p_info, ws);
        }

        //reload the hash table
        map.load_data_from_file(ht_file);
        ws.remove_file("ht_data");

        {
            std::cout<<"    Assigning ranks to the new phrases"<<std::flush;
            start = std::chrono::steady_clock::now();
            bv_t new_phrase_desc(tot_phrases, false);
            key_wrapper key_w{sym_width(p_info.tot_phrases), map.description_bits(), map.get_data()};

            size_t j=0;
            vector_t ranks;
            std::string ranks_file = ws.get_file("phr_ranks");
            sdsl::load_from_file(ranks, ranks_file);
            for(auto const& ptr : map){
                phrase_map_t::val_type val=0;
                map.get_value_from(ptr, val);
                val = ranks[j++];
                map.insert_value_at(ptr, val);
                //the first bit marks if the phrase is repeated or not. We need to shift it to get the real id
                new_phrase_desc[(val>>1UL)] = phrase_desc[key_w.read(ptr, 0)];

                //todo testing
                /*if(p_info.p_round>5){
                    std::cout<<" phrase : "<<(val>>1UL)<<" -> ";
                    for(size_t i=key_w.size(ptr);i-->0;){
                        std::cout<<key_w.read(ptr, i)<<" ";
                    }
                    std::cout<<""<<std::endl;
                }*/
                //
            }
            ws.remove_file("phr_ranks");
            std::string suffix_file = ws.get_file("suffix_file");
            sdsl::store_to_file(new_phrase_desc, suffix_file);
            end = std::chrono::steady_clock::now();
            report_time(start, end, 15);
        }

        std::cout<<"    Creating the parse of the text"<<std::flush;
        start = std::chrono::steady_clock::now();
        load_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
        psize = p_strategy.parse_text();
        end = std::chrono::steady_clock::now();
        report_time(start, end, 19);

        {
            //keep track of the phrases that have to be rephrased
            std::string suffix_file = ws.get_file("suffix_file");
            bv_t new_phrase_desc;
            sdsl::load_from_file(new_phrase_desc, suffix_file);
            p_info.prev_alph = phrase_desc.size();
            phrase_desc.swap(new_phrase_desc);
        }

        p_info.max_sym_freq = max_freq;
        p_info.lms_phrases = map.size();
        p_info.tot_phrases = tot_phrases;
        p_info.p_round++;

        std::cout<<"    Stats:"<<std::endl;
        std::cout<<"      Parsing phrases:                  "<<p_info.lms_phrases<<std::endl;
        std::cout<<"      Number of symbols in the phrases: "<<dict_sym<<std::endl;
        std::cout<<"      Number of BWT blocks:             "<<p_info.tot_phrases<<std::endl;
        std::cout<<"      Parse size:                       "<<psize<<std::endl;

        if(psize==0){
            return 0;
        }else{
            return p_info.lms_phrases;
        }
    } else{ //just copy the input

        std::ifstream in(i_file, std::ios_base::binary);
        std::ofstream out(o_file, std::ios_base::binary);
        auto buffer = reinterpret_cast<char*>(malloc(BUFFER_SIZE));
        do {
            in.read(&buffer[0], BUFFER_SIZE);
            out.write(&buffer[0], in.gcount());
            psize+=in.gcount();
        } while (in.gcount() > 0);
        free(buffer);
        psize/=sizeof(typename parse_strategy_t::sym_type);
        p_strategy.remove_files();
        std::cout<<"    No new phrases found"<<std::endl;
        return 0;
    }
}

template size_t build_lc_gram_int<st_byte_parse_strategy>(std::string &i_file, std::string &o_file, st_byte_parse_strategy& p_strat, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
template size_t build_lc_gram_int<st_int_parse_strategy>(std::string &i_file, std::string &o_file, st_int_parse_strategy& p_strat, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
template size_t build_lc_gram_int<mt_byte_parse_strategy>(std::string &i_file, std::string &o_file, mt_byte_parse_strategy& p_strat, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
template size_t build_lc_gram_int<mt_int_parse_strategy>(std::string &i_file, std::string &o_file, mt_int_parse_strategy& p_start, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
