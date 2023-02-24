//
// Created by Diaz, Diego on 26.1.2023.
//

#include "opt_par_phase.hpp"

#include "bwt_io.h"
#include "opt_LMS_induction.h"

//#include "malloc_count.h"

namespace opt_algo {

    void produce_grammar(tmp_workspace& ws, dictionary& dict, parsing_info& p_info) {

        i_file_stream<size_t> s_sa(ws.get_file("s_sa"), BUFFER_SIZE);
        dict.alphabet+=3;
        dict.metasym_dummy= dict.alphabet+s_sa.size()+1;

        size_t alphabet = dict.alphabet;
        size_t meta_sym_dummy = dict.metasym_dummy;

        dict.phrases_has_hocc.resize(s_sa.size()+1);
        sdsl::util::set_to_value(dict.phrases_has_hocc, false);

        buffered_hash_table<size_t> phrases_ht;
        i_file_stream<size_t> nested_phrases(ws.get_file("nested_phrases"), BUFFER_SIZE);
        size_t i=0, sym, meta_sym, offset = s_sa.size()-1;
        bool is_last;
        string_t phrase(2, sym_width(alphabet));
        while(i<nested_phrases.size()){
            meta_sym = offset - nested_phrases.read(i++);
            dict.phrases_has_hocc[meta_sym] = true;
            do{
                sym = nested_phrases.read(i++);
                is_last = sym & 1UL;
                phrase.push_back(sym>>1UL);
            }while(!is_last);

            phrase.mask_tail();
            phrases_ht.insert(phrase.data(), phrase.n_bits(), (alphabet + meta_sym));
            phrase.clear();
        }
        nested_phrases.close(true);

        size_t pos, new_size=0, l_sym, r_sym;
        vector_t new_dict((s_sa.size()+1)*2, 0, sym_width(meta_sym_dummy));
        bv_t meta_sym_marks;
        sdsl::load_from_file(meta_sym_marks, ws.get_file("meta_sym_marks"));

        for(size_t u=s_sa.size();u-->0;) {
            pos = s_sa.read(u);
            assert(pos<dict.dict.size() && !dict.d_lim[pos]);
            pos++;
            while(!meta_sym_marks[pos] && !dict.d_lim[pos]) pos++;

            l_sym = dict.dict.read(pos-1);
            assert(l_sym<alphabet);

            if(meta_sym_marks[pos]){
                phrase.clear();
                do{
                    phrase.push_back(dict.dict[pos]);
                }while(!dict.d_lim[pos++]);
                phrase.mask_tail();
                assert(phrase.size()>1);
                auto res = phrases_ht.find(phrase.data(), phrase.n_bits());
                assert(res.second);

                r_sym = 0;
                phrases_ht.get_value_from(res.first, r_sym);

                new_dict.write(new_size++, l_sym);
                new_dict.write(new_size++, r_sym);
            }else{
                r_sym = dict.dict.read(pos);
                new_dict.write(new_size++, meta_sym_dummy);
                new_dict.write(new_size++, dict.is_suffix(r_sym) ? meta_sym_dummy : l_sym);
            }
        }

        new_dict.write(new_size++, meta_sym_dummy);
        if(p_info.p_round>0){
            new_dict.write(new_size, dict.end_str_dummy);
        }else{
            new_dict.write(new_size, 10);
        }

        dict.dict.swap(new_dict);
        dict.n_phrases = s_sa.size();

        s_sa.close(true);
        ws.remove_file("meta_sym_marks");

        sdsl::util::clear(dict.d_lim);
        std::string dict_file = ws.get_file("dict_lev_"+std::to_string(p_info.p_round));
        sdsl::store_to_file(dict, dict_file);
    }

/*
template<class value_type>
size_t get_pre_bwt2(dictionary &dict, value_type * sa, size_t sa_size, parsing_info& p_info, tmp_workspace& ws) {

    //TODO testing
    produce_grammar(ws, dict, p_info);
    //

    bool is_maximal, exist_as_phrase;
    size_t u=0, d_pos, pl_sym, bg_pos, acc_freq, freq, rank=0, l_sym, f_sa_pos, d_rank, n_breaks, dummy_run=0, longest_dummy_run=0;

    //now the dummy symbols and the end-string symbol are part of the alphabet
    dict.alphabet+=3;
    bv_rs_t d_lim_rs(&dict.d_lim);

    std::string pre_bwt_file = ws.get_file("pre_bwt_lev_"+std::to_string(p_info.p_round));

    size_t sb = INT_CEIL(sym_width(std::max(dict.alphabet, dict.prev_alphabet)), 8);
    size_t fb = INT_CEIL(sym_width(dict.t_size),8);
    bwt_buff_writer pre_bwt(pre_bwt_file, std::ios::out, sb, fb);

    size_t width = sym_width(dict.dict.size())+1;
    vector_t new_metasymbols(dict.n_phrases, 0, width);
    vector_t first_symbol(size_t(double(dict.n_phrases)*1.2), sym_width(dict.alphabet));

    while(u<sa_size) {
        d_pos = (sa[u]>>2UL) - 1;
        if((sa[u] & 2UL)==0) {
            if(dict.d_lim[d_pos]) {
                u++;
                while(u<sa_size && sa[u] & 1UL){
                    u++;
                }
            }else if((sa[u] & 2UL)==0) {
                f_sa_pos = d_pos;
                bg_pos = u;
                exist_as_phrase = false;
                n_breaks = 0;
                acc_freq = 0;
                pl_sym = std::numeric_limits<size_t>::max();
                do{
                    d_pos = (sa[u]>>2UL) - 1;
                    d_rank = d_lim_rs(d_pos);
                    freq = dict.freqs.read(d_rank);
                    acc_freq += freq;

                    if(d_pos==0 || dict.d_lim[d_pos-1]){
                        l_sym = dict.bwt_dummy;
                        exist_as_phrase = true;
                        new_metasymbols.write(d_rank, (rank<<1UL) | (freq>1));
                    }else{
                        l_sym = dict.dict[d_pos-1];
                        if(l_sym>=dict.alphabet) l_sym = first_symbol[l_sym-dict.alphabet];
                        assert(l_sym<dict.end_str_dummy);
                    }
                    n_breaks += (pl_sym!=l_sym);
                    pl_sym = l_sym;
                }while(++u<sa_size && sa[u] & 1UL);
                is_maximal = n_breaks>1;

                if(is_maximal || exist_as_phrase) {

                    //std::cout<<bg_pos<<" "<<u-1<<" "<<l_sym<<" "<<((sa[bg_pos]>>2)-1)<<std::endl;
                    l_sym = dict.bwt_dummy+((u-bg_pos)>1);

                    //pre_bwt.push_back(dict.bwt_dummy+((u-bg_pos)>1), acc_freq);
                    dummy_run+=acc_freq;
                    if(dummy_run>longest_dummy_run) longest_dummy_run = dummy_run;

                    size_t f_sym = dict.dict.read(f_sa_pos);
                    assert(f_sym<dict.end_str_dummy);
                    first_symbol.push_back(f_sym);
                    dict.phrases_has_hocc[rank] = (u-bg_pos>1);

                    for(size_t j=bg_pos;j<u;j++){
                        size_t pos = (sa[j]>>2UL)-1;
                        //assert(dict.dict.read(pos)==f_sym);
                        dict.dict.write(pos, dict.alphabet + rank);
                    }

                    sa[rank++] = f_sa_pos;
                }else{
                    //size_t test = dict.dict.read(((sa[bg_pos]>>2)-1)-1);
                    //if(test>=dict.alphabet) test = first_symbol[test-dict.alphabet];
                    //std::cout<<bg_pos<<" "<<u-1<<" "<<l_sym<<" "<<((sa[bg_pos]>>2)-1)<<" "<<test<<std::endl;
                    assert(l_sym<dict.end_str_dummy);
                    dummy_run=0;
                }

                if(pre_bwt.size()>0 && pre_bwt.last_sym() == l_sym){
                    pre_bwt.inc_freq_last(acc_freq);
                }else{
                    //std::cout<<"dif "<<bg_pos<<" ["<<bg_pos<<", "<<u-1<<"]"<<std::endl;
                    pre_bwt.push_back(l_sym, acc_freq);
                }
            }
        }else {
            f_sa_pos = d_pos;
            size_t f_sym = dict.dict.read(f_sa_pos);
            assert(f_sym<dict.end_str_dummy);


            do{
                d_pos = (sa[u]>>2UL) - 1;
                d_rank = d_lim_rs(d_pos);
                freq = dict.freqs.read(d_rank);
                size_t phrase_flag = freq & 3UL;
                assert(phrase_flag!=0 && phrase_flag<=3);
                freq>>=2UL;

                if(d_pos==0 || dict.d_lim[d_pos-1]) {//a full phrase expanding to a string suffix
                    l_sym = dict.bwt_dummy;

                    if(phrase_flag==2){
                        l_sym = dict.end_str_dummy;
                    }else if(phrase_flag==3){
                        //TODO this is the case when a phrase occurs as a suffix but also as an entire string
                        exit(0);
                    }

                    dummy_run +=freq;
                    if(dummy_run>longest_dummy_run) longest_dummy_run = dummy_run;

                     // ***
                     // * We create a metasymbol for a suffix A if it is not proper.
                     // * Still, this is only because we need to access its left-context
                     // * symbols and because it serves to solve the relative order of the
                     // * phrases preceding A in the text.
                     // * Caveat: If the A expands to a full string, then we con append a dummy.
                     //
                    new_metasymbols.write(d_rank, (rank<<1UL) | (freq>1));
                    first_symbol.push_back(f_sym);

                    sa[rank++] = d_pos;
                }else{
                    //we don't compress any phrase that ends with a metasymbol
                    // expanding to a string suffix
                    l_sym = dict.dict[d_pos-1];
                    assert(l_sym<dict.end_str_dummy);
                    dummy_run=0;
                }

                // * NOTE: here it is unnecessary to check if the suffix is left-maximal as we can decide
                // * an arbitrary order for the symbols. For the same reason, we don't compress the suffix
                // * when it is not proper.
                //
                //TODO testing
                //if(p_info.p_round==6){
                //    if(l_sym==dict.sym_end_string){
                //        std::cout<<"$ | "<<freq<<" "<<(rank-1)<<std::endl;
                //    }else if(l_sym==dict.sym_dummy){
                //        std::cout<<"* | "<<freq<<" "<<(rank-1)<<std::endl;
                //    }else{
                //
                //        std::cout<<l_sym<<" | "<<freq<<std::endl;
                //}
                //}
                if(pre_bwt.size()>0 && pre_bwt.last_sym() == l_sym){
                    pre_bwt.inc_freq_last(freq);
                }else{
                    //std::cout<<"dif "<<u<<std::endl;
                    pre_bwt.push_back(l_sym, freq);
                }
            }while(++u<sa_size && (sa[u] & 1UL));
        }
    }

    //TODO some boundaries
    //std::cout<<"\nmax number of runs: "<<max_n_runs<<" max run length: "<<max_run_len<<std::endl;
    size_t tmp=sa_size*sizeof(value_type);
    tmp+=INT_CEIL(dict.dict.size()*dict.dict.width(), 8);
    tmp+=INT_CEIL(first_symbol.size()*first_symbol.width(), 8);
    tmp+=INT_CEIL(dict.freqs.size()*dict.freqs.width(), 8);
    tmp+=INT_CEIL(new_metasymbols.size()*new_metasymbols.width(), 8);
    tmp+=INT_CEIL(dict.phrases_has_hocc.size(), 8);
    tmp+=INT_CEIL(dict.d_lim.size(), 8);
    tmp+=INT_CEIL(dict.desc_bv->size(), 8);
    std::cout<<"\nspace reduction: "<<double((tmp-cont*sizeof(value_type)))/double(tmp)<<" "<<tmp<<std::endl;

    std::cout<<"unusued cells "<<cont<<" "<<sa_size<<std::endl;
    std::cout<<"bytes SA "<<sa_size*sizeof(value_type)<<std::endl;
    std::cout<<"bytes dict "<<INT_CEIL(dict.dict.size()*dict.dict.width(), 8)<<std::endl;
    std::cout<<"bytes first_symbol "<<INT_CEIL(first_symbol.size()*first_symbol.width(), 8)<<std::endl;
    std::cout<<"bytes freqs"<<INT_CEIL(dict.freqs.size()*dict.freqs.width(), 8)<<std::endl;
    std::cout<<"bytes ranks"<<INT_CEIL(new_metasymbols.size()*new_metasymbols.width(), 8)<<std::endl;
    std::cout<<"bytes has_hocc "<<INT_CEIL(dict.phrases_has_hocc.size(), 8)<<std::endl;
    std::cout<<"bytes d_lim "<<INT_CEIL(dict.d_lim.size(), 8)<<std::endl;
    std::cout<<"bytes is_suffix "<<INT_CEIL(dict.desc_bv->size(), 8)<<std::endl;

    assert(first_symbol.size()>=rank);
    assert(rank<dict.dict.size());

    //TODO test
    size_t r_sym, r_len, r_sym2, r_len2;
    bwt_buff_reader pre_bwt2(ws.get_file("pre_bwt"));
    //std::cout<<"test "<<pre_bwt.size()<<" "<<pre_bwt2.size()<<" there are "<<hits<<" hits "<<std::endl;
    assert(pre_bwt2.size()==pre_bwt.size());
    for(size_t i=pre_bwt.size();i-->0;){
        pre_bwt.read_run(i, r_sym, r_len);
        pre_bwt2.read_run(i, r_sym2, r_len2);
        //std::cout<<i<<"-> ("<<r_sym<<", "<<r_len<<") "<<" ("<<r_sym2<<", "<<r_len2<<") "<<std::endl;
        assert(r_sym==r_sym2 && r_len==r_len2);
    }
    pre_bwt2.close();
    vector_t alternative_ranks;
    load_from_file(ws.get_file("new_phrases"), alternative_ranks);
    for(size_t i=0;i<new_metasymbols.size();i++){
        assert(alternative_ranks[i]==new_metasymbols[i]);
    }
    //

    pre_bwt.close();
    sa = (value_type *)realloc(sa, rank*sizeof(value_type));
    first_symbol.resize(rank);

    //I add one to the rank because the dummy is now part of the alphabet
    dict.phrases_has_hocc.resize(rank+1);
    sdsl::util::clear(d_lim_rs);
    dict.freqs.erase();
    store_to_file(ws.get_file("phr_ranks"), new_metasymbols);
    new_metasymbols.erase();

    dict.max_sym_freq = std::max<size_t>(dict.max_sym_freq, longest_dummy_run);
#ifdef __linux__
    malloc_trim(0);
#endif


    dict2gram2<value_type>(dict, sa, rank, first_symbol, p_info, ws);

    return rank;
}
*/

    void invert_data(tmp_workspace& ws, size_t n_meta_syms,  parsing_info& p_info) {

        //invert the preliminary BWT
        bwt_buff_reader pre_bwt_inv(ws.get_file("inv_pre_bwt"));

        std::string pre_bwt_file = ws.get_file("pre_bwt_lev_"+std::to_string(p_info.p_round));
        bwt_buff_writer pre_bwt(pre_bwt_file,
                                std::ios::out,
                                pre_bwt_inv.bytes_per_rsym(),
                                pre_bwt_inv.bytes_per_rlen());
        size_t sym, len, bwt_len=0;
        for(size_t i=pre_bwt_inv.size();i-->0;){
            pre_bwt_inv.read_run(i, sym, len);
            bwt_len+=len;
            pre_bwt.push_back(sym, len);
        }

        pre_bwt_inv.close(true);
        pre_bwt.close();

        //invert the metasymbols for the new phrases
        vector_t meta_sym_list;
        load_from_file(ws.get_file("phr_ranks"), meta_sym_list);
        size_t meta_sym, is_rep;
        n_meta_syms--;
        for(size_t i=0;i<meta_sym_list.size();i++){
            meta_sym = meta_sym_list.read(i);
            is_rep = meta_sym & 1UL;
            meta_sym>>=1UL;
            meta_sym = ((n_meta_syms-meta_sym) <<1UL) | is_rep;
            meta_sym_list.write(i, meta_sym);
        }
        store_to_file(ws.get_file("phr_ranks"), meta_sym_list);

        //invert the bit vector marking the suffix ranges
        std::vector<std::pair<size_t, size_t>> same_suffix_ranges;
        load_pl_vector(ws.get_file("same_suffix_ranges"), same_suffix_ranges);

        for(auto & range : same_suffix_ranges){
            range.first = bwt_len - range.first;
            range.second = bwt_len - (range.second+1);
            //std::cout<<range.first<<" "<<range.second<<" "<<bwt_len<<std::endl;
            //std::cout<<bwt_len-range.first<<" "<<bwt_len-(range.second+1)<<"\n"<<std::endl;
        }
        std::reverse(same_suffix_ranges.begin(), same_suffix_ranges.end());
        //TODO just testing
        for(auto & range : same_suffix_ranges){
            std::cout<<range.first<<" "<<range.second<<" "<<bwt_len<<std::endl;
        }
        //
        store_pl_vector(ws.get_file("same_suffix_"+std::to_string(p_info.p_round)), same_suffix_ranges);
        /*for(auto const& range : same_suffix_ranges){
            std::cout<<range.first<<" "<<range.second<<std::endl;
        }*/
    }

    template<class vector_type, class sa_type>
    size_t process_dictionary_int(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

        std::cout<<"    Sorting the dictionary and constructing the preliminary BWT"<<std::flush;
        auto start = std::chrono::steady_clock::now();
        size_t n_phrases = suffix_induction<vector_type, sa_type>(dict, ws);
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 2);

        //malloc_count_print_status();
        //malloc_count_reset_peak();

        std::cout<<"    Finishing the preliminary BWT"<<std::flush;
        start = std::chrono::steady_clock::now();
        invert_data(ws, n_phrases, p_info);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 32);

        std::cout<<"    Compressing the dictionary"<<std::flush;
        start = std::chrono::steady_clock::now();
        //size_t n_phrases = get_pre_bwt2<value_type>(dict, sa, dict.dict.size(), p_info, ws);
        produce_grammar(ws, dict, p_info);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 35);

#ifdef __linux__
        malloc_trim(0);
#endif
        return n_phrases;
    }

    size_t process_dictionary(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

        uint8_t width = sym_width(dict.dict.size())+2;
        size_t n_phrases;

        if(width<=8){
            using uint8_vector_t = std::vector<uint8_t, mallocator<uint8_t>>;
            n_phrases = process_dictionary_int<uint8_vector_t, uint8_vector_t>(dict, p_info, ws);
        }else if(width<=16){
            using uint16_vector_t = std::vector<uint16_t, mallocator<uint16_t>>;
            n_phrases = process_dictionary_int<uint16_vector_t, uint16_vector_t>(dict, p_info, ws);
        }else if(width<=32){
            using uint32_vector_t = std::vector<uint32_t, mallocator<uint32_t>>;
            n_phrases = process_dictionary_int<uint32_vector_t, uint32_vector_t>(dict, p_info, ws);
        }else{
            using uint64_vector_t = std::vector<uint64_t, mallocator<uint64_t>>;
            n_phrases = process_dictionary_int<uint64_vector_t, vector_t>(dict, p_info, ws);
        }
        return n_phrases;
    }

    template<class par_type>
    size_t par_phase_int(std::string& i_file, std::string& o_file, parsing_info& p_info,
                         size_t hbuff_size, size_t n_threads, bv_t& symbol_desc, tmp_workspace& ws){
        size_t alph_size;
        if (n_threads > 1) {
            {
                using par_strat_t = mt_parse_strat_t<par_type, opt_hash_functor, opt_parse_functor>;
                par_strat_t p_strat(i_file, o_file, p_info, hbuff_size, n_threads);
                alph_size = opt_algo::par_round<par_strat_t>(p_strat, p_info, symbol_desc, ws);
            }
        } else {
            {
                using par_strat_t = st_parse_strat_t<par_type, opt_hash_functor, opt_parse_functor>;
                par_strat_t p_strat(i_file, o_file, p_info);
                alph_size = opt_algo::par_round<par_strat_t>(p_strat, p_info, symbol_desc, ws);
            }
        }
        return alph_size;
    }

    size_t par_phase(std::string &i_file, size_t n_threads, size_t hbuff_size, str_collection& str_coll, tmp_workspace &ws) {

        std::cout<<"Parsing the text:    "<<std::endl;

        if(n_threads>1){
            std::cout<<"  Running with up to "<<n_threads<<" working threads "<<std::endl;
            std::cout<<"  Using "<<hbuff_size<<" bytes for the thread hash tables ("<< hbuff_size/n_threads<<" bytes each)"<<std::endl;
        }

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
        p_info.active_strings = str_coll.n_strings;

        size_t iter=1;
        size_t n_syms;

        std::cout<<"  Parsing round "<<iter++<<std::endl;
        auto start = std::chrono::steady_clock::now();
        n_syms = par_phase_int<uint8t_parser_t>(i_file, tmp_i_file, p_info,
                                                hbuff_size, n_threads, symbol_desc, ws);
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 4);

#ifdef __linux__
        empty_page_cache(i_file)
#endif
        //malloc_count_print_status();
        //malloc_count_reset_peak();

        while (n_syms > 0) {
            start = std::chrono::steady_clock::now();
            std::cout<<"  Parsing round "<<iter++<<std::endl;
            size_t bps = sym_width(n_syms)+1;
            if(bps<=8){
                n_syms = par_phase_int<uint8t_parser_t>(tmp_i_file, output_file, p_info,
                                                        hbuff_size, n_threads, symbol_desc, ws);
            }else if(bps<=16){
                n_syms = par_phase_int<uint16t_parser_t>(tmp_i_file, output_file, p_info,
                                                         hbuff_size, n_threads, symbol_desc, ws);
            } else if(bps<=32){
                n_syms = par_phase_int<uint32t_parser_t>(tmp_i_file, output_file, p_info,
                                                         hbuff_size, n_threads, symbol_desc, ws);
            } else{
                n_syms = par_phase_int<uint64t_parser_t>(tmp_i_file, output_file, p_info,
                                                         hbuff_size, n_threads, symbol_desc, ws);
            }
            end = std::chrono::steady_clock::now();
            report_time(start, end,4);
            remove(tmp_i_file.c_str());
            rename(output_file.c_str(), tmp_i_file.c_str());
            //malloc_count_print_status();
            //malloc_count_reset_peak();
        }

        sdsl::util::clear(symbol_desc);
        ws.remove_file("suffix_file");

        return iter-2;
    }

    template<class parse_strategy_t>
    size_t par_round(parse_strategy_t& p_strategy, parsing_info &p_info, bv_t &phrase_desc, tmp_workspace &ws) {

#ifdef __linux__
        malloc_trim(0);
#endif

        std::cout<<"    Hashing the phrases"<<std::flush;
        auto start = std::chrono::steady_clock::now();
        auto res = p_strategy.get_phrases();
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 42);

        store_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
        std::vector<long>().swap(p_info.str_ptrs);

#ifdef __linux__
        malloc_trim(0);
#endif

        phrase_map_t & map = p_strategy.map;
        size_t psize;//<- for the iter stats
        assert(map.size()>0);

        size_t dict_sym = res.first;
        size_t max_freq = res.second;

        //save a copy of the hash table into a file
        std::string ht_file = ws.get_file("ht_data");
        map.store_data_to_file(ht_file);

        //temporal unload of the hash table (not the data)
        map.destroy_table();
        size_t tot_phrases;
        {
            //create a dictionary from where the ids will be computed
            std::cout<<"    Creating the dictionary from the hash table"<<std::flush;
            start = std::chrono::steady_clock::now();
            dictionary dict(map, dict_sym, max_freq, phrase_desc,
                            p_strategy.text_size, p_info.prev_alph,
                            std::max<size_t>(p_info.max_sym_freq, p_info.active_strings));
            end = std::chrono::steady_clock::now();
            map.destroy_data();
            report_time(start, end, 18);

            //process the dictionary
            tot_phrases = process_dictionary(dict, p_info, ws);
            p_info.prev_alph = dict.alphabet;
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
            load_from_file(ranks_file, ranks);
            for(auto const& ptr : map){
                phrase_map_t::val_type val=0;
                map.get_value_from(ptr, val);
                val = ranks[j++];
                map.insert_value_at(ptr, val);
                //the first bit marks if the phrase is repeated or not.
                // We need to shift it to get the real id
                new_phrase_desc[(val>>1UL)] = phrase_desc[key_w.read(ptr, 0)];
            }
            ws.remove_file("phr_ranks");
            std::string suffix_file = ws.get_file("suffix_file");
            sdsl::store_to_file(new_phrase_desc, suffix_file);
            end = std::chrono::steady_clock::now();
            report_time(start, end, 27);
        }

        std::cout<<"    Creating the parse of the text"<<std::flush;
        start = std::chrono::steady_clock::now();
        load_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);

        size_t bps = sym_width(tot_phrases)+1;
        if(bps<=8){
            psize = p_strategy.template parse_text<uint8_t>();
        }else if(bps<=16){
            psize = p_strategy.template parse_text<uint16_t>();
        }else if(bps<=32){
            psize = p_strategy.template parse_text<uint32_t>();
        }else{
            psize = p_strategy.template parse_text<uint64_t>();
        }

        end = std::chrono::steady_clock::now();
        report_time(start, end, 31);

        {
            //keep track of the phrases that have to be rephrased
            std::string suffix_file = ws.get_file("suffix_file");
            bv_t new_phrase_desc;
            sdsl::load_from_file(new_phrase_desc, suffix_file);
            phrase_desc.swap(new_phrase_desc);
        }

        p_info.max_sym_freq = max_freq;
        p_info.lms_phrases = map.size();
        p_info.tot_phrases = tot_phrases;
        p_info.p_round++;

        std::cout<<"    Stats:"<<std::endl;
        std::cout<<"      Parsing phrases:                  "<<p_info.lms_phrases<<std::endl;
        std::cout<<"      Number of symbols in the phrases: "<<dict_sym<<std::endl;
        std::cout<<"      Number of unsolved BWT blocks:    "<<p_info.tot_phrases<<std::endl;
        std::cout<<"      Processed strings:                "<<p_info.active_strings<<std::endl;
        std::cout<<"      Parse size:                       "<<psize<<std::endl;

        map.destroy_data();
        map.destroy_table();

        if(psize==0){
            return 0;
        }else{
            return p_info.tot_phrases;
        }
        //}

        /*else{ //just copy the input

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
        }*/
    }

    //template size_t par_round<opt_st_byte_parse_strategy>(opt_st_byte_parse_strategy& p_strat, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
    //template size_t par_round<opt_st_int_parse_strategy>(opt_st_int_parse_strategy& p_strat, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
    //template size_t par_round<opt_mt_byte_parse_strategy>(opt_mt_byte_parse_strategy& p_strat, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
    //template size_t par_round<opt_mt_int_parse_strategy>(opt_mt_int_parse_strategy& p_start, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
}