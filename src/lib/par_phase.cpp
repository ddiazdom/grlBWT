//
// Created by Diaz, Diego on 26.1.2023.
//

#include "par_phase.hpp"
#include "LMS_sorting.h"
#include "bwt_io.h"

#ifdef USE_MALLOC_COUNT
#include "malloc_count.h"
#endif

template<class sa_type>
void produce_grammar(dictionary& dict, sa_type& s_sa,  phrase_map_t& new_phrases_ht, bv_t& phr_marks, parsing_info& p_info, tmp_workspace& ws) {

    string_t phrase(2, sym_width(dict.alphabet));

    dict.alphabet+=3;
    dict.metasym_dummy= dict.alphabet+s_sa.size()+1;

    size_t alphabet = dict.alphabet;
    size_t meta_sym_dummy = dict.metasym_dummy;

    size_t pos, new_size=0, l_sym, r_sym;
    vector_t new_dict((s_sa.size()+1)*2, 0, sym_width(meta_sym_dummy));

    //TODO testing
    //std::cout<<"\ndist stats"<<std::endl;
    //new_phrases_ht.ht_stats(10);
    //

    for(size_t u=0;u<s_sa.size();u++) {

        pos = s_sa[u];
        assert(pos<dict.dict.size() && (!dict.d_lim[pos] || dict.is_suffix(dict.dict[pos])));

        if(dict.d_lim[pos]){
            assert(dict.is_suffix(dict.dict[pos]));
            new_dict.write(new_size++, meta_sym_dummy);
            new_dict.write(new_size++, dict.dict[pos]);
        }else{
            pos++;
            while(!phr_marks[pos] && !dict.d_lim[pos]) pos++;

            l_sym = dict.dict.read(pos-1);
            assert(l_sym<alphabet);

            if(phr_marks[pos]){

                //TODO
                /*if(p_info.p_round==11){
                    std::cout<<"\n"<<pos<<std::endl;
                }*/
                //

                phrase.clear();
                do{
                    phrase.push_back(dict.dict[pos]);
                }while(!dict.d_lim[pos++]);
                phrase.mask_tail();

                //TODO
                /*if(p_info.p_round==11){
                    for(size_t i=0;i<phrase.size();i++){
                        std::cout<<phrase[i]<<" ";
                    }
                    std::cout<<""<<std::endl;
                }*/
                //

                auto res = new_phrases_ht.find(phrase.data(), phrase.n_bits());
                assert(res.second);

                r_sym = 0;
                new_phrases_ht.get_value_from(res.first, r_sym);
                r_sym+=dict.alphabet;

                new_dict.write(new_size++, l_sym);
                new_dict.write(new_size++, r_sym);
            }else{
                r_sym = dict.dict.read(pos);
                new_dict.write(new_size++, meta_sym_dummy);
                new_dict.write(new_size++, dict.is_suffix(r_sym) ? r_sym : l_sym);
            }
        }
    }

    dict.dict.swap(new_dict);
    dict.n_phrases = s_sa.size();

    sdsl::util::clear(dict.d_lim);
    std::string dict_file = ws.get_file("dict_lev_"+std::to_string(p_info.p_round));
    sdsl::store_to_file(dict, dict_file);
}

template<class vector_type, class sa_type>
size_t process_dictionary_int(dictionary &dict, parsing_info &p_info, tmp_workspace &ws, logger& log) {

    sa_type sa;
    if constexpr (std::is_same<sa_type, vector_t>::value) {
        size_t width = sym_width(dict.dict.size())+2;
        sa.set_width(width);
        sa.resize(dict.dict.size());
        sa.initialize(0, sa.size());
    }else{
        sa.resize(dict.dict.size());
        memset(sa.data(), 0, sa.size()*sizeof(typename sa_type::value_type));
    }

    log.info_start("Sorting the dictionary and constructing the preliminary BWT  ");
    auto start = std::chrono::steady_clock::now();
    suffix_induction<vector_type, sa_type>(dict, sa);

    phrase_map_t new_phrases_ht;
    bv_t phr_marks(dict.dict.size(), false);
    size_t n_phrases = produce_pre_bwt<sa_type>(dict, sa, new_phrases_ht, phr_marks, p_info, ws);
    auto end = std::chrono::steady_clock::now();

    std::string e_time = time2str(start, end);
    //report_time(start, end, 2);
    log.info_end("elap. time: "+e_time);

    //malloc_count_print_status();
    //malloc_count_reset_peak();

    log.info_start("Compressing the dictionary"+std::string(35, ' '));
    start = std::chrono::steady_clock::now();
    produce_grammar<sa_type>(dict, sa, new_phrases_ht, phr_marks, p_info, ws);
    end = std::chrono::steady_clock::now();
    e_time = time2str(start, end);
    log.info_end("elap. time: "+e_time);
    //report_time(start, end, 35);

#ifdef __linux__
    malloc_trim(0);
#endif
    return n_phrases;
}

template<class sa_type>
size_t
produce_pre_bwt(dictionary &dict, sa_type &sa, phrase_map_t &new_phrases_ht, bv_t &phr_marks, parsing_info &p_info,
                tmp_workspace &ws) {


    string_t phrase(2, sym_width(dict.alphabet));

    size_t n_breaks=0, n_phrases=0, acc_freq=0, null_sym=std::numeric_limits<size_t>::max();
    size_t pos, pl_sym=null_sym, bg_range=0, rank=0, l_sym;
    bool is_full_phrase, is_valid_suffix;

    bv_rs_t d_lim_rs(&dict.d_lim);

    std::string pre_bwt_file = ws.get_file("pre_bwt_lev_"+std::to_string(p_info.p_round));
    size_t sb = INT_CEIL(sym_width(dict.hocc_dummy), 8);
    size_t fb = INT_CEIL(sym_width(dict.t_size),8);
    bwt_buff_writer pre_bwt(pre_bwt_file, std::ios::out, sb, fb);

    vector_t ranks(dict.n_phrases, 0, sym_width(dict.dict.size())+1);
    dict.phrases_has_hocc.resize(dict.dict.size());
    sdsl::util::set_to_value(dict.phrases_has_hocc, false);

    size_t u = 0 ;
    while(u<sa.size()) {

        pos = (sa[u]>>1UL) - 1;
        is_valid_suffix = !dict.d_lim[pos] || dict.is_suffix(dict.dict[pos]);

        if(!is_valid_suffix) {
            u++;
            while(sa[u] & 1UL) u++;
            bg_range = u;
            continue;
        }

        is_full_phrase = (pos == 0 || dict.d_lim[pos - 1]);
        size_t freq = dict.freqs[d_lim_rs(pos)];
        if(is_full_phrase){
            ranks[d_lim_rs(pos)] = rank<<1UL | (freq>1);
            l_sym = dict.bwt_dummy;
        }else{
            l_sym =  dict.dict[pos-1];
        }

        n_breaks +=pl_sym!=l_sym;
        n_phrases += is_full_phrase;
        acc_freq += freq;

        if((u+1)==sa.size() || !(sa[u+1] & 1UL)) {

            if(n_breaks>1 || n_phrases==1) {

                l_sym = dict.bwt_dummy;
                if((u-bg_range+1) > 1) {

                    size_t samp_pos = pos;
                    phrase.clear();
                    do{
                        phrase.push_back(dict.dict[samp_pos]);
                    } while(!dict.d_lim[samp_pos++]);

                    phrase.mask_tail();
                    auto res = new_phrases_ht.insert(phrase.data(), phrase.n_bits(), rank);
                    assert(res.second);
                    dict.phrases_has_hocc[rank] = true;

                    for (size_t j = bg_range; j <= u; j++) {
                        phr_marks[(sa[j] >> 1UL) - 1] = true;
                    }
                    l_sym = dict.hocc_dummy;
                }
                sa[rank] = pos;
                rank++;
            }

            if (pre_bwt.size() > 1 && pre_bwt.last_sym() == l_sym) {
                pre_bwt.inc_freq_last(acc_freq);
            } else {
                pre_bwt.push_back(l_sym, acc_freq);
            }

            bg_range = u+1;
            l_sym = null_sym;
            n_breaks = 0;
            n_phrases = 0;
            acc_freq = 0;
        }
        pl_sym = l_sym;
        u++;
    }

    pre_bwt.close();
    sa.resize(rank);
    dict.phrases_has_hocc.resize(rank);

    sdsl::util::clear(d_lim_rs);
    dict.freqs.erase();
    store_to_file(ws.get_file("phr_ranks"), ranks);
    ranks.erase();
    //new_phrases_ht.shrink_databuff();

#ifdef __linux__
    malloc_trim(0);
#endif
    return rank;
}

size_t process_dictionary(dictionary &dict, parsing_info &p_info, tmp_workspace &ws, logger& log) {

    uint8_t width = sym_width(dict.dict.size()) + 1;
    size_t n_phrases;

    if (width <= 8) {
        using uint8_vector_t = std::vector<uint8_t, mallocator<uint8_t>>;
        n_phrases = process_dictionary_int<uint8_vector_t, uint8_vector_t>(dict, p_info, ws, log);
    } else if (width <= 16) {
        using uint16_vector_t = std::vector<uint16_t, mallocator<uint16_t>>;
        n_phrases = process_dictionary_int<uint16_vector_t, uint16_vector_t>(dict, p_info, ws, log);
    } else if (width <= 32) {
        using uint32_vector_t = std::vector<uint32_t, mallocator<uint32_t>>;
        n_phrases = process_dictionary_int<uint32_vector_t, uint32_vector_t>(dict, p_info, ws, log);
    } else {
        using uint64_vector_t = std::vector<uint64_t, mallocator<uint64_t>>;
        n_phrases = process_dictionary_int<uint64_vector_t, vector_t>(dict, p_info, ws, log);
    }
    return n_phrases;
}

template<class par_type>
size_t par_phase_int(std::string& i_file, std::string& o_file, parsing_info& p_info,
                     size_t hbuff_size, size_t n_threads, bv_t& symbol_desc, tmp_workspace& ws, logger& log){
    size_t alph_size;
    if (n_threads > 1) {
        using par_strat_t = mt_parse_strat_t<par_type, ext_hash_functor, ext_parse_functor>;
        par_strat_t p_strat(i_file, o_file, p_info, hbuff_size, n_threads);
        alph_size = par_round<par_strat_t>(p_strat, p_info, symbol_desc, ws, log);
    } else {
        using par_strat_t = st_parse_strat_t<par_type, ext_hash_functor, ext_parse_functor>;
        par_strat_t p_strat(i_file, o_file, p_info);
        alph_size = par_round<par_strat_t>(p_strat, p_info, symbol_desc, ws, log);
    }
    return alph_size;
}

template<class sym_type>
size_t par_phase(std::string &i_file, size_t n_threads, tmp_workspace &ws, logger& log) {

    log.info("Reading the file");
    str_collection str_coll = collection_stats<sym_type>(i_file);
    log.inc_pad();
    std::string stats = str_coll.to_string(log.get_padding());
    log.info(stats);
    log.dec_pad();

    float hbuff_frac = 0.15;
    auto hbuff_size = std::max<size_t>(64 * n_threads, size_t(ceil(float(str_coll.n_syms) * hbuff_frac)));

    log.info("Parsing the text:");
    if(n_threads>1){
        log.info("Running with up to "+std::to_string(n_threads)+" working threads ");
        log.info("Using "+std::to_string(float(hbuff_size)/1000000)+" megabytes for the thread hash tables ("+std::to_string((float(hbuff_size)/1000000)/float(n_threads))+" megabytes each)");
    }

    std::string output_file = ws.get_file("tmp_output");
    std::string tmp_i_file = ws.get_file("tmp_input");

    // mark which symbols represent string boundaries
    bv_t symbol_desc(str_coll.max_sym + 1, false);
    symbol_desc[str_coll.min_sym] = true;

    parsing_info p_info;
    p_info.max_sym_freq = str_coll.max_sym_freq;
    p_info.tot_phrases = str_coll.max_sym + 1;
    p_info.str_ptrs.swap(str_coll.str_ptrs);
    p_info.str_ptrs.push_back((long) str_coll.n_syms);
    p_info.str_ptrs.shrink_to_fit();
    p_info.longest_str = str_coll.longest_string;
    p_info.active_strings = str_coll.n_strings;

    size_t iter = 1;
    size_t n_syms;
    using f_parser_t = lms_parsing<i_file_stream<sym_type>, string_t, true>;

    log.inc_pad();
    log.info("Parsing round "+std::to_string(iter++));
    auto start = std::chrono::steady_clock::now();
    n_syms = par_phase_int<f_parser_t>(i_file, tmp_i_file, p_info, hbuff_size, n_threads, symbol_desc, ws, log);
    auto end = std::chrono::steady_clock::now();

    std::string e_time = time2str(start, end);
    //report_time(start, end, 4);
    log.info("Elap. time: "+e_time);

#ifdef USE_MALLOC_COUNT
    malloc_count_print_status();
    malloc_count_reset_peak();
#endif

    while (n_syms > 0) {
        start = std::chrono::steady_clock::now();
        log.info("Parsing round "+std::to_string(iter++));
        size_t bps = sym_width(n_syms)+1;
        if(bps<=8){
            n_syms = par_phase_int<uint8t_parser_t>(tmp_i_file, output_file, p_info,
                                                    hbuff_size, n_threads, symbol_desc, ws, log);
        }else if(bps<=16){
            n_syms = par_phase_int<uint16t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, symbol_desc, ws, log);
        } else if(bps<=32){
            n_syms = par_phase_int<uint32t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, symbol_desc, ws, log);
        } else{
            n_syms = par_phase_int<uint64t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, symbol_desc, ws, log);
        }
        end = std::chrono::steady_clock::now();
        e_time = time2str(start, end);
        //report_time(start, end, 4);
        log.info("Elap. time: "+e_time);

        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
#ifdef USE_MALLOC_COUNT
        malloc_count_print_status();
        malloc_count_reset_peak();
#endif
    }

    sdsl::util::clear(symbol_desc);
    ws.remove_file("suffix_file");
    log.dec_pad();

    return iter - 2;
}

template<class parse_strategy_t>
size_t par_round(parse_strategy_t &p_strategy, parsing_info &p_info, bv_t &phrase_desc, tmp_workspace &ws, logger& log) {

#ifdef __linux__
    malloc_trim(0);
#endif
    log.inc_pad();
    log.info_start("Computing the dictionary of LMS phrases"+std::string(22,' '));
    auto start = std::chrono::steady_clock::now();
    auto res = p_strategy.get_phrases();
    auto end = std::chrono::steady_clock::now();
    std::string e_time = time2str(start, end);
    log.info_end("elap. time: "+e_time);
    //report_time(start, end, 22);

    store_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
    std::vector<long>().swap(p_info.str_ptrs);

#ifdef __linux__
    malloc_trim(0);
#endif

    phrase_map_t &map = p_strategy.map;
    size_t psize;//<- for the iter stats
    assert(map.size() > 0);

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
        log.info_start("Compacting the dictionary"+std::string(36,' '));
        start = std::chrono::steady_clock::now();
        dictionary dict(map, dict_sym, max_freq, phrase_desc,
                        p_strategy.text_size, p_info.prev_alph, p_info.max_sym_freq);
        end = std::chrono::steady_clock::now();
        //report_time(start, end, 36);
        e_time = time2str(start, end);
        log.info_end("elap. time: "+e_time);

        map.destroy_data();
        //process the dictionary
        tot_phrases = process_dictionary(dict, p_info, ws, log);
        p_info.prev_alph = dict.alphabet;
    }

    //reload the hash table
    map.load_data_from_file(ht_file);
    ws.remove_file("ht_data");

    {
        log.info_start("Assigning metasymbols to the LMS phrases"+std::string(21,' '));
        start = std::chrono::steady_clock::now();
        bv_t new_phrase_desc(tot_phrases, false);
        key_wrapper key_w{sym_width(p_info.tot_phrases), map.description_bits(), map.get_data()};

        size_t j = 0;
        vector_t ranks;
        std::string ranks_file = ws.get_file("phr_ranks");
        load_from_file(ranks_file, ranks);
        for (auto const &ptr: map) {
            size_t val = ranks[j++];
            map.insert_value_at(ptr, val);

            //the first bit marks if the phrase is repeated or not.
            // We need to shift it to get the real id
            new_phrase_desc[(val >> 1UL)] = phrase_desc[key_w.read(ptr, 0)];
        }
        ws.remove_file("phr_ranks");
        std::string suffix_file = ws.get_file("suffix_file");
        sdsl::store_to_file(new_phrase_desc, suffix_file);
        end = std::chrono::steady_clock::now();

        e_time = time2str(start, end);
        log.info_end("elap. time: "+e_time);
        //report_time(start, end, 21);
    }

    log.info_start("Creating the parse of the text"+std::string(31, ' '));
    start = std::chrono::steady_clock::now();

    load_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
    size_t bps = sym_width(tot_phrases)+1;
    if(bps<=8){
        psize = p_strategy.template parse_text<uint8_t>();
    }else if(bps<=16){
        psize = p_strategy.template parse_text<uint16_t>();
    } else if(bps<=32){
        psize = p_strategy.template parse_text<uint32_t>();
    }else{
        psize = p_strategy.template parse_text<uint64_t>();
    }
    assert(psize>=map.size());//the parse can't be smaller than the number of phrases

    end = std::chrono::steady_clock::now();
    //report_time(start, end, 31);
    e_time = time2str(start, end);
    log.info_end("elap. time: "+e_time);

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

    log.info("Stats:");
    log.info("  Parsing phrases:                  "+std::to_string(p_info.lms_phrases));
    log.info("  Number of symbols in the phrases: "+std::to_string(dict_sym));
    log.info("  Number of unsolved BWT blocks:    "+std::to_string(p_info.tot_phrases));
    log.info("  Parse size:                       "+std::to_string(psize));

    map.destroy_data();
    map.destroy_table();

#ifdef __linux__
    malloc_trim(0);
#endif

    log.dec_pad();
    return (p_info.str_ptrs.size()-1) == psize ? 0 : p_info.tot_phrases;
}

template unsigned long par_phase<uint8_t>(std::string &i_file, size_t n_threads,  tmp_workspace &ws, logger& log);
template unsigned long par_phase<uint16_t>(std::string &i_file, size_t n_threads,  tmp_workspace &ws, logger& log);
template unsigned long par_phase<uint32_t>(std::string &i_file, size_t n_threads, tmp_workspace &ws, logger& log);
template unsigned long par_phase<uint64_t>(std::string &i_file, size_t n_threads, tmp_workspace &ws, logger& log);