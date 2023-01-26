//
// Created by Diaz, Diego on 26.1.2023.
//

#include "exact_par_phase.hpp"
#include "exact_LMS_induction.h"

#include "malloc_count.h"

namespace exact_algo {

    void produce_grammar(tmp_workspace& ws, dictionary& dict, parsing_info& p_info) {

        i_file_stream<size_t> s_sa(ws.get_file("s_sa"), BUFFER_SIZE);
        dict.alphabet+=3;
        dict.metasym_dummy= dict.alphabet+s_sa.size()+1;

        size_t alphabet = dict.alphabet;
        size_t meta_sym_dummy = dict.metasym_dummy;

        dict.phrases_has_hocc.resize(s_sa.size()+1);
        sdsl::util::set_to_value(dict.phrases_has_hocc, false);

        bit_hash_table<size_t> phrases_ht;
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

    template<class vector_type, class sa_type>
    size_t process_dictionary_int(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

        std::cout << "    Sorting the dictionary and constructing the preliminary BWT" << std::flush;
        auto start = std::chrono::steady_clock::now();
        size_t n_phrases = suffix_induction<vector_type, sa_type>(dict, ws);
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 2);

        malloc_count_print_status();
        malloc_count_reset_peak();

        std::cout << "    Compressing the dictionary" << std::flush;
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

        uint8_t width = sym_width(dict.dict.size()) + 2;
        size_t n_phrases;

        if (width <= 8) {
            using uint8_vector_t = std::vector<uint8_t, mallocator<uint8_t>>;
            n_phrases = process_dictionary_int<uint8_vector_t, uint8_vector_t>(dict, p_info, ws);
        } else if (width <= 16) {
            using uint16_vector_t = std::vector<uint16_t, mallocator<uint16_t>>;
            n_phrases = process_dictionary_int<uint16_vector_t, uint16_vector_t>(dict, p_info, ws);
        } else if (width <= 32) {
            using uint32_vector_t = std::vector<uint32_t, mallocator<uint32_t>>;
            n_phrases = process_dictionary_int<uint32_vector_t, uint32_vector_t>(dict, p_info, ws);
        } else {
            using uint64_vector_t = std::vector<uint64_t, mallocator<uint64_t>>;
            n_phrases = process_dictionary_int<uint64_vector_t, vector_t>(dict, p_info, ws);
        }
        return n_phrases;
    }

    size_t par_phase(std::string &i_file, size_t n_threads, size_t hbuff_size, str_collection &str_coll,
                           tmp_workspace &ws) {

        std::cout << "Parsing the text:    " << std::endl;
        std::string output_file = ws.get_file("tmp_output");
        std::string tmp_i_file = ws.get_file("tmp_input");

        // mark which symbols represent string boundaries
        bv_t symbol_desc(str_coll.alphabet.back() + 1, false);
        symbol_desc[str_coll.alphabet[0]] = true;

        parsing_info p_info;
        for (unsigned long sym_freq: str_coll.sym_freqs) {
            if (p_info.max_sym_freq < sym_freq) {
                p_info.max_sym_freq = sym_freq;
            }
        }
        p_info.tot_phrases = str_coll.alphabet.back() + 1;
        p_info.str_ptrs.swap(str_coll.str_ptrs);
        p_info.str_ptrs.push_back((long) str_coll.n_char);
        p_info.str_ptrs.shrink_to_fit();
        p_info.longest_str = str_coll.longest_string;
        p_info.active_strings = str_coll.n_strings;

        size_t iter = 1;
        size_t rem_phrases;

        std::cout << "  Parsing round " << iter++ << std::endl;
        auto start = std::chrono::steady_clock::now();
        if (n_threads > 1) {
            {
                ext_mt_byte_parse_strategy p_strat(i_file, tmp_i_file, p_info, hbuff_size, n_threads);
                rem_phrases = exact_algo::par_round<ext_mt_byte_parse_strategy>(p_strat, p_info, symbol_desc, ws);
            }
        } else {
            {
                ext_st_byte_parse_strategy p_strat(i_file, tmp_i_file, p_info, symbol_desc);
                rem_phrases = exact_algo::par_round<ext_st_byte_parse_strategy>(p_strat, p_info, symbol_desc, ws);
            }
        }
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 4);
        malloc_count_print_status();
        malloc_count_reset_peak();

        while (rem_phrases > 0) {
            start = std::chrono::steady_clock::now();
            std::cout << "  Parsing round " << iter++ << std::endl;
            if (n_threads > 1) {
                ext_mt_int_parse_strategy p_strat(tmp_i_file, output_file, p_info, hbuff_size, n_threads);
                rem_phrases = exact_algo::par_round<ext_mt_int_parse_strategy>(p_strat, p_info, symbol_desc, ws);
            } else {
                ext_st_int_parse_strategy p_strat(tmp_i_file, output_file, p_info, symbol_desc);
                rem_phrases = exact_algo::par_round<ext_st_int_parse_strategy>(p_strat, p_info, symbol_desc, ws);
            }
            end = std::chrono::steady_clock::now();
            report_time(start, end, 4);
            remove(tmp_i_file.c_str());
            rename(output_file.c_str(), tmp_i_file.c_str());
            malloc_count_print_status();
            malloc_count_reset_peak();
        }

        sdsl::util::clear(symbol_desc);
        ws.remove_file("suffix_file");

        return iter - 2;
    }

    template<class parse_strategy_t>
    size_t par_round(parse_strategy_t &p_strategy, parsing_info &p_info, bv_t &phrase_desc, tmp_workspace &ws) {

#ifdef __linux__
        malloc_trim(0);
#endif

        std::cout << "    Hashing the phrases" << std::flush;
        auto start = std::chrono::steady_clock::now();
        auto res = p_strategy.get_phrases();
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 42);

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
            std::cout << "    Creating the dictionary from the hash table" << std::flush;
            start = std::chrono::steady_clock::now();
            dictionary dict(map, dict_sym, max_freq, phrase_desc,
                            p_strategy.text_size, p_info.prev_alph,
                            std::max<size_t>(p_info.max_sym_freq, p_info.active_strings));
            end = std::chrono::steady_clock::now();
            map.destroy_data();
            report_time(start, end, 18);

            //process the dictionary
            tot_phrases = exact_algo::process_dictionary(dict, p_info, ws);
            p_info.prev_alph = dict.alphabet;
        }

        //reload the hash table
        map.load_data_from_file(ht_file);
        ws.remove_file("ht_data");

        {
            std::cout << "    Assigning ranks to the new phrases" << std::flush;
            start = std::chrono::steady_clock::now();
            bv_t new_phrase_desc(tot_phrases, false);
            key_wrapper key_w{sym_width(p_info.tot_phrases), map.description_bits(), map.get_data()};

            size_t j = 0;
            vector_t ranks;
            std::string ranks_file = ws.get_file("phr_ranks");
            load_from_file(ranks_file, ranks);
            for (auto const &ptr: map) {
                phrase_map_t::val_type val = 0;
                map.get_value_from(ptr, val);
                val = ranks[j++];
                map.insert_value_at(ptr, val);
                //the first bit marks if the phrase is repeated or not.
                // We need to shift it to get the real id
                new_phrase_desc[(val >> 1UL)] = phrase_desc[key_w.read(ptr, 0)];
            }
            ws.remove_file("phr_ranks");
            std::string suffix_file = ws.get_file("suffix_file");
            sdsl::store_to_file(new_phrase_desc, suffix_file);
            end = std::chrono::steady_clock::now();
            report_time(start, end, 27);
        }

        std::cout << "    Creating the parse of the text" << std::flush;
        start = std::chrono::steady_clock::now();
        load_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
        psize = p_strategy.parse_text();
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

        std::cout << "    Stats:" << std::endl;
        std::cout << "      Parsing phrases:                  " << p_info.lms_phrases << std::endl;
        std::cout << "      Number of symbols in the phrases: " << dict_sym << std::endl;
        std::cout << "      Number of unsolved BWT blocks:    " << p_info.tot_phrases << std::endl;
        std::cout << "      Processed strings:                " << p_info.active_strings << std::endl;
        std::cout << "      Parse size:                       " << psize << std::endl;

        map.destroy_data();
        map.destroy_table();

        if (psize == 0) {
            return 0;
        } else {
            return p_info.lms_phrases;
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

    template size_t
    par_round<ext_st_byte_parse_strategy>(ext_st_byte_parse_strategy &p_strat, parsing_info &p_gram, bv_t &phrase_desc,
                                            tmp_workspace &ws);

    template size_t
    par_round<ext_st_int_parse_strategy>(ext_st_int_parse_strategy &p_strat, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);

    template size_t
    par_round<ext_mt_byte_parse_strategy>(ext_mt_byte_parse_strategy &p_strat, parsing_info &p_gram, bv_t &phrase_desc,
                                            tmp_workspace &ws);

    template size_t
    par_round<ext_mt_int_parse_strategy>(ext_mt_int_parse_strategy &p_start, parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
}