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

void comp_dict_int_v2(dictionary &dict, phrase_map_t& phrases_ht, vector_t& s_sa, bv_t& phr_marks) {
    dict.dict_dummy = dict.alphabet+s_sa.size()+1;
    size_t pos, em_nt, new_size=0, sym;
    string_t phrase(2, sdsl::bits::hi(dict.alphabet)+1);
    vector_t new_dict(s_sa.size()*2, 0, sdsl::bits::hi(dict.dict_dummy)+1);

    size_t rank=0;
    for(auto && u : s_sa) {
        pos = u;

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

                auto res = phrases_ht.new_find(phrase.data(), phrase.n_bits());
                assert(res.second);
                //auto res2 = phrases_ht.new_find(phrase.data(), phrase.n_bits());
                //assert(res.second==res2.second);

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
}

size_t compress_dictionary_v2(dictionary &dict, vector_t &sa, parsing_info& p_info, tmp_workspace& ws) {

    size_t p_round = p_info.p_round;
    string_t phrase(2, sdsl::bits::hi(dict.alphabet)+1);
    bool is_maximal, exist_as_phrase;
    size_t u=0, d_pos, pl_sym, bg_pos, freq, rank=0, l_sym, dummy_sym = dict.alphabet+1, f_sa_pos;

    bv_rs_t d_lim_rs(&dict.d_lim);

    std::string pre_bwt_file = ws.get_file("pre_bwt_lev_"+std::to_string(p_round));
    size_t sb = INT_CEIL(sdsl::bits::hi(dummy_sym)+1, 8);
    size_t fb = INT_CEIL(sdsl::bits::hi(dict.t_size)+1,8);
    bwt_buff_writer pre_bwt(pre_bwt_file, std::ios::out, sb, fb);

    phrase_map_t new_phrases_ht;

    size_t width = sdsl::bits::hi(dict.alphabet+dict.dict.size()-dict.n_phrases)+1;
    vector_t ranks(dict.n_phrases, 0, width);

    bv_t phr_marks(dict.dict.size(), false);

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
                ranks[d_lim_rs(d_pos)] = rank;
            }else{
                exist_as_phrase = false;
            }

            u++;
            while(u<sa.size() && sa[u] & 1UL){
                d_pos = (sa[u]>>1UL) - 1;
                freq += dict.freqs[d_lim_rs(d_pos)];
                l_sym = d_pos==0 || dict.d_lim[d_pos-1] ? dummy_sym : dict.dict[d_pos-1];
                if(!is_maximal && l_sym!=pl_sym) is_maximal = true;
                if(!exist_as_phrase && l_sym==dummy_sym){
                    ranks[d_lim_rs(d_pos)] = rank;
                    exist_as_phrase = true;
                }
                pl_sym = l_sym;
                u++;
            }

            if(is_maximal || exist_as_phrase){
                if(u-bg_pos>1){
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
                    dict.phrases_has_hocc[rank] = true;

                    for(size_t j=bg_pos;j<u;j++){
                        phr_marks[(sa[j]>>1UL)-1] = true;
                    }
                }

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
    pre_bwt.close();

    sa.resize(rank);
    dict.phrases_has_hocc.resize(rank);
    sdsl::util::clear(d_lim_rs);
    sdsl::util::clear(dict.freqs);
    sdsl::store_to_file(ranks, ws.get_file("phr_ranks"));
    sdsl::util::clear(ranks);
    new_phrases_ht.shrink_databuff();

#ifdef __linux__
    malloc_trim(0);
#endif

    auto start = std::chrono::steady_clock::now();
    comp_dict_int_v2(dict, new_phrases_ht, sa, phr_marks);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 4);

    std::string dict_file = ws.get_file("dict_lev_"+std::to_string(p_round));
    sdsl::store_to_file(dict, dict_file);
    return sa.size();
}

size_t assign_ids(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

    std::cout<<"Suffix induction"<<std::endl;
    auto start = std::chrono::steady_clock::now();
    vector_t sa(dict.dict.size(), 0, sdsl::bits::hi(dict.dict.size())+2);

    uint8_t width = sdsl::bits::hi(dict.dict.size())+1;

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

    std::cout<<"Number of phrases : "<<dict.n_phrases<<std::endl;
    std::cout<<"Alphabet :          "<<dict.alphabet<<std::endl;

    std::cout<<"dict  "<<double(sdsl::size_in_bytes(dict.dict))/1000000<<std::endl;
    std::cout<<"freqs "<<double(sdsl::size_in_bytes(dict.freqs))/1000000<<std::endl;
    std::cout<<"sa    "<<double(sdsl::size_in_bytes(sa))/1000000<<std::endl;
    std::cout<<"has_hocc    "<<double(sdsl::size_in_bytes(dict.phrases_has_hocc))/1000000<<std::endl;

    start = std::chrono::steady_clock::now();
    size_t new_number_of_phrases = compress_dictionary_v2(dict, sa, p_info, ws);
    end = std::chrono::steady_clock::now();
    std::cout<<"Sort and compress dictionary"<<std::endl;
    report_time(start, end, 4);
    return new_number_of_phrases;
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
size_t build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size, alpha_t& alphabet, tmp_workspace &ws) {

    std::cout<<"Parsing the text:    "<<std::endl;
    std::string output_file = ws.get_file("tmp_output");
    std::string tmp_i_file = ws.get_file("tmp_input");

    // mark which symbols represent string boundaries
    bv_t symbol_desc(alphabet.back().first+1,false);
    symbol_desc[alphabet[0].first] = true;

    parsing_info p_info;
    for(auto const& pair : alphabet){
        if(p_info.max_sym_freq<pair.second){
            p_info.max_sym_freq = pair.second;
        }
    }

    size_t iter=1;
    size_t rem_phrases;

    typedef lc_parser_t<i_file_stream<uint8_t>, string_t> byte_parser_t;
    typedef lc_parser_t<i_file_stream<size_t>, string_t> int_parser_t;

    std::cout<<"  Parsing round "<<iter++<<std::endl;
    auto start = std::chrono::steady_clock::now();
    rem_phrases = build_lc_gram_int<byte_parser_t>(i_file, tmp_i_file, n_threads, hbuff_size,
                                                   p_info, symbol_desc, ws);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 4);

    while (rem_phrases > 0) {
        start = std::chrono::steady_clock::now();
        std::cout<<"  Parsing round "<<iter++<<std::endl;
        rem_phrases = build_lc_gram_int<int_parser_t>(tmp_i_file, output_file, n_threads,
                                                      hbuff_size, p_info, symbol_desc, ws);
        end = std::chrono::steady_clock::now();

        report_time(start, end,4);
        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
    }

    sdsl::util::clear(symbol_desc);
    return iter-2;
}

template<class parser_t, class out_sym_t>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                       parsing_info &p_info, bv_t &phrase_desc, tmp_workspace &ws) {

#ifdef __linux__
    malloc_trim(0);
#endif

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

    //join the different phrase files
    std::vector<std::string> phrases_files;
    for(size_t i=0;i<threads_data.size();i++){
        phrases_files.push_back(threads_data[i].thread_dict.dump_file());
    }
    auto join_res = join_thread_phrases(mp_table, phrases_files);

    size_t psize=0;//<- for the iter stats
    if(mp_table.size()!=p_info.lms_phrases) {

        size_t width = sdsl::bits::hi(phrase_desc.size())+1;
        size_t dict_syms = join_res.first/width;
        size_t max_freq = join_res.second;

        //save a copy of the hash table into a file
        std::string ht_file = ws.get_file("ht_data");
        mp_table.store_data_to_file(ht_file);

        //temporal unload of the hash table (not the data)
        mp_table.destroy_table();
        size_t tot_phrases=0;
        {
            //create a dictionary from where the ids will be computed
            std::cout<<"    Creating the dictionary from the hash table"<<std::endl;
            dictionary dict(mp_table, dict_syms, max_freq, phrase_desc,
                            threads_data[0].ifs.size(), p_info.prev_alph,
                            p_info.max_sym_freq);
            mp_table.destroy_data();

            //rename phrases according to their lexicographical ranks
            std::cout<<"    Assigning identifiers to the phrases"<<std::endl;
            tot_phrases = assign_ids(dict, p_info, ws);
        }

        //reload the hash table
        mp_table.load_data_from_file(ht_file);
        if(remove(ht_file.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }

        {
            bv_t new_phrase_desc(tot_phrases, false);
            key_wrapper key_w{width, mp_table.description_bits(), mp_table.get_data()};

            size_t j=0;
            vector_t ranks;
            std::string ranks_file = ws.get_file("phr_ranks");
            sdsl::load_from_file(ranks, ranks_file);
            for(auto const& ptr : mp_table){
                phrase_map_t::val_type val=0;
                mp_table.get_value_from(ptr, val);
                val = ranks[j++];
                mp_table.insert_value_at(ptr, val);
                new_phrase_desc[val] = phrase_desc[key_w.read(ptr, 0)];
            }
            std::string suffix_file = ws.get_file("suffix_file");
            sdsl::store_to_file(new_phrase_desc, suffix_file);
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
            //key_wrapper key_w{width, mp_table.description_bits(), mp_table.get_data()};

            std::string suffix_file = ws.get_file("suffix_file");
            bv_t new_phrase_desc;
            sdsl::load_from_file(new_phrase_desc, suffix_file);
            /*auto it = mp_table.begin();
            auto it_end = mp_table.end();
            size_t sym;
            while (it != it_end) {
                auto val = it.value();
                //read the (reversed) last symbol
                sym = key_w.read(*it, 0);
                new_phrase_desc[val] = phrase_desc[sym];
                ++it;
            }*/
            p_info.prev_alph = phrase_desc.size();
            phrase_desc.swap(new_phrase_desc);
        }

        p_info.max_sym_freq = max_freq;
        p_info.lms_phrases = mp_table.size();
        p_info.tot_phrases = tot_phrases;
        p_info.p_round++;

        std::cout<<"    Stats:"<<std::endl;
        std::cout<<"      LMS phrases:    "<<p_info.lms_phrases<<std::endl;
        std::cout<<"      Total phrases:  "<<p_info.tot_phrases<<std::endl;
        std::cout<<"      Parse size:     "<<psize<<std::endl;

        return p_info.lms_phrases;

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

template size_t build_lc_gram<lms_parsing>(std::string &i_file, size_t n_threads, size_t hbuff_size, alpha_t& alphabet, tmp_workspace &ws);

template size_t build_lc_gram_int<lms_parsing<i_file_stream<uint8_t>, string_t>>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                                                                 parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);

template size_t build_lc_gram_int<lms_parsing<i_file_stream<size_t>, string_t>>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                                                                parsing_info &p_gram, bv_t &phrase_desc, tmp_workspace &ws);
