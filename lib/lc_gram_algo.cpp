//
// Created by Diego Diaz on 4/7/20.
//

#include "lc_gram_algo.hpp"
#include <cmath>
#include <thread>

void assign_ids(phrase_map_t &mp_map, size_t max_sym, size_t min_sym,
                key_wrapper &key_w, bvb_t &r_lim, sdsl::cache_config &config, ivb_t &r) {

    std::string dict_file = sdsl::cache_file_name("dict_file", config);
    std::string dict_lim_file = sdsl::cache_file_name("dict_lim_file", config);
    std::string sa_file = sdsl::cache_file_name("sa_file", config);
    {
        ivb_t dict(dict_file, std::ios::out);
        sdsl::int_vector_buffer<1> d_lim(dict_lim_file, std::ios::out);
        for (auto const &ptr : mp_map) {
            for(size_t i=key_w.size(ptr);i-->0;){
                dict.push_back(key_w.read(ptr, i)-min_sym);
                d_lim.push_back(false);
            }
            d_lim[d_lim.size()-1] = true;
        }
        dict.close();
        d_lim.close();
    }

    vector_t ranks(mp_map.size(), 0, sdsl::bits::hi(mp_map.size())+1);
    {
        suffix_induction(sa_file, dict_file, dict_lim_file, max_sym-min_sym+1);
        bv_t d_lim;
        vector_t dict, sa;
        sdsl::load_from_file(d_lim, dict_lim_file);
        sdsl::load_from_file(dict, dict_file);
        sdsl::load_from_file(sa, sa_file);
        bv_rs_t d_lim_rs(&d_lim);
        size_t rank = 0;

        for(auto pos : sa){

            if(pos==0) continue;

            pos--;
            if(pos==0 || d_lim[pos-1]){
                ranks[d_lim_rs(pos)] = rank++;
                size_t j=pos;
                while(!d_lim[j]){
                    r.push_back(min_sym+dict[j++]);
                    r_lim.push_back(false);
                }
                assert(d_lim[j]);
                r.push_back(min_sym+dict[j]);
                r_lim.push_back(true);
            }
        }

        /*for(auto const& pair : tmp_map){
            if(pair.second>1){
                std::cout<<pair.first<<" "<<pair.second<<" "<<dict.size()<<" "<<d_lim[pair.first-1]<<std::endl;
            }
        }
        std::cout<<rank<<" "<<mp_map.size()<<std::endl;
        assert(rank==mp_map.size());*/
    }

    //assign the ranks
    size_t j=0;
    for(auto const& ptr : mp_map){
        //modify the key value
        phrase_map_t::val_type val=0;
        mp_map.get_value_from(ptr, val);
        val |= (max_sym+ranks[j++]+1)<<1UL;
        mp_map.insert_value_at(ptr, val);
        //
    }
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


void join_thread_phrases(phrase_map_t& map, std::vector<std::string> &files) {

    bool rep;

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
            rep = bits.read(next_bit, next_bit);
            next_bit+=33;

            auto res = map.insert(key, key_bits, rep);
            if(!res.second){
                map.insert_value_at(res.first, 1UL);
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
}

template<class parser_t>
std::vector<std::pair<size_t, size_t>> compute_thread_ranges(size_t n_threads,
                                                             std::string& i_file,
                                                             parser_t& parser) {
    std::vector<std::pair<size_t, size_t>> thread_ranges;

    typename parser_t::stream_type is(i_file, BUFFER_SIZE);
    size_t n_chars = is.tot_cells;
    assert(n_chars>0);
    size_t sym_per_thread = INT_CEIL(n_chars, n_threads);
    size_t start, end;
    size_t eff_threads = INT_CEIL(n_chars, sym_per_thread);

    for(size_t i=0;i<eff_threads;i++){
        start = (i * sym_per_thread);
        end = std::min<size_t>(((i + 1) * sym_per_thread), n_chars-1);

        start = start==0? 0 : size_t(parser.prev_break(start, is)+1);
        long long tmp_end = parser.prev_break(end, is);

        end = tmp_end<0?  0 : size_t(tmp_end);
        if(start<end){
            thread_ranges.emplace_back(start, end);
        }
    }
    is.close();
    return thread_ranges;
}

template<template<class, class> class lc_parser_t>
void build_lc_gram(std::string &i_file, size_t n_threads, size_t hbuff_size,
                      gram_info_t &p_gram, alpha_t alphabet, sdsl::cache_config &config){

    std::cout<<"  Generating a locally consistent grammar:    "<<std::endl;
    std::string output_file = sdsl::cache_file_name("tmp_output", config);
    std::string tmp_i_file = sdsl::cache_file_name("tmp_input", config);

    // given an index i in symbol_desc
    //0 symbol i is in alphabet is unique
    //1 symbol i is repeated
    //>2 symbol i is sep symbol
    sdsl::int_vector<2> symbol_desc(alphabet.back().first+1,0);

    for(auto & sym : alphabet){
        symbol_desc[sym.first] = sym.second > 1;
    }
    symbol_desc[alphabet[0].first]+=2;

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

    std::cout<<"    Parsing round "<<iter++<<std::endl;
    rem_phrases = build_lc_gram_int<byte_parser_t>(i_file, tmp_i_file,
                                             n_threads, hbuff_size,
                                             p_gram, rules, rules_lim,
                                             symbol_desc, config);

    while (rem_phrases > 0) {
        std::cout<<"    Parsing round "<<iter++<<std::endl;
        rem_phrases = build_lc_gram_int<int_parser_t>(tmp_i_file, output_file,
                                                      n_threads, hbuff_size,
                                                      p_gram, rules, rules_lim,
                                                      symbol_desc, config);
        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
    }

    sdsl::util::clear(symbol_desc);

    {//put the compressed string at end
        std::ifstream c_vec(tmp_i_file, std::ifstream::binary);
        c_vec.seekg(0, std::ifstream::end);
        size_t tot_bytes = c_vec.tellg();
        c_vec.seekg(0, std::ifstream::beg);
        auto *buffer = reinterpret_cast<size_t*>(malloc(BUFFER_SIZE));
        size_t read_bytes =0;
        p_gram.c=0;
        while(read_bytes<tot_bytes){
            c_vec.read((char *) buffer, BUFFER_SIZE);
            read_bytes+=c_vec.gcount();
            assert((c_vec.gcount() % sizeof(size_t))==0);
            for(size_t i=0;i<c_vec.gcount()/sizeof(size_t);i++){
                rules.push_back(buffer[i]);
                rules_lim.push_back(false);
                p_gram.c++;
            }
        }
        rules_lim[rules_lim.size() - 1] = true;
        p_gram.r++;
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
    std::swap(p_gram.rules_breaks, rule_breaks);

    rules.close();
    rules_lim.close();

    std::cout<<"  Locally consistent grammar finished"<<std::endl;
    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Number of terimnals:    "<<(int)p_gram.sigma<<std::endl;
    std::cout<<"      Number of nonterminals: "<<p_gram.rules_breaks.back()-(p_gram.max_tsym+1)+1<<std::endl;
    std::cout<<"      Grammar size:           "<<p_gram.g<<std::endl;
    std::cout<<"      Compressed string:      "<<p_gram.c<<std::endl;

    if(remove(tmp_i_file.c_str())){
        std::cout<<"Error trying to delete file "<<tmp_i_file<<std::endl;
    }
}

template<class parser_t, class out_sym_t>
size_t build_lc_gram_int(std::string &i_file, std::string &o_file,
                         size_t n_threads, size_t hbuff_size,
                         gram_info_t &p_gram, ivb_t &rules,
                         bvb_t &rules_lim, sdsl::int_vector<2> &phrase_desc,
                         sdsl::cache_config &config) {

    typedef typename parser_t::sym_type          sym_type;
    typedef typename parser_t::stream_type       stream_type;
    typedef parse_data_t<stream_type, out_sym_t> parse_data_type;

    parser_t parser(phrase_desc);

    phrase_map_t mp_table(0, "", 0.85);

    auto thread_ranges = compute_thread_ranges<parser_t>(n_threads, i_file, parser);

    std::vector<parse_data_type> threads_data;
    threads_data.reserve(thread_ranges.size());

    //how many size_t cells we can fit in the buffer
    size_t buff_cells = hbuff_size/sizeof(size_t);

    //number of bytes per thread
    size_t hb_bytes = (buff_cells / thread_ranges.size()) * sizeof(size_t);

    void *buff_addr = malloc(hbuff_size);
    auto tmp_addr = reinterpret_cast<char*>(buff_addr);

    size_t k=0;
    for(auto &range : thread_ranges) {
        std::string tmp_o_file = o_file.substr(0, o_file.size() - 5);
        tmp_o_file.append("_range_"+std::to_string(range.first)+"_"+std::to_string(range.second));
        threads_data.emplace_back(i_file, tmp_o_file, mp_table, range.first, range.second, hb_bytes,
                                  tmp_addr + (k*hb_bytes));
        k++;
    }

    std::cout<<"      Computing the phrases in the text"<<std::endl;
    {
        std::vector<std::thread> threads(threads_data.size());
        hash_functor<parse_data_type, parser_t> hf;
        for(size_t i=0;i<threads_data.size();i++){
            threads[i] = std::thread(hf, std::ref(threads_data[i]), std::ref(parser));
        }

        for(size_t i=0;i<threads_data.size();i++) {
            threads[i].join();
        }
    }
    free(buff_addr);

    //join the different phrase files
    std::vector<std::string> phrases_files;
    for(size_t i=0;i<threads_data.size();i++){
        phrases_files.push_back(threads_data[i].thread_dict.dump_file());
    }
    join_thread_phrases(mp_table, phrases_files);

    size_t psize=0;//<- for the iter stats
    if(mp_table.size()>0){

        size_t min_sym = p_gram.rules_breaks.empty() ? 0 : p_gram.r - p_gram.rules_breaks.back();
        size_t max_sym = p_gram.r-1;

        p_gram.rules_breaks.push_back(mp_table.size());
        size_t width = sdsl::bits::hi(p_gram.r+1)+1;
        const bitstream<ht_buff_t>& stream = mp_table.get_data();
        key_wrapper key_w{width, mp_table.description_bits(), stream};

        //temporal unload of the hash table (not the data)
        std::string st_table = sdsl::cache_file_name("ht_data", config);
        mp_table.unload_table(st_table);

        //rename phrases according to their lexicographical ranks
        std::cout<<"      Assigning identifiers to the phrases"<<std::endl;
        assign_ids(mp_table, max_sym, min_sym, key_w, rules_lim, config, rules);

        //reload the hash table
        mp_table.load_table(st_table);
        if(remove(st_table.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }

        std::cout<<"      Creating the parse of the text"<<std::endl;
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
            phrase_desc.resize(p_gram.r+mp_table.size());
            std::cout << "      Updating symbols status" << std::endl;
            auto it = mp_table.begin();
            auto it_end = mp_table.end();
            size_t tmp_value, sym;

            while (it != it_end) {

                tmp_value = 0;
                auto val = it.value();

                //more than one occurrence of the phrase
                if (val & 1UL) {
                    tmp_value += 1;
                }

                //read the (reversed) last symbol
                sym = key_w.read(*it, 0);
                if (phrase_desc[sym] & 2U) {//phrase is suffix of some string
                    tmp_value += 2;
                }

                phrase_desc[val >> 1UL] = tmp_value;
                ++it;
            }
        }
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
    }

    p_gram.r +=mp_table.size();
    std::cout<<"      Stats:"<<std::endl;
    std::cout<<"        Parse size:          "<<psize<<std::endl;
    std::cout<<"        New nonterminals:    "<<mp_table.size()<<std::endl;

    if(psize>1){
        return mp_table.size();
    }else{
        return 0;
    }
}

template void build_lc_gram<lms_parsing>(std::string &i_file, size_t n_threads, size_t hbuff_size, gram_info_t &p_gram, alpha_t alphabet, sdsl::cache_config &config);

template size_t build_lc_gram_int<lms_parsing<i_file_stream<uint8_t>, string_t>>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                           gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                                           sdsl::int_vector<2> &phrase_desc, sdsl::cache_config &config);

template size_t build_lc_gram_int<lms_parsing<i_file_stream<size_t>, string_t>>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                          gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                                          sdsl::int_vector<2> &phrase_desc, sdsl::cache_config &config);
