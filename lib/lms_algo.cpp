//
// Created by Diego Diaz on 4/7/20.
//

#include "lms_algo.hpp"
#include <cmath>
#include "cdt/parallel_string_sort.hpp"

void assign_ids(phrase_map_t &mp_map, size_t max_sym, key_wrapper &key_w, ivb_t &r,
                      bvb_t &r_lim, size_t n_threads, sdsl::cache_config &config) {

    std::string syms_file = sdsl::cache_file_name("syms_file", config);
    {
        sdsl::int_vector_buffer<> syms_buff(syms_file, std::ios::out);
        for (auto const &phrase : mp_map) {
            syms_buff.push_back(phrase);
        }
        syms_buff.close();
    }
    auto compare = [&](const size_t &l, const size_t &r) -> bool {
        return key_w.compare(l, r);
    };
    auto access = [&](const size_t &val, size_t idx) -> size_t {
        return key_w.read(val, key_w.size(val)-1-idx);
    };
    parallel_str_sort(syms_file, compare, access, max_sym+1, n_threads, config);

    sdsl::int_vector<> k_list;
    sdsl::load_from_file(k_list, syms_file);

    if(remove(syms_file.c_str())){
        std::cout<<"Error trying to remove file "<<syms_file<<std::endl;
    }

    for(size_t m_pos=0; m_pos < mp_map.size(); m_pos++){

        size_t len = key_w.size(k_list[m_pos]);
        for(size_t i=len; i-- > 1;){
            r.push_back(key_w.read(k_list[m_pos], i));
            r_lim.push_back(false);
        }
        r.push_back( key_w.read(k_list[m_pos], 0));
        r_lim.push_back(true);

        //modify the key value
        phrase_map_t::val_type val=0;
        mp_map.get_value_from(k_list[m_pos], val);
        val |= (max_sym+m_pos+1)<<1UL;
        mp_map.insert_value_at(k_list[m_pos], val);
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
                map.insert_value_at(*res.first, 1UL);
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

template<class sym_t>
void * hash_phrases(void * data) {


    auto lms_data = (parsing_info<sym_t> *) data;

    bool s_type, prev_s_type = S_TYPE;
    sym_t curr_sym, prev_sym;

    string_t curr_lms(2, lms_data->sym_width);
    prev_sym = lms_data->ifs.read(lms_data->end);
    curr_lms.push_back(prev_sym);

    for(size_t i = lms_data->end; i-- > lms_data->start;){

        curr_sym = lms_data->ifs.read(i);

        //                                     L_TYPE   S_TYPE*
        //                                        ---- ----
        //this is a junction between two strings = ...$ $...
        if(lms_data->is_suffix(curr_sym)){
            bool full_str = curr_lms.size()==1 && lms_data->is_suffix(curr_lms[0]);
            if(!curr_lms.empty() && !full_str){
                lms_data->hash_phrase(curr_lms);
            }
            curr_lms.clear();
            s_type = S_TYPE;
        } else {
            if (curr_sym < prev_sym) {//S_TYPE type
                s_type = S_TYPE;
            } else if (curr_sym == prev_sym) {
                s_type = prev_s_type;
            } else {//L_TYPE type
                s_type = L_TYPE;

                if(prev_s_type == S_TYPE) {//Leftmost S-type suffix
                    curr_lms.pop_back();
                    if(!curr_lms.empty()){
                        lms_data->hash_phrase(curr_lms);
                        curr_lms.clear();
                    }
                    curr_lms.push_back(prev_sym);
                }
            }
        }
        curr_lms.push_back(curr_sym);
        prev_sym = curr_sym;
        prev_s_type = s_type;
    }

    assert(curr_lms[0]!=1);
    bool full_str = curr_lms.size()==1 &&
                    lms_data->is_suffix(curr_lms[0]) &&
                    (lms_data->start==0 || lms_data->is_suffix(lms_data->ifs.read(lms_data->start-1)));

    if(!curr_lms.empty() && !full_str){
        lms_data->hash_phrase(curr_lms);
    }

    lms_data->thread_dict.flush();
    pthread_exit(nullptr);
}
template<class sym_t>
void * record_phrases(void *data) {

    auto lms_data = (parsing_info<sym_t> *) data;

    bool s_type, prev_s_type = S_TYPE;
    sym_t curr_sym, prev_sym;

    string_t curr_lms(2, lms_data->sym_width);
    prev_sym = lms_data->ifs.read(lms_data->end);
    curr_lms.push_back(prev_sym);

    for(size_t i = lms_data->end; i-- > lms_data->start;){

        curr_sym = lms_data->ifs.read(i);

        //                                     L_TYPE   S_TYPE*
        //                                        ---- ----
        //this is a junction between two strings = ...$ $...
        if(lms_data->is_suffix(curr_sym)){
            if(!curr_lms.empty()){
                lms_data->store_phrase(curr_lms);
            }
            curr_lms.clear();
            s_type = S_TYPE;
        } else {
            if (curr_sym < prev_sym) {//S_TYPE type
                s_type = S_TYPE;
            } else if (curr_sym == prev_sym) {
                s_type = prev_s_type;
            } else {//L_TYPE type
                s_type = L_TYPE;

                if(prev_s_type == S_TYPE) {//Left-most suffix
                    curr_lms.pop_back();
                    if(!curr_lms.empty()){
                        lms_data->store_phrase(curr_lms);
                        curr_lms.clear();
                    }
                    curr_lms.push_back(prev_sym);
                }
            }
        }
        curr_lms.push_back(curr_sym);
        prev_sym = curr_sym;
        prev_s_type = s_type;
    }

    assert(curr_lms[0]!=1);
    if(!curr_lms.empty()){
        lms_data->store_phrase(curr_lms);
    }

    lms_data->ofs.close();
    lms_data->ifs.close();
    pthread_exit(nullptr);
}
template<class sym_type>
std::vector<std::pair<size_t, size_t>> compute_thread_ranges(size_t n_threads, std::string& i_file,
                                                             sdsl::int_vector<2>& phrase_desc) {
    std::vector<std::pair<size_t, size_t>> thread_ranges;

    i_file_stream<sym_type> is(i_file, BUFFER_SIZE);
    size_t n_chars = is.tot_cells;
    assert(n_chars>0);
    size_t sym_per_thread = INT_CEIL(n_chars, n_threads);
    size_t start, end;
    size_t eff_threads = INT_CEIL(n_chars, sym_per_thread);

    for(size_t i=0;i<eff_threads;i++){
        start = (i * sym_per_thread);
        end = std::min<size_t>(((i + 1) * sym_per_thread), n_chars-1);

        start = start==0? 0 : size_t(prev_lms_sym(start, is, phrase_desc)+1);
        long long tmp_end = prev_lms_sym(end, is, phrase_desc);
        end = tmp_end<0?  0 : size_t(tmp_end);
        if(start<end){
            thread_ranges.emplace_back(start, end);
        }
    }
    is.close();
    return thread_ranges;
}

template<class sym_type>
size_t lc_algo(std::string &i_file, std::string &o_file,
               size_t n_threads, size_t hbuff_size,
               gram_info_t &p_gram, ivb_t &rules,
               bvb_t &rules_lim, sdsl::int_vector<2> &phrase_desc,
               sdsl::cache_config &config) {

    phrase_map_t mp_table(0, "", 0.85);

    auto thread_ranges = compute_thread_ranges<sym_type>(n_threads, i_file, phrase_desc);

    std::vector<parsing_info<sym_type>> threads_data;
    threads_data.reserve(thread_ranges.size());

    //how many size_t cells we can fit in the buffer
    size_t buff_cells = hbuff_size/sizeof(size_t);

    //number of bytes per thread
    size_t hb_bytes = (buff_cells / thread_ranges.size()) * sizeof(size_t);

    void *buff_addr = malloc(hbuff_size);
    auto tmp_addr = reinterpret_cast<char*>(buff_addr);

    size_t k=0;
    for(auto &range : thread_ranges) {

        std::stringstream ss;
        ss << o_file.substr(0, o_file.size() - 5) << "_range_" << range.first << "_" << range.second;
        std::string tmp_o_file = ss.str();
        threads_data.emplace_back(i_file, tmp_o_file, mp_table, range.first, range.second,
                                  p_gram.r, hb_bytes, tmp_addr + (k*hb_bytes), phrase_desc);
        k++;
    }

    std::cout<<"      Computing the LMS phrases in the text"<<std::endl;
    {
        std::vector<pthread_t> threads(threads_data.size());
        for(size_t i=0;i<threads_data.size();i++){
            int ret =  pthread_create(&threads[i],
                                      nullptr,
                                      &hash_phrases<sym_type>,
                                      (void*)&threads_data[i]);
            if(ret != 0) {
                printf("Error: pthread_create() failed\n");
                exit(EXIT_FAILURE);
            }
        }

        for(size_t i=0;i<threads_data.size();i++) {
            pthread_join(threads[i], nullptr);
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

        p_gram.rules_breaks.push_back(mp_table.size());
        size_t width = sdsl::bits::hi(p_gram.r+1)+1;
        const bitstream<ht_buff_t>& stream = mp_table.get_data();
        key_wrapper key_w{width, mp_table.description_bits(), stream};

        //temporal unload of the hash table (not the data)
        std::string st_table = sdsl::cache_file_name("ht_data", config);
        mp_table.unload_table(st_table);

        //rename phrases according to their lexicographical ranks
        std::cout<<"      Assigning identifiers to the phrases"<<std::endl;
        assign_ids(mp_table, p_gram.r-1,  key_w, rules, rules_lim, n_threads, config);

        //reload the hash table
        mp_table.load_table(st_table);
        if(remove(st_table.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }

        std::cout<<"      Creating the parse of the text"<<std::endl;
        {//store the phrases into a new file
            std::vector<pthread_t> threads(threads_data.size());
            for(size_t i=0;i<threads_data.size();i++){
                int ret =  pthread_create(&threads[i],
                                          nullptr,
                                          &record_phrases<sym_type>,
                                          (void*)&threads_data[i]);
                if(ret != 0) {
                    printf("Error: pthread_create() failed\n");
                    exit(EXIT_FAILURE);
                }
            }

            for(size_t i=0;i<threads_data.size();i++) {
                pthread_join(threads[i], nullptr);
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
            //keep track of the lms phrases that have to be rephrased
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

template size_t lc_algo<uint8_t>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                 gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                                 sdsl::int_vector<2> &phrase_desc, sdsl::cache_config &config);

template size_t lc_algo<size_t>(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                gram_info_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                                sdsl::int_vector<2> &phrase_desc, sdsl::cache_config &config);
