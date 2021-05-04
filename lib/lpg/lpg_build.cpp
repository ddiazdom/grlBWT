//
// Created by Diego Diaz on 4/7/20.
//

#include "lpg/lpg_build.hpp"
#include <cmath>
#include <sdsl/select_support_mcl.hpp>
#include "cdt/parallel_string_sort.hpp"

pthread_mutex_t thread_mutex=PTHREAD_MUTEX_INITIALIZER;

void lpg_build::check_plain_grammar(std::string& g_file, std::string& uncomp_file) {

    plain_grammar_t p_gram;
    p_gram.load_from_file(g_file);

    sdsl::int_vector<> r;

    bv_t r_lim;
    sdsl::load_from_file(r, p_gram.rules_file);
    sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
    bv_t::select_1_type r_lim_ss;
    sdsl::util::init_support(r_lim_ss, &r_lim);

    std::cout<<"Checking the grammar produces the exact input string"<<std::endl;
    std::cout<<"  Terminals:              "<<(size_t)p_gram.sigma<<std::endl;
    std::cout<<"  Number of nonterminals: "<<p_gram.r-p_gram.sigma<<std::endl;
    std::cout<<"  Compressed string:      "<<p_gram.c<<std::endl;

    std::vector<size_t> tmp_decomp;

    i_file_stream<uint8_t> if_stream(uncomp_file, BUFFER_SIZE);
    uint8_t buff_symbol;

    size_t pos=0, curr_sym, start, end;
    std::stack<size_t> stack;

    size_t f = r_lim_ss(p_gram.r - 1) + 1;
    size_t l = r_lim_ss(p_gram.r);

    for(size_t i=f; i <= l; i++){
        tmp_decomp.clear();
        stack.push(r[i]);

        while(!stack.empty()){

            curr_sym = stack.top() ;
            stack.pop();

            if(curr_sym==0){
                start = 0;
            }else{
                start = r_lim_ss(curr_sym)+1;
            }

            end = r_lim_ss(curr_sym+1);

            if(r[start] == curr_sym){
                assert((end-start+1)==1);
                tmp_decomp.push_back(curr_sym);
            }else{
                for(size_t j=end+1; j-->start;){
                    stack.push(r[j]);
                }
            }
        }

        for(auto const& tmp_sym : tmp_decomp){
            buff_symbol = if_stream.read(pos++);
            //assert(tmp_sym == buff_symbol);
            assert(p_gram.symbols_map[tmp_sym] == buff_symbol);
        }
    }
    std::cout<<"\tGrammar is correct!!"<<std::endl;
}

void lpg_build::compute_LPG(std::string &i_file, std::string &p_gram_file, size_t n_threads, sdsl::cache_config &config,
                            size_t hbuff_size, alpha_t &alphabet, bool simp, bool rl_comp) {

    std::string rules_file = sdsl::cache_file_name("rules", config);
    std::string rules_len_file = sdsl::cache_file_name("rules_len", config);
    std::string lmsg_as_sp_file = sdsl::cache_file_name("lmsg_as_sp", config);

    plain_grammar_t p_gram(rules_file, rules_len_file, lmsg_as_sp_file);
    p_gram.sigma = alphabet.size();

    // given an index i in symbol_desc
    //0 symbol i is in alphabet is unique
    //1 symbol i is repeated
    //>2 symbol i is sep symbol
    sdsl::int_vector<2> symbol_desc(alphabet.back().first+1,0);

    for(auto & sym : alphabet){
        p_gram.symbols_map.push_back(sym.first);
        symbol_desc[sym.first] = sym.second > 1;
    }
    p_gram.r = alphabet.back().first + 1;
    symbol_desc[alphabet[0].first]+=2;

    ivb_t rules(p_gram.rules_file, std::ios::out, BUFFER_SIZE);
    bvb_t rules_lim(p_gram.rules_lim_file, std::ios::out);
    for(size_t i=0;i<p_gram.r; i++){
        rules.push_back(i);
        rules_lim.push_back(true);
    }
    for(unsigned char i : p_gram.symbols_map){
        rules[i] = i;
    }

    std::string output_file = sdsl::cache_file_name("tmp_output", config);
    std::string tmp_i_file = sdsl::cache_file_name("tmp_input", config);

    size_t iter=1;
    size_t rem_phrases;


    std::cout<<"  Iteration "<<iter++<<std::endl;
    rem_phrases = compute_LPG_int<uint8_t>(i_file, tmp_i_file,
                                           n_threads, hbuff_size,
                                           p_gram, rules, rules_lim,
                                           symbol_desc, config);

    while (rem_phrases > 0) {
        std::cout<<"  Iteration "<<iter++<<std::endl;
        rem_phrases = compute_LPG_int<size_t>(tmp_i_file, output_file,
                                              n_threads, hbuff_size,
                                              p_gram, rules, rules_lim,
                                              symbol_desc, config);
        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());
    }
    sdsl::util::clear(symbol_desc);

    {//put the compressed string at position rules[0]
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
    rules.close();
    rules_lim.close();

    simplify_grammar(p_gram, config);
    p_gram.save_to_file(p_gram_file);

    //TODO testing
    check_plain_grammar(p_gram_file, i_file);
    exit(1);
    //

    std::cout<<"  Collapsing nonterminal rules"<<std::endl;
    collapse_grammar(p_gram, iter, config);

    /*if(simp){
        std::cout<<"  Simplifying the grammar"<<std::endl;
        simplify_grammar(p_gram, config);
    }*/

    if(rl_comp){
        std::cout<<"  Creating run-length compressed nonterminals"<<std::endl;
    }

    p_gram.save_to_file(p_gram_file);

    //TODO testing
    check_plain_grammar(p_gram_file, i_file);
    //

    if(remove(tmp_i_file.c_str())){
        std::cout<<"Error trying to delete file "<<tmp_i_file<<std::endl;
    }
}

template<class sym_type>
size_t lpg_build::compute_LPG_int(std::string &i_file, std::string &o_file,
                                  size_t n_threads, size_t hbuff_size,
                                  plain_grammar_t &p_gram, ivb_t &rules,
                                  bvb_t &rules_lim, sdsl::int_vector<2> &phrase_desc,
                                  sdsl::cache_config &config) {

    phrase_map_t mp_table(0, "", 0.85);
    size_t rem_phrases=0;

    auto thread_ranges = compute_thread_ranges<sym_type>(n_threads, i_file, phrase_desc);

    std::vector<lms_info<sym_type>> threads_data;
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

    std::cout<<"    Computing the phrases of the text"<<std::endl;
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
        phrases_files.push_back(threads_data[i].thread_map.dump_file());
    }
    join_thread_phrases(mp_table, phrases_files);

    size_t psize=0;//<- for the iter stats
    if(mp_table.size()>0){

        size_t width = sdsl::bits::hi(p_gram.r+1)+1;
        const bitstream<buff_t>& stream = mp_table.get_data();
        key_wrapper key_w{width, mp_table.description_bits(), stream};

        //temporal unload of the hash table (not the data)
        std::string st_table = sdsl::cache_file_name("ht_data", config);
        mp_table.unload_table(st_table);

        //rename phrases according to their lexicographical ranks
        std::cout<<"    Assigning identifiers to the phrases"<<std::endl;
        assign_ids(mp_table, p_gram.r-1,  key_w, rules, rules_lim, n_threads, config);

        //reload the hash table
        mp_table.load_table(st_table);
        if(remove(st_table.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }

        std::cout<<"    Creating the parse of the text"<<std::endl;
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
            std::cout << "    Updating symbols status" << std::endl;
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
        rem_phrases=0;
    }

    p_gram.r +=mp_table.size();
    std::cout<<"    Iter. stats:"<<std::endl;
    std::cout<<"      Parse size:          "<<psize<<std::endl;
    std::cout<<"      New nonterminals:    "<<mp_table.size()<<std::endl;
    return mp_table.size();
}

void
lpg_build::assign_ids(phrase_map_t &mp_map, size_t max_sym, key_wrapper &key_w, ivb_t &r,
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

void lpg_build::join_parse_chunks(const std::string &output_file, std::vector<std::string> &chunk_files) {

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

template<class sym_t>
void * lpg_build::hash_phrases(void * data) {

    auto lms_data = (lms_info<sym_t> *) data;

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
    lms_data->thread_map.flush();

    pthread_exit(nullptr);
}

void lpg_build::join_thread_phrases(phrase_map_t& map, std::vector<std::string> &files) {

    bool rep;

    for(auto const& file : files){

        std::ifstream text_i(file, std::ios_base::binary);

        text_i.seekg (0, std::ifstream::end);
        size_t tot_bytes = text_i.tellg();
        text_i.seekg (0, std::ifstream::beg);

        auto buffer = reinterpret_cast<char *>(malloc(tot_bytes));

        text_i.read(buffer, tot_bytes);

        bitstream<buff_t> bits;
        bits.stream = reinterpret_cast<buff_t*>(buffer);

        size_t next_bit = 32;
        size_t tot_bits = tot_bytes*8;
        size_t key_bits;
        void* key=nullptr;
        size_t max_key_bits=0;

        while(next_bit<tot_bits){

            key_bits = bits.read(next_bit-32, next_bit-1);

            size_t n_bytes = INT_CEIL(key_bits, bitstream<buff_t>::word_bits)*sizeof(buff_t);
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
void * lpg_build::record_phrases(void *data) {

    auto lms_data = (lms_info<sym_t> *) data;

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
                    curr_lms.clear();
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

void lpg_build::collapse_grammar(plain_grammar_t &r_data, size_t &n_iter, sdsl::cache_config &config) {

    //load the rule limits
    bv_t r_lim;
    sdsl::load_from_file(r_lim, r_data.rules_lim_file);
    bv_t::select_1_type r_lim_ss;
    sdsl::util::init_support(r_lim_ss, &r_lim);
    size_t desc = 0, sym;

    {
        //load the rule symbols
        size_t pos = 0;
        ivb_t r_buff(r_data.rules_file, std::ios::in);

        //adjust the width of the symbol rules
        sdsl::int_vector<> r(r_buff.size(), 0, sdsl::bits::hi(r_data.r) + 1);
        while (pos < r_buff.size()) {
            r[pos] = r_buff[pos];
            pos++;
        }
        r_buff.close();

        //fix symbol references
        std::vector<std::pair<size_t, size_t>> pred(n_iter);
        size_t first_rule, start, l;//, len, end;
        bv_t marked_rules(r_data.r, false);

        for (size_t j = r_data.r; j-- > 1;) {
            if (marked_rules[j]) continue;
            start = r_lim_ss(j) + 1;
            if (r_lim[start] && (r[start] != j && r[start] != 0)) {
                l = 0;
                first_rule = j;

                while (r_lim[start]) {
                    pred[l++] = {start, first_rule};
                    first_rule = r[start];
                    start = r_lim_ss(first_rule) + 1;
                }
                for (size_t i = 0; i < l; i++) {
                    r[pred[i].first] = first_rule;
                    marked_rules[pred[i].second] = true;
                }
            }
        }
        r_data.c = r_lim_ss(r_data.r) - r_lim_ss(r_data.r - 1);

        //update references within the rules
        for (size_t j = 0; j < r.size(); j++) {
            //skip rules of length one
            desc = (desc << 1UL) | r_lim[j];
            if ((desc & 3UL) == 3UL) continue;

            //skip empty rules
            sym = r[j];
            if (sym == 0) continue;

            start = r_lim_ss(sym) + 1;

            if (r_lim[start] && sym != r[start]) {
                r[j] = r[start];
            }
        }
        sdsl::util::clear(r_lim_ss);
        sdsl::store_to_file(r, r_data.rules_file);
    }

    //collapse the grammar
    std::string tmp_r_file = sdsl::cache_file_name("tmp_r_file", config);
    std::string map_file = sdsl::cache_file_name("new_r_ids", config);
    sdsl::int_vector_buffer<> r_old(r_data.rules_file, std::ios::in, 1024 * 1024, sdsl::bits::hi(r_data.r) + 1);//old version of the grammar
    sdsl::int_vector_buffer<> r_new(tmp_r_file, std::ios::out, 1024*1024, sdsl::bits::hi(r_data.r) + 1); //new version of the grammar

    {
        sdsl::int_vector_buffer<> r_map(map_file, std::ios::out, 1024*1024, sdsl::bits::hi(r_data.r) + 1);
        size_t new_r = 0, next_pos_av = 0, curr_r = 0;
        desc = 0;
        for (size_t j = 0; j < r_old.size(); j++) {
            sym = r_old[j];
            desc = (desc << 1UL) | r_lim[j];

            if (sym != 0 && (sym == j || (desc & 3UL) != 3UL)) {
                r_new.push_back(sym);
                if ((desc & 1UL)) {
                    r_map[curr_r] = new_r;
                    r_lim[next_pos_av] = true;
                    new_r++;
                } else {
                    r_lim[next_pos_av] = false;
                }
                next_pos_av++;
            }
            if (r_lim[j]) curr_r++;
        }
        r_lim.resize(next_pos_av);
        r_data.r = new_r;
        r_map.close();
    }

    sdsl::int_vector<> r_map;
    sdsl::load_from_file(r_map, map_file);
    for(auto && j : r_new) j = r_map[j];
    r_data.g = r_new.size();
    r_new.close();

    if(remove(r_data.rules_file.c_str())){
        std::cout<<"Error trying to remove temporal file"<<std::endl;
        std::cout<<"Aborting"<<std::endl;
        exit(1);
    }
    if(rename(tmp_r_file.c_str(), r_data.rules_file.c_str())){
        std::cout<<"Error trying to rename temporal file"<<std::endl;
        std::cout<<"Aborting"<<std::endl;
        exit(1);
    }
    if(remove(map_file.c_str())){
        std::cout<<"Error trying to remove temporal file"<<std::endl;
        std::cout<<"Aborting"<<std::endl;
        exit(1);
    }
    sdsl::store_to_file(r_lim, r_data.rules_lim_file);
}

template<class sym_type>
std::vector<std::pair<size_t, size_t>> lpg_build::compute_thread_ranges(size_t n_threads, std::string& i_file,
                                                                        sdsl::int_vector<2>& phrase_desc) {
    std::vector<std::pair<size_t, size_t>> thread_ranges;

    i_file_stream<sym_type> is(i_file, BUFFER_SIZE);
    size_t n_chars = is.tot_cells;
    assert(n_chars>0);
    size_t sym_per_thread = INT_CEIL(n_chars, n_threads);
    size_t start, end;

    for(size_t i=0;i<n_threads;i++){
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

void lpg_build::decomp(size_t nt, sdsl::int_vector<>& rules,
                       bv_t& r_lim, bv_t::select_1_type& rlim_ss,
                       bv_t& rem_nt, bv_t::rank_1_type& rm_nt_rs,
                       ivb_t & buff){

    std::stack<size_t> stack;
    stack.push(nt);

    size_t start, end, tmp;
    while(!stack.empty()){

        tmp = stack.top();
        stack.pop();

        if(rem_nt[tmp]){
            start = rlim_ss(tmp)+1;
            end = rlim_ss(tmp+1);
            for(size_t j=end+1;j-->start;){
                stack.push(rules[j]);
            }
        }else{
            buff.push_back(tmp-rm_nt_rs(tmp));
        }
    }
}

void lpg_build::simplify_grammar(lpg_build::plain_grammar_t &p_gram, sdsl::cache_config &config) {

    std::cout<<"  Simplifying the grammar"<<std::endl;

    bv_t r_lim;
    sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
    bv_t::select_1_type r_lim_ss(&r_lim);

    sdsl::int_vector<> rules;
    sdsl::load_from_file(rules, p_gram.rules_file);

    size_t max_tsym = p_gram.symbols_map.back();

    ivb_t new_rules(p_gram.rules_file, std::ios::out);
    sdsl::int_vector_buffer<1> new_r_lim(p_gram.rules_lim_file, std::ios::out);

    bv_t rem_nts(p_gram.r + 1, false);
    {//mark duplicated nonterminals
        sdsl::int_vector<2> rep_nts(p_gram.r + 1, 0);
        for(size_t i=max_tsym+1;i<rules.size();i++){
            if(rep_nts[rules[i]]<2) rep_nts[rules[i]]++;
        }
        //mark the grammar symbols to remove
        for(size_t i=0;i<p_gram.r;i++){
            rem_nts[i] = rep_nts[i]==0 || (rep_nts[i]==1 && i > max_tsym);
        }
    }

    //unmark unique nonterminals that
    // appear in the compressed string
    size_t pos;
    for(size_t i=r_lim_ss(p_gram.r-1)+1; i<rules.size();i++){
        rem_nts[rules[i]] = false;
        //unmark the children of a nonterimnal
        // that appear in the compressed string
        pos = r_lim_ss(rules[i])+1;
        while(!r_lim[pos]){
            rem_nts[rules[pos++]]=false;
        }
        rem_nts[rules[pos]]=false;
    }
    rem_nts[p_gram.r-1] = false;//unmark the compressed string

    bv_t::rank_1_type rem_nts_rs(&rem_nts);
    float rm_per = float(rem_nts_rs(rem_nts.size()))/float(p_gram.r)*100;
    std::cout<<"    "<<rm_per<<"% of the rules will be removed"<<std::endl;

    for(size_t k=0;k<p_gram.sigma;k++){
        new_r_lim.push_back(true);
        new_rules.push_back(k);
    }

    //rules are one-based
    for(size_t i=max_tsym+1,curr_rule=max_tsym+1;i<rules.size();curr_rule++){
        assert(r_lim[i-1]);
        pos = i;
        while(!r_lim[i]) i++;
        i++;

        if(!rem_nts[curr_rule]){
            for(size_t j=pos;j<i;j++){
                if(rem_nts[rules[j]]){
                    decomp(rules[j], rules, r_lim, r_lim_ss, rem_nts, rem_nts_rs, new_rules);
                }else{
                    new_rules.push_back(rules[j]-rem_nts_rs(rules[j]));
                }
            }
            new_r_lim[new_rules.size()-1]=true;
        }
    }

    new_rules.close();
    new_r_lim.close();
    p_gram.r -= rem_nts_rs(rem_nts.size());
    p_gram.g = new_rules.size();
}

void lpg_build::plain_grammar_t::save_to_file(std::string& output_file){

    size_t buffer[255];

    std::ofstream of_stream(output_file, std::ofstream::binary);

    //write number of rules and alphabet size
    buffer[0] = sigma;
    buffer[1] = r;
    buffer[2] = c;
    buffer[3] = g;
    of_stream.write((char *) buffer, sizeof(size_t)*4);

    assert(symbols_map.size()==sigma);
    //write symbols m_map
    for(size_t i=0;i<symbols_map.size();i++){
        buffer[i] = symbols_map[i];
    }
    of_stream.write((char *) buffer, sizeof(size_t)*symbols_map.size());

    buffer[0] = rules_file.size();
    of_stream.write((char *) buffer, sizeof(size_t));
    of_stream.write(rules_file.c_str(), rules_file.size());

    buffer[0] = rules_lim_file.size();
    of_stream.write((char *) buffer, sizeof(size_t));
    of_stream.write(rules_lim_file.c_str(), rules_lim_file.size());

    buffer[0] = lms_as_sp_file.size();
    of_stream.write((char *) buffer, sizeof(size_t));
    of_stream.write(lms_as_sp_file.c_str(), lms_as_sp_file.size());

    of_stream.close();
}

void lpg_build::plain_grammar_t::load_from_file(std::string &g_file){

    size_t buffer[255];
    std::ifstream fp(g_file, std::ifstream::binary);

    fp.read((char *)buffer, sizeof(size_t)*4);

    sigma = buffer[0];
    r = buffer[1];
    c = buffer[2];
    g = buffer[3];

    fp.read((char *)buffer, sizeof(size_t)*sigma);
    symbols_map.resize(sigma);
    for(size_t i=0;i<sigma;i++){
        symbols_map[i] = buffer[i];
    }

    fp.read((char *)buffer, sizeof(size_t));
    auto tmp_file = reinterpret_cast<char *>(malloc(buffer[0]+1));
    fp.read(tmp_file, buffer[0]);
    tmp_file[buffer[0]]='\0';
    rules_file = std::string(tmp_file);

    fp.read((char *)buffer, sizeof(size_t));
    tmp_file = reinterpret_cast<char *>(realloc(tmp_file, buffer[0]+1));
    fp.read(tmp_file, buffer[0]);
    tmp_file[buffer[0]] = '\0';
    rules_lim_file = std::string(tmp_file);

    fp.read((char *)buffer, sizeof(size_t));
    tmp_file = reinterpret_cast<char *>(realloc(tmp_file, buffer[0]+1));
    fp.read(tmp_file, buffer[0]);
    tmp_file[buffer[0]] = '\0';
    lms_as_sp_file = std::string(tmp_file);

    fp.close();
    free(tmp_file);
}
