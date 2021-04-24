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

            //std::cout<<curr_sym<<" - ";
            if(curr_sym==0){
                start = 0;
            }else{
                start = r_lim_ss(curr_sym)+1;
            }

            end = r_lim_ss(curr_sym+1);

            if(r[start] == curr_sym){
                //std::cout<<curr_sym<<std::endl;
                assert((end-start+1)==1);
                tmp_decomp.push_back(curr_sym);
            }else{
                assert((end-start+1)>1);
                for(size_t j=end+1; j-->start;){
                    stack.push(r[j]);
                }
                /*std::cout<<start<<"-"<<end<<" = ";
                for(size_t j=start;j<=end;j++){
                    std::cout<<r[j]<<" ";
                }
                std::cout<<" "<<std::endl;*/
            }
        }

        for(auto const& tmp_sym : tmp_decomp){
            buff_symbol = if_stream.read(pos++);
            //std::cout<<r_data.symbols_map.size()<<std::endl;
            //std::cout<<tmp_sym<<" "<<(int)r_data.symbols_map[tmp_sym]<<" "<<(int)buff_symbol<<std::endl;
            assert(p_gram.symbols_map[tmp_sym] == buff_symbol);
        }
    }

    /*bool greater=false;
    std::h_vec<uint32_t > prev_seq;
    for(auto const& sym : top_rules){
        if(sym!=1) {
            tmp_decomp.clear();
            prev_dc_step.clear();
            decomp=false;

            for (auto const &smaller_sym : rules[sym]) {
                prev_dc_step.push_back(smaller_sym);
            }

            while(!decomp){
                decomp = true;
                for(auto const& tmp_sym : prev_dc_step){
                    if(tmp_sym>=6){
                        for(auto const& smaller_sym : rules[tmp_sym]){
                            if(smaller_sym>=6) decomp = false;
                            tmp_decomp.push_back(smaller_sym);
                        }
                    }else{
                        tmp_decomp.push_back(tmp_sym);
                    }
                }

                if(!decomp){
                    std::swap(tmp_decomp, prev_dc_step);
                    tmp_decomp.clear();
                }
            }

            greater = false;
            if(!prev_seq.empty()){
                bool is_prefix=false;
                for(size_t i=0;i<std::min(tmp_decomp.size(), prev_seq.size());i++){
                    if(tmp_decomp[i]!=prev_seq[i]){
                        is_prefix=true;
                        greater = tmp_decomp[i]>prev_seq[i];
                        break;
                    }
                }
                if(!is_prefix) greater = true;
                assert(greater);
            }
            std::swap(tmp_decomp, prev_seq);
        }
    }*/
    std::cout<<"\tGrammar is correct!!"<<std::endl;
}

//0 symbol doesn't exists
//1 symbol exists but it is unique
//2 symbol appears more than once
sdsl::int_vector<2> lpg_build::compute_alphabet(std::string &i_file){
    //TODO this can be done in parallel if the input is too big
    sdsl::int_vector<2> alphabet(256,0);
    for(size_t i=0;i<256;i++) alphabet[i]=0;

    i_file_stream<uint8_t> if_stream(i_file, BUFFER_SIZE);
    uint8_t tmp_symbol;

    for(size_t i=0;i<if_stream.tot_cells;i++){
        tmp_symbol = if_stream.read(i);
        //READ_SYMBOL(if_stream, i, tmp_symbol);
        if(alphabet[tmp_symbol]<2){
            alphabet[tmp_symbol]++;
        }
    }
    return alphabet;
}

void lpg_build::compute_LPG(std::string &i_file, std::string &p_gram_file, size_t n_threads,
                            sdsl::cache_config &config, size_t hbuff_size, uint8_t sep_symbol) {

    std::cout<<"Computing the alphabet"<<std::endl;
    auto symbol_desc = compute_alphabet(i_file);

    std::string rules_file = sdsl::cache_file_name("rules", config);
    std::string rules_len_file = sdsl::cache_file_name("rules_len", config);
    std::string lmsg_as_sp_file = sdsl::cache_file_name("lmsg_as_sp", config);

    plain_grammar_t p_gram(rules_file, rules_len_file, lmsg_as_sp_file);

    size_t max_symbol=0, min_symbol=0;
    for(size_t i=0;i<symbol_desc.size();i++){
        if(symbol_desc[i]!=0){
            p_gram.symbols_map.push_back(i);
            p_gram.sigma++;
            max_symbol = i;
            symbol_desc[i]--;
        }
    }
    p_gram.r = max_symbol + 1;
    symbol_desc[sep_symbol]+=2;

    ivb_t rules(p_gram.rules_file, std::ios::out, BUFFER_SIZE);
    bvb_t rules_lim(p_gram.rules_lim_file, std::ios::out);
    for(size_t i=0;i<p_gram.r; i++){
        rules.push_back(0);
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
                                           min_symbol, max_symbol,
                                           p_gram, rules, rules_lim,
                                           symbol_desc, config);

    while (rem_phrases > 0) {
        std::cout<<"  Iteration "<<iter++<<std::endl;
        rem_phrases = compute_LPG_int<size_t>(tmp_i_file, output_file,
                                              n_threads, hbuff_size,
                                              min_symbol, max_symbol,
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
        while(read_bytes<tot_bytes){
            c_vec.read((char *) buffer, BUFFER_SIZE);
            read_bytes+=c_vec.gcount();
            assert((c_vec.gcount() % sizeof(size_t))==0);
            for(size_t i=0;i<c_vec.gcount()/sizeof(size_t);i++){
                rules.push_back(min_symbol + buffer[i]);
                rules_lim.push_back(false);
            }
        }
        rules_lim[rules_lim.size() - 1] = true;
        p_gram.r++;
        c_vec.close();
        free(buffer);
    }
    rules.close();
    rules_lim.close();

    std::cout<<"  Collapsing nonterminal rules"<<std::endl;
    collapse_grammar(p_gram, iter, config);

    p_gram.save_to_file(p_gram_file);

    /*
    //std::cout<<"Checking the grammar";
    //check_plain_grammar(r_data, i_file);*/

    if(remove(tmp_i_file.c_str())){
        std::cout<<"Error trying to delete file "<<tmp_i_file<<std::endl;
    }
}

template<class sym_type>
size_t lpg_build::compute_LPG_int(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                                  size_t &min_symbol, size_t &max_symbol, plain_grammar_t &p_gram, ivb_t &rules,
                                  bvb_t &rules_lim, sdsl::int_vector<2> &phrase_desc,
                                  sdsl::cache_config &config) {

    static_map tr_table;
    tr_table.bv = bv_t(max_symbol - min_symbol + 1, false);

    phrase_map_t mp_table(0, "", 0.85);
    size_t rem_phrases=0;

    auto thread_ranges = compute_thread_ranges<sym_type>(n_threads, i_file, phrase_desc);
    size_t loop_alph = max_symbol-min_symbol+1;

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
        threads_data.emplace_back(i_file, tmp_o_file,
                                  mp_table, tr_table,
                                  range.first, range.second,
                                  loop_alph, hb_bytes, tmp_addr + (k*hb_bytes), phrase_desc);
        //std::cout<<uintptr_t(tmp_addr)<<" "<<k*hb_bytes<<" "<<"-"<<((k+1)*hb_bytes)-1<<" "<<hbuff_size<<std::endl;
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

    size_t nnt=0, tr=0, psize=0;//<- for the iter stats
    if(mp_table.size()>0){

        sdsl::util::init_support(tr_table.bv_rs, &tr_table.bv);

        size_t width = sdsl::bits::hi(max_symbol-min_symbol+1)+1;
        const bitstream<buff_t>& stream = mp_table.get_data();
        key_wrapper key_w{width, mp_table.description_bits(), stream};

        //temporal unload of the hash table (not the data)
        std::string st_table = sdsl::cache_file_name("ht_data", config);
        mp_table.unload_table(st_table);

        //rename phrases according to their lexicographical ranks
        std::cout<<"    Assigning identifiers to the phrases"<<std::endl;
        assign_ids(mp_table, tr_table, min_symbol, max_symbol, p_gram, key_w, rules, rules_lim,
                   n_threads, config);

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
            sdsl::int_vector<2> new_phrase_desc(max_symbol-min_symbol+1, false);

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
                    rem_phrases++;
                }

                //read the (reversed) last symbol
                sym = key_w.read(*it, 0);
                if (phrase_desc[sym] & 2U) {//phrase is suffix of some string
                    tmp_value += 2;
                }

                new_phrase_desc[val >> 1UL] = tmp_value;
                ++it;
            }

            for (size_t i = 0; i < tr_table.bv.size(); i++) {
                if (tr_table.bv[i]) {
                    new_phrase_desc[tr_table[i]] = phrase_desc[i];
                }
            }
            phrase_desc.swap(new_phrase_desc);
        }
        nnt = mp_table.size();
        tr = tr_table.size();
    }else{ //just copy the input

        sdsl::util::clear(phrase_desc);
        sdsl::util::clear(tr_table.bv);
        sdsl::util::clear(tr_table.map_vector);

        std::ifstream in(i_file, std::ios_base::binary);
        std::ofstream out(o_file, std::ios_base::binary);

        auto buffer = reinterpret_cast<char*>(malloc(BUFFER_SIZE));
        do {
            in.read(&buffer[0], BUFFER_SIZE);
            out.write(&buffer[0], in.gcount());
        } while (in.gcount() > 0);
        free(buffer);

        //remove remaining files
        for(size_t i=0;i<threads_data.size();i++){
            std::string tmp_file =  threads_data[i].ofs.file;
            if(remove(tmp_file.c_str())){
                std::cout<<"Error trying to delete file "<<tmp_file<<std::endl;
            }
        }
    }

    std::cout<<"    Iter. stats:"<<std::endl;
    std::cout<<"      Parse size:          "<<psize<<std::endl;
    std::cout<<"      New nonterminals:    "<<nnt<<std::endl;
    std::cout<<"      Transferred symbols: "<<tr<<std::endl;

    return rem_phrases;
}

void
lpg_build::assign_ids(phrase_map_t &mp_map, static_map &tr_table, size_t &min_gsym, size_t &max_gsym, plain_grammar_t &r_data,
                      key_wrapper &key_w, ivb_t &r, bvb_t &r_lim, size_t n_threads, sdsl::cache_config &config) {

    //some aliases to make the code more readable
    size_t &n_rules = r_data.r;

    //sort phrases
    /*std::vector<size_t> k_list; //list of keys (hash_table_offset, key bits)
    k_list.reserve(mp_map.size());

    for (auto const &phrase : mp_map) {
        k_list.emplace_back(phrase);
    }
    std::sort(k_list.begin(), k_list.end(),
              [key_w](auto &left, auto &right) {
        return key_w.compare(left, right);
    });*/

    //TODO this is new
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
    parallel_str_sort(syms_file, compare, access, max_gsym-min_gsym+1, n_threads, config);

    sdsl::int_vector<> k_list;
    sdsl::load_from_file(k_list, syms_file);

    if(remove(syms_file.c_str())){
        std::cout<<"Error trying to remove file "<<syms_file<<std::endl;
    }
    //

    size_t idx_uniq=0;
    for(size_t i=0;i<tr_table.bv.size();i++){
        if(tr_table.bv[i]){
            idx_uniq=i;
            break;
        }
    }

    tr_table.map_vector.width(sdsl::bits::hi(mp_map.size()+ tr_table.size())+1);
    tr_table.map_vector.resize(tr_table.size());

    n_rules += mp_map.size() + tr_table.size();

    size_t m_pos=0, s_pos=0, rank=0;
    while(s_pos < tr_table.size() && m_pos < mp_map.size()){

        size_t len = key_w.size(k_list[m_pos]);
        size_t first_sym = key_w.read(k_list[m_pos], len-1);

        if(idx_uniq < first_sym){
            r.push_back(min_gsym+idx_uniq);
            r_lim.push_back(true);

            tr_table.map_vector[tr_table.bv_rs(idx_uniq)] = rank;

            idx_uniq++, s_pos++;
            while(idx_uniq < tr_table.bv.size() && !tr_table.bv[idx_uniq]) idx_uniq++;
        }else{

            for(size_t i=len; i-- > 1;){
                r.push_back(min_gsym + key_w.read(k_list[m_pos], i));
                r_lim.push_back(false);
            }
            r.push_back(min_gsym + key_w.read(k_list[m_pos], 0));
            r_lim.push_back(true);

            //modify the key value
            phrase_map_t::val_type val=0;
            mp_map.get_value_from(k_list[m_pos], val);
            val |= rank<<1UL;
            mp_map.insert_value_at(k_list[m_pos], val);
            //

            m_pos++;
        }
        rank++;
    }

    while(s_pos < tr_table.size()){
        r.push_back(min_gsym+idx_uniq);
        r_lim.push_back(true);

        tr_table.map_vector[tr_table.bv_rs(idx_uniq)] = rank;

        idx_uniq++; rank++; s_pos++;
        while(idx_uniq < tr_table.bv.size() && !tr_table.bv[idx_uniq]) idx_uniq++;
    }

    while(m_pos < mp_map.size()){

        size_t len = key_w.size(k_list[m_pos]);

        for(size_t i=len; i-- > 1;){
            r.push_back(min_gsym + key_w.read(k_list[m_pos], i));
            r_lim.push_back(false);
        }
        r.push_back(min_gsym + key_w.read(k_list[m_pos], 0));
        r_lim.push_back(true);

        //modify the key value
        phrase_map_t::val_type val=0;
        mp_map.get_value_from(k_list[m_pos], val);
        val |= rank<<1UL;
        mp_map.insert_value_at(k_list[m_pos], val);
        //

        rank++; m_pos++;
    }

    min_gsym = max_gsym+1;
    max_gsym += mp_map.size() + tr_table.size();
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

    //temporal data structures used by the thread
    //std::string dump_file = lms_data->ofs.file + "_phrases";
    //string_map_t tmp_m_map(lms_data->buff_size, dump_file);
    //

    bool s_type, prev_s_type = S_TYPE;
    sym_t curr_sym, prev_sym;

    string_t curr_lms(2, lms_data->sym_width);

    //assert(lms_data->start==0 || lms_data->is_suffix(lms_data->ifs.read(lms_data->start - 1)));

    prev_sym = lms_data->ifs.read(lms_data->end);
    //READ_SYMBOL(lms_data->ifs, lms_data->end, prev_sym)

    //assert(lms_data->is_suffix(prev_sym));
    curr_lms.push_back(prev_sym);

    for(size_t i = lms_data->end; i-- > lms_data->start;){

        curr_sym = lms_data->ifs.read(i);
        //READ_SYMBOL(lms_data->ifs, i, curr_sym)

        //                                     L_TYPE   S_TYPE*
        //                                        ---- ----
        //this is a junction between two strings = ...$ $...
        if(lms_data->is_suffix(curr_sym)){
            if(!curr_lms.empty()){
                lms_data->hash_phrase(curr_lms);
                curr_lms.clear();
            }
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
    if(!curr_lms.empty()){
        lms_data->hash_phrase(curr_lms);
    }


    pthread_mutex_lock(&thread_mutex);
    for(size_t i=0;i<lms_data->thread_tr_bv.size();i++){
        if(lms_data->thread_tr_bv[i]){
            lms_data->tr_map.bv[i] = true;
        }
    }
    pthread_mutex_unlock(&thread_mutex);

    lms_data->thread_map.flush();
    sdsl::util::clear(lms_data->thread_tr_bv);

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

    //assert(lms_data->start==0 || lms_data->is_suffix(lms_data->ifs.read(lms_data->start - 1)));

    prev_sym = lms_data->ifs.read(lms_data->end);
    //READ_SYMBOL(lms_data->ifs, lms_data->end, prev_sym)

    //assert(lms_data->is_suffix(prev_sym));

    curr_lms.push_back(prev_sym);

    for(size_t i = lms_data->end; i-- > lms_data->start;){

        curr_sym = lms_data->ifs.read(i);
        //READ_SYMBOL(lms_data->ifs, i, curr_sym)

        //                                     L_TYPE   S_TYPE*
        //                                        ---- ----
        //this is a junction between two strings = ...$ $...
        if(lms_data->is_suffix(curr_sym)){
            if(!curr_lms.empty()){
                lms_data->store_phrase(curr_lms);
                curr_lms.clear();
            }
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

    //TODO : checking the grammar
    //sdsl::util::assign(r_lim_ss, &r_lim);
    //check_plain_grammar(r, r_lim_ss, r_data, "../tests/plain_reads.txt");
    //
    //sdsl::store_to_file(r, r_data.rules_file);

    sdsl::store_to_file(r_lim, r_data.rules_lim_file);
}

template<class sym_type>
std::vector<std::pair<size_t, size_t>> lpg_build::compute_thread_ranges(size_t n_threads, std::string& i_file, sdsl::int_vector<2>& phrase_desc) {
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
