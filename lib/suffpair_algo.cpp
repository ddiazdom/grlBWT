//
// Created by Diaz, Diego on 24.11.2021.
//
#include "suffpair_algo.hpp"

//global variables for multi threading
std::mutex m;
std::condition_variable cv;
bool proceed;
bool merged;
size_t processed;
size_t rem_threads;
//

//create a set of new m_rules from the hash table
void create_new_rules(phrase_map_t& ht, pairing_data& p_data){
    key_wrapper key_w{p_data.s_width, ht.description_bits(), ht.get_data()};
    phrase_map_t::val_type val=0;
    //size_t first_pos = p_data.new_rules.size();
    for(auto elm : ht){
        ht.get_value_from(elm, val);
        if(val==1){//repeated pair but new
            val |= p_data.next_av_rule<<1UL;
            ht.insert_value_at(elm, val);
            p_data.new_rules.push_back(key_w.read(elm, 0));
            p_data.new_rules.push_back(key_w.read(elm, 1));
            p_data.new_rules.push_back(key_w.read(elm, 2));
            p_data.next_av_rule++;
        }
    }
}

//insert the suffix pair m_rules into the grammar
void update_grammar(pairing_data& p_data, gram_info_t& gram){

    i_file_stream<size_t> rules(p_data.r_file, BUFFER_SIZE);

    std::string tr_file = sdsl::cache_file_name("col_rules", p_data.config);
    ivb_t col_rules(tr_file, std::ios::out);

    //collapse compressed symbols
    size_t n_av=0, tmp_sym;
    for(size_t i=0; i < p_data.gsyms; i++){
        if(i<=gram.max_tsym){
            p_data.r_lim[n_av++] = true;
            col_rules.push_back(i);
        }else{
            tmp_sym = rules.read(i);
            if(tmp_sym != p_data.lim_id) {
                col_rules.push_back(tmp_sym);
                p_data.r_lim[n_av++] = false;
            }

            if(p_data.r_lim[i]){
                assert(!p_data.r_lim[n_av-1]);
                p_data.r_lim[n_av-1] = true;
            }
        }
    }

    //assert we actually compress the grammar
    //assert((n_av+p_data.new_rules.size()/2)<p_data.r_lim.size());

    size_t lsym, rsym, freq;
    for(size_t i=0;i<p_data.new_rules.size();i+=3){
        lsym = p_data.new_rules[i];
        rsym = p_data.new_rules[i+1];
        freq = p_data.new_rules[i+2];
        for(size_t j=0;j<freq;j++){
            col_rules.push_back(lsym);
            p_data.r_lim[n_av++] = false;
        }
        col_rules.push_back(rsym);
        p_data.r_lim[n_av++] = true;
    }

    //put array C at the end of the new m_rules
    for(size_t i=p_data.gsyms; i < rules.size(); i++){
        col_rules.push_back(rules.read(i));
        p_data.r_lim[n_av++] = false;
    }
    p_data.r_lim[n_av-1] = true;
    p_data.r_lim.resize(n_av);
    sdsl::store_to_file(p_data.r_lim, gram.rules_lim_file);

    gram.sp_rules.first = gram.r-1;
    gram.sp_rules.second = p_data.new_rules.size() / 3;

    gram.r += p_data.new_rules.size()/3;

    col_rules.close();
    rules.close();
    rename(col_rules.filename().c_str(), gram.rules_file.c_str());

    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Grammar size before:        "<<gram.g<< std::endl;
    std::cout<<"      Grammar size after:         "<<n_av<<std::endl;
    std::cout<<"      Number of new nonterminals: "<<p_data.new_rules.size()/3<<std::endl;
    std::cout<<"      Compression ratio:          "<<double(n_av)/double(gram.g) << std::endl;
    gram.g = n_av;
}

//change the width of R and compute the repeated symbols
void prepare_input(gram_info_t& gram_info, bv_t& rep_syms, sdsl::cache_config& config){

    ivb_t rules(gram_info.rules_file, std::ios::in);
    sdsl::int_vector_buffer<1> rules_lim(gram_info.rules_lim_file, std::ios::in);

    std::string tmp_string = sdsl::cache_file_name("tmp_rules", config);
    o_file_stream<size_t> tmp_r(tmp_string,  BUFFER_SIZE, std::ios::out);
    sdsl::int_vector<2> tmp_rep(gram_info.r - 1, 0);
    size_t rule_id=gram_info.max_tsym+1, sym;

    for(size_t i=0;i<=gram_info.max_tsym;i++) tmp_r.push_back(rules[i]);
    for (size_t i=gram_info.max_tsym+1;i<rules.size();i++){
        sym = rules[i];
        tmp_r.push_back(sym);
        if(rule_id<gram_info.rules_breaks[gram_info.n_p_rounds] && tmp_rep[sym]<2) tmp_rep[sym]++;
        if(rules_lim[i]) rule_id++;
    }
    rules.close(true);
    tmp_r.close();

    size_t j=0;
    for(auto const& i : tmp_rep){
        if(i>1) rep_syms[j] = true;
        j++;
    }
    rename(tmp_string.c_str(), gram_info.rules_file.c_str());
}

void collect_pairs(thread_data* d, i_file_stream<size_t>& p_list, o_file_stream<size_t>& r){

    size_t id, val, pos, lsym, rsym, freq;
    string_t tmp_pair(3, d->p_data.s_width);
    phrase_map_t ht(d->hb_size, d->ht_file, 0.8, d->hb_addr);

    for(size_t i=0;i<p_list.size();i++) {

        pos = p_list.read(i);
        lsym = r.read(pos);
        rsym = r.read(pos+1);

        assert(pos>=d->start && pos+1<=d->end);

        pos--;
        freq=1;
        while(!d->p_data.r_lim[pos] && r.read(pos)==lsym){
            pos--;
            freq++;
        }

        tmp_pair.write(0, lsym);
        tmp_pair.write(1, rsym);
        tmp_pair.write(2, freq);

        //check if the pair covers an entire rule
        id = d->p_data.r_lim[pos] ? (d->p_data.r_lim_rs(pos+1)<<1UL) : 0;
        id |= freq>1;

        auto res = ht.insert(tmp_pair.data(), tmp_pair.n_bits(), id);
        if(!res.second){//the pair already exists
            val = 0;
            ht.get_value_from(res.first, val);
            val |= (id | 1UL);
            ht.insert_value_at(res.first, val);
        }
    }
    ht.flush();
}

//replace the new pairs in R with their identifiers
void replace_pairs(const phrase_map_t& ht, const pairing_data& p_data, std::string& pl_file, o_file_stream<size_t>& r){

    i_file_stream<size_t> p_list(pl_file, BUFFER_SIZE);

    //temporal file to store the new suffix positions
    std::string npl_file = pl_file+"_tmp";
    o_file_stream<size_t> np_list(npl_file,  BUFFER_SIZE, std::ios::out);

    //temporal variables
    string_t tmp_pair(3, p_data.s_width);
    phrase_map_t::val_type  val=0;

    //update the suffixes with the new symbols
    size_t lsym, rsym;
    long long pos, freq;
    for(size_t i=0; i<p_list.size();i++) {

        pos = (long long)p_list.read(i);
        lsym = r.read(pos);
        rsym = r.read(pos+1);

        freq=1;
        size_t tmp_pos=i-1;
        tmp_pair.write(0, lsym);
        tmp_pair.write(1, rsym);
        while(!p_data.r_lim[tmp_pos] && r.read(tmp_pos)==lsym){
            tmp_pair.write(2, freq);
            auto res = ht.find(tmp_pair.data(), tmp_pair.n_bits());
            val = 0;
            ht.get_value_from(res.first, val);
            if(val & 1UL){
            }

            freq++;
            tmp_pos--;
        }

        /*if(!p_data.r_lim[tmp_pos]){ //the suffix can't span a complete rule

            tmp_pair.write(0, lsym);
            tmp_pair.write(1, rsym);
            tmp_pair.write(2, freq);

            auto res = ht.find(tmp_pair.data(), tmp_pair.n_bits());
            ht.get_value_from(res.first, val);

            if(val & 1UL){
                r.write(pos+1, val>>1UL);
                for(long long j=1;j<=freq;j++){
                    r.write(pos+1+j, p_data.lim_id); //<-dummies indicating a replaced suffix
                }

                if(p_data.rep_sym[r.read(pos)]){
                    np_list.push_back(pos);
                }
            }
        }*/
    }
    //update the list with the suffixes
    p_list.close();
    np_list.close();
    if(remove(p_list.filename().c_str())){
        std::cout<<"Error trying to remove file "<<p_list.filename()<<std::endl;
    }
    rename(npl_file.c_str(), p_list.filename().c_str());
}

void * suffpair_thread(void * data) {

    {
        std::unique_lock<std::mutex> lk(m);
        cv.wait(lk, [] { return proceed; });
    }

    auto d = reinterpret_cast<thread_data *>(data);
    o_file_stream<size_t> r(d->p_data.r_file, BUFFER_SIZE, std::ios::in | std::ios::out);
    {
        phrase_map_t ht(d->hb_size, d->ht_file, 0.8, d->hb_addr);
        string_t tmp_pair(3, d->p_data.s_width);
        size_t val, freq;

        //position in R of the suffixes of length two
        o_file_stream<size_t> p_list(d->pl_chunk_file, BUFFER_SIZE, std::ios::out);

        size_t start = d->start;
        size_t i = d->end;
        size_t lsym, rsym, id, tmp_pos;
        bool valid_suffix, rep_suffix;

        //skip the terminal symbols
        if (start == 0) while (d->p_data.r_lim[start]) start++;

        while (i > start) {

            assert(d->p_data.r_lim[i]);

            lsym = r.read(i-1);
            rsym = r.read(i);

            //new repeated suffix -> pattern 01 with both symbols appearing more than once in R
            valid_suffix = !d->p_data.r_lim[i-1] && d->p_data.r_lim[i];
            rep_suffix = d->p_data.rep_sym[lsym] && d->p_data.rep_sym[rsym];

            if(valid_suffix && rep_suffix) {

                freq=1;
                tmp_pos=i-1;
                tmp_pair.write(0, lsym);
                tmp_pair.write(1, rsym);

                while(!d->p_data.r_lim[tmp_pos] && r.read(tmp_pos)==lsym){
                    tmp_pair.write(2, freq);
                    id = d->p_data.r_lim[tmp_pos-1] ? (d->p_data.r_lim_rs(tmp_pos)<< 1UL) : 0;
                    auto res = ht.insert(tmp_pair.data(), tmp_pair.n_bits(), id);
                    if (!res.second) {
                        val = 0;
                        ht.get_value_from(res.first, val);
                        val |= (id | 1UL);
                        ht.insert_value_at(res.first, val);
                    }
                    freq++;
                    tmp_pos--;
                }
                p_list.push_back(i-1);
            } else {
                i--;
            }
            while(!d->p_data.r_lim[i]) i--;
        }
        ht.flush();
        p_list.close();
    }

    {//notify this thread finished getting the pairs
        // in its range, and it is waiting for the global merge
        // of the pairs
        std::unique_lock<std::mutex> lck(m);
        processed++;
        lck.unlock();
        cv.notify_all();

        lck.lock();
        cv.wait(lck, []{return merged;});
    }

    replace_pairs(d->ht, d->p_data, d->pl_chunk_file, r);
    r.flush();
    i_file_stream<size_t> p_list(d->pl_chunk_file, BUFFER_SIZE);

    {//notify the main thread we finish the pair replacement and wait
        // for the signal to start a new loop
        std::unique_lock<std::mutex> lck(m);
        processed++;
        if(p_list.size()==0) rem_threads--;
        lck.unlock();
        cv.notify_all();

        lck.lock();
        cv.wait(lck, []{return proceed;});
    }

    while(p_list.size()>0){

        //collect the pairs in the thread range
        collect_pairs(d, p_list, r);
        p_list.close();

        {//wait for the main thread to gather the new pairs
            std::unique_lock<std::mutex> lck(m);
            processed++;
            lck.unlock();
            cv.notify_all();

            lck.lock();
            cv.wait(lck, []{return merged;});
        }

        //replace the new pairs in the thread range
        replace_pairs(d->ht, d->p_data, d->pl_chunk_file, r);
        r.flush();
        p_list = i_file_stream<size_t>(d->pl_chunk_file, BUFFER_SIZE);

        {//notify the main thread we finish the loop, and we are waiting
            // for the signal to start a new loop
            std::unique_lock<std::mutex> lck(m);
            processed++;
            if(p_list.size()==0){
                rem_threads--;
                d->finished = true;
            }
            lck.unlock();
            cv.notify_all();

            lck.lock();
            cv.wait(lck, []{return proceed;});
        }
    }
    p_list.close();
    r.close();

    if(remove(p_list.filename().c_str())){
        std::cout<<"Error trying to remove "<<p_list.filename()<<std::endl;
    }
    pthread_exit(nullptr);
}

void merge_ht_data(std::vector<thread_data>& t_data) {

    //collect the pairs
    size_t id;
    size_t k=0;
    for(auto & data : t_data){

        if(data.finished) continue;

        std::ifstream text_i(data.ht_file, std::ios_base::binary);

        text_i.seekg (0, std::ifstream::end);
        size_t tot_bytes = text_i.tellg();
        text_i.seekg (0, std::ifstream::beg);

        auto k_buffer = reinterpret_cast<char *>(malloc(tot_bytes));

        text_i.read(k_buffer, tot_bytes);
        assert(text_i.good());

        bitstream<phrase_map_t::buff_t> bits;
        bits.stream = reinterpret_cast<phrase_map_t::buff_t*>(k_buffer);

        size_t next_bit = 32;
        size_t tot_bits = tot_bytes*8;
        size_t key_bits;
        void* key=malloc(sizeof(size_t)*2);
        memset(key, 0, sizeof(size_t)*2);

        while(next_bit<tot_bits){

            key_bits = bits.read(next_bit-32, next_bit-1);

            char *tmp = reinterpret_cast<char*>(key);
            tmp[INT_CEIL(key_bits, 8)-1] = 0;

            bits.read_chunk(key, next_bit, next_bit+key_bits-1);
            next_bit+=key_bits;
            id = bits.read(next_bit, next_bit+data.ht.value_bits()-1);
            next_bit+=data.ht.value_bits()+32;

            auto res = data.ht.insert(key, key_bits, id);

            if(!res.second){
                size_t val = 0;
                data.ht.get_value_from(res.first, val);
                val |= (id | 1UL);
                data.ht.insert_value_at(res.first, val);
            }
        }
        text_i.close();

        if(remove(data.ht_file.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }
        k++;
        free(key);
        free(k_buffer);
    }
    t_data[0].ht.shrink_databuff();
}

void suffixpair_int(pairing_data& p_data) {

    std::vector<thread_data> threads_data;
    std::vector<pthread_t> threads(p_data.n_threads);
    size_t rules_per_thread = INT_CEIL(p_data.tot_lms_rules, p_data.n_threads);
    size_t pos=0, tmp=0, tmp_rules, start, end;
    phrase_map_t ht;
    rem_threads = p_data.n_threads;

    //adjust the hash table buffer for every thread
    size_t buff_cells = p_data.hbuff_size/sizeof(size_t);
    size_t hb_bytes = (buff_cells / p_data.n_threads) * sizeof(size_t);
    void *buff_addr = malloc(p_data.hbuff_size);
    auto tmp_addr = reinterpret_cast<char*>(buff_addr);

    for(size_t i=0;i<p_data.n_threads;i++) {
        start = pos;
        tmp_rules = std::min<size_t>(((i + 1) * rules_per_thread), p_data.tot_lms_rules);
        while(tmp<tmp_rules){
            if(p_data.r_lim[pos++]){
                tmp++;
            }
        }
        end = pos-1;
        threads_data.emplace_back(p_data, start, end, ht, hb_bytes, tmp_addr + (i*hb_bytes));
    }

    for(size_t i=0;i<p_data.n_threads;i++){
        int ret =  pthread_create(&threads[i],
                                  nullptr,
                                  &suffpair_thread,
                                  (void*)&threads_data[i]);
        if(ret != 0) {
            printf("Error: pthread_create() failed\n");
            exit(EXIT_FAILURE);
        }
    }

    size_t act_threads;
    {//indicate the workers to start
        std::unique_lock<std::mutex> lk(m);
        proceed = true;
        processed = 0;
        merged = false;
        rem_threads = p_data.n_threads;
        act_threads = rem_threads;
    }
    cv.notify_all();

    while(act_threads>0) {

        {// wait for the workers
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [act_threads] { return processed == act_threads;});
            free(buff_addr);
        }

        merge_ht_data(threads_data);
        create_new_rules(ht, p_data);

        {// tell the worker they can replace the pairs by
            // their values
            std::unique_lock<std::mutex> lk(m);
            merged = true;
            processed = 0;
            proceed = false;
        }
        cv.notify_all();

        {// wait for the workers to replace the pairs with their symbols
            std::unique_lock<std::mutex> lk(m);
            cv.wait(lk, [act_threads] { return processed == act_threads;});
        }

        {// update the global variables and indicate the workers
            // they can continue with the next loop
            std::unique_lock<std::mutex> lk(m);
            processed = 0;
            merged = false;
            proceed = true;
            act_threads = rem_threads;
            ht.reset();
            buff_addr = malloc(p_data.hbuff_size);
            tmp_addr = reinterpret_cast<char*>(buff_addr);
            for(size_t i=0;i<p_data.n_threads;i++){
                threads_data[i].hb_addr = tmp_addr+(i*threads_data[i].hb_size);
            }
        }
        cv.notify_all();
    }

    for(size_t i=0;i<threads_data.size();i++) {
        pthread_join(threads[i], nullptr);
    }
}


size_t create_new_rules(phrase_map_t& ht, gram_info_t& p_gram, bv_t& rep_sym, sdsl::cache_config& config){

    std::string new_rules_file = sdsl::cache_file_name("new_rules", config);
    std::string new_rules_lim_file = sdsl::cache_file_name("new_rules_lim_file", config);
    ivb_t new_rules(new_rules_file, std::ios::out);
    sdsl::int_vector_buffer<1> new_r_lim(new_rules_lim_file, std::ios::out);

    key_wrapper key_w{sdsl::bits::hi(p_gram.r)+1, ht.description_bits(), ht.get_data()};
    string_t phrase(2,sdsl::bits::hi(p_gram.r)+1);
    string_t suffix(2,sdsl::bits::hi(p_gram.r)+1);

    size_t next_av_rule = p_gram.r-1, val;
    for(auto elm : ht){
        val = 0;
        ht.get_value_from(elm, val);
        if(val==1) { //repeated pair but new
            val |= next_av_rule<<1UL;
            ht.insert_value_at(elm, val);
            next_av_rule++;
        }
    }

    bool found;
    for(auto elm : ht){
        val = 0;
        ht.get_value_from(elm, val);
        if((val>>1UL)>=(p_gram.r-1)){ //repeated pair but new

            //create the new rules
            phrase.clear();
            for(size_t i=key_w.size(elm);i-->0;){
                phrase.push_back(key_w.read(elm, i));
            }

            suffix.clear();
            for(size_t i=phrase.size()-1; rep_sym[phrase[i]] && i>0;i--){
                suffix.push_back(phrase[i]);
            }

            found = false;
            while(suffix.size()>1) {
                suffix.mask_tail();
                auto res = ht.find(suffix.data(), suffix.n_bits());
                assert(res.second);
                val = 0;
                ht.get_value_from(res.first, val);
                if(val & 1UL){
                    size_t prefix = phrase.size()-suffix.size();
                    for(size_t i=0;i<prefix;i++){
                        new_rules.push_back(phrase[i]);
                        new_r_lim.push_back(false);
                    }
                    new_rules.push_back(val>>1UL);
                    new_r_lim.push_back(true);

                    found = true;
                    break;
                }
                suffix.pop_back();
            }

            if(!found){
                for(size_t i=0;i<phrase.size();i++){
                    new_rules.push_back(phrase[i]);
                    new_r_lim.push_back(false);
                }
                new_r_lim[new_r_lim.size()-1] = true;
            }
        }
    }
    new_rules.close();
    new_r_lim.close();
    return next_av_rule - (p_gram.r-1);
}

void compress_grammar(gram_info_t& p_gram, vector_t& rules, bv_t& r_lim, bv_t& rep_sym, phrase_map_t& ht, sdsl::cache_config &config) {

    size_t n_new_rules = create_new_rules(ht, p_gram, rep_sym, config);

    /*std::string col_file = sdsl::cache_file_name("col_rules", config);
    std::string col_lim_file = sdsl::cache_file_name("col_lim_file", config);

    ivb_t col_rules(col_file, std::ios::out);
    sdsl::int_vector_buffer<1> col_r_lim(col_lim_file, std::ios::out);

    for(size_t i=0;i<=p_gram.max_tsym;i++){
        col_rules.push_back(i);
        col_r_lim.push_back(r_lim[i]);
    }*/

    string_t rule(2, sdsl::bits::hi(p_gram.r)+1);
    size_t nt=p_gram.max_tsym+1, pos=p_gram.max_tsym+1, new_pos=p_gram.max_tsym+1, start, end, id, val, prefix;
    bool found;
    while(nt<p_gram.r-1) {

        /*if(p_gram.is_rl(nt)){
            col_rules.push_back(rules[pos]);
            col_rules.push_back(rules[pos+1]);

            col_r_lim.push_back(r_lim[pos]);
            col_r_lim.push_back(r_lim[pos+1]);
            pos+=2;
            nt++;
        }*/

        assert(r_lim[pos-1]);

        start = pos;
        while(!r_lim[pos]) pos++;
        end = pos;

        rule.clear();
        for(size_t i=end; rep_sym[rules[i]] && i>start;i--){
            rule.push_back(rules[i]);
        }

        found = false;
        while(rule.size()>1){

            rule.mask_tail();
            auto res = ht.find(rule.data(), rule.n_bits());
            if(res.second){
                val = 0;
                ht.get_value_from(res.first, val);
                if(val & 1UL){

                    id = val>>1UL;
                    prefix = (end-start+1) - rule.size();
                    for(size_t i=0;i<prefix;i++){
                        rules[new_pos] = rules[start+i];
                        r_lim[new_pos++] = false;
                    }
                    rules[new_pos] = id;
                    r_lim[new_pos++] = true;

                    /*std::cout<<"rule :";
                    for(size_t i=start;i<=end;i++){
                        std::cout<<rules[i]<<" ";
                    }
                    std::cout<<"\nsuffix :";
                    for(size_t i=rule.size();i-->0;){
                        std::cout<<rule[i]<<" ";
                    }
                    std::cout<<" -> "<<id<<std::endl;
                    std::cout<<"replacement: ";
                    for(size_t i=0;i<=prefix;i++){
                        std::cout<<col_rules[start+i]<<" ";
                    }
                    std::cout<<" "<<std::endl;*/
                    found = true;
                    break;
                }
            }
            rule.pop_back();
            /*for(size_t i=rule.size();i-->0;){
                std::cout<<rule[i]<<" ";
            }
            std::cout<<" "<<id<<std::endl;*/
        }

        if(!found){
            for(size_t i=start;i<=end;i++){
                rules[new_pos] = rules[i];
                r_lim[new_pos++] = r_lim[i];
            }
        }

        assert(r_lim[pos]);
        nt++;
        pos++;
    }

    std::string new_rules_file = sdsl::cache_file_name("new_rules", config);
    std::string new_rules_lim_file = sdsl::cache_file_name("new_rules_lim_file", config);
    ivb_t new_rules(new_rules_file, std::ios::in);
    sdsl::int_vector_buffer<1> new_r_lim(new_rules_lim_file, std::ios::in);
    for(size_t i=0;i<new_rules.size();i++){
        rules[new_pos] = new_rules[i];
        r_lim[new_pos++] = new_r_lim[i];
    }
    new_rules.close(true);
    new_r_lim.close(true);

    assert(r_lim[pos-1]);
    for(size_t i=pos;i<rules.size();i++){
        rules[new_pos] = rules[i];
        r_lim[new_pos++] = false;
    }
    r_lim[new_pos-1] = true;

    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Grammar size before:        "<<p_gram.g<< std::endl;
    std::cout<<"      Grammar size after:         "<<new_pos<<std::endl;
    std::cout<<"      Number of new nonterminals: "<<n_new_rules<<std::endl;
    std::cout<<"      Compression ratio:          "<<double(new_pos)/double(rules.size()) << std::endl;

    rules.resize(new_pos);
    r_lim.resize(new_pos);
    p_gram.sp_rules.first = p_gram.r-1;
    p_gram.sp_rules.second = n_new_rules;
    p_gram.r += n_new_rules;
    p_gram.g = rules.size();

    sdsl::store_to_file(rules, p_gram.rules_file);
    sdsl::store_to_file(r_lim, p_gram.rules_lim_file);
}

void mark_repeated_nt(bv_t& rep_sym, vector_t& rules, bv_t& r_lim, gram_info_t& p_gram){

    sdsl::int_vector<2> tmp_rep(p_gram.r, 0);
    size_t rule_id=p_gram.max_tsym+1, sym;
    for (size_t i=p_gram.max_tsym+1;i<rules.size();i++){
        sym = rules[i];
        if(rule_id<p_gram.rules_breaks[p_gram.n_p_rounds] && tmp_rep[sym]<2) tmp_rep[sym]++;
        if(r_lim[i]) rule_id++;
    }

    size_t j=0;
    for(auto const& i : tmp_rep){
        if(i>1){
            rep_sym[j] = true;
        }
        j++;
    }
}

void suffpair(gram_info_t& p_gram, sdsl::cache_config &config, size_t n_threads, size_t hbuff_size) {

    std::cout<<"  Suffpair algorithm"<<std::endl;

    string_t rule(2, sdsl::bits::hi(p_gram.r)+1);
    sdsl::int_vector<> rules;
    sdsl::load_from_file(rules, p_gram.rules_file);
    bv_t r_lim;
    sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
    bv_t rep_sym(p_gram.r, false);
    mark_repeated_nt(rep_sym, rules, r_lim, p_gram);

    phrase_map_t ht;
    size_t end=rules.size(), j, id, full_size, nt = p_gram.rules_breaks.back()-1, tmp_nt=p_gram.r, val;
    do{
        if(r_lim[--end]) tmp_nt--;
    } while(tmp_nt>nt);

    j=end;
    while(nt>p_gram.max_tsym){

        assert(r_lim[j]);
        if(!rep_sym[rules[j]]){
            j--;
            goto next_rule;
        }

        rule.clear();
        rule.push_back(rules[j--]);
        while(!r_lim[j] && rep_sym[rules[j]]){
            rule.push_back(rules[j--]);
        }
        full_size = r_lim[j];

        //size_t k=0;
        //while(rep_sym[rule[k]]) k++;
        //while(rule.size()>k) rule.pop_back();

        /*for(size_t i=rule.size();i-->0;){
            std::cout<<rule[i]<<" ";
        }
        std::cout<<" -> "<<std::endl;*/
        while(rule.size()>1){

            id = full_size ? nt<<1UL : 0;
            rule.mask_tail();
            auto res = ht.insert(rule.data(), rule.n_bits(), id);

            if(!res.second){
                val = 0;
                ht.get_value_from(res.first, val);
                val |= (id | 1UL);

                /*if(rule.size()==3 && rule[0]==65 && rule[1]==84 && rule[2]==84){
                    std::cout<<"asdasd: ";
                    for(size_t i=rule.size();i-->0;){
                        std::cout<<rule[i]<<" ";
                    }
                    std::cout<<" -> "<<id<<" "<<(val>>1UL)<<std::endl;
                }*/
                ht.insert_value_at(res.first, val);
            }
            /*for(size_t i=rule.size();i-->0;){
                std::cout<<rule[i]<<" ";
            }
            std::cout<<" "<<id<<std::endl;*/
            rule.pop_back();
            full_size=false;
        }

        next_rule:
        while(!r_lim[j]) j--;
        nt--;
    }
    compress_grammar(p_gram, rules, r_lim, rep_sym, ht, config);
}

