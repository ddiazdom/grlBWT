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

//create a set of new rules from the hash table
void create_new_rules(phrase_map_t& ht, pairing_data& p_data){
    key_wrapper key_w{p_data.s_width, ht.description_bits(), ht.get_data()};
    phrase_map_t::val_type val=0;
    for(auto elm : ht){
        ht.get_value_from(elm, val);

        if(val==1){//repeated pair but new
            val |= p_data.next_av_rule<<1UL;
            ht.insert_value_at(elm, val);
            p_data.new_rules.push_back(key_w.read(elm, 0));
            p_data.new_rules.push_back(key_w.read(elm, 1));
            p_data.next_av_rule++;
        }
    }
}

//insert the suffix pair rules into the grammar
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
    assert((n_av+p_data.new_rules.size()/2)<p_data.r_lim.size());

    for(size_t i=0;i<p_data.new_rules.size();i++){
        col_rules.push_back(p_data.new_rules[i]);
        p_data.r_lim[n_av++] = (i & 1UL);
    }

    //put array C at the end of the new rules
    for(size_t i=p_data.gsyms; i < rules.size(); i++){
        col_rules.push_back(rules.read(i));
        p_data.r_lim[n_av++] = false;
    }
    p_data.r_lim[n_av-1] = true;
    p_data.r_lim.resize(n_av);
    sdsl::store_to_file(p_data.r_lim, gram.rules_lim_file);

    gram.r += p_data.new_rules.size()/2;
    gram.rules_breaks.push_back(gram.rules_breaks.back() + p_data.new_rules.size() / 2);

    col_rules.close();
    rules.close();
    rename(col_rules.filename().c_str(), gram.rules_file.c_str());

    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Grammar size before:        "<<gram.g<< std::endl;
    std::cout<<"      Grammar size after:         "<<n_av<<std::endl;
    std::cout<<"      Number of new nonterminals: "<<p_data.new_rules.size()/2<<std::endl;
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

    size_t id, val;
    string_t tmp_pair(2, d->p_data.s_width);
    phrase_map_t ht(d->hb_size, d->ht_file, 0.8, d->hb_addr);

    for(size_t i=0;i<p_list.size();i++){
        size_t pos = p_list.read(i);

        assert(pos>=d->start && pos+1<=d->end);

        tmp_pair.write(0, r.read(pos));
        tmp_pair.write(1, r.read(pos+1));

        //check if the pair covers an entire rule
        id = d->p_data.r_lim[pos-1]? (d->p_data.r_lim_rs(pos)<<1UL) : 0;

        auto res = ht.insert(tmp_pair.data(), tmp_pair.n_bits(), id);
        if(!res.second){//the pair already exists
            val = res.first.value();
            val |= (id | 1UL);
            ht.insert_value_at(*res.first, val);
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
    string_t tmp_pair(2, p_data.s_width);
    phrase_map_t::val_type  val=0;

    //update the suffixes with the new symbols
    size_t pos;
    for(size_t i=0; i<p_list.size();i++){
        pos = p_list.read(i);
        if(!p_data.r_lim[pos-1]){ //the suffix can't span a complete rule
            tmp_pair.write(0, r.read(pos));
            tmp_pair.write(1, r.read(pos+1));

            auto res = ht.find(tmp_pair.data(), tmp_pair.n_bits());
            val = res.first.value();

            if(val & 1UL){
                r.write(pos, val>>1UL);
                r.write(pos+1, p_data.lim_id); //<-mask
                if(p_data.rep_sym[r.read(pos-1)]){
                    np_list.push_back(pos-1);
                }
            }
        }
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
        string_t tmp_pair(2, d->p_data.s_width);
        size_t val;
        //position in R of the suffixes of length two
        o_file_stream<size_t> p_list(d->pl_chunk_file, BUFFER_SIZE, std::ios::out);

        size_t start = d->start;
        size_t i = d->end;
        size_t sym, psym, desc, id;

        //skip the terminal symbols
        if (start == 0) while (d->p_data.r_lim[start]) start++;

        psym = 0;
        desc = 0;

        while (i >= start) {

            sym = r.read(i);
            desc = (desc << 1UL) | d->p_data.r_lim[i];

            //new repeated suffix -> pattern 10 with both symbols appearing more than once in R
            if ((desc & 3UL) == 2UL &&
                d->p_data.rep_sym[sym] &&
                d->p_data.rep_sym[psym]) {

                tmp_pair.write(0, sym);
                tmp_pair.write(1, psym);


                //pattern x101: pair covers an entire preexisting rule
                id = (d->p_data.r_lim[i-1]) ? ( d->p_data.r_lim_rs(i)<< 1UL) : 0;

                auto res = ht.insert(tmp_pair.data(), tmp_pair.n_bits(), id);

                if (!res.second) {
                    val = res.first.value();
                    val |= (id | 1UL);
                    ht.insert_value_at(*res.first, val);
                }
                p_list.push_back(i);
            }
            psym = sym;
            i--;
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
    //std::cout<<r.modified<<" "<<r.size()<<std::endl;
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

void merge_ht_data(std::vector<thread_data>& t_data){

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
                size_t val = res.first.value();
                val |= (id | 1UL);
                data.ht.insert_value_at(*res.first, val);
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

void suffpair(gram_info_t& p_gram, sdsl::cache_config &config, size_t n_threads, size_t hbuff_size) {

    std::cout<<"  Suffpair algorithm"<<std::endl;
    //prepare the grammar symbols for suffix pairing
    bv_t rep_syms(p_gram.r - 1, false);
    prepare_input(p_gram, rep_syms, config);

    //create an object with all the necessary
    // information for the suffix pair algorithm
    pairing_data p_data(p_gram, rep_syms, n_threads, hbuff_size, config);

    suffixpair_int(p_data);
    update_grammar(p_data, p_gram);
}

