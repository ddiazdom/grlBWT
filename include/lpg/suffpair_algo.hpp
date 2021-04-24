//
// Created by diego on 17-09-20.
//

#ifndef LPG_COMPRESSOR_SUFFPAIR_ALGO_H
#define LPG_COMPRESSOR_SUFFPAIR_ALGO_H
#include <iostream>
#include "lpg_build.hpp"

#include <mutex>
#include <condition_variable>

//global variables for multi threading
std::mutex m;
std::condition_variable cv;
bool proceed;
bool merged;
size_t processed;
size_t rem_threads;
//

using bv_t = sdsl::bit_vector;
using pair_t = int_array<size_t>;
using ivb = sdsl::int_vector_buffer<>;
using ht_t =  bit_hash_table<size_t, 44>;
using key_wrapper = lpg_build::key_wrapper;
using plain_gram_t = lpg_build::plain_grammar_t;

struct pairing_data{
    size_t                    next_av_rule;
    size_t                    lim_id;
    size_t                    n_threads;
    size_t                    hbuff_size;
    size_t                    tot_lms_rules;
    size_t                    gsyms;
    uint8_t                   s_width;
    uint8_t                   sigma;

    std::string               r_file;
    std::string               nr_file;
    std::string               lmsg_as_sp_file;
    ivb                       new_rules;

    bv_t                      r_lim;
    bv_t::rank_1_type         r_lim_rs;
    const bv_t&               rep_sym;
    const sdsl::cache_config& config;

    bool                      first_run=true;
    size_t elms_frun          =0;

    pairing_data(plain_gram_t& gram_info, bv_t& rep_syms_, size_t n_treads_, size_t hbuff_size_, sdsl::cache_config& config_):
            rep_sym(rep_syms_),
            config(config_){

        sigma = gram_info.sigma;
        hbuff_size = hbuff_size_;
        n_threads = n_treads_;
        next_av_rule = gram_info.r - 1; //I subtract one as the last rule is the compressed string
        tot_lms_rules = gram_info.r - 1; // the same
        lim_id = 2*gram_info.g - gram_info.r;
        s_width = sdsl::bits::hi(lim_id)+1;
        r_file = gram_info.rules_file;
        lmsg_as_sp_file = gram_info.lms_as_sp_file;
        nr_file = sdsl::cache_file_name("nr_file", config);
        new_rules = ivb(nr_file, std::ios::out, BUFFER_SIZE);
        gsyms = gram_info.g - gram_info.c;

        sdsl::load_from_file(r_lim, gram_info.rules_lim_file);
        sdsl::util::init_support(r_lim_rs, &r_lim);
    }

    ~pairing_data(){
        if(new_rules.is_open()) new_rules.close(true);
    }
};

struct thread_data{
    const pairing_data& p_data;
    const size_t        start;
    const size_t        end;

    bool                finished;
    std::string         pl_chunk_file;
    std::string         ht_file;
    ht_t&               ht;
    void *              hb_addr;
    size_t              hb_size;

    thread_data(const pairing_data& p_data_, const size_t start_, const size_t end_, ht_t& ht_, size_t hb_size_, void * hb_addr_):
                p_data(p_data_),
                start(start_),
                end(end_),
                ht(ht_),
                hb_addr(hb_addr_),
                hb_size(hb_size_){
        std::stringstream ss;
        ss << p_data.config.dir<< "/range_" << start << "_" << end;
        std::string prefix = ss.str();

        ht_file = prefix + "_pairs";
        pl_chunk_file = prefix+"_pt_chunk_file";
        finished = false;
    }
};

//create a set of new rules from the hash table
void create_new_rules(ht_t& ht, pairing_data& p_data){
    key_wrapper key_w{p_data.s_width, ht.description_bits(), ht.get_data()};
    ht_t::val_type val=0;
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

void compute_lms_as_sp(i_file_stream<size_t>& r, pairing_data& p_data){
    //mark which nonterminals have LMSg symbols as SP
    sdsl::int_vector_buffer<1> lmsg_as_sp(p_data.lmsg_as_sp_file, std::ios::out, BUFFER_SIZE);
    size_t j=0,k;
    //skip terminal symbols
    while(p_data.r_lim[j]) j++;
    size_t curr_rule=j;
    while(j<p_data.gsyms){
        while(!p_data.r_lim[j])j++;
        k=j;
        while(r.read(k)==p_data.lim_id) k--;
        if(k<j && r.read(k)<p_data.tot_lms_rules){
            lmsg_as_sp[curr_rule] = true;
        }
        curr_rule++;
        j++;
    }
    curr_rule+=(p_data.elms_frun/2);
    for(size_t i=p_data.elms_frun+1;i<p_data.new_rules.size();i+=2){
        if(p_data.new_rules[i]<p_data.tot_lms_rules){
            lmsg_as_sp[curr_rule] = true;
        }else{
            lmsg_as_sp[curr_rule] = false;
        }
        curr_rule++;
    }
    lmsg_as_sp[curr_rule] = false; //we are not considering the compressed string!!
    //
}

void get_rep_suff_rules(i_file_stream<size_t>& rules, pairing_data& p_data, bv_t& uniq_nr){
    //mark the new suff. phrases with one occurrence
    sdsl::int_vector<2> uniq_nr_tmp(p_data.new_rules.size() / 2, 0);
    uniq_nr.resize(p_data.new_rules.size()/2);

    for (size_t i = 0; i < p_data.new_rules.size(); i++) {
        if (i & 1UL &&
            p_data.new_rules[i] >= p_data.tot_lms_rules &&
            uniq_nr_tmp[p_data.new_rules[i] - p_data.tot_lms_rules] < 2) {
            uniq_nr_tmp[p_data.new_rules[i] - p_data.tot_lms_rules]++;
        }
    }
    for (size_t i = 0; i < rules.size(); i++) {
        if (rules.read(i) == p_data.lim_id) continue;
        if (rules.read(i) >= p_data.tot_lms_rules &&
            uniq_nr_tmp[rules.read(i) - p_data.tot_lms_rules] < 2) {
            uniq_nr_tmp[rules.read(i) - p_data.tot_lms_rules]++;
        }
    }
    for(size_t i=0;i<uniq_nr_tmp.size();i++) uniq_nr[i] = uniq_nr_tmp[i]==1UL;
}

//insert the suff. pair rules into the grammar
void update_grammar(pairing_data& p_data, plain_gram_t& gram){

    std::cout<<"  Updating the grammar"<<std::endl;
    i_file_stream<size_t> rules(p_data.r_file, BUFFER_SIZE);
    compute_lms_as_sp(rules, p_data);

    std::cout<<"    Computing unique SuffPair"<<std::endl;
    //mark which suff. pairs are unique to delete them
    bv_t uniq_sr;
    bv_t::rank_1_type uniq_sr_rs;
    get_rep_suff_rules(rules, p_data, uniq_sr);
    sdsl::util::init_support(uniq_sr_rs, &uniq_sr);
    p_data.new_rules.close();
    sdsl::int_vector<> new_rules;
    sdsl::load_from_file(new_rules, p_data.nr_file);
    //

    std::cout<<"    Collapsing LMS rules"<<std::endl;
    //collapsed rules
    std::string tr_file = sdsl::cache_file_name("col_rules", p_data.config);
    ivb col_rules(tr_file, std::ios::out, BUFFER_SIZE, p_data.s_width);

    //collapse compressed symbols
    size_t n_av=0, tmp_sym, sr_pos;
    for(size_t i=0; i < p_data.gsyms; i++){

        if(i>p_data.sigma && p_data.r_lim[i-1]){
            p_data.r_lim[n_av-1] = true;
        }

        tmp_sym = rules.read(i);
        if(tmp_sym != p_data.lim_id) {

            //I have to delete unique suff. nonterminals in the LMS phrases
            //Example:
            // suff phrase 4080
            // LMS phrases: A -> 64 4080, B-> 65 65 4080
            // Then we have:
            //  A -> 64 4080, B-> 65 A
            // Now the SP phrase 4080 has only one occurrence
            if(tmp_sym>=p_data.tot_lms_rules &&
               uniq_sr[tmp_sym-p_data.tot_lms_rules]){ //decompress unique suffix phrases

                sr_pos = 2*(tmp_sym-p_data.tot_lms_rules);

                while(true){
                    //decompress as long as the last symbol is also unique
                    col_rules.push_back(new_rules[sr_pos]);
                    p_data.r_lim[n_av++] = false;

                    tmp_sym = new_rules[sr_pos+1];
                    if(tmp_sym>=p_data.tot_lms_rules &&
                       uniq_sr[tmp_sym-p_data.tot_lms_rules]){
                        sr_pos = 2*(tmp_sym-p_data.tot_lms_rules);
                    }else{
                        if(tmp_sym>=p_data.tot_lms_rules) tmp_sym -= uniq_sr_rs(tmp_sym-p_data.tot_lms_rules);
                        col_rules.push_back(tmp_sym);
                        p_data.r_lim[n_av++] = true;
                        break;
                    }
                }
            }else{
                if(tmp_sym>=tmp_sym-p_data.tot_lms_rules) tmp_sym -= uniq_sr_rs(tmp_sym-p_data.tot_lms_rules);
                col_rules.push_back(tmp_sym);
                p_data.r_lim[n_av++] = i<p_data.sigma;
            }
        }
    }
    p_data.r_lim[n_av-1]=true;

    std::cout<<"    Collapsing new SuffPair rules"<<std::endl;
    //insert the symbols of the new_rules
    for(size_t i=0;i<new_rules.size();i+=2){

        if(!uniq_sr[i>>1UL]){
            col_rules.push_back(new_rules[i]);
            p_data.r_lim[n_av++] = false;
            tmp_sym = new_rules[i+1];

            if(tmp_sym>=p_data.tot_lms_rules &&
               uniq_sr[tmp_sym-p_data.tot_lms_rules]){//decompress unique suffix phrases

                sr_pos = 2*(tmp_sym-p_data.tot_lms_rules);

                while(true){
                    //decompress as long as the
                    // last symbol is also unique
                    assert(new_rules[sr_pos]<p_data.tot_lms_rules);
                    col_rules.push_back(new_rules[sr_pos]);
                    p_data.r_lim[n_av++] = false;

                    tmp_sym = new_rules[sr_pos+1];
                    if(tmp_sym>=p_data.tot_lms_rules &&
                       uniq_sr[tmp_sym-p_data.tot_lms_rules]){
                        sr_pos = 2*(tmp_sym-p_data.tot_lms_rules);
                    }else{
                        if(tmp_sym>=p_data.tot_lms_rules) tmp_sym -= uniq_sr_rs(tmp_sym-p_data.tot_lms_rules);
                        col_rules.push_back(tmp_sym);
                        p_data.r_lim[n_av++] = true;
                        break;
                    }
                }
            }else{
                if(tmp_sym>=p_data.tot_lms_rules) tmp_sym -= uniq_sr_rs(tmp_sym-p_data.tot_lms_rules);
                col_rules.push_back(tmp_sym);
                p_data.r_lim[n_av++] = true;
            }
        }
    }

    std::cout<<"    Inserting the compressed string"<<std::endl;
    //put array C at the end of the new rules
    for(size_t i=p_data.gsyms; i < rules.size(); i++){
        col_rules.push_back(rules.read(i));
        p_data.r_lim[n_av++] = false;
    }
    p_data.r_lim[n_av-1] = true;
    p_data.r_lim.resize(n_av);

    col_rules.close();
    rules.close();

    size_t eff_new_rules = uniq_sr.size()-uniq_sr_rs(uniq_sr.size());

    rename(col_rules.filename().c_str(), gram.rules_file.c_str());
    sdsl::store_to_file(p_data.r_lim, gram.rules_lim_file);

    std::cout<<"  SuffPair stats:"<<std::endl;
    std::cout<<"    Grammar size before:         "<<gram.g - gram.sigma << std::endl;
    std::cout<<"    Grammar size after:          "<<n_av-gram.sigma<<std::endl;
    std::cout<<"    Number of new nonterminals:  "<<eff_new_rules<<std::endl;
    std::cout<<"    Compression ratio:           "<<double(n_av)/double(gram.g) << std::endl;

    gram.g = n_av;
    gram.r += eff_new_rules;
}

//change the width of R and compute the repeated symbols
void prepare_input(plain_gram_t& gram_info, bv_t& rep_syms, sdsl::cache_config& config){

    ivb r(gram_info.rules_file, std::ios::in);

    std::string tmp_string = sdsl::cache_file_name("tmp_rules", config);
    o_file_stream<size_t> tmp_r(tmp_string,  BUFFER_SIZE, std::ios::out);
    sdsl::int_vector<2> tmp_rep(gram_info.r - 1, false);
    for (auto const &sym : r){
        tmp_r.push_back(sym);
        if(tmp_rep[sym]<2) tmp_rep[sym]++;
    }

    r.close(true);
    tmp_r.close();

    size_t j=0;
    for(auto && i : tmp_rep){
        assert(i==1 || i==2);
        rep_syms[j++] = i-1;
    }
    rename(tmp_string.c_str(), gram_info.rules_file.c_str());
}

void collect_pairs(thread_data* d, i_file_stream<size_t>& p_list, o_file_stream<size_t>& r){

    size_t id, val;
    pair_t tmp_pair(2, d->p_data.s_width);
    ht_t ht(d->hb_size, d->ht_file, 0.8, d->hb_addr);

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
void replace_pairs(const ht_t& ht, const pairing_data& p_data, std::string& pl_file, o_file_stream<size_t>& r){

    i_file_stream<size_t> p_list(pl_file, BUFFER_SIZE);

    //temporal file to store the new suffix positions
    std::string npl_file = pl_file+"_tmp";
    o_file_stream<size_t> np_list(npl_file,  BUFFER_SIZE, std::ios::out);

    //temporal variables
    pair_t tmp_pair(2, p_data.s_width);
    ht_t::val_type  val=0;

    //update the suffixes with the new symbols
    size_t pos;
    for(size_t i=0; i<p_list.size();i++){
        pos = p_list.read(i);
        if(!p_data.r_lim[pos-1]){ //the suffix cannot span a complete rule
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
        ht_t ht(d->hb_size, d->ht_file, 0.8, d->hb_addr);
        pair_t tmp_pair(2, d->p_data.s_width);
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
            if ((desc & 2UL) == 2UL &&
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
     // in its range and it is waiting for the global merge
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

        {//notify the main thread we finish the loop and we are waiting
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

void merge_pointer_data(std::vector<thread_data>& t_data, std::string& pl_file){

    //concatenate the files
    std::ofstream of(pl_file, std::ofstream::binary);
    size_t buff_size = BUFFER_SIZE/sizeof(size_t);
    size_t len, rem, to_read;
    auto *buffer = new size_t[buff_size];

    //collect the pointers to the pairs
    for(auto & data : t_data){

        std::ifstream i_file(data.pl_chunk_file, std::ifstream::binary);

        i_file.seekg (0, std::ifstream::end);
        len = i_file.tellg()/sizeof(size_t);
        i_file.seekg (0, std::ifstream::beg);

        rem=len;
        to_read = std::min<size_t>(buff_size, len);

        while(true){

            i_file.seekg( (rem - to_read) * sizeof(size_t));
            i_file.read((char *)buffer, sizeof(size_t)*to_read);
            assert(i_file.good());

            of.write((char *)buffer, sizeof(size_t)*to_read);
            assert(of.good());

            rem -= i_file.gcount()/sizeof(size_t);
            to_read = std::min<size_t>(buff_size, rem);
            if(to_read == 0) break;
        }
        i_file.close();

        if(remove(data.pl_chunk_file.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }
    }
    delete[] buffer;
    of.close();
}

void merge_ht_data(std::vector<thread_data>& t_data){

    //collect the pairs
    using buff_t = ht_t::buff_t;
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

        bitstream<buff_t> bits;
        bits.stream = reinterpret_cast<buff_t*>(k_buffer);

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
    ht_t ht;
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

        //store how many elements were created in the first run
        if(p_data.first_run){
            p_data.elms_frun = p_data.new_rules.size();
            p_data.first_run = false;
        }

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

void suffpair(std::string &g_file, sdsl::cache_config &config, size_t n_threads, size_t hbuff_size) {

    //Load the grammar information from file
    plain_gram_t gram_info;
    gram_info.load_from_file(g_file);

    //prepare the grammar symbols for suffix pairing
    bv_t rep_syms(gram_info.r - 1, false);
    prepare_input(gram_info, rep_syms, config);

    //create an object with all the necessary
    // information for the suffixpair algorithm
    pairing_data p_data(gram_info, rep_syms, n_threads, hbuff_size, config);

    suffixpair_int(p_data);
    update_grammar(p_data, gram_info);
    gram_info.save_to_file(g_file);
}
#endif //LPG_COMPRESSOR_SUFFPAIR_ALGO_H
