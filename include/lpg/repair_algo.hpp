//
// Created by diego on 22-02-20.
//

#ifndef LMS_COMPRESSOR_PAIRING_ALGORITHMS_H
#define LMS_COMPRESSOR_PAIRING_ALGORITHMS_H
#include <cassert>
#include <iostream>
#include <vector>

#include <sdsl/bit_vectors.hpp>
#include "cdt/hash_table.hpp"
#include "cdt/int_array.h"
#include "cdt/si_int_array.h"
#include "lpg_build.hpp"

typedef sdsl::bit_vector bv_t;
typedef sdsl::int_vector_buffer<> ivb;
typedef std::vector<size_t> pt_vector_t;
typedef lpg_build::plain_grammar_t plain_gram_t;
typedef lpg_build::key_wrapper key_wrapper;
typedef bit_hash_table<size_t, 44> ht_t;

//information of a pair
struct pair_data_t {
    size_t freq;      // frequency of the pair
    size_t q_locus;   // locus of the pair in the priority queue
    size_t id;        // ID assigned to the pair
    pt_vector_t *occ; // positions in list of rules where the pair occurs
    pair_data_t(size_t freq_, size_t q_locus_,
                size_t id_, pt_vector_t *occ_) : freq(freq_),
                                                 q_locus(q_locus_),
                                                 id(id_),
                                                 occ(occ_) {};
    pair_data_t() : freq(0),
                    q_locus(0),
                    id(0),
                    occ(nullptr) {};
};

typedef bit_hash_table<pair_data_t> repair_ht_t;

struct repair_data{

    plain_gram_t &       gram;         // plain grammar where the RePair rules are built upon
    sdsl::cache_config & config;
    size_t               s_width;
    int_array<size_t>    seq;          // list of rules
    bv_t                 seq_lim;
    bv_t::rank_1_type    seq_lim_rs;
    bv_t                 rep_syms;     //which symbols are repeated in seq?
    bv_t::rank_1_type    rep_syms_rs;  //rank support of the repeated symbols
    int_array<size_t>    sym_dists;
    repair_ht_t          ht;           //hash table with the pairs
    size_t               next_av_rule; //next available id for a rule
    ivb                  new_rules;    //buffer with the new RePair rules
    size_t               dummy_sym;    //dummy symbol

    explicit repair_data(plain_gram_t & gram_info_, sdsl::cache_config& config_): gram(gram_info_),
                                                                                  config(config_){

        ivb r(gram.rules_file, std::ios::in);
        s_width = sdsl::bits::hi(r.size()*2) + 1;
        seq = int_array<size_t>(r.size(), s_width);
        dummy_sym = int_array<size_t>::stream_t::masks[s_width];
        sdsl::load_from_file(seq_lim, gram.rules_lim_file);

        sdsl::int_vector<2> rep_record(r.size(), 0);
        size_t sym;
        for(auto && i : r){
            sym = i;
            if(rep_record[sym] < 2){
                rep_record[sym]++;
            }
            seq.push_back(sym);
        }
        r.close();

        rep_syms = bv_t(r.size(), false);
        size_t max_rep_len=0;
        size_t tmp_rep_len=0;
        for(size_t i=gram.max_tsym+1; i < seq.size(); i++){
            if(rep_record[seq[i]] > 1){
                rep_syms[i] = true;
                if(!seq_lim[i-1]){
                    tmp_rep_len++;
                }else{
                    if(tmp_rep_len>max_rep_len){
                        max_rep_len = tmp_rep_len;
                    }
                    tmp_rep_len = 1;
                }
            }else{
                if(tmp_rep_len>max_rep_len){
                    max_rep_len = tmp_rep_len;
                }
                tmp_rep_len = 0;
            }
        }
        if(tmp_rep_len>max_rep_len){
            max_rep_len = tmp_rep_len;
        }
        sdsl::util::init_support(rep_syms_rs, &rep_syms);

        size_t d_width = sdsl::bits::hi(max_rep_len)+1;
        sym_dists = int_array<size_t>(rep_syms_rs(rep_syms.size())*2, d_width);

        size_t mask = (1UL<<d_width)-1UL;
        size_t dist;
        bool prev_is_rep=false;
        for(size_t i=gram.max_tsym+1; i < seq.size(); i++){

            if(rep_syms[i]){
                dist=0;
                if(!seq_lim[i-1]){
                    dist = prev_is_rep? 1 : mask;
                }
                sym_dists.write(rep_syms_rs(i)*2, dist);

                dist=0;
                if(!seq_lim[i]){
                    dist = rep_syms[i+1] ? 1 : mask;
                }
                sym_dists.write(rep_syms_rs(i)*2+1, dist);
            }
            prev_is_rep = rep_syms[i];
        }
        sym_dists.write(sym_dists.size()-1, mask);

        next_av_rule= gram.r - 1;
        std::string nr_file=sdsl::cache_file_name("new_rp_rules", config);
        new_rules = ivb(nr_file, std::ios::out);
        sdsl::util::init_support(seq_lim_rs, &seq_lim);
    };

    ~repair_data(){
        if(new_rules.is_open()){
            new_rules.close(true);
        }
    }
};

struct pair_greater {
    const bitstream<ht_t::buff_t>& ht_data;
    size_t dist_bits;//distance from the pointer to the position where the pair_value starts
    size_t pt_bits=sizeof(size_t)*8;//number of bits used for the frequency

    pair_greater(const bitstream<ht_t::buff_t>& ht_str, size_t db): ht_data(ht_str), dist_bits(db) {};

    bool operator()(const size_t& pt1, const size_t& pt2) const {
        size_t freq1 = ht_data.read(pt1+dist_bits, pt1+dist_bits+pt_bits-1);
        size_t freq2 = ht_data.read(pt2+dist_bits, pt2+dist_bits+pt_bits-1);
        return freq1>freq2;
    }

    inline size_t get_value(size_t pt)const{
        return ht_data.read(pt+dist_bits, pt+dist_bits+pt_bits-1);
    }
};

struct priority_queue{

    repair_ht_t& ht;
    pair_greater comp;
    std::vector<size_t> arr;
    typedef long long int size_type;

    inline void swap(size_t l, size_t r){
        pair_data_t pd_l, pd_r;
        ht.get_value_from(arr[l], pd_l);
        ht.get_value_from(arr[r], pd_r);

        std::swap(arr[l], arr[r]);

        pd_l.q_locus = r;
        pd_r.q_locus = l;
        ht.insert_value_at(arr[r], pd_l);
        ht.insert_value_at(arr[l], pd_r);
    }

    inline static size_type parent(size_type idx) {
        return (idx-1)/2;
    }

    inline static size_type l_child(size_type idx) {
        return (2*idx + 1);
    }

    inline static size_type r_child(size_type idx) {
        return (2*idx + 2);
    }

    inline void trim_leaves() {
        auto pt = arr.back();
        while(!arr.empty() && comp.get_value(pt)<2){

            pair_data_t pd;
            ht.get_value_from(pt, pd);
            delete pd.occ;

            arr.pop_back();
            if(!arr.empty()) pt = arr.back();
        }
    }

    void sift_down(size_t idx){

        size_t largest, left, right;

        while(true){

            largest = idx;
            left = l_child((size_type)idx);
            right = r_child((size_type)idx);

            if(left<arr.size() && comp(arr[left], arr[idx])){
                largest = left;
            }

            if(right<arr.size() && comp(arr[right], arr[largest])){
                largest = right;
            }

            if(largest!=idx){
                swap(idx, largest);
                idx = largest;
            }else{
                break;
            }
        }
    }

    void sift_up(size_t idx){
        while (idx != 0 && comp(arr[idx], arr[parent((size_type)idx)]) ){
            swap(idx, parent((size_type)idx));
            idx = parent((size_type)idx);
        }
    }

    explicit priority_queue(pair_greater& comp_, repair_ht_t& ht_): ht(ht_), comp(comp_) {};

    inline bool empty() const{
        return arr.empty();
    }

    inline size_t size() const{
        return arr.size();
    }

    size_type extract_max(){
        assert(!empty());
        // get the maximum value, and remove it from heap
        size_t root = arr[0];
        swap(0, arr.size()-1);
        arr.pop_back();
        sift_down(0);
        return (size_type)root;
    }

    void insert(size_t elm) {
        arr.push_back(elm);

        //insert the queue locus in
        // the hash table
        pair_data_t pd;
        ht.get_value_from(elm, pd);
        pd.q_locus = arr.size()-1;
        ht.insert_value_at(elm, pd);
        //

        sift_up(arr.size()-1);
    }
};

void update_grammar(repair_data& rp_data){

    size_t new_g_size=0;
    for(size_t i=0;i<rp_data.seq.size();i++){
        if(rp_data.seq[i]!=rp_data.dummy_sym) new_g_size++;
    }
    size_t new_gsyms = rp_data.new_rules.size();
    new_g_size += new_gsyms;

    sdsl::int_vector_buffer<1> r_lim(rp_data.gram.rules_lim_file, std::ios::out, BUFFER_SIZE);
    ivb r(rp_data.gram.rules_file, std::ios::out);

    for(size_t i=0;i<=rp_data.gram.max_tsym; i++){
        r.push_back(rp_data.seq[i]);
        r_lim.push_back(true);
    }

    size_t new_c_size=0;
    size_t gsyms = rp_data.gram.g - rp_data.gram.c;
    for(size_t i=rp_data.gram.max_tsym+1; i < gsyms; i++){
        if(rp_data.seq[i]!=rp_data.dummy_sym){
            r.push_back(rp_data.seq[i]);
        }
        if(rp_data.seq_lim[i]){
            r_lim[r.size()-1] = true;
        }
    }

    for(size_t i=0;i<rp_data.new_rules.size();i++){
        r.push_back(rp_data.new_rules[i]);
        r_lim.push_back(i & 1UL);
    }

    for(size_t i=gsyms;i<rp_data.seq.size();i++){
        if(rp_data.seq[i]!=rp_data.dummy_sym){
            r.push_back(rp_data.seq[i]);
            new_c_size++;
        }
    }
    r_lim[r.size()-1] = true;

    r.close();
    r_lim.close();

    //mark the repair rules as non run-length compressed
    sdsl::int_vector_buffer<1> is_rl(rp_data.gram.is_rl_file, std::ios::in);
    for(size_t j = rp_data.gram.r; j<rp_data.gram.r+new_gsyms-1;j++){
        is_rl[j] = false;
    }

    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Grammar size before:        "<<rp_data.gram.g << std::endl;
    std::cout<<"      Grammar size after:         "<<new_g_size<<std::endl;
    std::cout<<"      Number of new nonterminals: "<<new_gsyms/2<<std::endl;
    std::cout<<"      Compression ratio:          "<<double(new_g_size)/double(rp_data.gram.g) << std::endl;

    rp_data.gram.g = new_g_size;
    rp_data.gram.r += new_gsyms / 2;
    rp_data.gram.c = new_c_size;
}

static void repair_prim_pairs(repair_data& rp_data, priority_queue& p_queue){

    int_array<size_t> pair(2, rp_data.s_width);
    size_t rep, lim, curr_rule;

    bv_t is_rl;
    sdsl::load_from_file(is_rl, rp_data.gram.is_rl_file);

    {//get the frequency of the first pairs

        rep = rp_data.rep_syms[rp_data.gram.max_tsym+1] << 1UL;
        lim = 4;// 100

        ht_t rep_map;

        size_t val=0;
        curr_rule=rp_data.gram.max_tsym+1;
        for (size_t i = rp_data.gram.max_tsym+1; i < rp_data.seq.size() - 1; i++) {

            rep |= rp_data.rep_syms[i+1];
            lim |= rp_data.seq_lim[i+1];

            //11: both symbols are repeated
            //!10: pair is not a junction of two rules
            if ((rep & 3UL) == 3UL && (lim & 3UL) != 2U && !is_rl[curr_rule]) {
                pair.write(0, rp_data.seq[i]);
                pair.write(1, rp_data.seq[i + 1]);
                auto res = rep_map.insert(pair.data(), pair.n_bits(), 1);
                if (!res.second) {
                    rep_map.get_value_from(*res.first, val);
                    rep_map.insert_value_at(*res.first, val + 1);
                }
            }
            rep <<= 1U;
            lim <<= 1U;
            if(rp_data.seq_lim[i]) curr_rule++;
        }

        size_t freq=0;
        auto key = reinterpret_cast<char*>(malloc(sizeof(size_t)*2));
        memset(key, 0, sizeof(size_t)*2);

        for(auto pt : rep_map){
            rep_map.get_value_from(pt, freq);
            if(freq>1){
                size_t k_bits = rep_map.get_key(pt, key);
                assert(k_bits==rp_data.s_width*2);
                auto occ = new pt_vector_t();
                pair_data_t new_pair{freq,0,0,occ};
                rp_data.ht.insert(key, k_bits, new_pair);
            }
        }
        free(key);
    }

    //store the positions of the repeated pairs
    rep = rp_data.rep_syms[rp_data.gram.max_tsym+1] << 1UL;
    lim = 4;// 100
    curr_rule=rp_data.gram.max_tsym+1;
    for (size_t i = rp_data.gram.max_tsym+1; i < rp_data.seq.size() - 1; i++) {
        rep |= rp_data.rep_syms[i+1];
        lim |= rp_data.seq_lim[i+1];

        if((rep & 3UL)==3UL && (lim & 3UL) != 2U && !is_rl[curr_rule]){ //==11 repeated symbol

            pair.write(0, rp_data.seq[i]);
            pair.write(1, rp_data.seq[i+1]);

            auto res = rp_data.ht.find(pair.data(), pair.n_bits());

            if(res.second){
                pair_data_t tmp;
                rp_data.ht.get_value_from(*res.first, tmp);
                if(rp_data.seq_lim[i-1] && !rp_data.seq_lim[i] && rp_data.seq_lim[i+1]){
                    tmp.id = rp_data.seq_lim_rs(i);
                }
                tmp.occ->push_back(i);
                rp_data.ht.insert_value_at(*res.first, tmp);
            }
        }
        rep <<= 1U;
        lim <<= 1U;
        if(rp_data.seq_lim[i]) curr_rule++;
    }

    for(auto pt : rp_data.ht){
        p_queue.insert(pt);
    }
}

static void re_pair_int(repair_data& rp_data){

    pair_greater pg{rp_data.ht.get_data(), rp_data.ht.description_bits()+(2*rp_data.s_width)};
    priority_queue p_queue(pg, rp_data.ht);

    repair_prim_pairs(rp_data, p_queue);

    size_t l_dist, r_dist, n_r_dist;
    size_t t_pos, t_r_pos;

    size_t curr_pt;
    pair_data_t p_data;

    int_array<size_t> tmp_pair(2, rp_data.s_width);
    int_array<size_t> curr_pair(2, rp_data.s_width);

    key_wrapper key_w{rp_data.s_width,
                      rp_data.ht.description_bits(),
                      rp_data.ht.get_data()};

    size_t d_mask = (1UL<<rp_data.sym_dists.width())-1UL;
    bool is_rule;

    while(!p_queue.empty()) {

        curr_pt = p_queue.extract_max();

        rp_data.ht.get_value_from(curr_pt, p_data);
        curr_pair.write(0, key_w.read(curr_pt, 0));
        curr_pair.write(1, key_w.read(curr_pt, 1));

        if(p_data.freq>1){

            if(p_data.id==0){
                p_data.id = rp_data.next_av_rule++;
                rp_data.new_rules.push_back(curr_pair[0]);
                rp_data.new_rules.push_back(curr_pair[1]);
            }

            for (auto const &pos : *(p_data.occ)) {

                //do the pair replacement
                t_pos = rp_data.rep_syms_rs(pos);

                l_dist = rp_data.sym_dists[t_pos*2];
                r_dist = rp_data.sym_dists[t_pos*2+1];

                //the pair was modified and the neighbor is unique
                if (r_dist == d_mask) continue;

                //right distance of the right symbol of the pair
                t_r_pos = rp_data.rep_syms_rs(pos + r_dist);
                n_r_dist = rp_data.sym_dists[(t_r_pos*2)+1];

                //check if the current occurrence of
                // the pair covers an entire rule
                is_rule = (l_dist == 0) && (n_r_dist == 0);

                tmp_pair.write(0, rp_data.seq[pos]);
                tmp_pair.write(1, rp_data.seq[pos+r_dist]);

                //We need to check if the pairs are equals in case the
                // occurrence was invalidated by other pair replacement
                if (tmp_pair == curr_pair && !is_rule) {

                    //replace the pair with its id
                    rp_data.seq.write(pos, p_data.id);
                    rp_data.seq.write(pos+r_dist, rp_data.dummy_sym);

                    //invalidate the distances of the right symbol
                    rp_data.sym_dists.write((t_r_pos*2), d_mask);
                    rp_data.sym_dists.write((t_r_pos*2)+1, d_mask);

                    if (n_r_dist != d_mask && n_r_dist!=0) {
                        r_dist += n_r_dist;
                        t_r_pos = rp_data.rep_syms_rs(pos + r_dist);
                        //left distance of the symbol next to the right pair symbol
                        rp_data.sym_dists.write(t_r_pos*2, r_dist);
                    } else {
                        r_dist = n_r_dist;
                    }

                    rp_data.sym_dists.write((t_pos*2)+1, r_dist);

                    //hash the new pairs
                    if (l_dist!=0 && l_dist != d_mask) {

                        //decrease the frequency of the pair to the left
                        pair_data_t pd1;
                        tmp_pair.write(0, rp_data.seq[pos-l_dist]);
                        tmp_pair.write(1, curr_pair[0]);
                        auto res1 = rp_data.ht.find(tmp_pair.data(), tmp_pair.n_bits());
                        if(res1.second){
                            rp_data.ht.get_value_from(*res1.first, pd1);
                            if(pd1.freq>0){
                                pd1.freq--;
                            }
                            rp_data.ht.insert_value_at(*res1.first, pd1);
                            p_queue.sift_down(pd1.q_locus);
                        }
                        //

                        //increase the frequency of the new pair to the right
                        pair_data_t pd2;
                        tmp_pair.write(0, rp_data.seq[pos-l_dist]);
                        tmp_pair.write(1, rp_data.seq[pos]);

                        //check if the new pair covers an entire rule
                        size_t t_l_pos = rp_data.rep_syms_rs(pos - l_dist);
                        size_t n_l_dist = rp_data.sym_dists[(t_l_pos*2)];
                        if((n_l_dist == 0) && (r_dist == 0)){
                            pd2.id = rp_data.seq_lim_rs(pos-l_dist);
                        }
                        //

                        auto res2 = rp_data.ht.find(tmp_pair.data(), tmp_pair.n_bits());
                        if(!res2.second) {
                            pd2.occ = new pt_vector_t();
                            pd2.freq=1;
                            auto ires = rp_data.ht.insert(tmp_pair.data(), tmp_pair.n_bits(), pd2);
                            p_queue.insert(*ires.first);
                        }else{
                            rp_data.ht.get_value_from(*res2.first, pd2);
                            pd2.freq++;
                            rp_data.ht.insert_value_at(*res2.first, pd2);
                            p_queue.sift_up(pd2.q_locus);
                        }
                        pd2.occ->push_back(pos - l_dist);
                    }

                    if (r_dist!=0 && r_dist != d_mask) {

                        //decrease the frequency of the pair to the right
                        pair_data_t pd1;
                        tmp_pair.write(0, curr_pair[1]);
                        tmp_pair.write(1, rp_data.seq[pos+r_dist]);
                        auto res1 = rp_data.ht.find(tmp_pair.data(), tmp_pair.n_bits());
                        if(res1.second){
                            rp_data.ht.get_value_from(*res1.first, pd1);
                            if(pd1.freq>0){
                                pd1.freq--;
                            }
                            rp_data.ht.insert_value_at(*res1.first, pd1);
                            p_queue.sift_down(pd1.q_locus);
                        }
                        //

                        //increase the frequency of the new right pair
                        pair_data_t pd2;
                        tmp_pair.write(0, rp_data.seq[pos]);
                        tmp_pair.write(1, rp_data.seq[pos+r_dist]);

                        //check if the current occurrence covers an entire rule
                        t_r_pos = rp_data.rep_syms_rs(pos + r_dist);
                        n_r_dist = rp_data.sym_dists[(t_r_pos*2)+1];

                        if((l_dist==0) && (n_r_dist==0)){
                            pd2.id = rp_data.seq_lim_rs(pos);
                        }

                        auto res = rp_data.ht.find(tmp_pair.data(), tmp_pair.n_bits());
                        if(!res.second) {
                            pd2.occ = new pt_vector_t();
                            pd2.freq = 1;
                            auto ires = rp_data.ht.insert(tmp_pair.data(), tmp_pair.n_bits(), pd2);
                            p_queue.insert(*ires.first);
                        }else{
                            rp_data.ht.get_value_from(*res.first, pd2);
                            pd2.freq++;
                            rp_data.ht.insert_value_at(*res.first, pd2);
                            p_queue.sift_up(pd2.q_locus);
                        }
                        pd2.occ->push_back(pos);
                    }
                }
            }
            p_queue.trim_leaves();
        }
        //release memory
        delete p_data.occ;
    }

    rp_data.ht.reset();
    sdsl::util::clear(rp_data.rep_syms);
    sdsl::util::clear(rp_data.rep_syms_rs);
    update_grammar(rp_data);
}

void repair(plain_gram_t& p_gram, sdsl::cache_config &config) {
    std::cout<<"  Running RePair over the grammar's rules"<<std::endl;
    //Load the grammar from file
    repair_data rp_data{p_gram, config};
    re_pair_int(rp_data);
}
#endif //LMS_COMPRESSOR_PAIRING_ALGORITHMS_H
