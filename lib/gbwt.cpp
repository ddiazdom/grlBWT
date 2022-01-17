//
// Created by Diaz, Diego on 23.11.2021.
//
#include "gbwt.hpp"
#include "utils.hpp"
#ifdef __linux__
#include <malloc.h>
#endif

void check_plain_grammar(gram_info_t& p_gram, std::string& uncomp_file) {

    sdsl::int_vector<> r;
    bv_t r_lim;
    sdsl::load_from_file(r, p_gram.rules_file);
    sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
    bv_t::select_1_type r_lim_ss;
    sdsl::util::init_support(r_lim_ss, &r_lim);

    std::cout<<"  Checking the grammar produces the exact input string"<<std::endl;
    std::cout<<"    This step is optional and for debugging purposes"<<std::endl;
    std::cout<<"    Terminals:              "<<(size_t)p_gram.sigma<<std::endl;
    std::cout<<"    Number of nonterminals: "<<p_gram.r-p_gram.sigma<<std::endl;
    std::cout<<"    Compressed string:      "<<p_gram.c<<std::endl;

    std::vector<size_t> tmp_decomp;

    i_file_stream<uint8_t> if_stream(uncomp_file, BUFFER_SIZE);

    size_t start, end;
    std::stack<std::pair<size_t, uint8_t>> stack;

    size_t f = r_lim_ss(p_gram.r - 1) + 1;
    size_t l = r_lim_ss(p_gram.r);
    size_t idx = 0;


    for(size_t i=f; i <= l; i++){
        //std::cout<<i-f<<" "<<l-f<<std::endl;
        tmp_decomp.clear();
        stack.emplace(r[i], 1);
        assert(stack.size()<=if_stream.size());

        while(!stack.empty()){

            auto curr_sym = stack.top() ;
            stack.pop();

            if(curr_sym.first==0){
                start = 0;
            }else{
                start = r_lim_ss(curr_sym.first)+1;
            }

            end = r_lim_ss(curr_sym.first+1);

            if(r[start] == curr_sym.first){
                //std::cout<<idx<<" "<<p_gram.sym_map[curr_sym.first]<<" "<<(char)if_stream.read(idx)<<std::endl;
                //if(p_gram.sym_map[curr_sym.first]!=if_stream.read(idx)){
                //    std::cout<<f<<" "<<i<<" "<<i-f<<std::endl;
                //}
                assert(p_gram.sym_map[curr_sym.first]==if_stream.read(idx));
                idx++;
            }else{
                //std::cout<<curr_sym.first<<" -> ";
                 uint8_t is_sp;
                 uint8_t is_first;
                 if(curr_sym.second!=0){
                     for(size_t j=end+1; j-->start;){
                         is_sp =  j==end && p_gram.parsing_level(r[j])>p_gram.parsing_level(r[j-1]);
                         is_first = j==start && curr_sym.second & 1UL;
                         stack.emplace(r[j], is_first | (is_sp<<1UL));
                     }
                 }else{
                     for(size_t j=end; j>start;j--){
                         is_sp =  j==end && p_gram.parsing_level(r[j])>p_gram.parsing_level(r[j-1]);
                         stack.emplace(r[j], is_sp<<1UL);
                     }
                 }

                 //for(size_t j=start; j<=end;j++){
                 //    std::cout<<r[j]<<", "<<p_gram.parsing_level(r[j])<<" ";
                //}
                //std::cout<<""<<std::endl;
            }
        }
    }
    std::cout<<"\tGrammar is correct!!"<<std::endl;
}

alpha_t get_alphabet(std::string &i_file) {

    std::cout <<"Reading input file:" << std::endl;

    //TODO this can be done in parallel if the input is too big
    size_t alph_frq[256] = {0};
    alpha_t alphabet;

    i_file_stream<uint8_t> if_stream(i_file, BUFFER_SIZE);
    for (size_t i = 0; i < if_stream.tot_cells; i++) {
        alph_frq[if_stream.read(i)]++;
    }

    for (size_t i = 0; i < 256; i++) {
        if (alph_frq[i] > 0) alphabet.emplace_back(i, alph_frq[i]);
    }

    std::cout<<"  Number of characters: "<< if_stream.size() << std::endl;
    std::cout<<"  Number of strings:    "<< alphabet[0].second << std::endl;
    std::cout<<"  Alphabet:             "<< alphabet.size() << std::endl;
    std::cout<<"  Smallest symbol:      "<< (int) alphabet[0].first << std::endl;
    std::cout<<"  Greatest symbol:      "<< (int) alphabet.back().first << std::endl;

    if (if_stream.read(if_stream.size() - 1) != alphabet[0].first) {
        std::cout << "Error: sep. symbol " << alphabet[0].first << " differs from last symbol in file "
                  << if_stream.read(if_stream.size() - 1) << std::endl;
        exit(1);
    }
    return alphabet;
}

//extract freq symbols from bwt[j] onwards and put them in new_bwt
void extract_rl_syms(ivb_t& bwt, ivb_t& new_bwt, size_t& j, size_t freq, size_t& pos){
    size_t tmp_freq, sym;
    while(freq>0) {
        sym = bwt[j];
        tmp_freq = bwt[j+1];
        if(tmp_freq<=freq){
            freq-=tmp_freq;
            j+=2;
        }else{
            bwt[j+1]-=freq;
            tmp_freq = std::exchange(freq, 0);
        }
        if(pos>0 && new_bwt[pos-2]==sym){
            new_bwt[pos-1] += tmp_freq;
        } else {
            new_bwt.push_back(sym);
            new_bwt.push_back(tmp_freq);
            pos+=2;
        }
    }
}

void infer_lvl_bwt(sdsl::cache_config& config, size_t p_round) {

    dictionary dict;
    std::string dict_file = sdsl::cache_file_name("dict_lvl_"+std::to_string(p_round), config);
    sdsl::load_from_file(dict, dict_file);

    std::string prev_bwt_file = sdsl::cache_file_name("bwt_lvl_"+std::to_string(p_round+1), config);
    ivb_t bwt(prev_bwt_file, std::ios::in);

    bv_rs_t hocc_rs(&dict.phrases_has_hocc);
    vector_t hocc(dict.hocc_buckets[dict.hocc_buckets.size()-1]*2, 0, sdsl::bits::hi(std::max(dict.t_size, dict.alphabet+1))+1);
    size_t sym, left_sym, pos, freq, ptr, rank, dummy_sym = dict.alphabet+1;

    //TODO
    //vector_t bwt_buckets(dict.n_phrases+1, 0);
    //

    //std::cout<<dict.t_size+1<<" "<<dummy_sym+1<<std::endl;
    //assert(dict.t_size+1>=dummy_sym+1);
    //std::cout<<dict.phrases_has_hocc.size()<<" "<<dict.n_phrases<<std::endl;
    assert(dict.phrases_has_hocc.size() == dict.n_phrases);
    std::cout<<"    Inferring the order of symbols' hidden occurrences"<<std::endl;
    for(size_t i=0;i<bwt.size();i+=2) {

        sym = bwt[i];
        freq = bwt[i+1];
        //bwt_buckets[sym]+=freq;
        if(dict.phrases_has_hocc[sym]){
            rank = hocc_rs(sym);
            ptr = dict.hocc_buckets[rank]*2;

            if(ptr>2 && hocc[ptr-2]==(dummy_sym+1)){
                hocc[ptr-1] += freq;
            }else{
                //std::cout<<"sym: "<<sym<<" with freq: "<<freq<<" is transformed to dummy in position ("<<i<<", "<<i+1<<")"<<std::endl;
                hocc[ptr] = dummy_sym+1;
                hocc[ptr+1] = freq;
                dict.hocc_buckets[rank]++;
            }
        }

        pos = 2*sym+1;
        sym = dict.dict[pos];
        while(sym>=dict.alphabet){
            sym -= dict.alphabet;
            left_sym = dict.dict[pos-1];
            assert(left_sym<dict.alphabet && dict.phrases_has_hocc[sym]);

            rank = hocc_rs(sym);
            ptr = dict.hocc_buckets[rank]*2;

            if(ptr>0 && hocc[ptr-2]==(left_sym+1)){
                hocc[ptr-1] += freq;
            }else{
                hocc[ptr] = left_sym+1;
                hocc[ptr+1] = freq;
                dict.hocc_buckets[rank]++;
            }
            pos = 2*sym+1;
            sym = dict.dict[pos];
        }
        bwt[i] =  sym;
        assert(dict.dict[pos-1]==dict.dict_dummy);
    }

    sdsl::util::clear(dict.dict);
    sdsl::util::clear(dict.hocc_buckets);
    sdsl::util::clear(hocc_rs);
    pos=0;
    for(size_t i=0;i<hocc.size();i+=2){
        if(hocc[i]!=0){
            hocc[pos++] = hocc[i]-1;
            hocc[pos++] = hocc[i+1];
        }
    }
    hocc.resize(pos);

#ifdef __linux__
    malloc_trim(0);
#endif

    std::cout<<"    Merging the sorted lists into the new BWT"<<std::endl;
    std::string pre_bwt_file = sdsl::cache_file_name("pbwt_lvl_"+std::to_string(p_round), config);
    ivb_t pre_bwt(pre_bwt_file, std::ios::in);

    std::string new_bwt_file = sdsl::cache_file_name("bwt_lvl_"+std::to_string(p_round), config);
    ivb_t new_bwt(new_bwt_file, std::ios::out);

    //TODO testing
    //size_t acc=0, tmp;
    //for(size_t i=0;i<bwt_buckets.size()-1;i++){
    //    tmp = bwt_buckets[i];
    //    bwt_buckets[i] = acc;
    //    acc +=tmp;
    //}
    //bwt_buckets[bwt_buckets.size()-1] = acc;
    //


    size_t i=0, j=0, k=0, pbwt_freq, new_bwt_size=0;
    pos=0;
    rank = 0;
    //acc=0;
    while(i<pre_bwt.size()) {

        sym = pre_bwt[i];
        pbwt_freq = pre_bwt[i+1];
        //std::cout<<sym<<" "<<pbwt_freq<<std::endl;

        if(sym==dummy_sym) {

            if(dict.phrases_has_hocc[rank]) {
                //std::cout<<i<<" "<<pre_bwt[i]<<" "<<pre_bwt[i+1]<<std::endl;
                //std::cout<<pre_bwt[i+2]<<" "<<pre_bwt[i+3]<<std::endl;

                //copy from hocc+bwt
                while(pbwt_freq>0) {

                    sym = hocc[k];
                    freq = hocc[k+1];
                    //std::cout<<"pbwt freq: "<<pbwt_freq<<" hocc freq: "<<freq<<" sym: "<<sym<<" ("<<k<<","<<k+1<<") "<<rank<<std::endl;
                    if(freq<=pbwt_freq){
                        pbwt_freq-=freq;
                        k+=2;
                    }else{
                        assert(hocc[k+1]>pbwt_freq);
                        hocc[k+1]-=pbwt_freq;
                        //std::cout<<hocc[k+1]<<" "<<pbwt_freq<<std::endl;
                        freq = std::exchange(pbwt_freq, 0);
                    }

                    if(sym==dummy_sym) {//from bwt
                        //std::cout<<"dummy sym in hocc has freq "<<bwt[j+1]<<" in bwt ("<<j<<", "<<(j+1)<<") "<<rank<<" -> "<<acc<<" "<<bwt_buckets[rank]<<std::endl;
                        extract_rl_syms(bwt, new_bwt, j, freq, pos);
                        //acc+=freq;
                    }else{//from hocc
                        if(pos>0 && new_bwt[pos-2]==sym){
                            new_bwt[pos-1] += freq;
                        } else {
                            new_bwt.push_back(sym);
                            new_bwt.push_back(freq);
                            pos+=2;
                        }
                    }
                    new_bwt_size+=freq;
                }
            }else{//copy from bwt
                extract_rl_syms(bwt, new_bwt, j, pbwt_freq, pos);
                //acc+=pbwt_freq;
                new_bwt_size+=pbwt_freq;
            }
            rank++;
        }else{
            if(pos>0 && new_bwt[pos-2]==sym){
                new_bwt[pos-1] += pbwt_freq;
            }else{
                new_bwt.push_back(sym);
                new_bwt.push_back(pbwt_freq);
                pos+=2;
            }
            new_bwt_size+=pbwt_freq;
        }
        i+=2;
    }
    std::cout<<"    Stats:       "<<std::endl;
    std::cout<<"      BWT size (n):       "<<new_bwt_size<<std::endl;
    std::cout<<"      Number of runs (r): "<<new_bwt.size()/2<<std::endl;
    std::cout<<"      r/n:                "<<(double(new_bwt.size())/2)/double(new_bwt_size)<<std::endl;

    bwt.close(true);
    pre_bwt.close(true);
    new_bwt.close();
}

void parse2bwt(sdsl::cache_config& config, size_t p_round){

    std::string parse_file = sdsl::cache_file_name("tmp_input", config);
    std::ifstream c_vec(parse_file, std::ifstream::binary);
    c_vec.seekg(0, std::ifstream::end);
    size_t tot_bytes = c_vec.tellg();
    c_vec.seekg(0, std::ifstream::beg);

    auto *buffer = reinterpret_cast<size_t*>(malloc(BUFFER_SIZE));
    size_t read_bytes =0;

    size_t pos = 0, len = tot_bytes/sizeof(size_t);
    vector_t bwt(len*2, 0, sdsl::bits::hi(len)+1);
    while(read_bytes<tot_bytes){
        c_vec.read((char *) buffer, BUFFER_SIZE);
        read_bytes+=c_vec.gcount();
        assert((c_vec.gcount() % sizeof(size_t))==0);
        for(size_t i=0;i<c_vec.gcount()/sizeof(size_t);i++){

            if(pos==0 || buffer[i]!=bwt[pos-2]){
                bwt[pos++] = buffer[i];
                bwt[pos++] = 1;
            }else{
                assert(bwt[pos-2]==buffer[i]);
                bwt[pos-1]++;
            }
        }
    }
    bwt.resize(pos);
    sdsl::store_to_file(bwt, sdsl::cache_file_name("bwt_lvl_"+std::to_string(p_round), config));
    c_vec.close();
    free(buffer);
    if(remove(parse_file.c_str())){
        std::cout<<"Error trying to delete file "<<parse_file<<std::endl;
    }
    std::cout<<"    Size of the BWT:       "<<len<<std::endl;
    std::cout<<"    Size run-length rep.:  "<<pos/2<<std::endl;
}

void infer_bwt(sdsl::cache_config& config, size_t p_round){
    std::cout<<"Inferring the BWT"<<std::endl;

    std::cout<<"  Computing the BWT for parse "<<p_round<<std::endl;
    parse2bwt(config, p_round);

    while(p_round-->0){
        std::cout<<"  Inducing the BWT for parse "<<p_round<<std::endl;
        auto start = std::chrono::steady_clock::now();
        infer_lvl_bwt(config, p_round);
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 4);
        report_mem_peak();
    }
}

void g_bwt_algo(std::string &i_file, std::string& o_file, std::string& tmp_folder, size_t n_threads, float hbuff_frac) {

    auto alphabet = get_alphabet(i_file);
    size_t n_chars = 0;
    for (auto const &sym : alphabet) n_chars += sym.second;
    auto hbuff_size = std::max<size_t>(64 * n_threads, size_t(std::ceil(float(n_chars) * hbuff_frac)));

    sdsl::cache_config config(false, tmp_folder);

    std::string rules_file = sdsl::cache_file_name("m_rules", config);
    std::string rules_len_file = sdsl::cache_file_name("rules_len", config);
    gram_info_t p_gram(rules_file, rules_len_file);
    p_gram.sigma = alphabet.size();
    p_gram.last_dict_size = p_gram.sigma;
    for(auto & sym : alphabet){
        p_gram.sym_map[sym.first] = sym.first;
    }
    p_gram.max_tsym = alphabet.back().first;
    p_gram.r = p_gram.max_tsym + 1;

    size_t p_rounds = build_lc_gram<lms_parsing>(i_file, n_threads, hbuff_size, p_gram, alphabet, config);
    infer_bwt(config, p_rounds);
    //check_plain_grammar(p_gram, i_file);
}

