//
// Created by Diaz, Diego on 23.11.2021.
//
#include "gbwt.hpp"
#include "utils.hpp"
#include "bwt_io.h"
#ifdef __linux__
#include <malloc.h>
#endif

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
void extract_rl_syms(bwt_buff_writer& bwt_buff, bwt_buff_writer& new_bwt_buff, size_t& j, size_t freq){

    size_t tmp_freq, sym;

    while(freq>0) {

        bwt_buff.read_run(j, sym, tmp_freq);

        if(tmp_freq<=freq){
            freq-=tmp_freq;
            j++;
        }else{
            bwt_buff.dec_freq(j, freq);
            tmp_freq = std::exchange(freq, 0);
        }

        if(new_bwt_buff.size()>0 && new_bwt_buff.last_sym()==sym){
            new_bwt_buff.inc_freq_last(tmp_freq);
        } else {
            new_bwt_buff.push_back(sym, tmp_freq);
        }
    }
}

size_t compute_hocc_size(dictionary& dict, bv_rs_t& hocc_rs, vector_t& hocc_buckets,
                         size_t p_round, sdsl::cache_config & config){

    std::string prev_bwt_f = sdsl::cache_file_name("bwt_lev_"+std::to_string(p_round+1), config);
    bwt_buff_reader bwt_buff(prev_bwt_f);

    size_t sym, pos, dummy_sym=dict.alphabet+2, left_sym, freq=0, al_b, fr_b, bps;
    al_b = INT_CEIL(sdsl::bits::hi(dummy_sym)+1, 8);
    fr_b = INT_CEIL(sdsl::bits::hi(dict.max_sym_freq)+1,8);
    bps = al_b + fr_b;
    size_t tot_bytes = bps * hocc_rs(dict.phrases_has_hocc.size());
    auto * hocc_counts = (char *)malloc(tot_bytes);
    memset(hocc_counts, 0, tot_bytes);

    char * ptr;

    for(size_t i=0;i<bwt_buff.size();i++) {

        sym=bwt_buff.read_sym(i);
        if(dict.phrases_has_hocc[sym]){
            ptr = hocc_counts + hocc_rs(sym)*bps;
            if(memcmp(ptr, &dummy_sym, al_b)!=0){
                memcpy(ptr, &dummy_sym, al_b);
                memcpy(&freq,ptr+al_b, fr_b);
                freq++;
                memcpy(ptr+al_b, &freq, fr_b);
            }
        }

        pos = 2*sym+1;
        sym = dict.dict[pos];
        while(sym>=dict.alphabet){
            sym -= dict.alphabet;
            left_sym = dict.dict[pos-1]+1;
            ptr = hocc_counts + hocc_rs(sym)*bps;
            if(memcmp(ptr, &left_sym, al_b)!=0){
                memcpy(ptr, &left_sym, al_b);
                memcpy(&freq,ptr+al_b, fr_b);
                freq++;
                memcpy(ptr+al_b, &freq, fr_b);
            }
            pos = 2*sym+1;
            sym = dict.dict[pos];
        }
    }

    bwt_buff.close();
    size_t acc=0;
    ptr = hocc_counts;
    for(size_t i=0;i<hocc_buckets.size()-1;i++){
        hocc_buckets[i] = acc;
        memcpy(&freq,ptr+al_b, fr_b);
        acc+=freq;
        ptr+=bps;
    }

    assert(acc<=dict.t_size);
    hocc_buckets[hocc_buckets.size()-1] = acc;
    free(hocc_counts);
    return acc;
}

void infer_lvl_bwt(sdsl::cache_config& config, size_t p_round) {

    dictionary dict;
    std::string dict_file = sdsl::cache_file_name("dict_lev_"+std::to_string(p_round), config);
    sdsl::load_from_file(dict, dict_file);
    bv_rs_t hocc_rs(&dict.phrases_has_hocc);

    size_t sym, left_sym, pos, freq, rank, dummy_sym = dict.alphabet+2;

    std::cout<<"    Computing the number of induced symbols"<<std::endl;
    vector_t hocc_buckets(hocc_rs(dict.phrases_has_hocc.size())+1, 0, sdsl::bits::hi(dict.t_size)+1);
    size_t n_runs = compute_hocc_size(dict, hocc_rs, hocc_buckets, p_round, config);

    size_t al_b = INT_CEIL(sdsl::bits::hi(dummy_sym)+1, 8);
    size_t fr_b = INT_CEIL(sdsl::bits::hi(dict.max_sym_freq)+1,8);
    size_t bps = al_b + fr_b;

    auto * hocc = (char *)malloc(n_runs*bps);
    memset(hocc, 0, n_runs*bps);
    char *hocc_ptr;

    std::string prev_bwt_f = sdsl::cache_file_name("bwt_lev_"+std::to_string(p_round+1), config);
    bwt_buff_writer bwt_buff(prev_bwt_f, std::ios::in);

    std::cout<<"    Performing the induction from the previous BWT"<<std::endl;
    for(size_t i=0;i<bwt_buff.size();i++) {

        bwt_buff.read_run(i, sym, freq);

        if(dict.phrases_has_hocc[sym]){
            rank = hocc_rs(sym);
            hocc_ptr = hocc + hocc_buckets[rank]*bps;

            if(hocc_ptr!=hocc && memcmp(hocc_ptr-bps, &dummy_sym, al_b)==0){
                size_t new_freq=0;
                memcpy(&new_freq, hocc_ptr-fr_b, fr_b);
                new_freq+=freq;
                memcpy(hocc_ptr-fr_b, &new_freq, fr_b);
            }else{
                memcpy(hocc_ptr, &dummy_sym, al_b);
                memcpy(hocc_ptr+al_b, &freq, fr_b);
                hocc_buckets[rank]++;
            }
        }

        pos = 2*sym+1;
        sym = dict.dict[pos];
        while(sym>=dict.alphabet){
            sym -= dict.alphabet;
            left_sym = dict.dict[pos-1];
            assert(left_sym<dict.alphabet && dict.phrases_has_hocc[sym]);

            rank = hocc_rs(sym);
            hocc_ptr = hocc + hocc_buckets[rank]*bps;
            left_sym++;

            if(hocc_ptr!=hocc && memcmp(hocc_ptr-bps, &left_sym, al_b)==0){
                size_t new_freq=0;
                memcpy(&new_freq, hocc_ptr-fr_b, fr_b);
                new_freq+=freq;
                memcpy(hocc_ptr-fr_b, &new_freq, fr_b);
            }else{
                memcpy(hocc_ptr, &left_sym, al_b);
                memcpy(hocc_ptr+al_b, &freq, fr_b);
                hocc_buckets[rank]++;
            }
            pos = 2*sym+1;
            sym = dict.dict[pos];
        }
        assert(dict.dict[pos-1]==dict.dict_dummy);
        bwt_buff.write_sym(i, sym);
    }

    std::cout<<"dictionary: "<<double(sdsl::size_in_bytes(dict.dict))/1000000<<std::endl;
    std::cout<<"hocc_buckets: "<<double(sdsl::size_in_bytes(hocc_buckets))/1000000<<std::endl;
    std::cout<<"has_hocc: "<<double(sdsl::size_in_bytes(dict.phrases_has_hocc))/1000000<<std::endl;
    std::cout<<"hocc_rs: "<<double(sdsl::size_in_bytes(hocc_rs))/1000000<<std::endl;
    std::cout<<"hocc: "<<double(n_runs*bps)/1000000<<std::endl;

    sdsl::util::clear(dict.dict);
    sdsl::util::clear(hocc_buckets);
    sdsl::util::clear(hocc_rs);

#ifdef __linux__
    malloc_trim(0);
#endif
    std::cout<<"    Assembling the new BWT"<<std::endl;
    std::string new_bwt_f = sdsl::cache_file_name("bwt_lev_"+std::to_string(p_round), config);
    if(dict.alphabet<dict.p_alpha_size){
        al_b = INT_CEIL(sdsl::bits::hi(dict.p_alpha_size+2)+1, 8);
    }
    bwt_buff_writer new_bwt_buff(new_bwt_f, std::ios::out, al_b, fr_b);

    std::string p_bwt_file = sdsl::cache_file_name("pre_bwt_lev_"+std::to_string(p_round), config);
    bwt_buff_reader p_bwt(p_bwt_file);

    size_t i=0, j=0, pbwt_freq, new_bwt_size=0;
    rank = 0;
    hocc_ptr = hocc;
    char *last = hocc + n_runs*bps;
    freq=0;

    while(i<p_bwt.size()) {

        p_bwt.read_run(i, sym, pbwt_freq);

        if(sym==(dummy_sym-1)) {

            if(dict.phrases_has_hocc[rank]) {

                //copy from hocc+bwt
                while(pbwt_freq>0) {
                    memcpy(&sym, hocc_ptr, al_b);
                    memcpy(&freq, hocc_ptr+al_b, fr_b);

                    while(sym==0 && hocc_ptr!=last){
                        hocc_ptr+=bps;
                        memcpy(&sym, hocc_ptr, al_b);
                        memcpy(&freq, hocc_ptr+al_b, fr_b);
                    }

                    if(sym==0){
                        std::cout<<(hocc_ptr==last)<<" "<<bps<<" "<<fr_b<<" "<<al_b<<std::endl;
                    }
                    assert(sym>0);

                    if(freq<=pbwt_freq){
                        pbwt_freq-=freq;
                        hocc_ptr+=bps;
                    }else{
                        size_t tmp=0;
                        memcpy(&tmp, hocc_ptr+al_b, fr_b);
                        assert(tmp>pbwt_freq);
                        tmp -=pbwt_freq;
                        memcpy(hocc_ptr+al_b, &tmp, fr_b);
                        freq = std::exchange(pbwt_freq, 0);
                    }

                    sym--;
                    if(sym==(dummy_sym-1)) {//from bwt
                        extract_rl_syms(bwt_buff, new_bwt_buff, j, freq);
                    }else{//from hocc
                        if(new_bwt_buff.size()>0 && new_bwt_buff.last_sym()==sym){
                            new_bwt_buff.inc_freq_last(freq);
                        } else {
                            assert(sym<dict.alphabet);
                            new_bwt_buff.push_back(sym, freq);
                        }
                    }
                    new_bwt_size+=freq;
                }
            }else{//copy from bwt
                extract_rl_syms(bwt_buff, new_bwt_buff, j, pbwt_freq);
                new_bwt_size+=pbwt_freq;
            }
            rank++;
        }else{
            if(new_bwt_buff.size()>0 && new_bwt_buff.last_sym()==sym){
                new_bwt_buff.inc_freq_last(pbwt_freq);
            } else {
                assert(sym<dict.alphabet);
                new_bwt_buff.push_back(sym, pbwt_freq);
            }
            new_bwt_size+=pbwt_freq;
        }
        i++;
    }

    p_bwt.close();
    bwt_buff.close();
    new_bwt_buff.close();

    std::cout<<"    Stats:       "<<std::endl;
    std::cout<<"      BWT size (n):       "<<new_bwt_size<<std::endl;
    std::cout<<"      Number of runs (r): "<<new_bwt_buff.size()<<std::endl;
    std::cout<<"      n/r:                "<<double(new_bwt_size)/double(new_bwt_buff.size())<<std::endl;
    std::cout<<"      bytes per run:      "<<bps<<std::endl;
    std::cout<<"        sym bytes:        "<<al_b<<std::endl;
    std::cout<<"        freq bytes:       "<<fr_b<<std::endl;
    free(hocc);
}

void parse2bwt(sdsl::cache_config& config, size_t p_round) {

    std::string parse_file = sdsl::cache_file_name("tmp_input", config);
    std::ifstream c_vec(parse_file, std::ifstream::binary);
    c_vec.seekg(0, std::ifstream::end);
    size_t tot_bytes = c_vec.tellg();
    c_vec.seekg(0, std::ifstream::beg);

    auto *buffer = reinterpret_cast<size_t*>(malloc(BUFFER_SIZE));
    size_t read_bytes =0;

    size_t len = tot_bytes/sizeof(size_t);

    std::string bwt_lev_file = sdsl::cache_file_name("bwt_lev_"+std::to_string(p_round), config);
    size_t sb = INT_CEIL(sdsl::bits::hi(len)+1, 8);
    size_t fb = INT_CEIL(sdsl::bits::hi(len)+1, 8);
    bwt_buff_writer bwt_buff(bwt_lev_file, std::ios::out, sb, fb);

    while(read_bytes<tot_bytes){
        c_vec.read((char *) buffer, BUFFER_SIZE);
        read_bytes+=c_vec.gcount();
        assert((c_vec.gcount() % sizeof(size_t))==0);
        for(size_t i=0;i<c_vec.gcount()/sizeof(size_t);i++){
            if(bwt_buff.size()==0 || buffer[i]!=bwt_buff.last_sym()){
                bwt_buff.push_back(buffer[i], 1);
            }else{
                assert(buffer[i]==bwt_buff.last_sym());
                bwt_buff.inc_freq_last(1);
            }
        }
    }
    c_vec.close();
    free(buffer);
    bwt_buff.close();

    if(remove(parse_file.c_str())){
        std::cout<<"Error trying to delete file "<<parse_file<<std::endl;
    }
    std::cout<<"    Size of the BWT:       "<<len<<std::endl;
    std::cout<<"    Size run-length rep.:  "<<bwt_buff.size()<<std::endl;
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

    size_t p_rounds = build_lc_gram<lms_parsing>(i_file, n_threads, hbuff_size, alphabet, config);
    infer_bwt(config, p_rounds);
}

