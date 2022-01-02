//
// Created by Diaz, Diego on 23.11.2021.
//
#include "grammar_build.hpp"

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

    size_t curr_sym, start, end;
    std::stack<size_t> stack;

    size_t f = r_lim_ss(p_gram.r - 1) + 1;
    size_t l = r_lim_ss(p_gram.r);
    size_t idx = 0;

    for(size_t i=f; i <= l; i++){
        //std::cout<<i<<" "<<l<<std::endl;
        tmp_decomp.clear();
        stack.push(r[i]);
        assert(stack.size()<=if_stream.size());

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
                assert(p_gram.sym_map[curr_sym]==if_stream.read(idx));
                idx++;
            }else{
                if(p_gram.is_rl(curr_sym)){
                    assert(end-start+1==2);
                    for(size_t k=0;k<r[end];k++){
                        stack.push(r[start]);
                    }
                }else{
                    for(size_t j=end+1; j-->start;){
                        stack.push(r[j]);
                    }
                }
            }
        }
    }
    std::cout<<"\tGrammar is correct!!"<<std::endl;
}

void run_length_compress(gram_info_t &p_gram, sdsl::cache_config& config) {

    std::cout<<"  Run-length compressing the grammar"<<std::endl;

    ivb_t rules(p_gram.rules_file, std::ios::in);
    sdsl::int_vector_buffer<1> r_lim(p_gram.rules_lim_file, std::ios::in);

    std::string rl_rules_file = sdsl::cache_file_name("tmp_rl_file", config);
    ivb_t rl_rules(rl_rules_file, std::ios::out);
    std::string rl_r_lim_file = sdsl::cache_file_name("tmp_rl_lim_file", config);
    sdsl::int_vector_buffer<1> rl_r_lim(rl_r_lim_file, std::ios::out);

    size_t run_len=1;
    size_t new_id = p_gram.r-1;
    size_t tmp_sym;
    phrase_map_t ht;
    string_t pair(2,sdsl::bits::hi(rules.size())+1);

    for(size_t i=0;i<=p_gram.max_tsym;i++){
        rl_rules.push_back(i);
        rl_r_lim.push_back(true);
    }

    size_t i=p_gram.max_tsym+2;
    while(i<=rules.size()-p_gram.c){
        if(rules[i]!=rules[i-1] || r_lim[i-1]){
            if(run_len>1){
                pair.write(0, rules[i-1]);
                pair.write(1, run_len);
                auto res = ht.insert(pair.data(), pair.n_bits(), 0);
                if(res.second){
                    tmp_sym = new_id++;
                    ht.insert_value_at(res.first, tmp_sym);
                }else{
                    tmp_sym = 0;
                    ht.get_value_from(res.first, tmp_sym);
                    //tmp_sym = res.first.value();
                }
            }else{
                tmp_sym = rules[i-1];
            }
            run_len=0;

            rl_rules.push_back(tmp_sym);
            if(r_lim[i-1]){
                rl_r_lim[rl_rules.size()-1] = true;
            }
        }
        run_len++;
        i++;
    }

    const bitstream<phrase_map_t::buff_t>& stream = ht.get_data();
    key_wrapper key_w{pair.width(), ht.description_bits(), stream};
    for(auto const& phrase : ht){
        rl_rules.push_back(key_w.read(phrase, 0));
        rl_rules.push_back(key_w.read(phrase, 1));
        rl_r_lim[rl_rules.size()-1] = true;
    }

    for(size_t k=rules.size()-p_gram.c;k<rules.size();k++){
        rl_rules.push_back(rules[k]);
    }
    rl_r_lim[rl_rules.size()-1] = true;

    p_gram.rl_rules.first = p_gram.r-1;
    p_gram.rl_rules.second = ht.size();
    //p_gram.rules_breaks.push_back(p_gram.rules_breaks.back() + ht.size());

    p_gram.r+= ht.size();
    p_gram.g = rl_rules.size();

    std::cout<<"    Stats:"<<std::endl;
    std::cout<<"      Grammar size before:        "<<rules.size()<<std::endl;
    std::cout<<"      Grammar size after:         "<<rl_rules.size()<<std::endl;
    std::cout<<"      Number of new nonterminals: "<<ht.size()<<std::endl;
    std::cout<<"      Compression ratio:          "<<float(rl_rules.size())/float(rules.size())<<std::endl;

    rules.close();
    r_lim.close();
    rl_rules.close();
    rl_r_lim.close();

    rename(rl_rules_file.c_str(), p_gram.rules_file.c_str());
    rename(rl_r_lim_file.c_str(), p_gram.rules_lim_file.c_str());
}

bv_t mark_disposable_symbols(gram_info_t &p_gram) {

    bv_t r_lim;
    sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
    bv_ss_t r_lim_ss(&r_lim);

    sdsl::int_vector<> rules;
    sdsl::load_from_file(rules, p_gram.rules_file);

    bv_t rem_nts(p_gram.r, false);

    //compute which nonterminals are repeated and
    // which have a replacement of length 1
    sdsl::int_vector<2> rep_nts(p_gram.r + 1, 0);

    size_t max_tsym = p_gram.max_tsym;
    size_t r_len=1, curr_rule=max_tsym+1,k=max_tsym+1;
    while(k<rules.size()){
        if(p_gram.is_rl(curr_rule)){ //run-length compressed rules
            rep_nts[rules[k++]] = 2;
            assert(r_lim[k] && r_len==1);
            r_len=0;
            curr_rule++;
        }else{
            //get the frequency of every symbol
            if(rep_nts[rules[k]]<2){
                rep_nts[rules[k]]++;
            }
            if(r_lim[k]){
                r_len=0;
                curr_rule++;
            }
        }
        r_len++;
        k++;
    }

    //mark the rules to remove
    //1) m_rules whose left-hand side has length one
    //2) terminal symbols between [min_sym..max_sym] with
    // frequency zero: to compress the alphabet
    for(size_t i=0;i<p_gram.rules_breaks[p_gram.n_p_rounds];i++){
        //mark the m_rules with frequency one
        if(!rem_nts[i]){
            rem_nts[i] = rep_nts[i]==0 || (rep_nts[i]==1 && i > max_tsym);
        }
    }

    if(p_gram.sp_rules.second>0){
        //assert(p_gram.sp_rules.first+p_gram.sp_rules.second==p_gram.r-1);
        for(size_t i=p_gram.sp_rules.first;i<p_gram.sp_rules.first+p_gram.sp_rules.second;i++){
            //mark the rules with frequency one
            if(!rem_nts[i]){
                rem_nts[i] = rep_nts[i]==0 || (rep_nts[i]==1 && i > max_tsym);
            }
        }
    }
    return rem_nts;
}

alpha_t get_alphabet(std::string &i_file) {

    std::cout << "  Reading input file:" << std::endl;

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

    std::cout<<"    Number of characters: "<< if_stream.size() << std::endl;
    std::cout<<"    Number of strings:    "<< alphabet[0].second << std::endl;
    std::cout<<"    Alphabet:             "<< alphabet.size() << std::endl;
    std::cout<<"    Smallest symbol:      "<< (int) alphabet[0].first << std::endl;
    std::cout<<"    Greatest symbol:      "<< (int) alphabet.back().first << std::endl;

    if (if_stream.read(if_stream.size() - 1) != alphabet[0].first) {
        std::cout << "Error: sep. symbol " << alphabet[0].first << " differs from last symbol in file "
                  << if_stream.read(if_stream.size() - 1) << std::endl;
        exit(1);
    }
    return alphabet;
}
void decomp(size_t nt, sdsl::int_vector<> &rules, bv_ss_t &rlim_ss, bv_t &rem_nt,
            bv_rs_t &rem_nt_rs, ivb_t &dec) {

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
            dec.push_back(tmp - rem_nt_rs(tmp));
        }
    }
}

void simplify_grammar(gram_info_t &p_gram, bool full_simplification) {


    if(full_simplification) {

        std::cout<<"  Simplifying the grammar"<<std::endl;

        bv_t rem_nts = mark_disposable_symbols(p_gram);
        bv_rs_t rem_nts_rs(&rem_nts);

        bv_t r_lim;
        sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
        bv_ss_t r_lim_ss(&r_lim);
        sdsl::int_vector<> rules;
        sdsl::load_from_file(rules, p_gram.rules_file);
        size_t max_tsym = p_gram.max_tsym;

        ivb_t new_rules(p_gram.rules_file, std::ios::out);
        sdsl::int_vector_buffer<1> new_r_lim(p_gram.rules_lim_file, std::ios::out);

        //compress the alphabet
        size_t cont=0;
        uint8_t byte_sym;
        std::unordered_map<size_t, uint8_t> new_sym_map;
        for(size_t k=0;k<=p_gram.max_tsym;k++){
            auto res = p_gram.sym_map.find(k);
            if(res!=p_gram.sym_map.end()){
                byte_sym = res->second;
                new_sym_map.insert({cont, byte_sym});
                cont++;
            }
        }
        std::swap(new_sym_map, p_gram.sym_map);
        assert(p_gram.sym_map.size()==p_gram.sigma);

        for(size_t k=0;k<p_gram.sigma;k++){
            new_r_lim.push_back(true);
            new_rules.push_back(k);
        }

        size_t pos, tr_rule=p_gram.sigma, c_start;
        for(size_t i=max_tsym+1,curr_rule=max_tsym+1;i<rules.size();curr_rule++){
            assert(r_lim[i-1]);
            pos = i;
            while(!r_lim[i]) i++;
            i++;

            if((i-pos)==p_gram.c) c_start = new_rules.size();
            if(!rem_nts[curr_rule]){
                if(!p_gram.is_rl(curr_rule)){//regular rule
                    for(size_t j=pos;j<i;j++){
                        if(rem_nts[rules[j]]){
                            decomp(rules[j], rules, r_lim_ss, rem_nts, rem_nts_rs, new_rules);
                        }else{
                            new_rules.push_back(rules[j]-rem_nts_rs(rules[j]));
                        }
                    }
                }else{//run-length rule
                    assert((i-pos)==2);
                    new_rules.push_back(rules[pos]-rem_nts_rs(rules[pos]));
                    new_rules.push_back(rules[pos+1]);
                }
                new_r_lim[new_rules.size()-1]=true;
                tr_rule++;
            }
        }

        size_t rm_nt =rem_nts_rs(rem_nts.size());
        float rm_per = float(rm_nt)/float(p_gram.r)*100;
        float comp_rat = float(new_rules.size())/float(rules.size());

        std::cout<<"    Stats:"<<std::endl;
        std::cout<<"      Grammar size before:  "<<p_gram.g<<std::endl;
        std::cout<<"      Grammar size after:   "<<new_rules.size()<<std::endl;
        std::cout<<"      Deleted nonterminals: "<<rm_nt<<" ("<<rm_per<<"%)"<<std::endl;
        std::cout<<"      Compression ratio:    "<<comp_rat<<std::endl;

        p_gram.c = new_rules.size()-c_start;
        p_gram.r -= rm_nt;
        p_gram.g = new_rules.size();

        for(auto &sym : p_gram.rules_breaks){
            sym = sym - rem_nts_rs(sym);
        }

        if(p_gram.rl_rules.second>0){
            size_t last_rl = p_gram.rl_rules.first+p_gram.rl_rules.second-1;
            p_gram.rl_rules.first -= rem_nts_rs(p_gram.rl_rules.first);
            last_rl -= rem_nts_rs(last_rl);
            p_gram.rl_rules.second = last_rl - p_gram.rl_rules.first + 1;
        }

        if(p_gram.sp_rules.second>0){
            size_t last_sp = p_gram.sp_rules.first+p_gram.sp_rules.second-1;
            p_gram.sp_rules.first -= rem_nts_rs(p_gram.sp_rules.first);
            last_sp -= rem_nts_rs(last_sp);
            p_gram.sp_rules.second = last_sp - p_gram.sp_rules.first + 1;
        }
        new_rules.close();
        new_r_lim.close();
    }else{

        std::cout<<"  Compressing the alphabet of terminals"<<std::endl;

        //compress the alphabet
        size_t cont=0;
        uint8_t byte_sym;
        std::unordered_map<size_t, uint8_t> new_sym_map;
        uint8_t inv_sym_map[255]={0};
        for(size_t k=0;k<=p_gram.max_tsym;k++){
            auto res = p_gram.sym_map.find(k);
            if(res!=p_gram.sym_map.end()){
                byte_sym = res->second;
                new_sym_map.insert({cont, byte_sym});
                inv_sym_map[byte_sym] = cont;
                cont++;
            }
        }
        std::swap(new_sym_map, p_gram.sym_map);
        assert(p_gram.sym_map.size()==p_gram.sigma);

        sdsl::int_vector<> rules;
        sdsl::load_from_file(rules, p_gram.rules_file);
        sdsl::int_vector_buffer<> new_rules(p_gram.rules_file, std::ios::out);

        bv_t r_lim;
        sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
        sdsl::int_vector_buffer<1> new_r_lim(p_gram.rules_lim_file, std::ios::out);

        for(size_t k=0;k<p_gram.sigma;k++){
            new_rules.push_back(k);
            new_r_lim.push_back(true);
        }

        size_t delta = (p_gram.max_tsym+1)-p_gram.sigma, nt=p_gram.max_tsym+1, i=p_gram.sigma, j=p_gram.max_tsym+1;
        while(j<rules.size()){
            if(rules[j]<=p_gram.max_tsym){
                new_rules[i] = inv_sym_map[rules[j]];
            } else {
                new_rules[i] = rules[j]-delta;
            }
            if(p_gram.is_rl(nt)){
                new_r_lim[i] = r_lim[j];
                new_rules[++i] = rules[++j];
            }
            new_r_lim[i] = r_lim[j];
            if(r_lim[j]) nt++;
            i++;
            j++;
        }

        new_rules.close();
        new_r_lim.close();
        p_gram.r -= delta;
        p_gram.g = new_rules.size();

        for(auto &sym : p_gram.rules_breaks){
            sym = sym - delta;
        }

        if(p_gram.rl_rules.second>0){
            size_t last_rl = p_gram.rl_rules.first+p_gram.rl_rules.second-1;
            p_gram.rl_rules.first -= delta;
            last_rl -= delta;
            p_gram.rl_rules.second = last_rl - p_gram.rl_rules.first + 1;
        }

        if(p_gram.sp_rules.second>0){
            size_t last_sp = p_gram.sp_rules.first+p_gram.sp_rules.second-1;
            p_gram.sp_rules.first -= delta;
            last_sp -= delta;
            p_gram.sp_rules.second = last_sp - p_gram.sp_rules.first + 1;
        }
    }
}

void build_gram(std::string &i_file, std::string &p_gram_file,
                uint8_t comp_lvl, std::string& tmp_folder,
                size_t n_threads, float hbuff_frac) {

    auto alphabet = get_alphabet(i_file);
    size_t n_chars = 0;
    for (auto const &sym : alphabet) n_chars += sym.second;
    auto hbuff_size = std::max<size_t>(64 * n_threads, size_t(std::ceil(float(n_chars) * hbuff_frac)));

    sdsl::cache_config config(false, tmp_folder);
    std::string g_info_file = sdsl::cache_file_name("g_info_file", config);

    std::string rules_file = sdsl::cache_file_name("m_rules", config);
    std::string rules_len_file = sdsl::cache_file_name("rules_len", config);

    gram_info_t p_gram(rules_file, rules_len_file);
    p_gram.sigma = alphabet.size();
    for(auto & sym : alphabet){
        p_gram.sym_map[sym.first] = sym.first;
    }
    p_gram.max_tsym = alphabet.back().first;
    p_gram.r = p_gram.max_tsym + 1;

    build_lc_gram<lms_parsing>(i_file, n_threads, hbuff_size, p_gram, alphabet, config);
    //run_length_compress(p_gram, config);
    suffpair(p_gram, config, n_threads, hbuff_size);
    simplify_grammar(p_gram, true);
    //assert(p_gram.r-1==p_gram.rules_breaks[p_gram.n_p_rounds + 2]);
    check_plain_grammar(p_gram, i_file);

    std::cout<<"  Final grammar: " << std::endl;
    std::cout<<"    Number of terminals:            "<< (size_t) p_gram.sigma << std::endl;
    std::cout<<"    Number of nonterminals:         "<< p_gram.r - p_gram.sigma<<std::endl;
    std::cout<<"    Grammar size:                   "<< p_gram.g<<std::endl;
    std::cout<<"    Breakdown:"<<std::endl;

    for(size_t i=0;i<p_gram.n_p_rounds;i++){
        std::cout<<"      Rules of parsing round "<<(i+1);
        if(i<9) std::cout<<" ";
        std::cout<<":    "<<p_gram.rules_breaks[i+1]-p_gram.rules_breaks[i]<<std::endl;
    }
    std::cout<<"      Run-length rules:             "<<p_gram.rl_rules.second<<std::endl;
    std::cout<<"      SuffPair rules:               "<<p_gram.sp_rules.second<<std::endl;

    std::cout<<"  Compression stats: " << std::endl;
    std::cout<<"    Text size in MB:        " << double(n_chars)/1000000<<std::endl;
    std::cout<<"    Grammar size in MB:     " << INT_CEIL(p_gram.g*(sdsl::bits::hi(p_gram.r)+1),8)/double(1000000)<< std::endl;
    std::cout<<"    Compression ratio:      " << INT_CEIL(p_gram.g*(sdsl::bits::hi(p_gram.r)+1),8)/double(n_chars) << std::endl;

    std::cout<<"  Storing the final grammar in " << p_gram_file <<std::endl;

    if(comp_lvl==1){
        grammar<sdsl::int_vector<>> final_gram(p_gram, n_chars, alphabet[0].second);
        sdsl::store_to_file(final_gram, p_gram_file);
        final_gram.space_breakdown();
    }else if(comp_lvl==2){
        grammar<huff_vector<>> final_gram(p_gram, n_chars, alphabet[0].second);
        sdsl::store_to_file(final_gram, p_gram_file);
        final_gram.space_breakdown();
    }
}

