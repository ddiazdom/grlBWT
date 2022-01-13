//
// Created by Diaz, Diego on 23.11.2021.
//
#include "gbwt.hpp"

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

void sp_sa2bwt(std::string& sp_sa, std::string& output){

}

void induce_sym_order(gram_info_t& p_gram, size_t p_round){
    //TODO perform counting sort
}

void gram2sp_sa(gram_info_t& p_gram){
    for(size_t i=p_gram.n_p_rounds-1;i-->0;){
        induce_sym_order(p_gram, i);
    }
}

void g_bwt_algo(std::string &i_file, std::string& o_file, std::string& tmp_folder, size_t n_threads, float hbuff_frac) {

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
    p_gram.last_dict_size = p_gram.sigma;

    for(auto & sym : alphabet){
        p_gram.sym_map[sym.first] = sym.first;
    }
    p_gram.max_tsym = alphabet.back().first;
    p_gram.r = p_gram.max_tsym + 1;

    build_lc_gram<lms_parsing>(i_file, n_threads, hbuff_size, p_gram, alphabet, config);
    //gram2sp_sa(p_gram);
    //run_length_compress(p_gram, config);
    //check_plain_grammar(p_gram, i_file);
    //
    /*std::cout<<"  Final grammar: " << std::endl;
    std::cout<<"    Number of terminals:            "<< (size_t) p_gram.sigma << std::endl;
    std::cout<<"    Number of nonterminals:         "<< p_gram.r - p_gram.sigma<<std::endl;
    std::cout<<"    Grammar size:                   "<< p_gram.g<<std::endl;
    std::cout<<"    Breakdown:"<<std::endl;
    for(size_t i=0;i<p_gram.n_p_rounds;i++){
        std::cout<<"      Rules of parsing round "<<(i+1);
        if(i<9) std::cout<<" ";
        std::cout<<":    "<<p_gram.rules_breaks[i+1]-p_gram.rules_breaks[i]<<std::endl;
    }
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
    }*/
}

