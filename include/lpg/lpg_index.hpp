//
// Created by diediaz on 11-12-19.
//

#ifndef LMS_GRAMMAR_REP_HPP
#define LMS_GRAMMAR_REP_HPP

#include <pthread.h>
#include "lpg_build.hpp"
#include "grammar_tree.hpp"
#include "grid.hpp"

class lpg_index{

private:
    typedef typename lpg_build::plain_grammar_t         plain_grammar_t;
    typedef typename lpg_build::alpha_t                 alpha_t;

    grammar_tree_t                                      grammar_tree;
    grid_t                                              grid;
    sdsl::int_vector<8>                                 symbols_map; // map a compressed terminal to its original byte symbol
    uint8_t                                             m_sigma{}; //alphabet of terminal symbols
    uint8_t                                             parsing_rounds{}; //number of LMS parsing rounds during the grammar construction
    bool                                                rl_compressed{}; // is the grammar run-length compressed?

    void build_index(plain_grammar_t &p_gram, sdsl::cache_config& config){
        m_sigma = p_gram.sigma;
        parsing_rounds = p_gram.rules_per_level.size();
        symbols_map.resize(p_gram.symbols_map.size());
        for(size_t i=0; i < p_gram.symbols_map.size(); i++){
            symbols_map[i] = p_gram.symbols_map[i];
        }
        //build the grammar tree and grid from p_gram here!!
    };

    static alpha_t get_alphabet(std::string& i_file){

        std::cout<<"Reading input file"<<std::endl;

        //TODO this can be done in parallel if the input is too big
        size_t alph_frq[256] = {0};
        alpha_t alphabet;

        i_file_stream<uint8_t> if_stream(i_file, BUFFER_SIZE);
        for(size_t i=0;i<if_stream.tot_cells;i++){
            alph_frq[if_stream.read(i)]++;
        }

        for(size_t i=0;i<256;i++){
            if(alph_frq[i]>0) alphabet.emplace_back(i, alph_frq[i]);
        }

        std::cout<<"  Number of characters: "<<if_stream.size()<<std::endl;
        std::cout<<"  Alphabet:             "<<alphabet.size()<<std::endl;
        std::cout<<"  Smallest symbol:      "<<(int)alphabet[0].first<<std::endl;
        std::cout<<"  Greatest symbol:      "<<(int)alphabet.back().first<<std::endl;

        if(if_stream.read(if_stream.size()-1)!=alphabet[0].first){
            std::cout<<"Error: sep. symbol "<<alphabet[0].first <<" differs from last symbol in file "
                     <<if_stream.read(if_stream.size()-1)<<std::endl;
            exit(1);
        }
        return alphabet;
    }

    template<typename proc>
    static void lms_scan(proc& task, size_t pos, std::vector<size_t>& parse, size_t& n_lms){

        int_array<size_t> lms_phrase(2, 32);
        size_t curr_sym, prev_sym;
        bool s_type, prev_s_type;

        if(parse[pos]<parse[pos+1]){
            prev_s_type = S_TYPE;
        }else{
            prev_s_type = L_TYPE;
        }

        prev_sym = parse[pos];
        lms_phrase.push_back(prev_sym);
        n_lms=0;
        for(size_t i = pos; i-- > 0;){
            curr_sym = parse[i];
            if (curr_sym < prev_sym) {//S_TYPE type
                s_type = S_TYPE;
            } else if (curr_sym == prev_sym) {
                s_type = prev_s_type;
            } else {//L_TYPE type
                s_type = L_TYPE;
                if(prev_s_type == S_TYPE) {//LMS-type
                    lms_phrase.pop_back();
                    if(n_lms>0){
                        task(lms_phrase);
                    }
                    lms_phrase.clear();
                    lms_phrase.push_back(prev_sym);
                    n_lms++;
                }
            }
            lms_phrase.push_back(curr_sym);
            prev_sym = curr_sym;
            prev_s_type = s_type;
        }
    }

    std::pair<std::vector<size_t>, uint8_t> compute_pattern_cut(const std::string& pattern) const {

        std::vector<size_t> parse;
        std::vector<size_t> next_parse;

        //right elements are index in the first
        // level (terminals) where the S* string occur
        // left element is the index in the last level
        std::vector<std::pair<size_t, size_t>> lms_cuts;

        parse.reserve(pattern.size());
        next_parse.reserve(pattern.size()/2);
        for(auto const& sym : pattern) parse.push_back(sym);

        bit_hash_table<size_t, 44> ht;
        size_t n_lms, pos, rank=0;
        uint8_t lvl=0;

        auto hash_task = [&](auto& phrase){
            phrase.mask_tail();
            ht.insert(phrase.data(), phrase.n_bits(),0);
        };

        //hash the LMS phrases in the text
        while(true){
            assert(parse.size()>1);
            pos = parse.size()-2;
            while(pos>0 && parse[pos]==parse[pos+1]) pos--;
            if(pos==0) break;
            n_lms++;

            lms_scan(hash_task, pos, parse, n_lms);

            if(n_lms<=5){//report the cuts
                if(n_lms==0){
                    //go back to the previous parse
                }else{
                    std::vector<size_t> final_cuts;
                    for(auto const& cut: lms_cuts){
                        final_cuts.push_back(cut.second);
                    }
                    return std::make_pair(std::move(final_cuts), lvl);
                }
            }else if(lvl == parsing_rounds){
                //pattern is longer than the strings in the grammar
                return std::pair<std::vector<size_t>, uint8_t>();
            } else{
                //assign ranks to the lms phrases
                {
                    const bitstream<bit_hash_table<size_t,44>::buff_t>& stream = ht.get_data();
                    lpg_build::key_wrapper key_w{44, ht.description_bits(), stream};
                    std::vector<size_t> sorted_phrases;
                    sorted_phrases.reserve(ht.size());

                    for (auto const &phrase : ht) {
                        sorted_phrases.push_back(phrase);
                    }
                    std::sort(sorted_phrases.begin(), sorted_phrases.end(),
                              [&](const size_t& l, const size_t &r){
                                  return key_w.compare(l, r);
                              });
                    for (auto const &phrase_ptr : sorted_phrases) {
                        ht.insert_value_at(phrase_ptr, rank++);
                    }
                }

                // create the parse
                auto parse_task = [&](auto& phrase){
                    phrase.mask_tail();
                    auto res = ht.find(phrase.data(), phrase.n_bits());
                    next_parse.push_back(res.first.value());
                };
                lms_scan(parse_task, pos, parse, n_lms);
                std::reverse(next_parse.begin(), next_parse.end());
                std::swap(parse, next_parse);
            }
            next_parse.clear();
            ht.flush();
            lvl++;
        }
    }

public:
    typedef size_t                          size_type;

    lpg_index(std::string &input_file, std::string &tmp_folder, size_t n_threads, float hbuff_frac) {

        std::cout<<"Input file: "<<input_file<<std::endl;
        auto alphabet =  get_alphabet(input_file);

        size_t n_chars = 0;
        for(auto const &sym : alphabet) n_chars+=sym.second;

        //create a temporary folder
        std::string tmp_path = tmp_folder+"/lpg_index.XXXXXX";
        char temp[200]={0};
        tmp_path.copy(temp, tmp_path.size() + 1);
        temp[tmp_path.size()+1] = '\0';
        auto res = mkdtemp(temp);
        if(res==nullptr){
            std::cout<<"Error trying to create a temporal folder"<<std::endl;
        }
        std::cout<<"Temporal folder: "<<std::string(temp)<<std::endl;
        //

        sdsl::cache_config config(false, temp);
        std::string g_file = sdsl::cache_file_name("g_file", config);

        //maximum amount of RAM allowed to spend in parallel for the hashing step
        auto hbuff_size = std::max<size_t>(64*n_threads, size_t(std::ceil(float(n_chars)*hbuff_frac)));

        std::cout<<"Computing the LPG grammar"<<std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        lpg_build::compute_LPG(input_file, g_file, n_threads, config, hbuff_size, alphabet);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout<<"  Elap. time (secs): "<<elapsed.count()<<std::endl;

        //plain representation of the grammar
        plain_grammar_t plain_gram;
        plain_gram.load_from_file(g_file);
        std::cout <<"Resulting grammar: "<<std::endl;
        std::cout <<"  Terminals:                " << (size_t)plain_gram.sigma << std::endl;
        std::cout << "  Nonterminals:             " << plain_gram.r - plain_gram.sigma << std::endl;
        std::cout << "  Size of the comp. string: " << plain_gram.c << std::endl;
        std::cout << "  Grammar size:             " << plain_gram.g - plain_gram.sigma << std::endl;

        std::cout<<"Building the self-index"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        build_index(plain_gram, config);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout<<"  Elap. time (secs): "<<elapsed.count()<<std::endl;
    }

    lpg_index()= default;

    lpg_index(const lpg_index &other){
        grammar_tree = other.grammar_tree;
        grid = other.grid;
        symbols_map = other.symbols_map;
        m_sigma = other.m_sigma;
        parsing_rounds = other.parsing_rounds;
        rl_compressed = other.rl_compressed;
    }

    void swap(lpg_index&& other){
        std::swap(grammar_tree, other.grammar_tree);
        std::swap(grid, other.grid);
        std::swap(symbols_map, other.symbols_map);
        std::swap(m_sigma, other.m_sigma);
        std::swap(parsing_rounds, other.parsing_rounds);
        std::swap(rl_compressed, other.rl_compressed);
    }

    lpg_index(lpg_index &&other) noexcept {
        swap(std::forward<lpg_index>(other));
    }

    lpg_index& operator=(lpg_index&& other) noexcept {
        swap(std::forward<lpg_index>(other));
        return *this;
    }

    lpg_index& operator=(lpg_index const& other) noexcept = default;

    //statistics about the text: number of symbols, number of documents, etc
    void text_stats(std::string& list){
    }

    //extract text[start, end] from the index
    void extract(size_t start, size_t end){
    }

    //search for a list of patterns
    void search(std::vector<std::string>& list){
        for(auto const& pattern : list){
            auto cut = compute_pattern_cut(pattern);
        }
    }

    //search for a list of pattern (input is a file)
    void search(std::string& list){
    }

    void load(std::istream &in){
        grammar_tree.load(in);
        grid.load(in);
        symbols_map.load(in);
        sdsl::read_member(m_sigma, in);
        sdsl::read_member(parsing_rounds, in);
        sdsl::read_member(rl_compressed, in);
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += grammar_tree.serialize(out, child, "grammar_tree");
        written_bytes += grid.serialize(out, child, "grid");
        written_bytes += symbols_map.serialize(out, child, "symbols_map");
        written_bytes += sdsl::write_member(m_sigma, out, child, "sigma");
        written_bytes += sdsl::write_member(parsing_rounds, out, child, "sigma");
        written_bytes += sdsl::write_member(rl_compressed, out, child, "sigma");
        return written_bytes;
    }
};
#endif //LMS_GRAMMAR_REP_HPP
