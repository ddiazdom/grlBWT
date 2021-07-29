//
// Created by diediaz on 11-12-19.
//

#ifndef LMS_GRAMMAR_REP_HPP
#define LMS_GRAMMAR_REP_HPP

#include <pthread.h>
#include <iostream>
#include <cstdlib>
#include <mem_monitor/mem_monitor.hpp>
#include "lpg_build.hpp"
#include "grammar_tree.hpp"
#include "grid.hpp"
#include "macros.hpp"

class lpg_index {

public:
    typedef sdsl::sd_vector<> bv_y;
    typedef typename lpg_build::plain_grammar_t plain_grammar_t;
    typedef typename lpg_build::alpha_t alpha_t;

    grammar_tree_t grammar_tree;
    grid m_grid;

    bv_y Y;
    bv_y::rank_1_type rank_Y;

    sdsl::int_vector<8> symbols_map; // map a compressed terminal to its original byte symbol
    uint8_t m_sigma{}; //alphabet of terminal symbols
    uint8_t parsing_rounds{}; //number of LMS parsing rounds during the grammar construction
    bool rl_compressed{}; // is the grammar run-length compressed?

    void build_index(const std::string &i_file, plain_grammar_t &p_gram, const size_t &text_length,
                     sdsl::cache_config &config) {
        m_sigma = p_gram.sigma;
        parsing_rounds = p_gram.rules_per_level.size();

        {
            sdsl::bit_vector _y(p_gram.r, 0);
            symbols_map.resize(p_gram.sym_map.size());
            size_type ii = 0;
            for (const auto &item : p_gram.sym_map) {

                _y[item.first] = true;
                symbols_map[ii] = item.second;
                ii++;
            }
            std::sort(symbols_map.begin(), symbols_map.end());
            Y = bv_y(_y);
            rank_Y = bv_y::rank_1_type(&Y);
#ifdef DEBUG_PRINT
            utils::pretty_printer_v(symbols_map, "sym");
            utils::pretty_printer_bv(Y, "Y");
#endif
        }

        utils::lenght_rules lenghts;
        size_type S;
        utils::nav_grammar NG = build_nav_grammar(p_gram, S);
        grammar_tree.build(NG, p_gram, text_length, lenghts, S);
        std::vector<utils::sfx> grammar_sfx;
        const auto &T = grammar_tree.getT();
        compute_grammar_sfx(NG, p_gram, lenghts, grammar_sfx);
        NG.clear();
        lenghts.clear();
#ifdef DEBUG_INFO
        std::cout << "sort_suffixes[" << grammar_sfx.size() << "]\n";
#endif
        utils::sort_suffixes(i_file, grammar_sfx);
#ifdef DEBUG_PRINT
        int i = 0;
        for (const auto &sfx : grammar_sfx) {
            std::cout<<i++<<":\t"<<":O["<<sfx.off<<"]P["<<sfx.preorder<<"]"<<"R["<<sfx.rule<<"]\t";
            print_suffix_grammar(sfx.preorder,10);
            std::cout<<std::endl;
        }
#endif
#ifdef DEBUG_INFO
        std::cout << "compute_grid_points\n";
#endif
        std::vector<grid_point> points;
        compute_grid_points(grammar_sfx, points);

        grammar_sfx.clear();
//        grid = grid_t(points,p_gram.rules_per_level.size());
        m_grid = grid(points);
#ifdef DEBUG_INFO
        std::cout << "build grid\n";
        breakdown_space();
        std::cout << "end build index\n";
#endif
    };

    void breakdown_space() const {
        std::cout << "breakdown-space-index\n";
        std::cout << "Text-length," << grammar_tree.get_text_len() << std::endl;
        std::cout << "Number-rules," << grammar_tree.get_size_rules() << std::endl;
        std::cout << "Grammar-size," << grammar_tree.get_grammar_size() << std::endl;
        std::cout << "Grammar-Tree," << sdsl::size_in_bytes(grammar_tree) << std::endl;
        grammar_tree.breakdown_space();
        std::cout << "Grid," << sdsl::size_in_bytes(m_grid) << std::endl;
        m_grid.breakdown_space();
        std::cout << "symbols_map," << sdsl::size_in_bytes(symbols_map);
        std::cout << "m_sigma," << sizeof(m_sigma);
        std::cout << "parsing_rounds," << sizeof(parsing_rounds);
        std::cout << "rl_compressed," << sizeof(rl_compressed);
    }

    struct parsing_data {
        size_t idx{};
        size_t n_lms;
        size_t round=0;
        bool tail=true;
        std::vector<uint16_t> lms_pos;
        std::vector<uint16_t> new_lms_pos;
        std::vector<size_t> cuts;
        std::vector<size_t> parse;

        explicit parsing_data(const std::string& pattern){
            lms_pos.resize(pattern.size());
            std::iota(lms_pos.begin(), lms_pos.end(), 0);
            new_lms_pos.resize(lms_pos.size());
            parse.reserve(pattern.size());
            for (auto const &sym : pattern) parse.push_back(sym);
            n_lms = parse.size();
        }

        void extract_cuts(std::pair<uint8_t, size_t>& lms_data){

            std::cout<<"extract_cuts:before:new_lms_pos"<<std::endl;
            for (const auto &item : new_lms_pos) {
                std::cout<<item<<" ";
            }

            std::cout<<std::endl;
            // the first text position is either S-type or the first symbol of the original text
            if(lms_data.first || round==0) cuts.push_back(lms_pos[0]);

            // check if the prefix in the text is a run.
            // If so, create a cut in the after the rightmost symbol of that run
            if(round==0){
                size_t pos = 0;
                while(pos<parse.size()-1 && parse[pos]==parse[pos+1]) pos++;
                if(pos>0 && pos<parse.size()-1) cuts.push_back(lms_pos[pos]);
            }
            // cuts for the leftmost and rightmost
            // incomplete phrases

            if(n_lms>0){
                cuts.push_back(new_lms_pos[n_lms-1]);
                if(n_lms>1){
                    //TODO: we can discard this cut if we determine that the last
                    // phrase in the text is LMS

                    if( n_lms < 4 ){
                        cuts.push_back(new_lms_pos[1]);
                    }

                    cuts.push_back(new_lms_pos[0]);
                }
            }

            // when the suffix of the text is a run, and the
            // symbol before that run is L-type, then create
            // a cut after the leftmost symbol of that run
            if(lms_data.second<parse.size()-1 &&
               lms_data.second>0 &&
               parse[lms_data.second-1]>parse[lms_data.second] &&
               (cuts.empty() || cuts.back()!=lms_pos[lms_data.second])){
                cuts.push_back(lms_pos[lms_data.second]);
            }

            std::cout<<"extract_cuts:after:cuts"<<std::endl;
            for (const auto &item : cuts) {
                std::cout<<item<<" ";
            }
            std::cout<<std::endl;
        }

        void update_lms_pos(){

            std::cout<<"update_lms_pos:lms_pos"<<std::endl;
            for (const auto &item : lms_pos) {
                std::cout<<item<<" ";
            }
            std::cout<<std::endl;

            std::cout<<"update_lms_pos:new_lms_pos"<<std::endl;
            for (const auto &item : new_lms_pos) {
                std::cout<<item<<" ";
            }
            std::cout<<std::endl;

            new_lms_pos.resize(n_lms-1);
            std::reverse(new_lms_pos.begin(), new_lms_pos.end());
            std::swap(new_lms_pos, lms_pos);
        }
    };

    static alpha_t get_alphabet(std::string &i_file) {

        std::cout << "Reading input file" << std::endl;

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

        std::cout << "  Number of characters: " << if_stream.size() << std::endl;
        std::cout << "  Alphabet:             " << alphabet.size() << std::endl;
        std::cout << "  Smallest symbol:      " << (int) alphabet[0].first << std::endl;
        std::cout << "  Greatest symbol:      " << (int) alphabet.back().first << std::endl;

        if (if_stream.read(if_stream.size() - 1) != alphabet[0].first) {
            std::cout << "Error: sep. symbol " << alphabet[0].first << " differs from last symbol in file "
                      << if_stream.read(if_stream.size() - 1) << std::endl;
            exit(1);
        }
        return alphabet;
    }

    //returns a boolean that indicates if the first symbol could belong to another phrase
    // 0 : the first symbol belong to a phrase in parse
    // 1 : the first symbol belong (or might) to a phrase to the left of the parse
    template<typename proc>
    static std::pair<uint8_t, size_t> lms_scan(proc &task, std::vector<size_t> &parse) {

        /*for(size_t i=0;i<parse.size();i++){
            std::cout<<parse[i]<<" ";
        }
        std::cout<<""<<std::endl;*/

        int_array<size_t> lms_phrase(2, 32);
        size_t curr_sym, prev_sym,pos, rm_run;
        bool s_type, prev_s_type;

        lms_phrase.push_back(parse.back());
        pos = parse.size() - 2;
        while (parse[pos] == parse[pos + 1]) {
            if(pos==0) {
                return {true, 0};
            }
            if(pos==1){
                return {parse[pos-1]<=parse[pos], parse[pos-1]!=parse[pos]};
            }

            lms_phrase.push_back(parse[pos--]);
        }
        rm_run = pos+1;
        prev_s_type = parse[pos] < parse[pos+1];
        prev_sym = parse[pos];
        lms_phrase.push_back(prev_sym);

        for (size_t i = pos; i-- > 0;) {
//            std::cout<<parse[i]<<std::endl;
            curr_sym = parse[i];
            if (curr_sym < prev_sym) {//S_TYPE type
                s_type = S_TYPE;
            } else if (curr_sym == prev_sym) {
                s_type = prev_s_type;
            } else {//L_TYPE type
                s_type = L_TYPE;
                if (prev_s_type == S_TYPE) {//LMS-type
//                    std::cout<<prev_sym<<std::endl;
                    lms_phrase.pop_back();
                    task(lms_phrase);
                    lms_phrase.clear();
                    lms_phrase.push_back(prev_sym);
                }
            }
            lms_phrase.push_back(curr_sym);
            prev_sym = curr_sym;
            prev_s_type = s_type;
        }
        task(lms_phrase);

        return {prev_s_type == S_TYPE, rm_run};
    }


public:

    std::pair<std::vector<size_t>, uint8_t> compute_pattern_cuts(const std::string &pattern) const {

        parsing_data p_data(pattern);
        //hash table to hash the LMS phrases
        bit_hash_table<size_t, 44> ht;
        //lambda function to hash the LMS phrases
        auto hash_task = [&](auto &phrase) {
            if (p_data.tail) {
                p_data.tail = false;
            } else {
                phrase.mask_tail();
                if (p_data.idx > phrase.size()) {
                    ht.insert(phrase.data(), phrase.n_bits(), 0);
                }
            }

            if (p_data.idx > phrase.size()) {
                p_data.idx -= phrase.size();
                p_data.new_lms_pos[p_data.n_lms++] = p_data.lms_pos[p_data.idx];
            }
        };
        //lambda function to create the LMS parse
        auto parse_task = [&](auto &phrase) {
            if (p_data.tail) {
                p_data.tail = false;
            } else if (p_data.n_lms > 1) {
                phrase.mask_tail();
                auto res = ht.find(phrase.data(), phrase.n_bits());
                assert(res.second);
                p_data.parse[p_data.idx--] = res.first.value();
                p_data.n_lms--;
            }
        };
        //uint8_t p_round=0;
        //hash the LMS phrases in the text
        while (true) {
            //the pattern is longer than the text
            if(p_data.round > parsing_rounds){
                return {{},parsing_rounds};
            }
            assert(p_data.parse.size() > 1);
            p_data.idx = p_data.parse.size() - 1;
            p_data.n_lms = 0;
            p_data.tail = true;

            auto lms_data = lms_scan(hash_task, p_data.parse);


            if (p_data.n_lms < 4) {
                //report the cuts
                p_data.extract_cuts(lms_data);
                return {p_data.cuts, p_data.round};
            } else {
                //assign ranks to the LMS phrases
                {
                    size_t rank = 0;
                    const bitstream<bit_hash_table<size_t, 44>::buff_t> &stream = ht.get_data();
                    lpg_build::key_wrapper key_w{32, ht.description_bits(), stream};
                    std::vector<size_t> sorted_phrases;
                    sorted_phrases.reserve(ht.size());

                    for (auto const &phrase : ht) {
                        sorted_phrases.push_back(phrase);
                    }
                    std::sort(sorted_phrases.begin(), sorted_phrases.end(),
                              [&](const size_t &l, const size_t &r) {
                                  return key_w.compare(l, r);
                              });
                    for (auto const &phrase_ptr : sorted_phrases) {
                        ht.insert_value_at(phrase_ptr, rank++);
                    }
                }

                p_data.extract_cuts(lms_data);
                p_data.update_lms_pos();

                //create the new parse
                p_data.tail = true;
                p_data.idx = p_data.parse.size() - 1;

                lms_scan(parse_task, p_data.parse);

                p_data.parse.erase(p_data.parse.begin(), p_data.parse.begin() + p_data.idx + 1);
            }
            ht.flush();
            p_data.round++;
        }
    }

public:
    typedef size_t size_type;

    bool just_one_zero(std::string &input_file){
        std::string text;
        utils::readFile(input_file,text);
        int zero_count = 0;
        for (const auto &item : text) {
            if(item == 0 ) zero_count++;
            if(zero_count > 1)
                return false;
        }
        return true;
    }

    lpg_index(std::string &input_file, std::string &tmp_folder, size_t n_threads, float hbuff_frac) {

        mem_monitor mem(input_file + "-mem.csv");
        std::cout<<"measuring peak memory\n";
        mem.event("LPG-BUILD-GRAMMAR");
        std::cout << "Input file: " << input_file << std::endl;
        if (!just_one_zero(input_file)) {
            std::cout << "More than one zero error" << std::endl;
            return;
        }
        auto alphabet = get_alphabet(input_file);

        /*TODO este if puede reemplazar a just_one_zero
        if(alphabet[0].first==0 && alphabet[0].second>1){
            exit(1);
        }*/

        size_t n_chars = 0;
        for (auto const &sym : alphabet) n_chars += sym.second;

        //create a temporary folder
        std::string tmp_path = tmp_folder + "/lpg_index.XXXXXX";
        char temp[200] = {0};
        tmp_path.copy(temp, tmp_path.size() + 1);
        temp[tmp_path.size() + 1] = '\0';
        auto res = mkdtemp(temp);
        if (res == nullptr) {
            std::cout << "Error trying to create a temporal folder" << std::endl;
        }
        std::cout << "Temporal folder: " << std::string(temp) << std::endl;
        //

        sdsl::cache_config config(false, temp);
        std::string g_file = sdsl::cache_file_name("g_file", config);

        //maximum amount of RAM allowed to spend in parallel for the hashing step
        auto hbuff_size = std::max<size_t>(64 * n_threads, size_t(std::ceil(float(n_chars) * hbuff_frac)));


        std::cout << "Computing the grammar for the self-index" << std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        lpg_build::compute_LPG(input_file, g_file, n_threads, config, hbuff_size, alphabet);
        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed_grammar = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "  Elap. time (microsec): " << elapsed_grammar.count() << std::endl;

        //plain representation of the grammar
        plain_grammar_t plain_gram;
        plain_gram.load_from_file(g_file);

#ifdef DEBUG_PRINT
        plain_gram.print_grammar();
#endif

        mem.event("LPG-BUILD-INDEX");
        std::cout << "Building the self-index" << std::endl;
        start = std::chrono::high_resolution_clock::now();
        build_index(input_file, plain_gram, n_chars, config);
        end = std::chrono::high_resolution_clock::now();
        auto elapsed_index = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
        std::cout << "  Elap. time (microsec): " << elapsed_index.count() << std::endl;
        double text_size = grammar_tree.get_text_len();
        double index_size = sdsl::size_in_bytes(*this);
        std::cout << "  Index size(bytes) " << sdsl::size_in_bytes(*this) << std::endl;
        std::cout << "  Text size(bytes) " << grammar_tree.get_text_len() << std::endl;
        std::cout << "  (bits/sym) " << index_size * 8 /text_size << std::endl;

        std::ofstream csv_file;
        csv_file.open (input_file + "build-time.csv");
        csv_file << "name,time-grammar,time-index\n";
        csv_file << "LPG-INDEX-RRR,"+std::to_string(elapsed_grammar.count())+","+std::to_string(elapsed_index.count())+"\n";


    }

    lpg_index() = default;

    lpg_index(const lpg_index &other) {
        grammar_tree = other.grammar_tree;
        m_grid = other.m_grid;
        symbols_map = other.symbols_map;
        m_sigma = other.m_sigma;
        parsing_rounds = other.parsing_rounds;
        rl_compressed = other.rl_compressed;
    }

    void swap(lpg_index &&other) {
        std::swap(grammar_tree, other.grammar_tree);
        std::swap(m_grid, other.m_grid);
        std::swap(symbols_map, other.symbols_map);
        std::swap(m_sigma, other.m_sigma);
        std::swap(parsing_rounds, other.parsing_rounds);
        std::swap(rl_compressed, other.rl_compressed);
    }

    lpg_index(lpg_index &&other) noexcept {
        swap(std::forward<lpg_index>(other));
    }

    lpg_index &operator=(lpg_index &&other) noexcept {
        swap(std::forward<lpg_index>(other));
        return *this;
    }

    lpg_index &operator=(lpg_index const &other) noexcept = default;

    //statistics about the text: number of symbols, number of documents, etc
    void text_stats(std::string &list) {}

    void locate(const std::string &pattern, std::set<uint64_t> &pos) const;
    void locate_all_cuts(const std::string &pattern, std::set<uint64_t> &pos) const;
    void locate_split_time(const std::string &pattern, std::set<uint64_t> &pos, uint64_t&, uint64_t&) const;
    //extract text[start, end] from the index
    void extract(size_t start, size_t end) {}

    static void bt_search(const std::string &str,const std::string &sub, std::set<size_t> &positions){
        size_t pos = str.find(sub, 0);
        while(pos != std::string::npos)
        {
            positions.insert(pos);
            pos = str.find(sub,pos+1);
        }
    }
    //search for a list of patterns
    void search(std::vector<std::string> &list
#ifdef CHECK_OCC
            ,const std::string& file
#endif
    ) {
#ifdef CHECK_OCC
        std::string data;
        utils::readFile(file,data);
        int ii = 0;
#endif
        std::cout << "Locate pattern list["<<list.size()<<"]" << std::endl;
        size_t total_occ = 0,total_time = 0 , total_occ_bt = 0;
        for (auto const &pattern : list) {
#ifdef DEBUG_PRINT
            std::cout << pattern << ":";
#endif
//            std::cout<<++ii<<"--"<<pattern<<std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            std::set<size_type> occ;
            locate(pattern, occ);
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            total_occ += occ.size();
            total_time+=elapsed;


#ifdef CHECK_OCC
            //
            std::set<size_type> positions;
            bt_search(data,pattern,positions);
            total_occ_bt += positions.size();
            if(positions != occ){
                std::cout<<"Locate error\n";
                std::cout<<"pattern:"<<pattern<<std::endl;
                std::cout<<"occ:"<<occ.size()<<std::endl;
                std::cout<<"bt-occ:"<<positions.size()<<std::endl;
                if(positions.size() > occ.size())
                {
                    std::vector<size_type> X;
                    X.resize(positions.size(),0);
                    auto it = std::set_difference(positions.begin(),positions.end(),occ.begin(),occ.end(),X.begin());
                    X.resize(it - X.begin());
                    std::cout<<"missing positions["<<X.size()<<"]\n";
                }
////                return;
            }

#endif

        }

        auto text_size = (double)grammar_tree.get_text_len();
        auto index_size = (double)sdsl::size_in_bytes(*this);
        std::cout << "Elap. time (microsec): " << total_time << std::endl;
        std::cout << "Total occ: " << total_occ << std::endl;
        std::cout << "Real Total occ: " << total_occ_bt << std::endl;
        double time_per_occ = (double)total_time/(double)total_occ;
        std::cout << "Time/occ (microsec): " << time_per_occ << std::endl;
        std::cout << "Index size " << sdsl::size_in_bytes(*this) << std::endl;
        std::cout << "Text size " << grammar_tree.get_text_len() << std::endl;
        std::cout << "Bps " << index_size * 8 /text_size << std::endl;


    }

    void search_all_cuts(std::vector<std::string> &list
#ifdef CHECK_OCC
            ,const std::string& file
#endif
    ) {
#ifdef CHECK_OCC
        std::string data;
        utils::readFile(file,data);
        int ii = 0;
#endif
        std::cout << "Locate pattern list["<<list.size()<<"]" << std::endl;
        size_t total_occ = 0,total_time = 0 , total_occ_bt = 0;
        for (auto const &pattern : list) {
#ifdef DEBUG_PRINT
            std::cout << pattern << ":";
#endif
//            std::cout<<++ii<<"--"<<pattern<<std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            std::set<size_type> occ;
            locate_all_cuts(pattern, occ);
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            total_occ += occ.size();
            total_time+=elapsed;


#ifdef CHECK_OCC
            //
            std::set<size_type> positions;
            bt_search(data,pattern,positions);
            total_occ_bt += positions.size();
            if(positions != occ){
                std::cout<<"Locate error\n";
                std::cout<<"pattern:"<<pattern<<std::endl;
                std::cout<<"occ:"<<occ.size()<<std::endl;
                std::cout<<"bt-occ:"<<positions.size()<<std::endl;
                if(positions.size() > occ.size())
                {
                    std::vector<size_type> X;
                    X.resize(positions.size(),0);
                    auto it = std::set_difference(positions.begin(),positions.end(),occ.begin(),occ.end(),X.begin());
                    X.resize(it - X.begin());
                    std::cout<<"missing positions["<<X.size()<<"]\n";
                }
////                return;
            }

#endif

        }

        auto text_size = (double)grammar_tree.get_text_len();
        auto index_size = (double)sdsl::size_in_bytes(*this);
        std::cout << "Elap. time (microsec): " << total_time << std::endl;
        std::cout << "Total occ: " << total_occ << std::endl;
        std::cout << "Real Total occ: " << total_occ_bt << std::endl;
        double time_per_occ = (double)total_time/(double)total_occ;
        std::cout << "Time/occ (microsec): " << time_per_occ << std::endl;
        std::cout << "Index size " << sdsl::size_in_bytes(*this) << std::endl;
        std::cout << "Text size " << grammar_tree.get_text_len() << std::endl;
        std::cout << "Bps " << index_size * 8 /text_size << std::endl;


    }

    void search_split_time(std::vector<std::string> &list
#ifdef CHECK_OCC
    ,const std::string& file
#endif
    ) {
#ifdef CHECK_OCC
        std::string data;
        utils::readFile(file,data);
        int ii = 0;
#endif
        std::cout << "Locate pattern list["<<list.size()<<"]" << std::endl;
        size_t total_occ = 0,total_time = 0, total_time_p = 0,total_time_s = 0 , total_occ_bt = 0;
        for (auto const &pattern : list) {
#ifdef DEBUG_PRINT
            std::cout << pattern << ":";
#endif
//            std::cout<<++ii<<"--"<<pattern<<std::endl;
            auto start = std::chrono::high_resolution_clock::now();
            std::set<size_type> occ;
            locate_split_time(pattern, occ,total_time_p,total_time_s);
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
            total_occ += occ.size();
            total_time+=elapsed;


#ifdef CHECK_OCC
//
//            bt_search(data,pattern,positions);
//            total_occ_bt += positions.size();
//            if(positions != occ){
//                std::cout<<"Locate error\n";
//                std::cout<<"pattern:"<<pattern<<std::endl;
//                std::cout<<"occ:"<<occ.size()<<std::endl;
//                std::cout<<"bt-occ:"<<positions.size()<<std::endl;
//                if(positions.size() > occ.size())
//                {
//                    std::vector<size_type> X;
//                    X.resize(positions.size(),0);
//                    auto it = std::set_difference(positions.begin(),positions.end(),occ.begin(),occ.end(),X.begin());
//                    X.resize(it - X.begin());
//                    std::cout<<"missing positions["<<X.size()<<"]\n";
//                }
//////                return;
//            }

#endif

        }

        auto text_size = (double)grammar_tree.get_text_len();
        auto index_size = (double)sdsl::size_in_bytes(*this);
        std::cout << "Elap. time (microsec): " << total_time << std::endl;
        std::cout << "Elap. time primary occ (microsec): " << total_time_p << std::endl;
        std::cout << "Elap. time secondary occ (microsec): " << total_time_s << std::endl;
        std::cout << "Total occ: " << total_occ << std::endl;
        std::cout << "Real Total occ: " << total_occ_bt << std::endl;
        double time_per_occ = (double)total_time/(double)total_occ;
        std::cout << "Time/occ (microsec): " << time_per_occ << std::endl;
        std::cout << "Index size " << sdsl::size_in_bytes(*this) << std::endl;
        std::cout << "Text size " << grammar_tree.get_text_len() << std::endl;
        std::cout << "Bps " << index_size * 8 /text_size << std::endl;


    }

    //search for a list of pattern (input is a file each line is a pattern)
    void search(std::string &input_file) const {
        std::fstream in(input_file, std::ios::in);
        std::string pattern;
        long total_time = 0;
        while (in.good() && std::getline(in, pattern)) {
            auto start = std::chrono::high_resolution_clock::now();
                std::set<size_type> occ;
                locate(pattern, occ);
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = (end - start).count();
            total_time += elapsed;
        }
        std::cout << "  Elap. time (secs): " << total_time << std::endl;
    }

    void load(std::istream &in) {
        grammar_tree.load(in);
        m_grid.load(in);
        symbols_map.load(in);
        sdsl::read_member(m_sigma, in);
        sdsl::read_member(parsing_rounds, in);
        sdsl::read_member(rl_compressed, in);
        Y.load(in);
        rank_Y.load(in);


        rank_Y = bv_y::rank_1_type(&Y);
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += grammar_tree.serialize(out, child, "grammar_tree");
        written_bytes += m_grid.serialize(out, child, "m_grid");
        written_bytes += symbols_map.serialize(out, child, "symbols_map");
        written_bytes += sdsl::write_member(m_sigma, out, child, "sigma");
        written_bytes += sdsl::write_member(parsing_rounds, out, child, "sigma");
        written_bytes += sdsl::write_member(rl_compressed, out, child, "sigma");
        written_bytes += Y.serialize(out, child, "Y");
        written_bytes += rank_Y.serialize(out, child, "rank_Y");
        return written_bytes;
    }
    template<typename F>
    bool dfs_mirror_leaf_base_case(const uint64_t &preorder_node, const uint64_t &node, const F &f) const {
        size_type _x = grammar_tree.get_rule_from_preorder_node(preorder_node);//get the rule of the leaf
        if (is_terminal(_x)) {//check if the rule is a terminal symbol
            //terminal case
            bool keep = f(preorder_node, node, _x);//call process function
            return keep; // keep says if we must to stop the expansion of the rule
        }
        //second occ.
        auto fpre_node = grammar_tree.first_occ_from_rule(_x);//get preorder of first mention node
        const auto &T = grammar_tree.getT();
        auto fnode = T[fpre_node];//preorder select

        return dfs_mirror_leaf(fpre_node, fnode, f); // call recursive in the new node...
    }

    template<typename F>
    bool dfs_leaf_base_case(const uint64_t &preorder_node, const uint64_t &node, const F &f) const {
        size_type _x = grammar_tree.get_rule_from_preorder_node(preorder_node);//get the rule of the leaf
        if (is_terminal(_x)) {//check if the rule is a terminal symbol
            //terminal case
            bool keep = f(preorder_node, node, _x);//call process function
            return keep; // keep says if we must to stop the expansion of the rule
        }
        //second occ.
        auto fpre_node = grammar_tree.first_occ_from_rule(_x);//get preorder of first mention node
        const auto &T = grammar_tree.getT();
        auto fnode = T[fpre_node];//preorder select

        return dfs_leaf(fpre_node, fnode, f); // call recursive in the new node...
    }

    template<typename F>
    bool dfs_mirror_leaf(const uint64_t &preorder_node, const uint64_t &node, const F &f) const {
        if (grammar_tree.isLeaf(preorder_node)) { //leaf case
            return dfs_mirror_leaf_base_case(preorder_node, node, f);
        } else {
            auto len = grammar_tree.is_run(preorder_node);
            const auto &T = grammar_tree.getT();
            if (len) {
//                uint64_t X = grammar_tree.get_rule_from_preorder_node(grammar_tree.getT().pre_order(node));
//                std::cout<<"("<<node<<")NODO-LEN["<<len<<"] X:"<<X<<std::endl;
                for (uint32_t i = 0; i < len; ++i) {
                    auto chnode = T.child(node, 1);
                    bool keep = dfs_mirror_leaf(T.pre_order(chnode), chnode, f);
                    if (!keep) return keep;
                }
                return true;
            } else {
                uint32_t n = T.children(node); //number of children
//                uint64_t X = grammar_tree.get_rule_from_preorder_node(grammar_tree.getT().pre_order(node));
//                std::cout<<"("<<node<<")NODO-BLOCK["<<n<<"] X:"<<X<<std::endl;
                for (uint32_t i = n; i > 0; --i) { // from right to left
                    auto chnode = T.child(node, i); //get i-child node
                    bool keep = dfs_mirror_leaf(T.pre_order(chnode), chnode,
                                                f);// compute preorder and visit recursively
                    if (!keep) return keep; // keep says if we must to stop the expansion of the rule
                }
                return true;
            }

        }
    }

    template<typename F>
    bool dfs_leaf(const uint64_t &preorder_node, const uint64_t &node, const F &f) const {
        if (grammar_tree.isLeaf(preorder_node)) { //leaf case
            return dfs_leaf_base_case(preorder_node, node, f);
        } else {
            auto len = grammar_tree.is_run(preorder_node);
            const auto &T = grammar_tree.getT();
            if (len) {
                auto chnode = T.child(node, 1);
                for (uint32_t i = 0; i < len; ++i) {
                    bool keep = dfs_leaf(T.pre_order(chnode), chnode, f);
                    if (!keep) return keep;
                }
            } else {
                uint32_t n = T.children(node); //number of children
                for (uint32_t i = 0; i < n; ++i) { // from left to right
                    auto chnode = T.child(node, i + 1); //get i-child node
                    bool keep = dfs_leaf(T.pre_order(chnode), chnode, f);// compute preorder and visit recursively
                    if (!keep) return keep; // keep says if we must to stop the expansion of the rule
                }
            }
            return true;
        }
    }

    template<typename F>
    void process_prefix_rule(const size_type &preorder_node, const F & f) const{
        const auto &m_tree = grammar_tree.getT();
        auto node = m_tree[preorder_node];
        dfs_mirror_leaf(preorder_node, node,f);

    }

    void down_up_print_mirror(const uint64_t &preorder_node,const std::string& prefix) const {

        auto node = grammar_tree.getT().operator[](preorder_node);
        auto leaf = grammar_tree.getT().lastleaf(node);

        std::function<void(const uint64_t &node,uint l)> expand_node;
        expand_node = [this,prefix,expand_node](const uint64_t &node,uint l){
                uint64_t X = grammar_tree.get_rule_from_preorder_node(grammar_tree.getT().pre_order(node));
                    std::cout<<prefix<<" "<<X<<std::endl;
                    print_prefix_rule(X,10);
                    std::cout<<std::endl;
                if(l > 0)
                    expand_node(grammar_tree.getT().parent(node),l-1);
        };

        expand_node(leaf,3);
    }


    void down_up_print(const uint64_t &preorder_node,const std::string& prefix) const {

        auto node = grammar_tree.getT().operator[](preorder_node);
        auto rleaf = grammar_tree.getT().leafrank(node);
        auto leaf = grammar_tree.getT().leafselect(rleaf);

        std::function<void(const uint64_t &node,uint l)> expand_node;
        expand_node = [this,prefix,expand_node](const uint64_t &node,uint l){
            uint64_t pp = grammar_tree.getT().pre_order(node);
            uint64_t X = grammar_tree.get_rule_from_preorder_node(pp);
            std::cout<<prefix<<" "<<X<<std::endl;
            print_suffix_grammar(pp,10);
            std::cout<<std::endl;
            if(l > 0)
                expand_node(grammar_tree.getT().parent(node),l-1);
        };
        expand_node(leaf,3);
    }


    void print_prefix_rule(const size_type &preorder_node,const  size_type& l) const {
        size_type cont = l;
        auto cmp = [this,&cont](const uint64_t &prenode, const uint64_t &node,const uint64_t &X) {
            auto c1 = get_symbol(X);//symbols_map[X];
            std::cout<<(char)c1;
            if(cont == 0)
                return false;
            cont--;
            return true;
        };
        process_prefix_rule(preorder_node,cmp);
        std::cout<<"\n";
    }


    int cmp_prefix_rule(const size_type &preorder_node, const std::string &str, const uint32_t &i) const {
        long ii = i;
        int r = 0;
        bool match = false;
        const auto &m_tree = grammar_tree.getT();
        auto node = m_tree[preorder_node];
        auto cmp = [&str, &r, &match, &ii, this](const size_type &preorder_node, const size_type &node,
                                                 const size_type &X) {
            // X is a terminal rule
            auto c1 = get_symbol(X);//symbols_map[X];
            auto c2 = (uint8_t) str[ii];

            if (c1 > c2) {
                r = 1;
                return false;
            }
            if (c1 < c2) {
                r = -1;
                return false;
            }
            --ii;
            if (ii < 0) {
                match = true;
                return false;
            }
            return true;
        };
        dfs_mirror_leaf(preorder_node, node, cmp);
        // case : rev(grammar rule)  is longer than the rev(string prefix) and rev(string prefix)  is prefix of rev(grammar rule)
        if (r == 0 && match) return 0;
        // case : rev(string prefix) is longer than the rev(grammar rule) and rev(grammar rule) is prefix of rev(string prefix)
        if (r == 0 && !match) return -1;
        // in other case they have a at least a symbol diferent
        return r;

    }

    void print_suffix_grammar(const size_type &preorder_node, const  size_type& l) const {
        size_type cont = l;
        auto cmp = [this,&cont](const uint64_t &prenode, const uint64_t &node,const uint64_t &X) {
            auto c1 = get_symbol(X);//symbols_map[X];
            std::cout<<(char)c1;
            if(cont == 0)
                return false;
            cont--;
            return true;
        };
        process_suffix_grammar(preorder_node,cmp);
    }
    template<typename F>
    int process_suffix_grammar(const size_type &preorder_node, const F & f) const{

        const auto &T = grammar_tree.getT();
        auto node = T[preorder_node];//preorder select

        auto fopen = T.bps.find_open(node - 1); // this is to compute parent and childrank
        auto parent = T.pred0(fopen) + 1;//  compute parent
        auto pre_parent = T.pre_order(parent);//  compute parent
        size_type len = grammar_tree.is_run(pre_parent);
        if(len){
            for (size_type j = 1; j < len; ++j){
                auto cnode = T.child(parent, j);
                if(!dfs_leaf(pre_parent + 1, cnode, f)) return 0;
            }
        }else{
            uint32_t ch = T.children(parent);//  compute children
            uint32_t chr = T.succ0(fopen) - fopen;//  compute childrank
            for (uint32_t j = chr; j <= ch; ++j) {
                auto cnode = T.child(parent, j);
                auto pcnode = T.pre_order(cnode);
                if(!dfs_leaf(pcnode, cnode, f)) return 0;
            }
        }
        return 1;
    }
    int cmp_suffix_grammar(const size_type &preorder_node, const std::string &str, const uint32_t &i) const {
        uint32_t ii = i, sfx_len = str.size();
        int r = 0;
        bool match = false;
        const auto &T = grammar_tree.getT();
        auto node = T[preorder_node];//preorder select
        auto cmp = [&str, &r, &match, &ii, &sfx_len, this](const uint64_t &prenode, const uint64_t &node,
                                                           const uint64_t &X) {
            auto c1 = get_symbol(X);//symbols_map[X];
            auto c2 = (uint8_t) str[ii];
            if (c1 > c2) {
                r = 1;
                return false;
            }
            if (c1 < c2) {
                r = -1;
                return false;
            }
            ii++;
            if (ii == sfx_len) {
                match = true;
                return false;
            }
            return true;
        };

        auto fopen = T.bps.find_open(node - 1); // this is to compute parent and childrank
        auto parent = T.pred0(fopen) + 1;//  compute parent
        auto pre_parent = T.pre_order(parent);//  compute parent


        size_type len = grammar_tree.is_run(pre_parent);
        if(len){
            for (size_type j = 1; j < len; ++j){
                auto cnode = T.child(parent, j);
                dfs_leaf(pre_parent + 1, cnode, cmp);
                // if they have a at least a symbol diferent return r
                if (r != 0) return r;
                // case : grammar sfx  is longer than the string sfx and string sfx  is prefix of grammar sfx
                if (r == 0 && match) return 0;
                // if r == 0 and !match we process next sibling
            }
            // case : string sfx is longer than the grammar sfx and grammar sfx is prefix of string sfx
            if (r == 0 && !match) return -1;
            // in other case they have a at least a symbol diferent
            return r;

        }else{
            uint32_t ch = T.children(parent);//  compute children
            uint32_t chr = T.succ0(fopen) - fopen;//  compute childrank
            for (uint32_t j = chr; j <= ch; ++j) {
                auto cnode = T.child(parent, j);
                auto pcnode = T.pre_order(cnode);
                dfs_leaf(pcnode, cnode, cmp);
                // if they have a at least a symbol diferent return r
                if (r != 0) return r;
                // case : grammar sfx  is longer than the string sfx and string sfx  is prefix of grammar sfx
                if (r == 0 && match) return 0;
                // if r == 0 and !match we process next sibling
            }
            // case : string sfx is longer than the grammar sfx and grammar sfx is prefix of string sfx
            if (r == 0 && !match) return -1;
            // in other case they have a at least a symbol diferent
            return r;
        }
    }

    inline bool is_terminal(const size_type &X) const {
        return Y[X];
    }

    inline uint8_t get_symbol(const size_type &X) const {
        return symbols_map[rank_Y(X)];
    }

    /**
     * Find the m_grid range to search
     * @return grid_query
     * */

    bool search_grid_range(const char *pattern, const uint32_t &len, const uint32_t &p, const uint32_t &level,
                           grid_query &q) const {
        // search rules range....
        auto cmp_rev_prefix_rule = [&p, &pattern, this](const size_type &rule_id) {
            // compute node definiton preorder of the rule
            uint64_t prenode = grammar_tree.first_occ_from_rule(rule_id);
            return cmp_prefix_rule(prenode, pattern, p - 1);
        };
        uint64_t row_1 = 0, row_2 = grammar_tree.get_size_rules() - 2;
        //search lower
        if (!utils::lower_bound(row_1, row_2, cmp_rev_prefix_rule)) return false;
        q.row1 = row_1;
        //search upper
        row_2 = grammar_tree.get_size_rules() - 2;
        if (!utils::upper_bound(row_1, row_2, cmp_rev_prefix_rule)) return false;
        q.row2 = row_2;
//
        //search suffixes

        auto cmp_suffix_grammar_rule = [&p, &pattern, this](const size_type &suffix_id) {
            //val is just to use the std lower bound method
            // compute node definiton preorder of the rule
            uint64_t prenode = m_grid.first_label_col(suffix_id);
            return cmp_suffix_grammar(prenode, pattern, p);
        };

        uint64_t col_1 = 1, col_2 = m_grid.size_cols();
        //search lower
        if (!utils::lower_bound(col_1, col_2, cmp_suffix_grammar_rule)) return false;
        q.col1 = col_1;
        //search upper
        col_2 = m_grid.size_cols();
        if (!utils::upper_bound(col_1, col_2, cmp_suffix_grammar_rule)) return false;
        q.col2 = col_2;
        //search suffixes

        return true;
    }


    void grid_search(const grid_query &range, const uint64_t &pattern_off, const uint32_t &m, const uint32_t &level,
                     std::vector<utils::primaryOcc> &occ) const {
        std::vector<uint64_t> sfx;
//        m_grid.search(range,level,sfx);
        m_grid.search_2d(range, sfx);
        occ.reserve(sfx.size());
        const auto &T = grammar_tree.getT();


        for (size_type i = 0; i < sfx.size(); ++i) {
            size_type preorder_node = m_grid.first_label_col(sfx[i]);
//            std::cout<<"preorder_node:"<<preorder_node<<std::endl;
            size_type node = T[preorder_node];
            size_type leaf = 0;
            size_type off = grammar_tree.offset_node(node, leaf);

            size_type parent = T.parent(node);
            size_type parent_off = grammar_tree.offset_node(parent);
            size_type parent_preorder = T.pre_order(parent);
            size_type run_len = grammar_tree.is_run(parent_preorder);


            if (run_len) {
                // add run length primary occ
                size_type first_child_size = off - grammar_tree.offset_node(parent);
                size_type num_leaves = grammar_tree.num_leaves();
                size_type end_node = leaf == num_leaves ?
                        grammar_tree.get_text_len() : grammar_tree.selectL(leaf + 1) - 1;

                uint32_t pattern_tail = m - pattern_off;
                while (end_node >= off + pattern_tail - 1 ) {
                    occ.emplace_back(parent, parent_preorder, parent_off, (off - parent_off) - pattern_off, true);
                    off += first_child_size ;
                }
            } else {
                size_type prefix_size = off - grammar_tree.offset_node(parent);
                occ.emplace_back(parent, parent_preorder, parent_off, prefix_size - pattern_off, true);
            }

        }
    }


    void find_secondary_occ(const utils::primaryOcc &p_occ, std::set<size_type> &occ) const {

        //queue for node processing
        std::deque<utils::primaryOcc> Q;
        const auto &T = grammar_tree.getT();
        //auxiliar functions
        auto insert_second_mentions = [&Q, &T, this](const utils::primaryOcc &occ) {
            grammar_tree.visit_secondary_occ(occ.preorder, [&T, &Q, &occ, this](const size_type &preorder) {
                size_type node = T[preorder];
                size_type node_off = grammar_tree.offset_node(node);
                Q.emplace_back(node, preorder, node_off, occ.off_pattern);
            });
        };
        //initialize the queue
        if (p_occ.preorder == 1) {
            occ.insert(p_occ.off_pattern);
            return;
        }
        Q.emplace_back(p_occ); // insert a primary occ for the node and all its second mentions
        insert_second_mentions(p_occ);

        while (!Q.empty()) {
            auto top = Q.front(); //first element
            if (top.preorder == 1) { //base case
                occ.insert(top.off_pattern);
            } else {
                //check if parent is run length
                size_type parent = T.parent(top.node);
                size_type preorder_parent = T.pre_order(parent);
                auto rlen = grammar_tree.is_run(preorder_parent);
                // check if parent is run-length
                if( rlen > 0 && preorder_parent + 1 == top.preorder){
                    // if parent is a run-length node
                    size_type fchild_len = grammar_tree.offset_node(T.child(parent, 2)) - top.off_node ;
                    // insert n times, the parent and all its second occ
                    for (size_type i = 0; i < rlen ; ++i) {
                        utils::primaryOcc s_occ;
                        s_occ.preorder = preorder_parent;
                        s_occ.node = parent;
                        s_occ.off_node = top.off_node; //as it is first child is the same offset
                        s_occ.off_pattern = top.off_pattern + fchild_len*i;
                        Q.emplace_back(s_occ);
                        insert_second_mentions(s_occ);
                    }
                }else{
                    if(rlen == 0){
                        // if parent is not a run-length node
                        size_type off_parent = grammar_tree.offset_node(parent);
                        // insert the parent and all its second occ
                        utils::primaryOcc s_occ(parent, preorder_parent, off_parent,top.off_pattern + (top.off_node - off_parent));
                        Q.emplace_back(s_occ);
                        insert_second_mentions(s_occ);
                    }
                    // if parent is a run-length node but the current node is not the first child do not process it
                    // first child insert all the occ in the second child with the parent node
                }


            }
            Q.pop_front();
        }
    }


    void compute_grammar_sfx(utils::nav_grammar &grammar, lpg_build::plain_grammar_t &G,
                             utils::lenght_rules &len, std::vector<utils::sfx> &grammar_sfx) const;



    typedef sdsl::int_vector_buffer<1>                   bvb_t;
    typedef sdsl::int_vector_buffer<>                    ivb_t;
    typedef std::unordered_map<size_type,std::vector<size_type>>  nav_grammar;

    nav_grammar build_nav_grammar(const lpg_build::plain_grammar_t& G, size_type& S) const {

        sdsl::int_vector<> rules_buff;
        std::ifstream _buf(G.rules_file,std::ios::in);
        sdsl::load(rules_buff,_buf);
//        ivb_t rules_buff(G.rules_file);

        size_type zero_count = 0;
        for (size_t i = 0; i < rules_buff.size(); ++i) {
            if(rules_buff[i] == 0)
                zero_count++;
        }
        if(zero_count != 2) std::cout<<"ERROR[RULES FILE] 0 APPEARS MORE THAN 1 TIME IN THE GRAMMAR:"<<zero_count<<std::endl;

        bvb_t rules_lim_buff(G.rules_lim_file);
        size_type id = 0;
        nav_grammar NG;
        std::vector<size_type> right_hand;
        zero_count = 0;
        for (size_t i = 0; i < rules_lim_buff.size(); ++i) {
            right_hand.push_back(rules_buff[i]);
            if(rules_buff[i] == 0) {
                zero_count++;
//                std::cout<<"right hand"<<std::endl;
//                for (const auto &item : right_hand) {
//                    std::cout<<item<<" ";
//                }
//                std::cout<<std::endl;
            }
            if(rules_lim_buff[i] == 1){
                NG[id] = right_hand;
                right_hand.clear();
                id++;
            }
        }
        S = NG[id-1][0];
//        std::cout<<"plain-grammar"<<std::endl;
//        for (const auto &item : NG) {
//            std::cout<<item.first<<"->";
//            for (const auto &second : item.second) {
//                std::cout<<second<<" ";
//            }
//            std::cout<<std::endl;
//        }


        return NG;
    }

    void uncompress_grammar(const std::string & file_dir) const {

        size_type cont = grammar_tree.get_text_len();
        std::fstream fout_lpg(file_dir,std::ios::out|std::ios::binary);
        auto cmp = [this,&cont,&fout_lpg](const uint64_t &prenode, const uint64_t &node,const uint64_t &X) {
            auto c1 = get_symbol(X);//symbols_map[X];
            fout_lpg.write((const char *)&c1,1);
            if(cont == 0)
                return false;
            cont--;
            return true;
        };
        process_suffix_grammar(2,cmp);
    }


};

void lpg_index::compute_grammar_sfx(
        utils::nav_grammar & grammar,
        lpg_build::plain_grammar_t& G,
        utils::lenght_rules& len,
        std::vector<utils::sfx>& grammar_sfx
)const{

//    utils::cuts_rules cuts;
//    utils::build_nav_cuts_rules(G,cuts); //read the cuts of each rule

    const auto& m_tree = grammar_tree.getT();
    sdsl::int_vector_buffer<1> is_rules_len(G.is_rl_file);
    for (const auto &r : grammar) {
        if(G.sym_map.find(r.first) == G.sym_map.end() && r.second.size() > 1) {

            size_type preorder = grammar_tree.first_occ_from_rule(r.first);
            auto node = m_tree[preorder];//preorder select
//            auto run_len = grammar_tree.is_run(preorder);

            if( is_rules_len[r.first] == true ){

                size_type  off = grammar_tree.offset_node(node);
                size_type  _child = m_tree.child(node,2);
                size_type ch_off = grammar_tree.offset_node(_child);
                size_type l = (ch_off-off) * (r.second[1] - 1);
                size_type  _ch_pre = m_tree.pre_order(_child);
                utils::sfx s(ch_off,l,r.second[0],_ch_pre,0);
                grammar_sfx.push_back(s);

            }else{
                size_type acc_len = 0;

                for (size_type i = 1; i < r.second.size(); ++i) {
                    size_type  _child = m_tree.child(node, i + 1);
                    size_type  _ch_pre = m_tree.pre_order(_child);
                    size_type off  = grammar_tree.offset_node(_child);
                    acc_len += len[r.second[i-1]].second;
                    utils::sfx s(
                            off,//rule off
                            len[r.first].second - acc_len,// len of parent - len of prev-sibling
                            r.second[i-1], //prev-sibling id
                            _ch_pre,//preorder
                            0
                    );

                    grammar_sfx.push_back(s);

                }
        }
    }
}


//    utils::dfs(init_r,grammar,[&G,&grammar,&m_tree, &mark,&preorder,&cuts,&len,&grammar_sfx,this](const size_type& id){
//
//        auto it_rule = grammar.find(id);
//        ++preorder;
//        if(mark.find(id) != mark.end())
//            return false; // stop descending if it is second mention
//        mark.insert(id);
//
//        if(is_terminal(id)) return false;
//
//        size_type n_children = it_rule->second.size();
//
//        if(n_children <= 0) return false; // if leaf break and stop descending
//
//        auto node = m_tree[preorder];//preorder select
//        size_type acc_len = 0;
//
//        for(size_type j = 1; j < n_children; ++j){
//            auto _child = m_tree.child(node, j + 1);
//            acc_len += len[it_rule->second[j-1]].second;
//            utils::sfx s(
//                    grammar_tree.offset_node(_child),//rule off
//                    len[it_rule->first].second - acc_len,// len of parent - len of prev-sibling
//                    it_rule->second[j-1], //prev-sibling id
//                    m_tree.pre_order(_child),//preorder
//                    cuts[it_rule->first][j-1]//cut level
//            );
//
//            grammar_sfx.push_back(s);
//#ifdef DEBUG_PRINT
//            std::cout<<"id:"<<id<<":P["<<m_tree.pre_order(_child)<<"]"<<"R["<<it_rule->second[j-1]<<"]";
//            print_suffix_grammar(m_tree.pre_order(_child));
//            std::cout<<std::endl;
//#endif
//        }
//        return true;
//    });
#ifdef DEBUG_INFO
    std::cout<<"compute_grammar_sfx\n";
#endif
}



void lpg_index::locate(const std::string &pattern, std::set<uint64_t> &pos)  const {
        //find primary occ
        auto partitions  = compute_pattern_cuts(pattern);
//        std::cout<<"cortes:\n";

        uint32_t level = partitions.second;
        for (const auto &item : partitions.first) {
//            std::cout<<item<<" ";
            grid_query range{};
            //range search
            if(search_grid_range(pattern.c_str(),pattern.size(),item + 1,level, range)){
                std::vector<utils::primaryOcc> pOcc;
                // grid search
                grid_search(range,item + 1,pattern.size(),level,pOcc);
                // find secondary occ
//GTAGGTAAAGAGTTCAACCACCTGGAAAAAAGAATAGAGAATTTAAATAAAAAAGTTGATGATGGTTTCCTGGACATTTGGACTTACAATGCCGAACTGT"
                for (const auto &occ : pOcc) {
                    find_secondary_occ(occ,pos);
                }
            }
        }
//        std::cout<<std::endl;
}



void lpg_index::locate_all_cuts(const std::string &pattern, std::set<uint64_t> &pos)  const {
    //find primary occ
//    auto partitions  = compute_pattern_cuts(pattern);
    uint32_t level = 0;
//    std::cout<<"ptt:"<<pattern<<std::endl;
//    std::cout<<"cortes occ primarias:";
    for (uint64_t item = 0; item < pattern.size() - 1; ++item) {
        grid_query range{};
        //range search
        if(search_grid_range(pattern.c_str(),pattern.size(),item + 1,level, range)){
            std::vector<utils::primaryOcc> pOcc;
            // grid search
            grid_search(range,item + 1,pattern.size(),level,pOcc);
            // find secondary occ
//            std::cout<<item<<" ";

            if( item == 6 ){
//                std::cout<<"66666666666666666666666666666666666666\n";
                    uint64_t init_preorder  = 2256704;
                    uint64_t node_1  = grammar_tree.getT().operator[](init_preorder);
                    uint64_t init_parent = grammar_tree.getT().parent(node_1);
                    uint64_t chr  = grammar_tree.getT().childrank(node_1);
                    uint64_t prev_sib  = grammar_tree.getT().child(init_parent,chr - 1);
                    uint64_t prev_sib_pre  = grammar_tree.getT().pre_order(prev_sib);
                    uint64_t foccX = grammar_tree.first_occ_preorder_node(prev_sib_pre);
                    uint64_t foccY = grammar_tree.first_occ_preorder_node(init_preorder);

//                    auto f = [this](const size_type& node){
//
//                        size_type preorder = grammar_tree.getT().pre_order(node);
//                        size_type X = grammar_tree.get_rule_from_preorder_node(preorder);
//                        std::cout<<X<<"["<<preorder<<"]:->";
//                        size_type n = grammar_tree.getT().children(node);
//
//                        for (uint64_t i = 1; i <= n; ++i){
//
//                            size_type node_ch = grammar_tree.getT().child(node,i);
//                            size_type pre_node_ch = grammar_tree.getT().pre_order(node_ch);
//                            X = grammar_tree.get_rule_from_preorder_node(pre_node_ch);
//                            std::cout<<X<<"["<<pre_node_ch<<"]"<<" ";
//
//                        }
////                        std::cout<<std::endl;
//                        return true;
//                    };
//                print_prefix_rule(prev_sib_pre,1000);
//                print_suffix_grammar(init_preorder,1000);

//                std::cout<<"left"<<std::endl;
//                down_up_print(foccX,"");
//                auto nn = grammar_tree.getT().operator[](foccX);
//                uint64_t nnend = 0;
//                grammar_tree.getT().subtree(nn,nnend);
//                while(nn <= nnend)
//                    std::cout<<grammar_tree.getT().bit_vector[nn++];
//                std::cout<<std::endl;
//                std::cout<<"rigth"<<std::endl;
//                down_up_print_mirror(foccY,"");
            }


            for (const auto &occ : pOcc) {
                find_secondary_occ(occ,pos);
            }
        }
    }

//    std::cout<<"\n";
}

void lpg_index::locate_split_time(const std::string &pattern, std::set<uint64_t> &pos, uint64_t& p_time, uint64_t& s_time) const {

    auto start = std::chrono::high_resolution_clock::now();
    //find primary occ
    std::vector<utils::primaryOcc> prim_occ;
    auto partitions  = compute_pattern_cuts(pattern);
    uint32_t level = partitions.second;
    for (const auto &item : partitions.first) {
        grid_query range{};
        //range search
        if(search_grid_range(pattern.c_str(),pattern.size(),item + 1,level, range)){
            std::vector<utils::primaryOcc> pOcc;
            // grid search
            grid_search(range,item + 1,pattern.size(),level,pOcc);
            // find secondary occ
            size_t prev_size = prim_occ.size();
            prim_occ.resize(prev_size + pOcc.size());
            std::copy(pOcc.begin(),pOcc.end(),prim_occ.begin()+prev_size);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    p_time += elapsed;

    start= std::chrono::high_resolution_clock::now();
    for (const auto &occ : prim_occ) {
        find_secondary_occ(occ,pos);
    }
    end = std::chrono::high_resolution_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    s_time += elapsed;

}


#endif //LMS_GRAMMAR_REP_HPP
