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

    void build_index(const std::string& i_file,plain_grammar_t &p_gram,const size_t& text_length, sdsl::cache_config& config){
        m_sigma = p_gram.sigma;
        parsing_rounds = p_gram.rules_per_level.size();
        symbols_map.resize(p_gram.symbols_map.size());
        for(size_t i=0; i < p_gram.symbols_map.size(); i++){
            symbols_map[i] = p_gram.symbols_map[i];
        }
        //build the grammar tree and grid from p_gram here!!

        utils::lenght_rules lenghts;
        utils::nav_grammar NG = utils::build_nav_grammar(p_gram);
        grammar_tree.build(NG,p_gram,text_length,lenghts);

        std::vector<utils::sfx> grammar_sfx;
        const auto& T = grammar_tree.getT();
        utils::compute_grammar_sfx(NG,p_gram,T,lenghts,grammar_sfx);
        NG.clear();
        utils::sort_suffixes(i_file,grammar_sfx);
        std::cout<<"sort_suffixes\n";
        std::vector<grid_point> points;
        compute_grid_points(grammar_sfx,points);
        std::cout<<"compute_grid_points\n";
        grammar_sfx.clear();
        grid = grid_t(points,p_gram.rules_per_level.size());
        std::cout<<"build grid\n";
        breakdown_space();
        std::cout<<"end build index\n";
    };

    void breakdown_space() const {
        std::cout<<"breakdown-space-index\n";
        std::cout<<"Text-length,"<<grammar_tree.get_text_len()<<std::endl;
        std::cout<<"Number-rules,"<<grammar_tree.get_size_rules()<<std::endl;
        std::cout<<"Grammar-size,"<<grammar_tree.get_grammar_size()<<std::endl;
        std::cout<<"Grammar-Tree,"<<sdsl::size_in_bytes(grammar_tree)<<std::endl;
        grammar_tree.breakdown_space();
        std::cout<<"Grid,"<<sdsl::size_in_bytes(grid)<<std::endl;;
        grid.breakdown_space();
        std::cout<<"symbols_map,"<<sdsl::size_in_bytes(symbols_map);
        std::cout<<"m_sigma,"<<sizeof(m_sigma);
        std::cout<<"parsing_rounds,"<<sizeof(parsing_rounds);
        std::cout<<"rl_compressed,"<<sizeof(rl_compressed);
    }

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

    lpg_index(std::string &input_file, std::string &tmp_folder, size_t n_threads, float hbuff_frac, bool rl_comp) {

        std::cout<<"Input file: "<<input_file<<std::endl;
        auto alphabet =  get_alphabet(input_file);

        size_t n_chars = 0;
        for(auto const sym : alphabet) n_chars+=sym.second;

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
        auto hbuff_size = size_t(std::ceil(float(n_chars)*hbuff_frac));

        std::cout<<"Computing the LPG grammar"<<std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        lpg_build::compute_LPG(input_file, g_file, n_threads, config, hbuff_size, alphabet,  rl_comp);
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


//        plain_gram.print_grammar();


        std::cout<<"Building the self-index"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        build_index(input_file,plain_gram,n_chars, config);
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
            std::cout<<pattern<<":";
            std::set<size_type> occ;
            locate(pattern,occ);
            for (const auto &item : occ) {
                std::cout<<item<<" ";
            }
            std::cout<<std::endl;
        }
    }

    //search for a list of pattern (input is a file each line is a pattern)
    void search(std::string& input_file){
        std::fstream in(input_file,std::ios::in);
        std::string pattern;
        while(in.good() && std::getline(in,pattern)){
            std::cout<<pattern<<":";
            std::set<size_type> occ;
            locate(pattern,occ);
            for (const auto &item : occ) {
                std::cout<<item<<" ";
            }
            std::cout<<std::endl;
        }

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

    template<typename F>
    bool dfs_leaf_base_case(const uint64_t &preorder_node, const uint64_t &node, const F &f)const {
        size_type _x = grammar_tree.get_rule_from_preorder_node(preorder_node);//get the rule of the leaf
        if(isTerminal(_x)){//check if the rule is a terminal symbol
            //terminal case
            bool keep = f(preorder_node, node,_x);//call process function
            return keep; // keep says if we must to stop the expansion of the rule
        }
        //second occ.
        auto fpre_node = grammar_tree.first_occ_from_rule(_x);//get preorder of first mention node
        const auto&T = grammar_tree.getT();
        auto fnode = T[fpre_node];//preorder select

        return dfs_mirror_leaf(fpre_node, fnode, f); // call recursive in the new node...
    }

    template<typename F>
    bool dfs_mirror_leaf(const uint64_t &preorder_node, const uint64_t &node, const F &f) const {
            if (grammar_tree.isLeaf(preorder_node)) { //leaf case
                return dfs_leaf_base_case(preorder_node,node,f);
            } else {
                    const auto&T = grammar_tree.getT();
                    uint32_t n = T.children(node); //number of children
                    for (uint32_t i = n; i > 0; --i) { // from right to left
                        auto chnode = T.child(node, i); //get i-child node
                        bool keep = dfs_mirror_leaf(T.pre_order(chnode), chnode, f);// compute preorder and visit recursively
                        if (!keep) return keep; // keep says if we must to stop the expansion of the rule
                    }
                    return true;
            }
    }

    template<typename F>
    bool dfs_leaf(const uint64_t &preorder_node, const uint64_t &node, const F &f) const {
        if (grammar_tree.isLeaf(preorder_node)) { //leaf case
            return dfs_leaf_base_case(preorder_node,node,f);
        }  else {

            const auto&T = grammar_tree.getT();
            uint32_t n = T.children(node); //number of children

            for (uint32_t i = 0; i > n; ++i) { // from left to right
                auto chnode = T.child(node, i+1); //get i-child node
                bool keep = dfs_leaf(T.pre_order(chnode), chnode, f);// compute preorder and visit recursively
                if (!keep) return keep; // keep says if we must to stop the expansion of the rule
            }
            return true;
        }
    }

    int cmp_prefix_rule(const size_type &preorder_node, const std::string& str, const uint32_t &i) const  {
        long ii = i;
        int r = 0;
        bool match = false;
        const auto& m_tree = grammar_tree.getT();
        auto node = m_tree[preorder_node];
        auto cmp = [&str, &r, &match, &ii, this](const size_type &preorder_node, const size_type &node, const size_type &X) {
            // X is a terminal rule
            auto c1 = (unsigned char) symbols_map[X];
            auto c2 = (unsigned char) str[ii];
            if (c1 > c2) {r = 1;return false;}
            if (c1 < c2) {r = -1;return false;}
            --ii;
            if ( ii < 0 ) {match = true;return false;}
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
    int cmp_suffix_grammar(const size_type &preorder_node, const std::string& str, const uint32_t &i) const  {
        uint32_t ii = i, sfx_len = str.size();
        int r = 0;
        bool match = false;
        const auto& T = grammar_tree.getT();
        auto node = T[preorder_node];//preorder select
        auto cmp = [&str, &r, &match, &ii, &sfx_len, this](const uint64_t &prenode, const uint64_t &node, const uint64_t &X) {
            auto c1 = (unsigned char) symbols_map[X];
            auto c2 = (unsigned char) str[ii];
            if (c1 > c2) {r = 1;return false;}
            if (c1 < c2) {r = -1;return false;}
            ii++;
            if (ii == sfx_len) {match = true;return false;}
            return true;
        };

        auto fopen = T.bps.find_open(node - 1); // this is to compute parent and childrank
        auto parent = T.pred0(fopen) + 1;//  compute parent
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
    inline bool isTerminal(const size_type& X) const{ return X < m_sigma;}

    /**
     * Find the grid range to search
     * @return grid_query
     * */

    bool search_grid_range(const char* pattern ,const uint32_t & len,const uint32_t & p,const uint32_t & level,grid_query& q) const  {
        // search rules range....
        auto cmp_rev_prefix_rule = [&p, &pattern,this](const size_type &rule_id) {
            // compute node definiton preorder of the rule
            uint64_t prenode = grammar_tree.first_occ_from_rule(rule_id);
            return cmp_prefix_rule(prenode, pattern, p - 1);
        };
        uint64_t row_1 = 1, row_2 = grammar_tree.get_size_rules();
        //search lower
        if (!utils::lower_bound(row_1, row_2, cmp_rev_prefix_rule)) return false;
        q.row1 = row_1;
        //search upper
        row_2 = grammar_tree.get_size_rules();
        if (!utils::upper_bound(row_1, row_2, cmp_rev_prefix_rule)) return false;
        q.row2 = row_2;
        //search suffixes
        auto cmp_suffix_grammar_rule = [ &level, &p, &pattern,this](const size_type &suffix_id,const size_type &val) {
            //val is just to use the std lower bound method
            // compute node definiton preorder of the rule
            uint64_t prenode = grid.get_preorder_node_from_suffix(suffix_id,level);
            return cmp_suffix_grammar(prenode, pattern, p);
        };

        std::vector<size_type> sfx_by_level;
        this->grid.map_suffixes_levels(level,sfx_by_level);
        //search lower
        auto it_low = std::lower_bound(sfx_by_level.begin(),sfx_by_level.end(),0,cmp_suffix_grammar_rule);
        //check we found a match...
        if(it_low == sfx_by_level.end() || cmp_suffix_grammar_rule(*it_low,0) != 0) return false;
        q.col1 = *it_low;
        //search upper
        //remove elements less than lower bound
        sfx_by_level.erase(sfx_by_level.begin(),it_low);
        auto it_upper = std::lower_bound(sfx_by_level.begin(),sfx_by_level.end(),0,cmp_suffix_grammar_rule);
        //check we found a match...
        if(it_upper == sfx_by_level.end() || cmp_suffix_grammar_rule(*it_upper,0) != 0) return false;
        q.col2 = *it_upper;
        return true;
    }


    void locate(const std::string& pattern, std::set<uint64_t> &pos) const {
        std::cout<<"all posible partitions\n";
//        auto partitions  = compute_pattern_cut(pattern);
        size_type l = this->grid.get_levels();
//        std::cout<<partitions.first.size()<<std::endl;
//        uint32_t level = partitions.second;
//        for (const auto &item : partitions.first) {
        for (int i = 0; i < pattern.size() ; ++i) {
            for (int level = 1; level <= l ; ++level) {

                //find primary occ
                grid_query range{};
                //range search
                if(search_grid_range(pattern.c_str(),pattern.size(),i + 1,level, range)){
                    std::vector<utils::primaryOcc> pOcc;
                    // grid search
                    grid_search(range,i,level,pOcc);
                    // find secondary occ
                    for (const auto &occ : pOcc) {
                        find_secondary_occ(occ,pos);
                    }
                }

            }
        }
    }

    void grid_search(const grid_query& range, const uint64_t & pattern_off,const uint32_t &level,std::vector<utils::primaryOcc>& occ) const{
        std::vector<uint64_t> sfx;
        this->grid.search(range,level,sfx);
        occ.resize(sfx.size());
        const auto &T = grammar_tree.getT();
        for (size_type i = 0; i < sfx.size(); ++i) {
            size_type preorder_node = this->grid.get_preorder_node_from_suffix(sfx[i],level);
            size_type node = T[preorder_node];
            assert(grammar_tree.offset_node(node) - pattern_off >= 0);
            size_type parent = T.parent(node);
            size_type off = grammar_tree.offset_node(node);
            occ[i] = utils::primaryOcc(parent,T.pre_order(parent),off,off - pattern_off);
        }
    }

    void find_secondary_occ(const utils::primaryOcc& p_occ, std::set<size_type>& occ) const {
        //queue for node processing
        std::deque<utils::primaryOcc> Q;
        const auto& T = grammar_tree.getT();

        //auxiliar functions
        auto insert_second_mentions = [&Q,&T,this](const utils::primaryOcc& occ){
            grammar_tree.visit_secondary_occ(occ.preorder,[&T,&Q,&occ,this](const size_type& preorder){
                size_type node = T[preorder];
                size_type node_off = grammar_tree.offset_node(T[preorder]);
                Q.emplace_back(node,preorder,node_off,occ.off_pattern);
            });
        };

        //initialize the queue
        if(p_occ.preorder == 1){
            occ.insert(p_occ.off_pattern);
            return;
        }
        Q.emplace_back(p_occ); // insert a primary occ for the node and all its second mentions
        insert_second_mentions(p_occ);

        while(!Q.empty()){
            auto top = Q.front(); //first element
            if(top.preorder == 1){ //base case
                occ.insert(top.off_pattern);
            } else {
                // insert in Q parent
                size_type parent = T.parent(top.node);
                size_type preorder_parent = T.pre_order(parent);
                size_type off_parent = grammar_tree.offset_node(parent);
                utils::primaryOcc s_occ(parent,preorder_parent,off_parent,top.off_pattern + (top.off_node - off_parent ));
                Q.push_back(s_occ);
                // insert in Q second occ of the node
                insert_second_mentions(s_occ);
            }
            Q.pop_front();
        }
    }


};


#endif //LMS_GRAMMAR_REP_HPP
