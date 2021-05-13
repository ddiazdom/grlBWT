//
// Created by ale on 11-05-21.
//

#ifndef LPG_COMPRESSOR_UTILS_HPP
#define LPG_COMPRESSOR_UTILS_HPP

#include "lpg_build.hpp"
#include "dfuds_tree.hpp"
#include "grid.hpp"
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/rmq_succinct_sada.hpp>

namespace utils {

    typedef size_t          size_type;
    typedef sdsl::int_vector_buffer<1>                   bvb_t;
    typedef sdsl::int_vector_buffer<>                    ivb_t;
    typedef std::map<size_type,std::vector<size_type>>  nav_grammar;
    typedef std::unordered_map<size_type,std::pair<size_type,size_type>> lenght_rules; //off,len
    typedef std::unordered_map<size_type,std::vector<size_type>> cuts_rules; //off,len

    struct sfx{
        size_type off{0};
        size_type len{0};
        size_type rule{0};
        size_type preorder{0};
        size_type level{0};


        sfx() = default;
        sfx(const size_type& _off,const size_type& _l,const size_type& _r,const size_type& _p,const size_type& _lv):off(_off),len(_l),rule(_r),preorder(_p),level(_lv){}
        sfx(const sfx& s) = default;

        std::string expansion(const std::string& text) const
        {
            std::string ss;ss.resize(len);
            std::copy(text.begin()+off, text.begin()+off+len,ss.begin());
            return ss;
        }
        std::string expansion(const unsigned char* text) const
        {
            std::string ss;ss.resize(len);
            std::copy(text+off, text+off+len,ss.begin());
            return ss;
        }
    };

    void build_nav_cuts_rules(lpg_build::plain_grammar_t& ,cuts_rules& );
    nav_grammar build_nav_grammar(const lpg_build::plain_grammar_t&);
    void sort_suffixes(const std::string& i_file, std::vector<sfx> &grammar_sfx);
    template <typename F>void dfs_posorder(const size_type& , nav_grammar&, const F&);
    template <typename F>void dfs(const size_type& , nav_grammar&, const F&);
    template <typename F>void dfs_2v(const size_type& , nav_grammar&, const F&);
    template<typename O>void pretty_printer_v( O&bv, const std::string &header);
    template<typename O>void pretty_printer_bv( O&bv, const std::string &header);
    size_type readFile(const std::string& i_file,unsigned char ** text);
    void compute_grammar_sfx(nav_grammar &,lpg_build::plain_grammar_t&,const dfuds_tree&,lenght_rules&,std::vector<sfx>&);
    void compute_grid_points(const std::vector<sfx>&, std::vector<grid_point>&);


    /**
     * IMPLEMENTATIONS
     * */


    void compute_grid_points(const std::vector<sfx>& suffixes,std::vector<grid_point>&points){
        points.resize(suffixes.size());
        for (size_type i = 0; i < suffixes.size() ; ++i) {
            points[i] = grid_point(suffixes[i].rule,i+1,suffixes[i].preorder,suffixes[i].level);
        }
    }

    nav_grammar build_nav_grammar(const lpg_build::plain_grammar_t& G){
        ivb_t rules_buff(G.rules_file);
        bvb_t rules_lim_buff(G.rules_lim_file);
        size_type id = 0;
        nav_grammar NG;
        std::vector<size_type> right_hand;
        for (auto i = (size_type)G.sigma; i < rules_lim_buff.size(); ++i) {
            right_hand.push_back(rules_buff[i]);
            if(rules_lim_buff[i] == 1){
                NG[id + G.sigma ] = right_hand;
                right_hand.clear();
                id++;
            }
        }
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

    void build_nav_cuts_rules(lpg_build::plain_grammar_t& G,cuts_rules& cuts){
        ivb_t breaks_buff(G.lvl_breaks_file);
        size_type i = 0;
        while ( i < breaks_buff.size()) {
            size_type rule = breaks_buff[i];
            size_type l = breaks_buff[++i];
            std::cout<<rule<<"["<<l<<"]\n";
            cuts[rule] = std::vector<size_type>(l,0);
            for (int j = 0; j < l; ++j) {
                cuts[rule][j] = breaks_buff[j+i+1];
            }
            i += l + 1;
        }
    }

    void compute_grammar_sfx(
                             nav_grammar & grammar,
                             lpg_build::plain_grammar_t& G,
                             const dfuds_tree& m_tree,
                             lenght_rules& len,
                             std::vector<sfx>& grammar_sfx
                             ){
        std::cout<<"compute_grammar_sfx\n";
        cuts_rules cuts;
        build_nav_cuts_rules(G,cuts); //read the cuts of each rule
        uint64_t preorder = 0;
        std::set<uint64_t> mark;

        dfs(G.sigma,grammar,[&G,&grammar,&m_tree, &mark,&preorder,&cuts,&len,&grammar_sfx](const size_type& id){

            auto it_rule = grammar.find(id);
            ++preorder;
            if(mark.find(id) != mark.end())
                return false; // stop descending if it is second mention
            mark.insert(id);

            if(id < G.sigma) return false;

            size_type n_children = it_rule->second.size();
            if(n_children <= 0) return false; // if leaf break and stop descending
            auto node = m_tree[preorder];//preorder select
            for(size_type j = 1; j < n_children; ++j){
                auto _child = m_tree.child(node, j + 1);
                sfx s(
                        len[it_rule->second[j]].first, //rule off
                        len[it_rule->first].second - len[it_rule->second[j-1]].second,// len of parent - len of prev-sibling
                        it_rule->second[j-1], //prev-sibling id
                        m_tree.pre_order(_child),//preorder
                        cuts[it_rule->first][j-1]//cut level
                    );
                grammar_sfx.push_back(s);
            }
            return true;
        });
    }

    void sort_suffixes(const std::string& i_file, std::vector<sfx> &grammar_sfx) {
        unsigned char * text = nullptr;
        size_type len;
        len = readFile(i_file, &text);

        if(text == nullptr) throw "ERROR TEXT INVALID";

        sdsl::int_vector<> SA;
        sdsl::int_vector<> SA_1;
        sdsl::lcp_bitcompressed<> LCP;
        sdsl::rmq_succinct_sada<> rmq;

        sdsl::cache_config config(false, ".", "cache-normal");
        sdsl::store_to_file((const char *)text, sdsl::conf::KEY_TEXT);
        sdsl::register_cache_file(sdsl::conf::KEY_TEXT, config);
        sdsl::construct(LCP, sdsl::conf::KEY_TEXT, config, 1);
        if (sdsl::cache_file_exists(sdsl::conf::KEY_SA, config)) {
            sdsl::load_from_cache(SA, sdsl::conf::KEY_SA, config);
            SA_1 = SA;
            for (uint32_t i = 0; i < SA.size(); i++) {
                SA_1[SA[i]] = i;
            }
            sdsl::util::clear(SA);
        }

        // Builds the RMQ Support.
        rmq = sdsl::rmq_succinct_sada<>(&LCP);
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_SA, config));
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_TEXT, config));
        sdsl::remove(sdsl::cache_file_name(sdsl::conf::KEY_LCP, config));

        std::sort(grammar_sfx.begin(),grammar_sfx.end(),
                  [&SA_1,&text,&LCP,&rmq](const sfx & a,
                                          const sfx &b )->bool{
                      /*
                       * offset of suffix in the text
                       * */
                      uint64_t a_pos = a.off;
                      uint64_t b_pos = b.off;

                      uint64_t size_a = a.len;
                      uint64_t size_b = b.len;

                      /*
                       * is start at the same position return the shortest
                       * */

                      if(a_pos == b_pos)
                          return size_a < size_b;

                      uint32_t _rmq;
                      if (SA_1[a_pos] < SA_1[b_pos]) {
                          _rmq = rmq(SA_1[a_pos] + 1, SA_1[b_pos]);
                      } else {
                          _rmq = rmq(SA_1[b_pos] + 1, SA_1[a_pos]);
                      }
                      if (size_a <= LCP[_rmq] && size_b <= LCP[_rmq]) {
                          return size_a < size_b;
                      } else if (size_a <= LCP[_rmq]) {
                          return true;
                      } else if (size_b <= LCP[_rmq]) {
                          return false;
                      } else {
                          /***
                           * Neither is a prefix of the other. Use ISA to find
                           *the order
                           ***/
                          return SA_1[a_pos] < SA_1[b_pos];
                      }
                  });
        std::cout<<"sort_suffixes\n";

    }



    template <typename F>
    void dfs_posorder(const size_type& id, nav_grammar& G, const F& f)  {
        std::function< void (const size_type& )> dfs_aux;
        dfs_aux = [&dfs_aux,& G,& f](const size_type& id){

            auto it = G.find(id);
            if(it != G.end()){
                // non - terminal case
                for (const auto &item : it->second) dfs_aux(item);
            }
            f(id);
        };
        dfs_aux(id);
    }


    template <typename F>
    void dfs(const size_type& id, nav_grammar& G, const F& f)  {
        std::function< void (const size_type& )> dfs_aux;
        dfs_aux = [&dfs_aux,& G,& f](const size_type& id){
            if(!f(id)) return ;
            auto it = G.find(id);
            if(it != G.end()){
                // non - terminal case
                for (const auto &item : it->second) dfs_aux(item);
            }
        };
        dfs_aux(id);
    }

    template <typename F>
    void dfs_2v(const size_type& id, nav_grammar& G, const F& f)  {
        std::function< void (const size_type& )> dfs_aux;
        dfs_aux = [&dfs_aux,& G,& f](const size_type& id){
            if(!f(id,1)) return ;
            auto it = G.find(id);
            if(it != G.end()){
                // non - terminal case
                for (const auto &item : it->second){
                    dfs_aux(item);
                }
            }
            f(id,0);
        };

        dfs_aux(id);

    }

//    size_type compute_rule_offsets( const plain_grammar& _Gr,utils::nav_grammar& grammar, const size_t& text_length,
//                                    Roff& rules_off,
//                                    sdsl::bit_vector _l) {
//
//
//        _l(text_length + 1, 0); // text-len bitvector marking start position of parser phrases..
//        size_type c_nodes = 0; // number of nodes in the grammar tree.
//        std::set<size_type> M; // auxiliar set to mark first occ of rules in the tree
//        size_type l_pos = 0; // current off in the text.
//        /*
//         * Count the number of nodes of compact tree representation
//         * Build bitvector L for marking the pos of every no terminal in text
//         *
//         * */
//        std::cout << "* Count the number of nodes of compact tree representation\n";
//        utils::dfs_2v(_Gr.sigma,grammar,[&_Gr,&_l, &l_pos, &M,&rules_off, &c_nodes](const size_type& id, bool first_visit){
//            if(first_visit){ // first visit
//                ++c_nodes; // count the node
//                if (M.find(id) != M.end()) { // second mention of a non-terminal => leaf tree
//                    _l[l_pos] = 1; // mark text position for this leaf tree
//                    l_pos += rules_off[id].second; // increase the off-text with pre-computed len-rule
//                    return false; // return false stop descending in the tree
//                }
//                M.insert(id); // first mention of a rule tree
//                if(id < _Gr.sigma){
//                    //if is a terminal symbol
//                    rules_off[id] = std::make_pair(l_pos,1); // store the off in the text
//                    ++l_pos;
//                    return false; // stop descending in the tree
//                }
//                rules_off[id] = std::make_pair(l_pos,1); // store the off in the text
//                return true; //  keep descending in the tree
//            }
//            //second visit (only in non-terminal)
//            rules_off[id].second = l_pos - rules_off[id].first; // update rule length with the offset diferences
//            return false; // to avoid warning.... :)
//        });
//
//        return c_nodes;
//    }



    size_type readFile(const std::string& i_file,unsigned char ** text){
        i_file_stream<uint8_t> if_stream(i_file, BUFFER_SIZE);
        *text = new unsigned char [if_stream.tot_cells];
        for(size_t i=0;i<if_stream.tot_cells;i++){
            (*text)[i] = if_stream.read(i);
        }
        return if_stream.tot_cells;
    }
    template<typename O>
    void pretty_printer_v( O&bv, const std::string &header){

        std::cout<<header<<std::endl;

        auto dec = [](const uint64_t & t){
            uint64_t tt = t, c = 0;
            while (tt/10){
                tt /= 10;
                c++;
            }
            return c;
        };

        std::cout<<"idx:";

        for(int i = 0 ; i < bv.size(); i++){
            uint64_t c = dec(bv[i]);
            for (int j = 0; j < c  ; ++j) {
                std::cout<<" ";
            }
            std::cout<<i%10<<" ";
        }
        std::cout<<std::endl<<"arr:";
        for(int i =0 ; i < bv.size(); i++)
            std::cout<<bv[i]<<" ";
        std::cout<<std::endl;

    }
    template<typename O>
    void pretty_printer_bv( O&bv, const std::string &header){

        std::cout<<header<<std::endl;
        std::cout<<"idx:";
        for(int i = 0 ; i < bv.size(); i++)std::cout<<i%10;
        std::cout<<std::endl;
        uint64_t acc = 0; //rank
        std::cout<<"rnk:";
        for(int i = 0 ; i < bv.size(); i++){
            std::cout<<acc%10;
            acc += bv[i];
        }
        std::cout<<acc%10;
        std::cout<<std::endl<<"arr:";
        for(int i = 0 ; i < bv.size(); i++){
            std::cout<<bv[i];
        }
        std::cout<<std::endl;
    }

    template<typename K>
    bool lower_bound( uint64_t &lr,uint64_t  hr, const K &f)  {
        bool found = false;
        while(lr < hr){
            uint64_t mid = (lr+hr)/2;
            int c = f(mid); // return 1 if  str < x
            if(c > 0){
                hr = mid -1;
            }else{
                if(c < 0) // return -1 if x > str
                    lr = mid+1;
                else{
                    hr = mid;
                    found = true;
                }
            }
        }
        if(found) return true;

        if(!found && lr == hr && f(lr) == 0 ) {
            return true;
        }
        return false;
    }

    template<typename K>
    bool upper_bound(uint64_t  lr,uint64_t  &hr, const K &f)  {
        bool found = false;
        while(lr < hr){
            uint64_t mid = ceil((lr+hr)/2.0);
            int c = f(mid);
            if(c > 0){
                hr = mid -1;
            }else{
                if(c < 0)
                    lr = mid+1;
                else{
                    lr = mid;
                    found = true;
                }
            }
        }
        if(found) return true;

        if(!found && lr == hr && f(lr) == 0 ) {
            return true;
        }
        return false;

    }

}




#endif //LPG_COMPRESSOR_UTILS_HPP
