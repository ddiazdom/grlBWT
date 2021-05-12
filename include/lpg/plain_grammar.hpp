//
// Created by ale on 11-05-21.
//

#ifndef LPG_COMPRESSOR_PLAIN_GRAMMAR_HPP
#define LPG_COMPRESSOR_PLAIN_GRAMMAR_HPP

#include "lpg_build.hpp"

namespace utils {

    typedef size_t          size_type;
    typedef sdsl::int_vector_buffer<1>                   bvb_t;
    typedef sdsl::int_vector_buffer<>                    ivb_t;

    typedef std::map<size_type,std::vector<size_type>>  nav_grammar;

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
        return NG;
    }



    template <typename F>
    void dfs(const size_type& id, nav_grammar& G, const F& f) const {
        std::function< void (const size_type& )> dfs_aux;
        dfs_aux = [&dfs_aux,& G,& f](const size_type& id){
            f(id);
            auto it = G.find(id);
            if(it != G.end()){
                // non - terminal case
                for (const auto &item : it->second) dfs_aux(item)
            }
        };

    }

    template <typename F>
    void dfs_posorder(const size_type& id, nav_grammar& G, const F& f) const {
        std::function< void (const size_type& )> dfs_aux;
        dfs_aux = [&dfs_aux,& G,& f](const size_type& id){

            auto it = G.find(id);
            if(it != G.end()){
                // non - terminal case
                for (const auto &item : it->second) dfs_aux(item)
            }
            f(id);
        };

    }


    template <typename F>
    void dfs(const size_type& id, nav_grammar& G, const F& f) const {
        std::function< void (const size_type& )> dfs_aux;
        dfs_aux = [&dfs_aux,& G,& f](const size_type& id){
            if(!f(id)) return ;
            auto it = G.find(id);
            if(it != G.end()){
                // non - terminal case
                for (const auto &item : it->second) dfs_aux(item)
            }
        };

    }

    template <typename F>
    void dfs_2v(const size_type& id, nav_grammar& G, const F& f) const {
        std::function< void (const size_type& )> dfs_aux;
        dfs_aux = [&dfs_aux,& G,& f](const size_type& id){
            if(!f(id,0)) return ;

            auto it = G.find(id);
            if(it != G.end()){
                // non - terminal case
                for (const auto &item : it->second) dfs_aux(item);
            }
            f(id,1);
        };

    }

//    typedef std::unordered_map<size_type,std::pair<size_type,size_type>> Roff; //off,len
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


}




#endif //LPG_COMPRESSOR_PLAIN_GRAMMAR_HPP
