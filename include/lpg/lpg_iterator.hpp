//
// Created by Diego Diaz on 4/16/20.
//

#ifndef LPG_COMPRESSOR_LPG_ITERATOR_HPP
#define LPG_COMPRESSOR_LPG_ITERATOR_HPP

#include <cstdlib>
#include "cdt/huff_vector.hpp"
#include "cdt/louds_tree.hpp"

//same as before but this also includes the pointers of the leaves
template <class grammar_t>
struct grammar_iterator {

private:
    typedef std::forward_iterator_tag   iterator_category;
    typedef size_t                      value_type;
    typedef const value_type&           reference;
    typedef const value_type*           pointer;

    const grammar_t&                    grammar;
    size_t                              tree_shape_size;
    louds_tree_iterator<louds_tree>     bfs_it;
    huff_vector<>::const_iterator       g_lab_it;
    value_type                          m_label;
    const value_type&     int_rank    = bfs_it.int_rank;

public:

    const value_type&     id          = bfs_it.id;
    const std::ptrdiff_t& map         = bfs_it.map;
    const value_type&     children    = bfs_it.children;
    const value_type&     label       = m_label;

    grammar_iterator(const grammar_t& _grammar, size_t _curr_node):
                     grammar(_grammar),
                     tree_shape_size(_grammar.m_grammar_tree.shape.size()),
                     bfs_it(grammar.m_grammar_tree, _curr_node),
                     g_lab_it(grammar.m_l_labels, 0),
                     m_label(0){

        if(id<tree_shape_size){
            size_t n_labs = grammar.m_grammar_tree.n_leaves(id);
            if(n_labs>1){
                g_lab_it+=n_labs-1;
                if(children==0){
                    m_label = *g_lab_it;
                }else{
                    m_label =int_rank+grammar.sigma;
                }
                ++g_lab_it;
            }else{
                m_label =int_rank+grammar.sigma;
            }
        }else{
            g_lab_it += grammar.m_l_labels.size();
        }
    };

    inline grammar_iterator& operator+=(ptrdiff_t diff) {
        diff = std::min(diff, (ptrdiff_t)(grammar.m_grammar_tree.nodes()-map+1));
        size_t old_idx = grammar.m_grammar_tree.n_leaves(bfs_it.id);
        bfs_it+=diff;
        size_t new_idx = grammar.m_grammar_tree.n_leaves(bfs_it.id);
        g_lab_it+= new_idx - old_idx;

        if(children>0){
            m_label = int_rank+grammar.sigma;
        }else{
            m_label = *g_lab_it;
            ++g_lab_it;
        }
        return *this;
    }

    inline grammar_iterator& operator-=(ptrdiff_t diff) {
        diff = std::min(bfs_it.map, diff);
        size_t old_idx = grammar.m_grammar_tree.n_leaves(bfs_it.id);
        bfs_it-=diff;
        size_t new_idx = grammar.m_grammar_tree.n_leaves(bfs_it.id);
        g_lab_it-= old_idx - new_idx;

        if(children>0){
            m_label = int_rank+grammar.sigma;
        }else{
            m_label = *g_lab_it;
            ++g_lab_it;
        }
        return *this;
    }

    inline grammar_iterator& operator++() {
        ++bfs_it;
        if(id<tree_shape_size){
            if(children == 0){
                m_label = *g_lab_it;
                ++g_lab_it;
            }else{
                assert(int_rank!=0);
                m_label = int_rank+grammar.m_sigma;
            }
        }
        return *this;
    }

    grammar_iterator operator++(int) {
        grammar_iterator tmp(*this);
        ++*this;
        return tmp;
    }

    inline bool operator==(grammar_iterator& other) const{
        return id==other.id;
    }

    inline bool operator!=(grammar_iterator& other) const{
        return id!=other.id;
    }
};
#endif //LPG_COMPRESSOR_LPG_ITERATOR_HPP
