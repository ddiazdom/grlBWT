//
// Created by Diego Diaz on 4/17/20.
//

#ifndef LPG_COMPRESSOR_ITERATORS_HPP
#define LPG_COMPRESSOR_ITERATORS_HPP

template<class vector_t>
struct huff_vec_iterator {
    typedef std::forward_iterator_tag        iterator_category;
    typedef typename vector_t::symbol_type   value_type;
    typedef const value_type&                reference;
    typedef const value_type*                pointer;

    const vector_t&                          h_vec;
    value_type                               curr_sym;
    size_t                                   bit_pos;
    std::ptrdiff_t                           m_arr_idx;

    huff_vec_iterator(const vector_t& _vector, size_t _idx):
            h_vec(_vector), bit_pos(_vector.get_bit_pos(_idx)), m_arr_idx(_idx){
        ++*this;
    };

    reference operator*() const {
        return curr_sym;
    }

    pointer operator->() const {
        return &(curr_sym);
    }

    huff_vec_iterator& operator+=(ptrdiff_t diff) {
        if((m_arr_idx + diff) < (std::ptrdiff_t)h_vec.size()){
            m_arr_idx +=diff;
            bit_pos = h_vec.get_bit_pos(m_arr_idx-1);
            curr_sym = h_vec.decode_symbol(bit_pos);
        }else{
            make_end();
        }
        return *this;
    }

    huff_vec_iterator& operator-=(ptrdiff_t diff) {
        diff = std::min(diff+1, m_arr_idx);//the index is always one position ahead!!
        m_arr_idx -=diff;
        bit_pos = h_vec.get_bit_pos(m_arr_idx);
        curr_sym = h_vec.decode_symbol(bit_pos);
        m_arr_idx++;
        return *this;
    }

    huff_vec_iterator& operator++() {
        if(m_arr_idx < (std::ptrdiff_t)h_vec.size()){
            curr_sym =  h_vec.decode_symbol(bit_pos);
        }
        m_arr_idx++;
        if(m_arr_idx>(std::ptrdiff_t)h_vec.size()){
            make_end();
        }
        return *this;
    }

    void make_end(){
        m_arr_idx = h_vec.size()+1;
        bit_pos = h_vec.m_written_bits;
        curr_sym=0;
    }

    huff_vec_iterator operator++(int) {
        huff_vec_iterator tmp(*this);
        ++*this;
        return tmp;
    }

    bool operator==(huff_vec_iterator& other) const{
        return m_arr_idx == other.m_arr_idx;
    }

    bool operator!=(const huff_vec_iterator& other) const{
        return m_arr_idx != other.m_arr_idx;
    }
};

template <class tree_t>
struct louds_tree_iterator {

    typedef std::forward_iterator_tag  iterator_category;
    typedef size_t                     value_type;
    typedef const value_type&          reference;
    typedef const value_type*          pointer;

private:
    const tree_t&                      m_tree;
    value_type                         m_id;
    std::ptrdiff_t                     m_map;
    value_type                         m_children;
    value_type                         m_int_cnt;
    value_type                         m_int_rank;

    void make_end(){
        m_id = m_tree.shape.size();
        m_children=0;
        m_int_cnt = m_tree.int_nodes()+1;
        m_int_rank = m_int_cnt;
        m_map = m_tree.nodes()+1;
    }

public:
    const size_t& id                 = m_id;
    const std::ptrdiff_t& map        = m_map;
    const size_t& children           = m_children;
    const size_t& int_rank           = m_int_rank;

    louds_tree_iterator(const tree_t& _tree, size_t _curr_node): m_tree(_tree),
                                                                 m_id(_curr_node),
                                                                 m_map(_tree.nodes()),
                                                                 m_children(0),
                                                                 m_int_cnt(_tree.int_nodes()),
                                                                 m_int_rank(_tree.int_nodes()){
        if(m_id<m_tree.shape_size()){
            m_map = m_tree.node_map(m_id);
            m_children = m_tree.n_children(m_id);
            m_int_cnt = m_tree.n_ints(m_id);
            m_int_rank = m_children==0? 0 : m_int_cnt;
        }
    };

    inline louds_tree_iterator& operator++() {

        if((m_id + m_children + 1) < m_tree.shape.size()){
            m_id += m_children + 1;
            m_children = m_tree.shape[m_id] ? m_tree.n_children(m_id) : 0 ;
            m_int_rank = m_children == 0 ? 0 : ++m_int_cnt;
            m_map++;
        }else{
            make_end();
        }
        return *this;
    }

    louds_tree_iterator& operator+=(std::ptrdiff_t diff){
        if((m_map + diff) < m_tree.nodes()){
            m_map += diff;
            m_id = m_tree.node_select(m_map);
            m_children = m_tree.n_children(m_id);
            m_int_cnt = m_tree.n_ints(m_id);
            m_int_rank = m_children==0 ? 0 : m_int_cnt;
        }else{
            make_end();
        }
    }

    louds_tree_iterator& operator-=(std::ptrdiff_t diff){
        diff = std::min(diff, m_map-1);
        m_map -= diff;
        m_id = m_tree.node_select(m_map);
        m_children = m_tree.n_children(m_id);
        m_int_cnt = m_tree.n_ints(m_id);
        m_int_rank = m_children==0 ? 0 : m_int_cnt;
        return *this;
    }

    louds_tree_iterator operator++(int) {
        louds_tree_iterator tmp(*this);
        ++*this;
        return tmp;
    }

    inline bool operator==(louds_tree_iterator &other) const{
        return other.m_id == m_id;
    }

    inline bool operator!=(louds_tree_iterator &other) const{
        return m_id != other.m_id;
    }
};
#endif //LPG_COMPRESSOR_ITERATORS_HPP
