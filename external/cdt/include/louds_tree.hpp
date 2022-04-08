//
// Created by diediaz on 04-12-19.
//

#ifndef LOUDS_TREE_HPP
#define LOUDS_TREE_HPP

#include<sdsl/bit_vectors.hpp>

class louds_tree {

private:
    typedef sdsl::bit_vector                bit_vector_t;
    bit_vector_t                            m_shape;
    typename bit_vector_t::select_0_type    ss_0;
    typename bit_vector_t::select_1_type    ss_1;
    typename sdsl::select_support_mcl<10,2> ss_10;
    typename bit_vector_t::rank_0_type      rs_0;
    typename bit_vector_t::rank_1_type      rs_1;
    typename sdsl::rank_support_v5<10,2>    rs_10;
    size_t                                  fchild_locus{};//locus of the first child of the root

public:
    typedef size_t          size_type;
    const sdsl::bit_vector& shape = m_shape;

    inline size_t succ_0(size_t bit_pos) const{

        //return ss_0(rs_0(get_bit_pos+1)+1);
        size_t pos = bit_pos;
        size_t bits = m_shape.get_int(pos);

        /*size_t tmp=node_select(2);
        if(bit_pos<tmp){
            std::cout<<"holaa "<<bit_pos<<std::endl;
        }*/

        while(bits==0xFFFFFFFFFFFFFFFF){
            pos+=64;
            bits = m_shape.get_int(pos);
        }
        pos += sdsl::bits::lo(~bits);
        return pos;
    }

    inline size_t succ_1(size_t bit_pos) const{
        return ss_1(rs_1(bit_pos+1)+1);
    }

    inline size_t pred_0(size_t bit_pos) const{
        //return ss_0(rs_0(get_bit_pos));
        //!Note: it is faster doing a linear scan as the
        //! degree of the nodes is not too high
        size_t pos = bit_pos>=64 ? bit_pos-64 : 0;
        size_t bits = m_shape.get_int(pos);

        /*size_t tmp=node_select(2);
        if(bit_pos<tmp){
            std::cout<<"holaa2 "<<bit_pos<<std::endl;
        }*/

        while(bits==0xFFFFFFFFFFFFFFFF && pos >=64){
            //std::cout<<"holaaa"<<std::endl;
            pos -=64;
            bits = m_shape.get_int(pos);
        }

        if(pos>0 && pos<64){//border case
            bits =  ((1ULL<<(64U-(pos+1U)))-1U) << (pos+1U) | m_shape.get_int(0);
            pos = 0;
        }

        return pos + sdsl::bits::hi(~bits);
    }

    inline size_t pred_1(size_t bit_pos) const{
        return ss_1(rs_1(bit_pos));
    }

    void init_tree_support(){
        sdsl::util::init_support(ss_0, &m_shape);
        sdsl::util::init_support(ss_1, &m_shape);
        sdsl::util::init_support(rs_0, &m_shape);
        sdsl::util::init_support(rs_1, &m_shape);
        sdsl::util::init_support(rs_10, &m_shape);
        sdsl::util::init_support(ss_10, &m_shape);
    }

    void swap(louds_tree&& other){
        std::swap(other.m_shape, m_shape);
        std::swap(other.ss_0, ss_0);
        std::swap(other.ss_1, ss_1);
        std::swap(other.rs_0, rs_0);
        std::swap(other.rs_1, rs_1);
        std::swap(other.rs_10, rs_10);
        std::swap(other.ss_10, ss_10);
        std::swap(other.fchild_locus, fchild_locus);
        ss_0.set_vector(&m_shape);
        ss_1.set_vector(&m_shape);
        rs_0.set_vector(&m_shape);
        rs_1.set_vector(&m_shape);
        rs_10.set_vector(&m_shape);
        ss_10.set_vector(&m_shape);
    }


    louds_tree()= default;

    explicit louds_tree(bit_vector_t &input_bv){
        m_shape.resize(input_bv.size());
        for(size_t i=0;i<input_bv.size();i++){
            m_shape[i] = input_bv[i];
        }
        init_tree_support();
        fchild_locus = child(2,1);
    }

    explicit louds_tree(bit_vector_t &&input_bv) : m_shape(input_bv) {
        init_tree_support();
        fchild_locus = child(2,1);
    }

    louds_tree(louds_tree&& other) noexcept{
        swap(std::forward<louds_tree>(other));
    }

    louds_tree& operator=(louds_tree&& other) noexcept{
        swap(std::forward<louds_tree>(other));
        return *this;
    }

    louds_tree& operator=(const louds_tree& other) noexcept{

        if(this!=&other) {
            m_shape = other.m_shape;
            ss_0 = other.ss_0;
            ss_1 = other.ss_1;
            rs_0 = other.rs_0;
            rs_1 = other.rs_1;
            rs_10 = other.rs_10;
            ss_10 = other.ss_10;
            ss_0.set_vector(&m_shape);
            ss_1.set_vector(&m_shape);
            rs_0.set_vector(&m_shape);
            rs_1.set_vector(&m_shape);
            rs_10.set_vector(&m_shape);
            ss_10.set_vector(&m_shape);
            fchild_locus = other.fchild_locus;
        }
        return *this;
    }

    //previous sibling
    inline size_t p_sibling(size_t v) const{
        size_t n_z = rs_0(v);
        if(m_shape[ss_1(n_z)-1]){
            return ss_0(n_z-1)+1;
        }else{
            return 0;
        }
    }

    inline bool is_lm_child(size_t v) const{
        size_t n_z = rs_0(v);
        return !m_shape[ss_1(n_z)-1];
    }

    //next sibling
    inline size_t n_sibling(size_t v) const{
        if(m_shape[ss_1(rs_0(v))+1]) {
            return m_shape[v] ? succ_0(v)+1 : v+1;
        }else{
            return 0;
        }
    }

    //count the number of right siblings of v
    inline size_t n_rsibs(size_t v) const{
        size_t child_pos = ss_1(rs_0(v));
        return succ_0(child_pos)-child_pos-1;
    }

    //right-most sibling
    inline size_t rm_sibling(size_t v) const{
        size_t nz = rs_0(v);
        size_t child_pos = ss_1(nz);

        if(m_shape[child_pos+1]) {
            return ss_0(nz+(succ_0(child_pos)-child_pos)-1)+1;
        }else{
            return 0;
        }
    }

    inline size_t parent(size_t v) const {
        assert(v>2);
        size_t j = ss_1( rs_0(v) );
        //this is hardcoded for the grammar!
        if(j>2 && j<fchild_locus) return 2;
        return pred_0(j)+1;
    }

    inline size_t n_children(size_t v) const{
        //size_t j = ss_0(rs_0(v)+1);
        //return j-v;
        if(v==2) return ss_0(rs_0(v)+1)-v;
        return succ_0(v)-v;
    }

    //child rank is 1-based
    inline size_t child(size_t v, size_t r) const {
        assert(v>=2);
        return ss_0(rs_1(v+r))+1;
    }

    //result is 1-based (node_map(2) = 1)
    inline size_t node_map(size_t v) const {
        return rs_0(v);
    }

    //input is 1-based (node_select(1)=2 is the root)
    inline size_t node_select(size_t v) const {
        return ss_0(v) + 1;
    }

    inline bool is_leaf(size_t v) const{
        assert(v>=2);
        return m_shape[v] == 0;
    }

    //report the rank of a leaf
    //precondition : v must be a leaf
    inline size_t leaf_rank(size_t v) const{
        assert(!m_shape[v]);
        return rs_0(v)-rs_10(v)+1;
    }

    //how many leaves are up to node with id v
    // (the rank also includes v if it is an internal node)
    inline size_t n_leaves(size_t v) const{
        size_t n_leaves = rs_0(v)-rs_10(v);
        if(v<m_shape.size() && !m_shape[v]){
            n_leaves++;
        }
        return n_leaves;
    }

    //number of internal nodes up node v
    // (the rank also includes v if it is an internal node)
    inline size_t n_ints(size_t v) const{
        return rs_10(v);
    }

    //input inclusive [0..v]!!
    //int_rank(root=2)=1
    inline size_t int_rank(size_t v) const{
        assert(m_shape[v] && !m_shape[v-1]);
        return rs_10(v);
    }

    //input 1-based: int_select(1)=root=2
    inline size_t int_select(size_t r) const{
        return pred_0(ss_10(r + 1)) + 1;
    }

    inline size_t nodes() const{
        return rs_0(m_shape.size())-1;
    }

    inline size_t int_nodes() const{
        return rs_10(m_shape.size())-1;
    }

    inline size_t leaves() const{
        return nodes()-int_nodes();
    }

    inline size_t shape_size() const{
        return m_shape.size();
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += m_shape.serialize(out, child, "m_shape");
        return written_bytes;
    }

    void load(std::istream &in){
        m_shape.load(in);
        init_tree_support();
        fchild_locus = child(2,1);
    }
};


#endif //CDT_LOUDS_TREE_HPP
