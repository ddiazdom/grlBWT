//
// Created by inspironXV on 8/17/2018.
//

#ifndef IMPROVED_GRAMMAR_INDEX_DFUDS_TREE_H
#define IMPROVED_GRAMMAR_INDEX_DFUDS_TREE_H


#include <sdsl/int_vector.hpp>
#include <sdsl/bp_support_sada.hpp>
#include <sdsl/rrr_vector.hpp>


class dfuds_tree {

public:
    typedef size_t                                               size_type;
    typedef sdsl::bit_vector                                            bv;
    typedef sdsl::bp_support_sada<>                         parenthesis_seq;

    //protected:
    bv bit_vector;
    parenthesis_seq bps;
    sdsl::rank_support_v<00, 2> rank_00;
    sdsl::select_support_mcl<00, 2> select_00;
    bv::select_0_type select_0;


public:
    dfuds_tree() = default;

    ~dfuds_tree() = default;

    void build(const sdsl::bit_vector &v) {

        auto _bv = sdsl::bit_vector(v.size() + 3);
        _bv[0] = 1;
        _bv[1] = 1;
        _bv[2] = 0;

        for (int i = 0; i < v.size(); ++i) {
            _bv[3 + i] = v[i];
        }
        bit_vector = bv(_bv);
        compute_aux_st();


    }


    static inline short root() { return 3; }

    /*
     * return the number of nodes in the subtree of the node
     * */
    inline size_type subtree(const size_type &v) const {
        return ((bps.fwd_excess(v - 1, -1) - v) / 2 + 1);
    }

    inline bool isleaf(const size_type &i) const {
        return !bit_vector[i];
    }

    /*
     * (preoreder) return the preorder of a node
     * */
    inline size_type pre_order(const size_type &v) const {
        return v - bps.rank(v) + bit_vector[v];
    }

    /*
     * (preoreder select) return the i-th node in preorder
     * */
    inline size_type operator[](const size_type &i) const {
        return select_0(i) + 1;
    }

    inline size_type nsibling(const size_type &v) const {

        return bps.fwd_excess(v - 1, -1) + 1;
    }

    inline size_type lchild(const size_type &v) const {
        if (bit_vector[v] == false) return 0;
        return bps.find_close(v) + 1;
    }

    inline size_type fchild(const size_type &v) const {
        if (bit_vector[v] == false) return 0;
        return succ0(v) + 1;
    }


    /*
     *  return the i-th child of a node
     * */
    inline size_type child(const size_type &v, const size_type &t) const {
        if (!bit_vector[v] /*|| t > children(v)*/) return 0;

        size_t close = bps.find_close(succ0(v) - t);

        if (close == bps.size())return 0;

        return close + 1;
    }

    /*
     * return the rank of a node between its brothers
     * */
    inline size_type childrank(const size_type &v) const {
        //if(v == 3) return 0;

        size_t bp_parent = bps.find_open(v - 1);

        return succ0(bp_parent) - bp_parent;

    }

    /*
     * return the number of children of a node
     * */
    inline size_type children(const size_type &v) const {
        if (!bit_vector[v]) return 0;

        size_t zeros = v - bps.rank(v) + 1;
        return select_0(zeros + 1) - v;
    }

    /*
     * return the number of leaves in the subtree of a node
     * */
    inline size_type leafnum(const size_type &v) const {
        return leafrank(bps.fwd_excess(v - 1, -1) + 1) - leafrank(v);
    }

    /*
     * return the rank between the leaf of the left most leaf of a node
     * */
    inline size_type leafrank(const size_type &v) const {
        return rank_00(v) + 1;
    }

    inline size_type nextTreeNode(const size_type &v) const { return bps.fwd_excess(v - 1, -1) + 1; }

    inline size_type lastleaf(const size_type &v) const { return leafrank(bps.fwd_excess(v - 1, -1) + 1) - 1; }

    /*
     * return the i-th leaf
     * */
    inline size_type leafselect(const size_type &i) const { return select_00(i); }

    /*
     * return the parent of a node
     * */
    inline size_type parent(const size_type &v) const {
        if (v == 3) return 0;
        /*auto p = bps.find_open(v-1);
        while(bit_vector[p-1]!=0)--p;
        return p;*/
        return pred0(bps.find_open(v - 1)) + 1;

    }

    inline bool is_ancestor(const size_type &u, const size_type &v) const {
        return (u < v) && (v < bps.fwd_excess(u - 1, -1));
    }

    /*
     * dfs preorder over the tree
     * */
    template<typename K>
    void dfs_preorder(const size_type &node, const K &f) const {

        auto keep = f(node);
        if (!keep)
            return;
        size_type n = children(node);
        for (size_type i = 0; i < n; ++i) {
            dfs_preorder(child(node, i + 1), f);
        }
    }

    /*
     * dfs posorder over the tree
     *
     * */
    template<typename K>
    void dfs_posorder(const size_type &node, const K &f) const {

        size_type n = children(node);
        for (size_type i = 0; i < n; ++i) {
            dfs_posorder(child(node, i + 1), f);
        }

        f(node);
    }

    /*
     * this function walks a path from root to a leaf
     *
     * */
    template<typename K>
    void path(size_type &node, const K &f) const {
        size_type next_node = 0;

        if (f(node, next_node))
            path(next_node, f);
    }

    template<typename K>
    uint find_child_dbs(const uint &node, uint &ls, uint &hs, const K &f) const {

        uint p2 = 1;
        /*
         * f return true if p2 < value else return false             *
         * */
        while (f(p2)) {

            if (p2 * 2 > hs) {
                p2 = hs;
                break;
            }
            p2 *= 2;
        }

        ls = (p2 == 1) ? 1 : p2 / 2;
        hs = p2;

        while (ls + 1 < hs) {
            uint m = (ls + hs) / 2;
            f(m) ? (ls = m) : (hs = m - 1);
        }

        if (f(hs) == true)
            return hs;
        return ls;

    }

    template<typename K>
    uint find_child_dbs_mirror(const uint &node, uint &ls, uint &hs, const K &f) const {

        uint p2 = hs;
        /*
         * f return true if p2 < value else return false             *
         * */

        while (p2 && f(p2)) {
            p2 /= 2;
        }


        ls = p2 ? p2 : 1;
        if (hs / 2 != ls)
            hs = (p2 * 2 < hs) ? p2 * 2 : hs;


        while (ls + 1 < hs) {
            uint m = (ls + hs) / 2;
            f(m) ? (hs = m) : (ls = m + 1);
        }

        if (f(ls) == false)
            return hs;
        return ls;

    }


    template<typename K>
    size_type
    find_child(const size_type &node, size_type &ls, size_type &hs, const K &f) const {
        while (ls + 1 < hs) {
            size_type rank_ch = (ls + hs) / 2;
            f(rank_ch) ? (hs = rank_ch - 1) : (ls = rank_ch);
        }
        if (f(hs) == true)
            return ls;
        return hs;
    }


    void load(std::istream &in) {

        sdsl::load(bit_vector, in);
        bps.load(in, &bit_vector);
        sdsl::load(rank_00, in);
        sdsl::load(select_00, in);
        sdsl::load(select_0, in);

        compute_aux_st();


    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += sdsl::serialize(bit_vector, out);
        written_bytes += sdsl::serialize(bps, out);
        written_bytes += sdsl::serialize(rank_00, out);
        written_bytes += sdsl::serialize(select_00, out);
        written_bytes += sdsl::serialize(select_0, out);


        return written_bytes;
    }


    inline size_type rank_1(const size_type &i) const {
        return bps.rank(i);
    }

    void print() const {
        for (unsigned long i : bit_vector) {
            std::cout << i;
        }
        std::cout << std::endl;
        for (unsigned long i : bit_vector) {
            if (i)
                std::cout << '(';
            else
                std::cout << ')';
        }
        std::cout << "\n";

    }

    dfuds_tree &operator=(const dfuds_tree &T) {
        bit_vector = T.bit_vector;
        compute_aux_st();
        return *this;
    }

public:
    inline size_t pred0(const size_type &i) const {
        if (bit_vector[i] == 0)
            return i;
        return select_0(i - bps.rank(i) + 1);
    }

    inline size_t succ0(const size_type &i) const {
        return select_0(i - bps.rank(i) + 1 + bit_vector[i]);
    }

protected:
    void compute_aux_st(){
        bps = sdsl::bp_support_sada<>(&bit_vector);
        rank_00 = sdsl::rank_support_v<00, 2>(&bit_vector);
        select_00 = sdsl::select_support_mcl<00, 2>(&bit_vector);
        select_0 = bv::select_0_type(&bit_vector);
    }
};






#endif //IMPROVED_GRAMMAR_INDEX_DFUDS_TREE_H
