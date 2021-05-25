//
// Created by diego on 5/5/21.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_TREE_H
#define LPG_COMPRESSOR_GRAMMAR_TREE_H


#include <sdsl/wt_gmr.hpp>
#include <sdsl/inv_perm_support.hpp>
#include <sdsl/construct.hpp>
#include "utils.hpp"
#include "dfuds_tree.hpp"
#include "macros.hpp"

class grammar_tree_t{

public:
    typedef size_t                                          size_type;
    typedef lpg_build::plain_grammar_t                      plain_grammar;
    typedef dfuds_tree                                      top_tree;

    typedef sdsl::wt_gmr<
                sdsl::int_vector<>,
                sdsl::inv_multi_perm_support<INV_PI_WX>
            >                                                wt_x;

    typedef sdsl::sd_vector<>                                bv_z;
    typedef sdsl::sd_vector<>                                bv_l;
    typedef sdsl::sd_vector<>                                bv_r;
    typedef sdsl::int_vector<>                                 vi;
    typedef sdsl::inv_perm_support<INV_PI_X>               inv_vi;



protected:
    /**
     * @T is the topology of the grammar
     * */
    top_tree                        T;
    /**
     * @Z bitmap with the first position of symbols in preorder sequence
     * */
    bv_z                            Z;
    bv_z::rank_1_type               rank1_Z;
    bv_z::select_1_type             select1_Z;
    bv_z::select_0_type             select0_Z;
    /**
     * @X sequence of preorder grammar tree symbols removing first mention
     * */
    wt_x                            X;
    /**
     *  @F is a permutation of Xj such that F[i] = Xj only and only if Xj is the ith diferent symbol in preorder
     *  @F_inv inverse permutation of F
     * */
    vi                              F;
    inv_vi                          F_inv;
    /**
     *  @L is a bitvector of length |T| that the start position on the text of the phrases created by the grammar tree
     *  (or the start posititon on the text of the leaves of the grammar tree )
     * */
    bv_l                            L; // marks the init position of each Xi in T
    bv_l::select_1_type             select_L;
    bv_l::rank_1_type               rank_L;

    bv_r                            R; // mark the run lenght rules. g bits
    bv_r::rank_1_type          rank_r;
    vi                              RL;// store the length of the run lenght rules. |rank(R)|


public:


    grammar_tree_t() = default;
    grammar_tree_t(const grammar_tree_t& _g):T(_g.T),Z(_g.Z),X(_g.X),F(_g.F),L(_g.L){compute_aux_st();}
    virtual ~grammar_tree_t() = default;

    void load(std::istream &in){

        sdsl::load(T,in);

        sdsl::load(Z,in);
        sdsl::load(rank1_Z,in);
        sdsl::load(select1_Z,in);
        sdsl::load(select0_Z,in);

        sdsl::load(X,in);

        sdsl::load(F,in);
        sdsl::load(F_inv,in);


        sdsl::load(L,in);
        sdsl::load(select_L,in);

        sdsl::load(R,in);
        sdsl::load(rank_r,in);
        sdsl::load(RL,in);

        compute_aux_st();
    }
    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += sdsl::serialize(T,out);
        written_bytes += sdsl::serialize(Z,out);
        written_bytes += sdsl::serialize(rank1_Z,out);
        written_bytes += sdsl::serialize(select1_Z,out);
        written_bytes += sdsl::serialize(select0_Z,out);
        written_bytes += sdsl::serialize(X,out);
        written_bytes += sdsl::serialize(F,out);
        written_bytes += sdsl::serialize(F_inv,out);
        written_bytes += sdsl::serialize(L,out);
        written_bytes += sdsl::serialize(select_L,out);
        written_bytes += sdsl::serialize(R,out);
        written_bytes += sdsl::serialize(rank_r,out);
        written_bytes += sdsl::serialize(RL,out);
        return written_bytes;
    }


    void build(utils::nav_grammar& NG, const plain_grammar& Gr, const size_t& text_length,utils::lenght_rules& rules_off, const size_type& init_rule){
#ifdef DEBUG_INFO
        std::cout<<"build_tree(Gr,NG,text_length);\n";
#endif
        build_tree(Gr,NG,text_length,rules_off,init_rule);
#ifdef DEBUG_INFO
        std::cout<<"compute_aux_st();\n";
#endif
        compute_aux_st();

    }
    const top_tree& getT()const {return T;}



    size_type first_occ_preorder_node(const size_type& preorder)const {
        // preorder has to be a non-terminal leaf
        if(!Z[preorder - 1]){
            //check if it is not a first mention
            size_type _x = X[preorder - rank1_Z(preorder) - 1 ];
            return select1_Z(F[_x]) + 1;
        }
    }
    size_type first_occ_from_rule(const size_type& _x)const {
        // preorder has to be a non-terminal leaf
//        std::cout<<F[_x]<<std::endl;
//        std::cout<<select1_Z(F[_x]) + 1 <<std::endl;
//        utils::pretty_printer_bv(Z,"Z");
        return select1_Z(F[_x]) + 1;
    }
    inline bool isLeaf(const size_type& preorder)const {
        size_type node = T[preorder];
        return T.isleaf(node);
    }
    inline size_type get_rule_from_preorder_node(const size_type& preorder) const {

        if(!Z[preorder - 1]){ //second mention case
            return X[preorder - rank1_Z(preorder) - 1 ];
        }
        //first mention case
        return F_inv[rank1_Z(preorder)];
    }
    // return 0 if is not a run
    // in other case return the length of the run
    inline size_type is_run(const size_type& preorder)const{
        if(!R[preorder - 1]) return 0;
        return RL[rank_r(preorder - 1)];
    }
    inline size_type get_size_rules()const {return F.size();}
    inline size_type get_grammar_size()const {return Z.size();}
    inline size_type offset_node(const size_type& node) const {
        size_type leaf = T.leafrank(node);
        return select_L(leaf);
    }
    inline size_type offset_node(const size_type& node, size_type& leaf) const {
        leaf = T.leafrank(node);
        return select_L(leaf);
    }
    template<typename F>
    void visit_secondary_occ(const size_type& preorder_node, const F& f) const {
        size_type  x = get_rule_from_preorder_node(preorder_node);
        size_type num_occ = X.rank(X.size(),x);
        for (size_type i = 0; i < num_occ ; ++i) {
            size_type  j = X.select(i+1,x);
            size_type preorder = select0_Z( j + 1 ) + 1;
            f(preorder);
        }
    }
    void breakdown_space() const {

        std::cout<<"T,>"<<sdsl::size_in_bytes(T)<<std::endl;
        std::cout<<"Z,>"<<sdsl::size_in_bytes(Z)<<std::endl;
        std::cout<<"X,>"<<sdsl::size_in_bytes(X)<<std::endl;
        std::cout<<"F,>"<<sdsl::size_in_bytes(F)<<std::endl;
        std::cout<<"L,>"<<sdsl::size_in_bytes(L)<<std::endl;

    }
    inline size_type get_text_len()const { return L.size();}
    inline size_type selectL(const size_type &i)const{ return select_L(i);}
    inline size_type num_leaves()const{return rank_L(L.size());}

protected:



    /**
     * build algorithms
     * */

    void build_tree( const plain_grammar& Gr,utils::nav_grammar& grammar, const size_t& text_length,utils::lenght_rules& rules_off
                     ,const size_type& init_rule){

        size_type nnodes = (Gr.g - Gr.sigma);
#ifdef DEBUG_INFO
        std::cout<<"Building dfuds representation of parser tree S:"<<init_rule<<std::endl;
        std::cout<<"Number of nodes:"<<nnodes<<std::endl;
#endif
        sdsl::bit_vector _bv(2 * nnodes  - 1 , 1);
        size_type pos = 0; // offset in tree topology bitvector
        std::set<size_type> M; //
        sdsl::bit_vector _z(nnodes, 0); // mark if the node i in preorder is a first mention node

        sdsl::int_vector_buffer<1> is_rules_len(Gr.is_rl_file);
        size_type n_l = 0; // count the number of runs
        for (size_type i = 0; i < is_rules_len.size(); ++i) {
            if(is_rules_len[i]){
                auto it = grammar[i];
                if(it.size() == 2 && it[1] > 2) //(X-> Xi^l === Xi,l ) check if l > 2 else is not consider as a run length node...
                    n_l++;
            }
        }
        sdsl::bit_vector _r(nnodes, 0); // mark if the node i is a run length
        sdsl::int_vector<> _rl(n_l, 0); // mark if the node i is a run length

        size_t z_pos = 0,z_rank = 0, x_pos = 0,_rl_pos = 0;
        sdsl::int_vector<> _f(Gr.r, 0); // permutation of rules
        sdsl::int_vector<> _x(nnodes - Gr.r + 1); // label of second mention nodes in the tree
        sdsl::bit_vector _l(text_length, 0); // text-len bitvector marking start position of parser phrases..
        size_type l_pos = 0; // current off in the text.

#ifdef DEBUG_INFO
    std::cout<<"Start dfs\n";
#endif


        utils::dfs_2v(init_rule,grammar,is_rules_len,[
                &Gr,&is_rules_len,&grammar,&_bv,&pos,&_z,&z_pos
                ,&z_rank,&_f,&_x,&x_pos,&_l, &l_pos
                ,&_r,&_rl,&_rl_pos,&M,&rules_off
                ](const size_type& id, bool first_visit){ //len is > 0  if it is the second child of a run length
            if(first_visit){ // first visit

                z_pos ++ ; // pre-incress offset in vector z;
                _l[l_pos] = true; // mark text position for first time visit

                if (M.find(id) != M.end()) { // second mention of a non-terminal => leaf tree

                    _bv[pos] = false; // add to the tree as leaf
                    ++pos; // increase off in tree

                    _x[x_pos] = id;
                    x_pos++; // store second mention node

                    l_pos += rules_off[id].second; // increase the off-text with pre-computed len-rule
                    return false; // return false stop descending in the tree
                }

                M.insert(id); // store first mention of a rule tree
                _z[ z_pos - 1 ] = true; //mark the node as first mention (we index in z_pos -1 because we pre increase the off)
                _f[id] = ++z_rank; // map the rule with the number of 1s in z;

                if(is_rules_len[id] && grammar[id][1] > 2){
//                     check if it is a run length rule
                    _r[ z_pos - 1 ] = true;
                    _rl[_rl_pos++] = grammar[id][1]; //store the exponent

                    size_type rhl = 2; // right hand length
                    _bv[pos + rhl] = false; // add 0 before rhl children
                    pos += rhl + 1; // increase off in tree

                    rules_off[id] = std::make_pair(l_pos,1); // store the off in the text
                    return true; //  keep descending in the tree

                }else{

                    if(Gr.isTerminal(id)){
                        //if is a terminal symbol
                        _bv[pos] = false; // add to the tree as leaf
                        ++pos; // increase off in tree

                        rules_off[id] = std::make_pair(l_pos,1); // store the off in the text
                        ++l_pos;
                        return false; // stop descending in the tree
                    }

                    size_type rhl = grammar[id].size(); // right hand length
                    _bv[pos + rhl] = false; // add 0 before rhl children
                    pos += rhl + 1; // increase off in tree


                    rules_off[id] = std::make_pair(l_pos,1); // store the off in the text
                    return true; //  keep descending in the tree

                }
            }
            //second visit (only in non-terminal)


            if(is_rules_len[id] && grammar[id][1] > 2){
                //case run-length
                size_type _len = grammar[id][1]; //run len
                size_type one_child_len = rules_off[grammar[id][0]].second; //one child len

                rules_off[id].second = _len * one_child_len; // total len
                l_pos += (_len - 2) * one_child_len;
                return false;
            }

            rules_off[id].second = l_pos - rules_off[id].first; // update rule length with the offset diferences
            return false; // to avoid warning.... :)
        });


        //check residual sizes
#ifdef DEBUG_INFO
        std::cout<<"check residual sizes\n";
        std::cout<<"bv:"<<_bv.size()<<"-"<<pos<<std::endl;
        std::cout<<"f:"<<_f.size()<<"-"<<_f.size()<<std::endl;
        std::cout<<"x:"<<_x.size()<<"-"<<x_pos<<std::endl;
        std::cout<<"z:"<<_z.size()<<"-"<<z_pos<<std::endl;
        std::cout<<"l:"<<_l.size()<<"-"<<l_pos<<std::endl<<" root-off:"<<rules_off[init_rule].first<<","<<rules_off[init_rule].second<<std::endl;
        std::cout<<"r:"<<_r.size()<<"-"<<z_pos<<std::endl;
        std::cout<<"rl:"<<_rl.size()<<"-"<<_rl_pos<<std::endl;
#endif

        T.build(_bv);
        sdsl::util::bit_compress(_x);
        std::ofstream x_file("x_file", std::ios::binary);
        sdsl::serialize(_x, x_file);
        x_file.close();
        sdsl::construct(X, "x_file", 0);
        F = vi (_f);
        Z = bv_z(_z);
        L = bv_l (_l);
        R = bv_r(_r);
        RL = vi(_rl);

#ifdef DEBUG_PRINT
        T.print();
        utils::pretty_printer_bv(Z,"Z");
        utils::pretty_printer_v(F,"F");
        utils::pretty_printer_v(X,"X");
        utils::pretty_printer_bv(L,"L");
        utils::pretty_printer_bv(R,"R");
        utils::pretty_printer_bv(RL,"RL");
#endif
//        compute_aux_st();

    }
    /**
     * Initialize rank and select structures of all components and inverse permutation....
     */
    void compute_aux_st(){

        rank1_Z         = bv_z ::rank_1_type (&Z);
        select1_Z       = bv_z ::select_1_type (&Z);
        select0_Z       = bv_z ::select_0_type (&Z);

        F_inv           = inv_vi(&F);
        sdsl::util::bit_compress(F);

        select_L        = bv_l::select_1_type(&L);
        rank_L          = bv_l::rank_1_type(&L);

        rank_r          = bv_r::rank_1_type (&R);
        sdsl::util::bit_compress(RL);



    }


};
#endif //LPG_COMPRESSOR_GRAMMAR_TREE_H

