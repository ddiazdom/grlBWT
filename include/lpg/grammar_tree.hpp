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

#define INV_PI_WX 8
#define INV_PI_X 8

class grammar_tree_t{

public:
    typedef size_t                                          size_type;
    typedef lpg_build::plain_grammar_t                      plain_grammar;
    typedef dfuds_tree                                      top_tree;

    typedef sdsl::wt_gmr<
                sdsl::int_vector<>,
                sdsl::inv_multi_perm_support<INV_PI_WX>
            >                                               wt_x;

    typedef sdsl::sd_vector<>                                bv_z;
    typedef sdsl::sd_vector<>                                bv_l;
    typedef sdsl::sd_vector<>                                bv_y;
    typedef sdsl::int_vector<>                               vi;
    typedef sdsl::inv_perm_support<INV_PI_X>           inv_vi;
    typedef std::vector<uint8_t>                             vi_alp;




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
     *  @Y is a bitvector of length g that marks which rules are terminal ( Xi -> a)
     * */
    bv_y                            Y; // todo esto quizas no hace falta las reglas son X.id < sigma... problema para buscar patrones de
    bv_y::rank_1_type               rank_Y;
    bv_y::select_1_type             select_Y;
    /**
     *  @L is a bitvector of length |T| that the start position on the text of the phrases created by the grammar tree
     *  (or the start posititon on the text of the leaves of the grammar tree )
     * */
    bv_l                            L; // marks the init position of each Xi in T
    bv_l::rank_1_type               rank_L;
    bv_l::select_1_type             select_L;
    /**
     * @alp is a array with the alphabeth symbols.
     * */
    vi_alp alp;

public:


    grammar_tree_t() = default;
    grammar_tree_t(const grammar_tree_t& _g):T(_g.T),Z(_g.Z),X(_g.X),F(_g.F),L(_g.L),alp(_g.alp){
          compute_aux_st();
    }

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

        sdsl::load(Y,in);
        sdsl::load(rank_Y,in);
        sdsl::load(select_Y,in);

        sdsl::load(L,in);
        sdsl::load(rank_L,in);
        sdsl::load(select_L,in);

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

        written_bytes += sdsl::serialize(Y,out);
        written_bytes += sdsl::serialize(rank_Y,out);
        written_bytes += sdsl::serialize(select_Y,out);


        written_bytes += sdsl::serialize(L,out);
        written_bytes += sdsl::serialize(rank_L,out);
        written_bytes += sdsl::serialize(select_L,out);


        return written_bytes;
    }

    void build(utils::nav_grammar& NG, const plain_grammar& Gr, const size_t& text_length,utils::lenght_rules& rules_off){
        std::cout<<"utils::nav_grammar NG = utils::build_nav_grammar(Gr);\n";
        build_tree(Gr,NG,text_length,rules_off);
        std::cout<<"build_tree(Gr,NG,text_length);\n";
        compute_aux_st();
        std::cout<<"compute_aux_st();\n";
    }

    const top_tree& getT(){return T;}

protected:



    /**
     * build algorithms
     * */

    void build_tree( const plain_grammar& Gr,utils::nav_grammar& grammar, const size_t& text_length,utils::lenght_rules& rules_off ){

        std::cout<<"Building dfuds representation of parser tree"<<std::endl;
        sdsl::bit_vector _bv(2 * Gr.g - 1, 1);
        size_type pos = 0; // offset in tree topology bitvector
        std::set<size_type> M; //
        sdsl::bit_vector _z(Gr.g, 0); // mark if the node i in preorder is a first mention node
        size_t z_pos = 0,z_rank = 0, x_pos = 0;
        sdsl::int_vector<> _f(Gr.r + 1, 0); // permutation of rules
        sdsl::int_vector<> _x(Gr.g - Gr.r); // label of second mention nodes in the tree
        sdsl::bit_vector _l(text_length + 1, 0); // text-len bitvector marking start position of parser phrases..
        size_type l_pos = 0; // current off in the text.


        utils::dfs_2v(Gr.sigma,grammar,[&Gr,&grammar,&_bv,&pos,&_z,&z_pos,&z_rank,&_f,&_x,&x_pos,&_l, &l_pos, &M,&rules_off](const size_type& id, bool first_visit){
            if(first_visit){ // first visit
                z_pos ++ ; // pre-incress offset in vector z;
                if (M.find(id) != M.end()) { // second mention of a non-terminal => leaf tree

                    _bv[pos] = false; // add to the tree as leaf
                    ++pos; // increase off in tree

                    _x[x_pos] = id;x_pos++; // store second mention node

                    _l[l_pos] = true; // mark text position for this leaf tree
                    l_pos += rules_off[id].second; // increase the off-text with pre-computed len-rule

                    return false; // return false stop descending in the tree
                }

                M.insert(id); // store first mention of a rule tree

                _z[ z_pos - 1 ] = true; //mark the node as first mention (we index in z_pos -1 because we pre increase the off)
                _f[id] = ++z_rank; // map the rule with the number of 1s in z;

                if(id < Gr.sigma){
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
            //second visit (only in non-terminal)
            rules_off[id].second = l_pos - rules_off[id].first; // update rule length with the offset diferences
            return false; // to avoid warning.... :)
        });


        T.build(_bv);
        sdsl::util::bit_compress(_x);
        std::ofstream x_file("x_file", std::ios::binary);
        sdsl::serialize(_x, x_file);
        x_file.close();
        sdsl::construct(X, "x_file", 0);
        F = vi (_f);
        sdsl::util::bit_compress(F);
        Z = bv_z(_z);
        L = bv_l (_l);

//        T.print();
//        utils::pretty_printer_bv(Z,"Z");
//        utils::pretty_printer_v(F,"F");
//        utils::pretty_printer_v(X,"X");
//        utils::pretty_printer_bv(L,"L");

    }

//    void build_Y(const plain_grammar& _Gr){
//        sdsl::bit_vector _y(_Gr.g + 1,0);
//        Y = bv_y (_y);
//    }
//    size_type build_L(const plain_grammar& _Gr,utils::nav_grammar& grammar, const size_t& text_length){
//
//        utils::Roff rules_len;
//        sdsl::bit_vector _l;
//        size_type n_nodes = utils::compute_rule_offsets(_Gr,grammar,text_length,rules_len,_l);
//        L = bv_l (_l);
//        return n_nodes;
//    }

    /**
     * Initialize rank and select structures of all components and inverse permutation....
     */
    void compute_aux_st(){

        rank1_Z         = bv_z ::rank_1_type (&Z);
        select1_Z       = bv_z ::select_1_type (&Z);
        select0_Z       = bv_z ::select_0_type (&Z);

//        F_inv           = inv_vi(&F);

        rank_Y          = bv_y::rank_1_type(&Y);
        select_Y        = bv_y::select_1_type(&Y);

        select_L        = bv_l::select_1_type(&L);
        rank_L          = bv_l::rank_1_type(&L);

    }


};
#endif //LPG_COMPRESSOR_GRAMMAR_TREE_H

