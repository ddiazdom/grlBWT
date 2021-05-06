//
// Created by diego on 5/5/21.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_TREE_H
#define LPG_COMPRESSOR_GRAMMAR_TREE_H
class grammar_tree_t{

public:

    typedef size_t                        size_type;
    void load(std::istream &in){
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        return written_bytes;
    }
};
#endif //LPG_COMPRESSOR_GRAMMAR_TREE_H
