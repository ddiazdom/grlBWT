//
// Created by Diaz, Diego on 10.11.2021.
//

#include "grammar.hpp"

grammar::size_type grammar::serialize(std::ostream& out, sdsl::structure_tree_node * v, std::string name) const {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child( v, name, sdsl::util::class_name(*this));
    size_type written_bytes= sdsl::write_member(orig_size, out, child, "orig_size");
    written_bytes+= sdsl::write_member(n_seqs, out, child, "n_seqs");
    written_bytes+= sdsl::write_member(grammar_size, out, child, "grammar_size");
    written_bytes+= sdsl::write_member(sigma, out, child, "n_ter");
    written_bytes+= sdsl::write_member(n_nter, out, child, "n_nter");
    written_bytes+= sdsl::write_member(comp_string_size, out, child, "comp_string_size");
    written_bytes+= symbols_map.serialize(out, child, "symbols_map");
    written_bytes+= rules.serialize(out, child, "rules");
    written_bytes+= nter_ptr.serialize(out, child, "nter_pointers");
    written_bytes+= seq_pointers.serialize(out, child, "seq_pointers");
    written_bytes+= is_rl.serialize(out, child, "is_rl");
    return written_bytes;
}

void grammar::load(std::ifstream& in){
    sdsl::read_member(orig_size, in);
    sdsl::read_member(n_seqs, in);
    sdsl::read_member(grammar_size, in);
    sdsl::read_member(sigma, in);
    sdsl::read_member(n_nter, in);
    sdsl::read_member(comp_string_size, in);
    symbols_map.load(in);
    rules.load(in);
    nter_ptr.load(in);
    seq_pointers.load(in);
    is_rl.load(in);
}

void grammar::mark_str_boundaries() {

    std::stack<size_t> stack;
    size_t pos = nter_ptr[n_nter-1], sym, suff_sym, sym_state, seq=0;
    //0 : the symbol's rule has not been visited yet
    //1 : the symbol recursively expands to a string suffix
    //2 : the symbol does not recursively expand to a string suffix
    sdsl::int_vector<2> state(n_nter, 0);
    seq_pointers.width(sdsl::bits::hi(comp_string_size)+1);
    seq_pointers.resize(n_seqs);
    state[0] = 1;

    while(pos<rules.size()){
        sym =  rules[pos];
        if(sym>=sigma && state[sym]==0){

            suff_sym = is_rl[sym] ? rules[nter_ptr[sym]] : rules[nter_ptr[sym+1]-1];

            while(suff_sym!=sym && state[suff_sym]==0){

                stack.push(suff_sym);
                sym = suff_sym;
                suff_sym = is_rl[sym] ? rules[nter_ptr[sym]] : rules[nter_ptr[sym+1]-1];
            }

            sym_state = state[suff_sym]==1 ? 1: 2;

            while(!stack.empty()){
                suff_sym = stack.top();
                state[suff_sym] = sym_state;
                stack.pop();
            }
            state[rules[pos]] = sym_state;
            if(sym_state==1) seq_pointers[seq++] = pos;
        }else if(state[sym]==1){
            seq_pointers[seq++] = pos;
        }
        pos++;
    }
    std::cout<<seq<<" "<<n_seqs<<std::endl;
}

void gram_info_t::save_to_file(std::string& output_file){

    size_t buffer[255];

    std::ofstream of_stream(output_file, std::ofstream::binary);

    //write number of rules and alphabet size
    buffer[0] = sigma;
    buffer[1] = r;
    buffer[2] = c;
    buffer[3] = g;
    buffer[4] = max_tsym;
    of_stream.write((char *) buffer, sizeof(size_t)*5);

    assert(sym_map.size()==sigma);

    //write symbol mappings
    size_t left, right;
    for(auto const pair : sym_map){
        left = pair.first;
        right = pair.second;
        of_stream.write((char *) &left, sizeof(size_t));
        of_stream.write((char *) &right, sizeof(size_t));
    }

    buffer[0] = rules_file.size();
    of_stream.write((char *) buffer, sizeof(size_t));
    of_stream.write(rules_file.c_str(), std::streamsize(rules_file.size()));

    buffer[0] = rules_lim_file.size();
    of_stream.write((char *) buffer, sizeof(size_t));
    of_stream.write(rules_lim_file.c_str(), std::streamsize(rules_lim_file.size()));

    buffer[0] = is_rl_file.size();
    of_stream.write((char *) buffer, sizeof(size_t));
    of_stream.write(is_rl_file.c_str(), std::streamsize(is_rl_file.size()));

    buffer[0] = lvl_breaks_file.size();
    of_stream.write((char *) buffer, sizeof(size_t));
    of_stream.write(lvl_breaks_file.c_str(), std::streamsize(lvl_breaks_file.size()));

    size_t len = rules_per_level.size();
    of_stream.write((char *) &len, sizeof(size_t));
    of_stream.write((char *)rules_per_level.data(), std::streamsize(sizeof(size_t)*rules_per_level.size()));

    of_stream.close();
}

void gram_info_t::load_from_file(std::string &g_file){

    size_t buffer[255];
    std::ifstream fp(g_file, std::ifstream::binary);

    fp.read((char *)buffer, sizeof(size_t)*5);

    sigma = buffer[0];
    r = buffer[1];
    c = buffer[2];
    g = buffer[3];
    max_tsym = buffer[4];

    for(size_t i=0;i<sigma;i++){
        fp.read((char *)&buffer[0], sizeof(size_t));
        fp.read((char *)&buffer[1], sizeof(size_t));
        sym_map.insert({buffer[0], buffer[1]});
    }

    fp.read((char *)buffer, sizeof(size_t));
    auto tmp_file = reinterpret_cast<char *>(malloc(buffer[0]+1));
    fp.read(tmp_file, std::streamsize(buffer[0]));
    tmp_file[buffer[0]]='\0';
    rules_file = std::string(tmp_file);

    fp.read((char *)buffer, sizeof(size_t));
    tmp_file = reinterpret_cast<char *>(realloc(tmp_file, buffer[0]+1));
    fp.read(tmp_file, std::streamsize(buffer[0]));
    tmp_file[buffer[0]] = '\0';
    rules_lim_file = std::string(tmp_file);

    fp.read((char *)buffer, sizeof(size_t));
    tmp_file = reinterpret_cast<char *>(realloc(tmp_file, buffer[0]+1));
    fp.read(tmp_file, std::streamsize(buffer[0]));
    tmp_file[buffer[0]] = '\0';
    is_rl_file = std::string(tmp_file);

    fp.read((char *)buffer, sizeof(size_t));
    tmp_file = reinterpret_cast<char *>(realloc(tmp_file, buffer[0]+1));
    fp.read(tmp_file, std::streamsize(buffer[0]));
    tmp_file[buffer[0]] = '\0';
    lvl_breaks_file = std::string(tmp_file);

    size_t len=0;
    fp.read((char *)&len, sizeof(size_t));
    auto tmp_arr = reinterpret_cast<size_t*>(malloc(sizeof(size_t)*len));
    fp.read((char*)tmp_arr, std::streamsize(sizeof(size_t)*len));

    for(size_t i=0;i<len;i++){
        rules_per_level.push_back(tmp_arr[i]);
    }
    fp.close();
    free(tmp_file);
}
