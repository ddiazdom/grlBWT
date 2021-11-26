//
// Created by Diaz, Diego on 10.11.2021.
//

#include "grammar.hpp"

grammar::size_type grammar::serialize(std::ostream& out, sdsl::structure_tree_node * v, std::string name) const {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child( v, name, sdsl::util::class_name(*this));
    size_type written_bytes= sdsl::write_member(text_size, out, child, "text_size");
    written_bytes+= sdsl::write_member(n_strings, out, child, "n_strings");
    written_bytes+= sdsl::write_member(grammar_size, out, child, "grammar_size");
    written_bytes+= sdsl::write_member(text_alph, out, child, "n_ter");
    written_bytes+= sdsl::write_member(gram_alph, out, child, "gram_alph");
    written_bytes+= sdsl::write_member(comp_string_size, out, child, "comp_string_size");
    written_bytes+= sdsl::write_member(n_p_rounds, out, child, "n_p_rounds");
    written_bytes+= rules_breaks.serialize(out, child, "rules_breaks");
    written_bytes+= symbols_map.serialize(out, child, "symbols_map");
    written_bytes+= rules.serialize(out, child, "rules");
    written_bytes+= nter_ptr.serialize(out, child, "nter_pointers");
    written_bytes+= seq_pointers.serialize(out, child, "seq_pointers");
    return written_bytes;
}

void grammar::load(std::ifstream& in){
    sdsl::read_member(text_size, in);
    sdsl::read_member(n_strings, in);
    sdsl::read_member(grammar_size, in);
    sdsl::read_member(text_alph, in);
    sdsl::read_member(gram_alph, in);
    sdsl::read_member(comp_string_size, in);
    sdsl::read_member(n_p_rounds, in);
    rules_breaks.load(in);
    symbols_map.load(in);
    rules.load(in);
    nter_ptr.load(in);
    seq_pointers.load(in);
}

void grammar::mark_str_boundaries() {

    std::stack<size_t> stack;
    size_t pos = nter_ptr[gram_alph - 1], sym, suff_sym, sym_state, seq=0;
    //0 : the symbol's rule has not been visited yet
    //1 : the symbol recursively expands to a string suffix
    //2 : the symbol does not recursively expand to a string suffix
    sdsl::int_vector<2> state(gram_alph, 0);
    seq_pointers.width(sdsl::bits::hi(comp_string_size)+1);
    seq_pointers.resize(n_strings);
    state[0] = 1;
    size_t c_start = nter_ptr[nter_ptr.size()-1];

    while(pos<rules.size()){
        sym =  rules[pos];
        if(sym >= text_alph && state[sym] == 0){

            suff_sym = is_rl(sym) ? rules[nter_ptr[sym]] : rules[nter_ptr[sym+1]-1];

            while(suff_sym!=sym && state[suff_sym]==0){
                stack.push(suff_sym);
                sym = suff_sym;
                suff_sym = is_rl(sym) ? rules[nter_ptr[sym]] : rules[nter_ptr[sym+1]-1];
            }

            sym_state = state[suff_sym]==1 ? 1: 2;

            while(!stack.empty()){
                suff_sym = stack.top();
                state[suff_sym] = sym_state;
                stack.pop();
            }
            state[rules[pos]] = sym_state;
            if(sym_state==1) seq_pointers[seq++] = pos-c_start;
        }else if(state[sym]==1){
            seq_pointers[seq++] = pos-c_start;
        }
        pos++;
    }
}

std::string grammar::decomp_str(size_t idx) {
    assert(idx<seq_pointers.size());
    std::string exp;
    size_t c_start = nter_ptr[nter_ptr.size()-1];
    size_t str_start, str_end;
    str_start = seq_pointers[idx]+c_start;
    str_end = seq_pointers[idx+1]+c_start;
    for(size_t j=str_start;j<str_end;j++){
        buff_decomp_nt(rules[j], exp);
    }
    return exp;
}

//buffered decompression
template<class vector_t>
void grammar::buff_decomp_nt_int(size_t sym, vector_t& exp, hash_table_t& ht) {

    size_t pos1, pos2, len, val, freq, start, end;
    if(sym<text_alph){
        exp.push_back((char)symbols_map[sym]);
    }else{
        pos1 = exp.size();
        auto res = ht.find(&sym, 64);
        if(res.second){
            ht.get_value_from(*res.first, val);
            pos2 = val >>32UL;
            len = val & ((1UL<<32UL)-1UL);
            //TODO use memcpy
            for(size_t j=pos2;j<pos2+len;j++){
                exp.push_back(exp[j]);
            }
            val = (pos1<<32UL) |  len;
            ht.insert_value_at(*res.first, val);
        }else{
            auto res2 = ht.insert(&sym, 64, 0);
            start = nter_ptr[sym];
            end = nter_ptr[sym+1]-1;
            if(is_rl(sym)){
                sym = rules[start];
                freq = rules[end];
                buff_decomp_nt_int<vector_t>(sym, exp, ht);
                len = exp.size()-pos1;
                for(size_t i=1;i<freq;i++){
                    //TODO use memcpy
                    for(size_t j=pos1;j<(pos1+len);j++){
                        exp.push_back(exp[j]);
                    }
                }
            }else{
                for(size_t i=start;i<=end;i++){
                    buff_decomp_nt_int<vector_t>(rules[i], exp, ht);
                }
            }
            len = exp.length() - pos1;
            val = (pos1<<32UL) |  len;
            ht.insert_value_at(*res2.first, val);
        }
    }
}

template<class vector_t>
void grammar::buff_decomp_nt(size_t nt, vector_t& exp, size_t ht_buff_size){
    void * buff = malloc(ht_buff_size);
    hash_table_t ht(ht_buff_size, "", 0.8, buff);
    buff_decomp_nt_int<vector_t>(nt, exp, ht);
    free(buff);
}

std::string grammar::decomp_nt(size_t nt) {
    std::string exp;
    assert(nt<gram_alph);
    size_t start, end, sym, rl_sym, freq;
    std::stack<size_t> stack;
    stack.push(nt);

    while(!stack.empty()){
        sym = stack.top();
        stack.pop();

        if(sym<text_alph){
            exp.push_back((char)symbols_map[sym]);
        }else{
            start = nter_ptr[sym];
            if(is_rl(sym)){
                rl_sym = rules[start];
                freq = rules[start+1];
                assert(freq>1);
                std::string res = decomp_nt(rl_sym);
                for(size_t i=0;i<freq;i++){
                    exp+=res;
                }
            }else{
                end = nter_ptr[sym+1];
                for(size_t j=end;j-->start;){
                    stack.push(rules[j]);
                }
            }
        }
    }
    return exp;
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

    size_t len = rules_breaks.size();
    of_stream.write((char *) &len, sizeof(size_t));
    of_stream.write((char *)rules_breaks.data(), std::streamsize(sizeof(size_t) * rules_breaks.size()));

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

    size_t len=0;
    fp.read((char *)&len, sizeof(size_t));
    auto tmp_arr = reinterpret_cast<size_t*>(malloc(sizeof(size_t)*len));
    fp.read((char*)tmp_arr, std::streamsize(sizeof(size_t)*len));

    for(size_t i=0;i<len;i++){
        rules_breaks.push_back(tmp_arr[i]);
    }
    fp.close();
    free(tmp_file);
}

template void grammar::buff_decomp_nt_int<std::string>(size_t nt, std::string& exp, grammar::hash_table_t& ht);
template void grammar::buff_decomp_nt<std::string>(size_t nt, std::string& exp, size_t ht_buff_size);
