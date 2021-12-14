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

void grammar::load(std::istream& in){
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
    size_t c_start = pos;

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
    size_t c_start = nter_ptr[gram_alph-1];
    size_t str_start, str_end;

    str_start = idx == 0 ? c_start : c_start+seq_pointers[idx-1]+1;
    str_end = seq_pointers[idx]+c_start;

    for(size_t j=str_start;j<=str_end;j++){
        buff_decomp_nt(rules[j], exp);
    }
    return exp;
}

template<>
bool grammar::copy_to_front<std::string>(std::string& stream, size_t src, size_t len, size_t freq){

    size_t n_syms, rem_syms=len*freq, dst=stream.size();
    stream.resize(dst+(len*freq));
    char * data = stream.data();
    rem_syms-=len;

    memcpy(&data[dst], &data[src], len);
    for(size_t i=1;rem_syms>0;i*=2){
        n_syms=std::min<size_t>(len*i, rem_syms);
        memcpy(&data[dst+(len*i)], &data[dst], n_syms);
        rem_syms-=n_syms;
    }
    return true;
}

template<>
bool grammar::copy_to_front<o_file_stream<char>>(o_file_stream<char>& stream, size_t src, size_t len, size_t freq){
    return stream.copy_to_front((buff_s_type) src, (buff_s_type) len, (buff_s_type) freq);
}

template<class vector_t>
bool grammar::copy_to_front(vector_t &stream, size_t src, size_t len, size_t freq) {
    size_t n_syms, rem_syms=len*freq, dst=stream.size();
    stream.resize(dst+(len*freq));
    char * data = stream.data();
    rem_syms-=len;

    memcpy(&data[dst], &data[src], len);
    for(size_t i=1;rem_syms>0;i*=2){
        n_syms=std::min<size_t>(len*i, rem_syms);
        memcpy(&data[dst+(len*i)], &data[dst], n_syms);
        rem_syms-=n_syms;
    }
    return true;
}

//buffered decompression
template<class vector_t>
void grammar::buff_it_decomp_nt_int(size_t root, vector_t& exp, hash_table_t& ht) {

    std::stack<node_t> st;

    size_t len=0, freq, exp_size;
    locus_t locus = {0,0};
    node_t node = {root, 0, 0, false};

    bool copied;
    size_t ht_addr;

    do {
        while(true) {

            st.push(node);

            if(node.sym<text_alph){
                break;
            } else {
                //check if the nonterminal was already decompressed
                auto res = ht.find(&node.sym, sdsl::bits::hi(node.sym)+1);
                copied = false;
                if(res.second){
                    locus = {0,0};
                    ht.get_value_from(res.first, locus);
                    if(copy_to_front(exp, locus.src, locus.exp_len, 1)){
                        len += locus.exp_len;
                        copied = true;
                        ht_addr = res.first;
                        break;
                    }
                }

                node.g_pos = nter_ptr[node.sym];
                node.is_rm_child = (is_rl(node.sym) || node.g_pos==(nter_ptr[node.sym+1]-1));
                node.sym = rules[node.g_pos];
            }
        }

        auto temp = st.top();
        st.pop();

        if(temp.sym<text_alph){
            exp.push_back((char)symbols_map[temp.sym]);
            len++;
        } else {
            locus.exp_len = len - temp.l_exp;
            locus.src = exp.size()-locus.exp_len;

            if(copied){
                ht.insert_value_at(ht_addr, locus);
            }else{
                auto res = ht.insert(&temp.sym, sdsl::bits::hi(temp.sym)+1, locus);
                if(!res.second) ht.insert_value_at(res.first, locus);
            }
        }

        while(!st.empty() && temp.is_rm_child) {

            temp = st.top();
            st.pop();

            if(is_rl(temp.sym)){
                exp_size = len - temp.l_exp;
                freq = rules[nter_ptr[temp.sym+1]-1];
                copy_to_front(exp, exp.size()-exp_size, exp_size, freq-1);
                len+=exp_size*(freq-1);
            }

            locus.exp_len = len - temp.l_exp;
            locus.src = exp.size() - locus.exp_len;

            auto res = ht.insert(&temp.sym, sdsl::bits::hi(temp.sym)+1, locus);
            if(!res.second) ht.insert_value_at(res.first, locus);
        }

        if (!st.empty()) {
            assert(!temp.is_rm_child);
            node.g_pos = temp.g_pos+1;
            node.is_rm_child = node.g_pos == (nter_ptr[st.top().sym+1]-1);
            node.sym = rules[node.g_pos];
            node.l_exp = len;
        }
    }while(!st.empty());
}

//buffered decompression
template<class vector_t>
void grammar::buff_decomp_nt_int(size_t sym, vector_t& exp, hash_table_t& ht) {

    size_t pos1, start, end;
    if(sym<text_alph){
        exp.push_back((char)symbols_map[sym]);
    }else{
        pos1 = exp.size();
        auto res = ht.find(&sym, sdsl::bits::hi(sym)+1);
        if(res.second){
            locus_t locus = {0,0};
            ht.get_value_from(res.first, locus);

            //if the string couldn't be copied, then we do it the hard way
            if(!copy_to_front(exp, locus.src, locus.exp_len, 1)){
                start = nter_ptr[sym];
                end = nter_ptr[sym+1]-1;
                if(is_rl(sym)){
                    buff_decomp_nt_int<vector_t>(rules[start], exp, ht);
                    copy_to_front(exp, pos1, exp.size() - pos1, rules[end] - 1);
                }else{
                    for(size_t i=start;i<=end;i++){
                        buff_decomp_nt_int<vector_t>(rules[i], exp, ht);
                    }
                }
                locus.src = pos1;
                auto res2 =  ht.insert(&sym, sdsl::bits::hi(sym)+1, locus);
                if(!res2.second) ht.insert_value_at(res2.first, locus);
            }else{
                locus.src = pos1;
                ht.insert_value_at(res.first, locus);
            }
        }else{
            start = nter_ptr[sym];
            end = nter_ptr[sym+1]-1;
            if(is_rl(sym)){
                buff_decomp_nt_int<vector_t>(rules[start], exp, ht);
                copy_to_front(exp, pos1, exp.size() - pos1, rules[end] - 1);
            }else{
                for(size_t i=start;i<=end;i++){
                    buff_decomp_nt_int<vector_t>(rules[i], exp, ht);
                }
            }

            uint32_t len = exp.size() - pos1;
            locus_t locus = {pos1, len};
            ht.insert(&sym, sdsl::bits::hi(sym)+1, locus);
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

void grammar::se_decomp_str(size_t start, size_t end, std::string& output_file,
                            std::string& tmp_folder, size_t n_threads,
                            size_t ht_buff_size, size_t file_buff_size) {

    o_file_stream<char> ofs(output_file, file_buff_size, std::ios::out);

    size_t c_start = nter_ptr[gram_alph-1];
    size_t str_start, str_end;

    void * ht_buff = malloc(ht_buff_size);
    hash_table_t ht(ht_buff_size, "", 0.8, ht_buff);

    for(size_t i=start;i<=end;i++){

        str_start = i == 0 ? c_start : c_start+seq_pointers[i-1]+1;
        str_end = seq_pointers[i]+c_start;

        for(size_t j=str_start;j<=str_end;j++){
            buff_it_decomp_nt_int(rules[j], ofs,ht);
        }
    }
    free(ht_buff);
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
