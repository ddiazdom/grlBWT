//
// Created by Diaz, Diego on 10.11.2021.
//

#include "grammar.hpp"

template<class vector_type>
typename grammar<vector_type>::size_type
grammar<vector_type>::serialize(std::ostream& out, sdsl::structure_tree_node * v, std::string name) const {
    sdsl::structure_tree_node* child = sdsl::structure_tree::add_child( v, name, sdsl::util::class_name(*this));
    size_type written_bytes= sdsl::write_member(comp_level, out, child, "comp_level");
    written_bytes+= sdsl::write_member(text_size, out, child, "text_size");
    written_bytes+= sdsl::write_member(n_strings, out, child, "n_strings");
    written_bytes+= sdsl::write_member(grammar_size, out, child, "grammar_size");
    written_bytes+= sdsl::write_member(text_alph, out, child, "n_ter");
    written_bytes+= sdsl::write_member(gram_alph, out, child, "gram_alph");
    written_bytes+= sdsl::write_member(comp_string_size, out, child, "comp_string_size");
    written_bytes+= sdsl::write_member(n_p_rounds, out, child, "n_p_rounds");
    written_bytes+= rules_breaks.serialize(out, child, "rules_breaks");
    written_bytes+= symbols_map.serialize(out, child, "symbols_map");
    written_bytes+= m_rules.serialize(out, child, "m_rules");
    written_bytes+= m_nter_ptr.serialize(out, child, "nter_pointers");
    written_bytes+= m_seq_pointers.serialize(out, child, "m_seq_pointers");
    return written_bytes;
}

template<class vector_type>
void grammar<vector_type>::load(std::istream& in){
    sdsl::read_member(comp_level, in);
    sdsl::read_member(text_size, in);
    sdsl::read_member(n_strings, in);
    sdsl::read_member(grammar_size, in);
    sdsl::read_member(text_alph, in);
    sdsl::read_member(gram_alph, in);
    sdsl::read_member(comp_string_size, in);
    sdsl::read_member(n_p_rounds, in);
    rules_breaks.load(in);
    symbols_map.load(in);
    m_rules.load(in);
    m_nter_ptr.load(in);
    m_seq_pointers.load(in);
}

template<class vector_type>
void grammar<vector_type>::mark_str_boundaries(std::string& rules_file) {

    std::stack<size_t> stack;
    size_t pos = m_nter_ptr[gram_alph - 1], sym, suff_sym, sym_state, seq=0;
    //0 : the symbol's rule has not been visited yet
    //1 : the symbol recursively expands to a string suffix
    //2 : the symbol does not recursively expand to a string suffix
    sdsl::int_vector<2> state(gram_alph, 0);
    m_seq_pointers.width(sdsl::bits::hi(comp_string_size) + 1);
    m_seq_pointers.resize(n_strings);
    state[0] = 1;
    size_t c_start = pos;

    sdsl::int_vector<> rules_arr;
    sdsl::load_from_file(rules_arr, rules_file);

    while(pos<rules_arr.size()){
        sym =  rules_arr[pos];
        if(sym >= text_alph && state[sym] == 0){

            suff_sym = is_rl(sym) ? rules_arr[m_nter_ptr[sym]] : rules_arr[m_nter_ptr[sym + 1] - 1];

            while(suff_sym!=sym && state[suff_sym]==0){
                stack.push(suff_sym);
                sym = suff_sym;
                suff_sym = is_rl(sym) ? rules_arr[m_nter_ptr[sym]] : rules_arr[m_nter_ptr[sym + 1] - 1];
            }

            sym_state = state[suff_sym]==1 ? 1: 2;

            while(!stack.empty()){
                suff_sym = stack.top();
                state[suff_sym] = sym_state;
                stack.pop();
            }
            state[rules_arr[pos]] = sym_state;
            if(sym_state==1) m_seq_pointers[seq++] = pos - c_start;
        }else if(state[sym]==1){
            m_seq_pointers[seq++] = pos - c_start;
        }
        pos++;
    }
}

template<class vector_type>
std::string grammar<vector_type>::decomp_str(size_t idx) const {
    assert(idx < m_seq_pointers.size());
    std::string exp;
    size_t c_start = m_nter_ptr[gram_alph - 1];
    size_t str_start, str_end;

    str_start = idx == 0 ? c_start : c_start + m_seq_pointers[idx - 1] + 1;
    str_end = m_seq_pointers[idx] + c_start;

    for(size_t j=str_start;j<=str_end;j++){
        decomp_nt(m_rules[j], exp);
    }
    return exp;
}

/*
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
}*/

template<>
template<>
bool grammar<sdsl::int_vector<>>::copy_to_front<o_file_stream<char>>(o_file_stream<char>& stream, size_t src, size_t len, size_t freq) const {
    return stream.copy_to_front((buff_s_type) src, (buff_s_type) len, (buff_s_type) freq);
}

template<>
template<>
bool grammar<huff_vector<>>::copy_to_front<o_file_stream<char>>(o_file_stream<char>& stream, size_t src, size_t len, size_t freq) const {
    return stream.copy_to_front((buff_s_type) src, (buff_s_type) len, (buff_s_type) freq);
}

template<class vector_type>
template<class vector_t>
bool grammar<vector_type>::copy_to_front(vector_t &stream, size_t src, size_t len, size_t freq) const {
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
template<class vector_type>
template<class vector_t>
void grammar<vector_type>::buff_it_decomp_nt_int(size_t root, vector_t& exp, std::vector<size_t>& dc_info) {

    std::stack<node_t> st;

    size_t len=0, freq, exp_size;
    size_t locus, src, exp_len;
    node_t node = {root, 0, 0, false};

    do {
        while(true) {
            st.push(node);
            if(node.sym<text_alph){
                break;
            } else {

                //check if the nonterminal was already decompressed
                locus = dc_info[node.sym];
                if(locus!=0){
                    src = locus & ((1UL<<40)-1UL);
                    exp_len = locus >> 40;
                    if(copy_to_front(exp, src, exp_len, 1)){
                        len += exp_len;
                        break;
                    }
                }
                node.g_pos = m_nter_ptr[node.sym];
                node.is_rm_child = (is_rl(node.sym) || node.g_pos==(m_nter_ptr[node.sym + 1] - 1));
                node.sym = m_rules[node.g_pos];
            }
        }

        auto temp = st.top();
        st.pop();

        if(temp.sym<text_alph){
            exp.push_back((char)symbols_map[temp.sym]);
            len++;
        } else {
            exp_len = len - temp.l_exp;
            src = exp.size()-exp_len;
            dc_info[temp.sym] = (exp_len<<40) | src;
        }

        while(!st.empty() && temp.is_rm_child) {

            temp = st.top();
            st.pop();

            if(is_rl(temp.sym)){
                exp_size = len - temp.l_exp;
                freq = m_rules[m_nter_ptr[temp.sym + 1] - 1];
                copy_to_front(exp, exp.size()-exp_size, exp_size, freq-1);
                len+=exp_size*(freq-1);
            }

            exp_len = len - temp.l_exp;
            src = exp.size() - exp_len;
            dc_info[temp.sym] = (exp_len<<40) | src;
        }

        if (!st.empty()) {
            assert(!temp.is_rm_child);
            node.g_pos = temp.g_pos+1;
            node.is_rm_child = node.g_pos == (m_nter_ptr[st.top().sym + 1] - 1);
            node.sym = m_rules[node.g_pos];
            node.l_exp = len;
        }
    }while(!st.empty());
}

//buffered decompression
template<class vector_type>
template<class vector_t>
void grammar<vector_type>::buff_decomp_nt(size_t sym, vector_t& exp, std::vector<size_t>& dc_info) const {

    size_t pos1, start, end, src, exp_len, locus;
    if(sym<text_alph){
        exp.push_back((char)symbols_map[sym]);
    }else{
        pos1 = exp.size();
        locus = dc_info[sym];
        if(locus!=0){

            src = locus & ((1UL<<40)-1UL);
            exp_len = locus >> 40;

            //if the string couldn't be copied, then we do it the hard way
            if(!copy_to_front(exp, src, exp_len, 1)){
                start = m_nter_ptr[sym];
                end = m_nter_ptr[sym + 1] - 1;
                if(is_rl(sym)){
                    buff_decomp_nt<vector_t>(m_rules[start], exp, dc_info);
                    copy_to_front(exp, pos1, exp.size() - pos1, m_rules[end] - 1);
                }else{
                    for(size_t i=start;i<=end;i++){
                        buff_decomp_nt<vector_t>(m_rules[i], exp, dc_info);
                    }
                }
            }
            dc_info[sym] = (exp_len<<40) | pos1;
        }else{
            start = m_nter_ptr[sym];
            end = m_nter_ptr[sym + 1] - 1;
            if(is_rl(sym)){
                buff_decomp_nt<vector_t>(m_rules[start], exp, dc_info);
                copy_to_front(exp, pos1, exp.size() - pos1, m_rules[end] - 1);
            }else{
                for(size_t i=start;i<=end;i++){
                    buff_decomp_nt<vector_t>(m_rules[i], exp, dc_info);
                }
            }
            dc_info[sym] = ((exp.size() - pos1)<<40) | pos1;
        }
    }
}

template<class vector_type>
void grammar<vector_type>::decomp_nt(size_t sym, std::string& exp) const {
    assert(sym<gram_alph);
    size_t start, end;
    if(sym<text_alph){
        exp.push_back((char)symbols_map[sym]);
    }else{
        start = m_nter_ptr[sym];
        end = m_nter_ptr[sym + 1] - 1;
        if(is_rl(sym)){
            size_t prev_size = exp.size();
            decomp_nt(m_rules[start], exp);
            copy_to_front(exp, prev_size, exp.size()-prev_size, m_rules[end] - 1);
        }else{
            for(size_t j=start;j<=end;j++){
                decomp_nt(m_rules[j], exp);
            }
        }
    }
}

template<class vector_type>
void grammar<vector_type>::se_decomp_str(size_t start, size_t end, std::string& output_file,
                            std::string& tmp_folder, size_t n_threads, size_t buff_size) const {

    /*
    //the buffer size cannot be greater than the text size
    buff_size = std::min<size_t>(buff_size, text_size);

    //each thread should cover at least 200 MiB
    size_t syms_per_threads = std::max<size_t>(1024*1024*200, INT_CEIL(text_size, n_threads));
    size_t t_start=0, t_end, i=1;
    while(t_start<text_size){
        t_end = std::min<size_t>((i*syms_per_threads), text_size)-1;
        std::cout<<t_start<<" "<<t_end<<std::endl;
        t_start=t_end+1;
        i++;
    }*/

    o_file_stream<char> ofs(output_file, buff_size, std::ios::out);
    size_t c_start = m_nter_ptr[gram_alph - 1];
    size_t str_start, str_end;

    std::vector<size_t> dc_info(gram_alph, 0);
    for(size_t i=start;i<=end;i++){
        str_start = i == 0 ? c_start : c_start + m_seq_pointers[i - 1] + 1;
        str_end = m_seq_pointers[i] + c_start;
        for(size_t j=str_start;j<=str_end;j++){
            buff_decomp_nt(m_rules[j], ofs, dc_info);
        }
    }

    /*auto worker = [&](const grammar& gram, size_t start, size_t end){
        std::vector<size_t> dc_info(gram.gram_alph, 0);
        size_t str_start, str_end;
        size_t c_start = gram.m_nter_ptr[gram.gram_alph-1];
        for(size_t i=start;i<=end;i++){
            str_start = i == 0 ? c_start : c_start+gram.m_seq_pointers[i-1]+1;
            str_end = gram.m_seq_pointers[i]+c_start;
            for(size_t j=str_start;j<=str_end;j++){
                gram.buff_decomp_nt(gram.m_rules[j], ofs, dc_info);
            }
        }
    };
    worker(*this, 0,1);*/
}



void gram_info_t::save_to_file(std::string& output_file){

    size_t buffer[255];

    std::ofstream of_stream(output_file, std::ofstream::binary);

    //write number of m_rules and alphabet size
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

template typename grammar<sdsl::int_vector<>>::size_type grammar<sdsl::int_vector<>>::serialize(std::ostream& out, sdsl::structure_tree_node * v, std::string name) const;
template typename grammar<huff_vector<>>::size_type grammar<huff_vector<>>::serialize(std::ostream& out, sdsl::structure_tree_node * v, std::string name) const;
template void grammar<huff_vector<>>::se_decomp_str(size_t start, size_t end, std::string& output_file, std::string& tmp_folder, size_t n_threads, size_t buff_size) const;
template void grammar<huff_vector<>>::load(std::istream& in);
template void grammar<sdsl::int_vector<>>::load(std::istream& in);
template void grammar<sdsl::int_vector<>>::decomp_nt(size_t idx, std::string& exp) const;
template void grammar<huff_vector<>>::decomp_nt(size_t idx, std::string& exp) const;
template void grammar<sdsl::int_vector<>>::se_decomp_str(size_t start, size_t end, std::string& output_file, std::string& tmp_folder, size_t n_threads, size_t buff_size) const;
template void grammar<sdsl::int_vector<>>::mark_str_boundaries(std::string& rules_file);
template void grammar<huff_vector<>>::mark_str_boundaries(std::string& rules_file);
template void grammar<sdsl::int_vector<>>::buff_decomp_nt<std::string>(size_t nt, std::string& exp, std::vector<size_t>& dc_info) const;
template void grammar<sdsl::int_vector<>>::buff_decomp_nt<o_file_stream<char>>(size_t nt, o_file_stream<char>& exp, std::vector<size_t>& dc_info) const;
