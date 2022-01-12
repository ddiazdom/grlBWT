//
// Created by Diaz, Diego on 10.11.2021.
//

#ifndef LPG_COMPRESSOR_GRAMMAR_HPP
#define LPG_COMPRESSOR_GRAMMAR_HPP

#include <iostream>
#include <cdt/huff_vector.hpp>
#include <sdsl/int_vector.hpp>
#include "cdt/hash_table.hpp"
#include "cdt/file_streams.hpp"

struct simple_stack {

    std::vector<size_t> container;
    size_t last_pos;

    simple_stack():container(100, 0), last_pos(0){}

    simple_stack(simple_stack&& other) noexcept : container(std::move(other.container)),
                                                  last_pos(other.last_pos){}

    inline void push(size_t new_elm){
        if(container.size()==last_pos) container.resize(container.size()*2);
        container[last_pos++] = new_elm;
    }

    [[nodiscard]] inline size_t top() const {
        assert(last_pos>0);
        return container[last_pos-1];
    }

    inline void pop (){
        if(last_pos>0) last_pos--;
    }

    [[nodiscard]] inline bool empty() const {
        return last_pos==0;
    }

    inline void reset(){
        last_pos=0;
    }
};

struct gram_info_t{
    uint8_t                            sigma{}; // terminal's alphabet
    size_t                             r{}; //r: number of grammar symbols (nter + ter)
    size_t                             c{}; //c: length of the right-hand of the start symbol
    size_t                             g{}; //g: sum of the m_rules' right-hand sides
    size_t                             max_tsym{}; //highest terminal symbol
    size_t                             last_dict_size{}; //size of the last dictionary
    std::unordered_map<size_t,uint8_t> sym_map; //map terminal symbols to their original byte symbols
    std::string                        rules_file; // m_rules are concatenated in this array
    std::string                        rules_lim_file; // bit vector marking the last symbol of every right-hand
    std::vector<size_t>                rules_breaks; // number of m_rules generated in every LMS parsing round
    size_t                             n_p_rounds{}; // number of parsing rounds
    std::pair<size_t, size_t>          rl_rules={0,0};

    gram_info_t() = default;
    gram_info_t(std::string& r_file_,
                std::string& r_lim_file_): rules_file(r_file_),
                                           rules_lim_file(r_lim_file_) {}

    void save_to_file(std::string& output_file);
    void load_from_file(std::string &g_file);

    [[nodiscard]] bool is_terminal(const size_t& id) const {
        return sym_map.find(id) != sym_map.end();
    }

    [[nodiscard]] inline bool is_rl(size_t symbol) const{
        if(rl_rules.second==0) return false;
        return symbol>=rl_rules.first && symbol < (rl_rules.first + rl_rules.second);
    }

    [[nodiscard]] inline long long int parsing_level(size_t symbol) const{
        if(symbol < sigma) return 0;
        if(symbol>=rules_breaks[n_p_rounds]){
            return -1;
        }

        for(long long int i=0;i<int(n_p_rounds);i++){
            if(rules_breaks[i]<=symbol && symbol<rules_breaks[i+1]) return i+1;
        }
        return -1;
    }
};

template<class grammar_t>
class grammar_iterator{

    const grammar_t& grammar;
    size_t sym{};
    size_t dummy{};
    simple_stack s_stack;
public:
    explicit grammar_iterator(const grammar_t& _grammar, size_t start, size_t end) : grammar(_grammar),
                                                                                     dummy(grammar.symbols()){

        if(end<grammar.rules.size()){
            for(size_t j=end;j>start;j--){
                s_stack.push(grammar.rules[j]);
            }
            sym = grammar.rules[start];
        }
    }

    explicit grammar_iterator(const grammar_t& _grammar, size_t _sym) : grammar(_grammar),
                                                                        sym(_sym),
                                                                        dummy(grammar.symbols()){
        if(sym>=grammar.ter() && sym<dummy){
            size_t start = grammar.nter_ptr[sym];
            if(grammar.is_rl(sym)){
                size_t rl_sym = grammar.rules[start];
                size_t freq = grammar.rules[start+1];
                for(size_t i=0;i<freq;i++){
                    s_stack.push(rl_sym);
                }
            }else{
                size_t end = grammar.nter_ptr[sym+1]-1;
                for(size_t i=end;i>=start;i--){
                    s_stack.push(grammar.rules[i]);
                }
            }
        }
    }

    explicit grammar_iterator(grammar_t&& other): grammar(other.grammar),
                                                  s_stack(std::move(other.stack)),
                                                  sym(other.sym),
                                                  dummy(other.dummy){}

    explicit grammar_iterator(const grammar_t& other): grammar(other.grammar),
                                                       s_stack(other.stack),
                                                       sym(other.sym),
                                                       dummy(other.dummy){}
    inline void operator++(){
        if(sym<grammar.ter()){
            if(!s_stack.empty()){
                sym = s_stack.top();
                s_stack.pop();
            }else{
                sym = dummy;
            }
        }else if(sym<dummy) {
            size_t start = grammar.nter_ptr[sym];
            if(grammar.is_rl(sym)){
                size_t rl_sym = grammar.rules[start];
                size_t freq = grammar.rules[start+1];
                for(size_t i=0;i<freq;i++){
                    s_stack.push(rl_sym);
                }
            }else{
                size_t end = grammar.nter_ptr[sym+1]-1;
                for(size_t i=end;i>=start;i--){
                    s_stack.push(grammar.rules[i]);
                }
            }
            sym = s_stack.top();
            s_stack.pop();
        }
    }

    inline void skip_subtree(){
        if(!s_stack.empty()){
            sym = s_stack.top();
            s_stack.pop();
        }else{
            sym = dummy;
        }
    }

    inline void update_it(size_t start, size_t end){
        s_stack.reset();
        if(end<grammar.rules.size()){
            for(size_t j=end;j>start;j--){
                s_stack.push(grammar.rules[j]);
            }
            sym = grammar.rules[start];
        }
    }

    inline size_t operator*() const{
        return sym;
    }

    inline bool operator==(grammar_iterator& other) const{
        return sym==other.sym;
    }

    inline bool operator!=(grammar_iterator& other) const {
        return sym!=other.sym;
    }

    inline bool operator<(grammar_iterator& other) {

        //move to the right as long as they have the same sequence
        while (sym != dummy &&
               other.sym != dummy &&
               sym == other.sym) {
            skip_subtree();
            other.skip_subtree();
        }

        if (sym == dummy || other.sym == dummy) {//one is a proper prefix of the other
            return !(sym == dummy && other.sym != dummy);
        } else {
            long long int lvl_a, lvl_b;
            while (sym != dummy && other.sym != dummy) {

                while (!grammar.is_lc(sym)) ++(*this);
                while (!grammar.is_lc(other.sym)) ++other;
                lvl_a = grammar.parsing_level(sym);
                lvl_b = grammar.parsing_level(other.sym);

                while (lvl_a != lvl_b) {
                    if (lvl_b < lvl_a) {
                        ++(*this);
                        while (!grammar.is_lc(sym)) ++(*this);
                        lvl_a = grammar.parsing_level(sym);
                    } else if (lvl_a < lvl_b) {
                        ++other;
                        while (!grammar.is_lc(other.sym)) ++other;
                        lvl_b = grammar.parsing_level(other.sym);
                    }
                }

                if (sym != other.sym) {
                    return sym < other.sym;
                } else {
                    skip_subtree();
                    other.skip_subtree();
                }
            }
            return !(sym == dummy && other.sym != dummy);
        }
    }
};

template<class vector_type=sdsl::int_vector<>>
class grammar {

private:

    struct node_t{
        size_t sym=0;
        size_t g_pos=0;
        size_t l_exp=0;
        bool is_rm_child=false;
    };

    //typedef bit_hash_table<size_t, 64, size_t, 6, true> hash_table_t;
    typedef typename o_file_stream<char>::size_type buff_s_type;

    size_t                     text_size{}; //original size of the text
    size_t                     n_strings{}; //original size of the text
    size_t                     grammar_size{}; //size of the grammar (i.e., sum of the lengths of the right-hand sides of the m_rules)
    size_t                     text_alph{}; //alphabet size of the text
    size_t                     gram_alph{}; //alphabet size of the grammar (ter + nters)
    size_t                     comp_string_size{}; //size of the compressed string
    size_t                     n_p_rounds{}; //number of parsing rounds
    uint8_t                    comp_level{}; //compression level
    sdsl::int_vector<>         rules_breaks; //mark the first rule of every parsing round
    std::pair<size_t, size_t>  rl_rules={0,0};
    std::pair<size_t, size_t>  sp_rules={0,0};
    sdsl::int_vector<8>        symbols_map; //map a terminal to its original byte symbol
    vector_type                m_rules; //list of m_rules
    sdsl::int_vector<>         m_nter_ptr; //pointers of the m_rules in the m_rules' array
    sdsl::int_vector<>         m_seq_pointers; //pointers to the boundaries of the strings in the text

    void mark_str_boundaries(std::string& rules_file);
    template<class vector_t>
    void buff_decomp_nt(size_t nt, vector_t& exp, std::vector<size_t>& dc_info) const;
    template<class vector_t>
    void buff_it_decomp_nt_int(size_t sym, vector_t& exp, std::vector<size_t>& dc_info);
    template<class vector_t>
    bool copy_to_front(vector_t& stream, size_t src, size_t len, size_t freq) const;
public:
    typedef size_t size_type;
    typedef grammar_iterator<grammar<vector_type>> iterator;
    const vector_type& rules = m_rules;
    const sdsl::int_vector<>& nter_ptr = m_nter_ptr;
    const sdsl::int_vector<>& seq_ptr = m_seq_pointers;

    grammar()=default;

    explicit grammar(gram_info_t& gram_info, size_t _orig_size, size_t _n_seqs): text_size(_orig_size),
                                                                                 n_strings(_n_seqs),
                                                                                 grammar_size(gram_info.g),
                                                                                 text_alph(gram_info.sigma),
                                                                                 gram_alph(gram_info.r),
                                                                                 comp_string_size(gram_info.c),
                                                                                 n_p_rounds(gram_info.n_p_rounds){

        //TODO check if the alphabet is compressed
        symbols_map.resize(gram_info.sigma);
        for(auto const pair: gram_info.sym_map){
            symbols_map[pair.first] = pair.second;
        }

        rules_breaks.resize(gram_info.rules_breaks.size());
        for(size_t i=0;i<gram_info.rules_breaks.size();i++){
            rules_breaks[i] = gram_info.rules_breaks[i];
        }

        rl_rules = gram_info.rl_rules;
        //sp_rules = gram_info.sp_rules;

        sdsl::int_vector_buffer<1> rules_lim_buff(gram_info.rules_lim_file, std::ios::in);
        m_nter_ptr.width(sdsl::bits::hi(grammar_size) + 1);
        m_nter_ptr.resize(gram_alph+1);
        m_nter_ptr[0] = 0;

        size_t i=1;
        for(size_t j=1;j<rules_lim_buff.size();j++){
            if(rules_lim_buff[j-1]){
                m_nter_ptr[i++] = j;
            }
        }
        m_nter_ptr[gram_alph] = rules_lim_buff.size();
        rules_lim_buff.close();
        mark_str_boundaries(gram_info.rules_file);

        sdsl::int_vector_buffer<> rules_buff(gram_info.rules_file, std::ios::in);
        if constexpr(std::is_same<vector_type, sdsl::int_vector<>>::value) {
            m_rules.width(sdsl::bits::hi(gram_alph) + 1);
            m_rules.resize(rules_buff.size());
            i = 0;
            for (auto const &sym: rules_buff) m_rules[i++] = sym;
            rules_buff.close();
            comp_level = 1;
        } else {
            m_rules = huff_vector<>(rules_buff);
            comp_level = 2;
        }
    }

    [[nodiscard]] inline bool is_lc(size_t symbol) const{
        return symbol<rules_breaks[n_p_rounds];
    }

    [[nodiscard]] inline bool is_rl(size_t symbol) const{
        if(rl_rules.second==0) return false;
        return symbol>=rl_rules.first && symbol < (rl_rules.first + rl_rules.second);
    }

    [[nodiscard]] inline bool is_sp(size_t symbol) const{
        if(sp_rules.second==0) return false;
        return symbol>=sp_rules.first && symbol < (sp_rules.first + sp_rules.second);
    }

    [[nodiscard]] inline bool in_rl_zone(size_t idx) const{
        size_t start = m_nter_ptr[rl_rules.first];
        size_t end = m_nter_ptr[rl_rules.first+rl_rules.second]-1;
        return start<=idx && idx<=end;
    }

    [[nodiscard]] inline std::pair<size_t, size_t> rl_zone() const{
        size_t start = m_nter_ptr[rl_rules.first];
        size_t end = m_nter_ptr[rl_rules.first+rl_rules.second]-1;
        return {start, end};
    }

    [[nodiscard]] inline long long int parsing_level(size_t symbol) const{
        if(symbol < text_alph) return 0;
        if(symbol>=rules_breaks[n_p_rounds]){
            assert(is_sp(symbol) || is_rl(symbol));
            return -1;
        }

        for(long long int i=0;i<int(n_p_rounds);i++){
            if(rules_breaks[i]<=symbol && symbol<rules_breaks[i+1]) return i+1;
        }
        return -1;
    }

    void space_breakdown() const {
        std::cout<<"Space breakdown for the grammar representation:"<<std::endl;
        std::cout<<"  Total size:         "<<sdsl::size_in_bytes(*this)/1000000.0<<" MB"<<std::endl;
        std::cout << "    Grammar rules:    " << sdsl::size_in_bytes(m_rules) / 1000000.0 << " MB" << std::endl;
        std::cout << "    Grammar pointers: " << sdsl::size_in_bytes(m_nter_ptr) / 1000000.0 << " MB" << std::endl;
        std::cout << "    String pointers:  " << sdsl::size_in_bytes(m_seq_pointers) / 1000000.0 << " MB" << std::endl;
        std::cout<<"    Symbols map:      "<<sdsl::size_in_bytes(symbols_map)/1000000.0<<" MB"<<std::endl;
        std::cout<<"    Rules breaks:     "<<sdsl::size_in_bytes(rules_breaks)/1000000.0<<" MB"<<std::endl;
    }

    inline grammar_iterator<grammar<vector_type>> begin() const {
        return grammar_iterator<grammar<vector_type>>(*this, gram_alph-1);
    }

    inline grammar_iterator<grammar<vector_type>> grammar_node(size_t sym) const {
        assert(sym<gram_alph);
        return grammar_iterator<grammar<vector_type>>(*this, sym);
    }

    inline grammar_iterator<grammar<vector_type>> range(size_t start, size_t end) const {
        assert(start<=end && end < rules.size());
        return grammar_iterator<grammar<vector_type>>(*this, start, end);
    }

    grammar_iterator<grammar<vector_type>> end(){
        return grammar_iterator<grammar<vector_type>>(*this, gram_alph);
    }

    [[nodiscard]] inline size_t n_rl_rules() const{
        return rl_rules.second;
    }

    [[nodiscard]] inline size_t n_sp_rules() const{
        return sp_rules.second;
    }

    [[nodiscard]] inline size_t n_lc_rules() const{
        return rules_breaks[rules_breaks.size()-1]-ter();
    }

    [[nodiscard]] inline size_t gram_size() const{
        return grammar_size;
    }

    [[nodiscard]] inline size_t nter() const{
        return gram_alph-text_alph;
    }

    [[nodiscard]] inline size_t ter() const{
        return text_alph;
    }

    [[nodiscard]] inline size_t symbols() const{
        return gram_alph;
    }

    [[nodiscard]] std::string decomp_str(size_t idx) const;

    [[nodiscard]] inline size_t strings() const {
        return m_seq_pointers.size()-1;
    }

    [[nodiscard]] inline size_t t_size() const {
        return text_size;
    }

    void se_decomp_str(size_t start, size_t end,
                       std::string& output,std::string& tmp_folder,
                       size_t n_threads, size_t buff_size) const;
    void decomp_nt(size_t idx, std::string& exp) const;
    size_type serialize(std::ostream& out, sdsl::structure_tree_node * v=nullptr, std::string name="") const;
    void load(std::istream& in);
};

#endif //LPG_COMPRESSOR_GRAMMAR_HPP
