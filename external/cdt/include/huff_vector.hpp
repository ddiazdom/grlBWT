//
// Created by diediaz on 02-12-19.
//

#ifndef CDT_HUFF_VECTOR_HPP
#define CDT_HUFF_VECTOR_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include <filesystem>
#include "inv_bitstream.hpp"
#include "iterators.hpp"
#include "macros.h"
#define BUFFER_SIZE 8388608

template<class arr_t=sdsl::int_vector_buffer<>,
         size_t smp_per=5> //sample the sample_per% of elements in the array
class huff_vector {

public:
    typedef typename arr_t::value_type                     symbol_type;
    typedef unsigned long                                  size_type;
    typedef huff_vec_iterator<huff_vector<arr_t, smp_per>> const_iterator;

private:
    using ivb    =                                         sdsl::int_vector_buffer<>;

    friend struct                                          huff_vec_iterator<huff_vector<arr_t, smp_per>>;
    uint64_t                                               *first=nullptr;
    uint64_t                                               *heights=nullptr;
    sdsl::int_vector<>                                     symbols;
    inv_bitstream<size_t>                                  bit_str;
    sdsl::int_vector<>                                     sampled_pointers;

    size_t                                                 max_depth;
    size_t                                                 m_size;
    size_t                                                 m_samp_window;
    size_t                                                 m_written_bits;

    struct h_int_node{
        bool is_leaf : 1;
        size_t freq : 63;
        h_int_node *next;

        h_int_node *left;
        h_int_node *right;
    };

    struct h_leaf{
        bool is_leaf : 1;
        size_t freq : 63;
        h_int_node * next;

        size_t symbol;
    };

    struct symbol_data{
        size_t symbol;
        uint8_t len : 7;
        size_t  code: 57;
    };

    void compute_code_lengths(ivb& sym_data, size_t& eff_alph, h_int_node* hn, size_t depth){
        if(hn->is_leaf){
            sym_data.push_back(reinterpret_cast<h_leaf*>(hn)->symbol);
            sym_data.push_back(depth);
        }else{
            compute_code_lengths(sym_data, eff_alph, hn->left, depth+1);
            compute_code_lengths(sym_data, eff_alph, hn->right, depth+1);
        }
    }

    h_int_node* build_tree(void * huff_mem, size_t& eff_alph){
        //connect the leaves as a linked list
        auto leaves = reinterpret_cast<h_leaf*>(huff_mem);
        for(size_t i=0; i < eff_alph-1; i++){
            leaves[i].next = reinterpret_cast<h_int_node*>(&leaves[i+1]);
        }

        auto huff_tree  = reinterpret_cast<h_int_node*>(huff_mem);
        auto tmp = huff_tree;

        //address where the internal nodes start
        auto av_buckets = reinterpret_cast<h_int_node*>(leaves+eff_alph);

        while(huff_tree->next!= nullptr){
            //create a new internal node
            h_int_node *new_node = av_buckets++;

            new_node->left = huff_tree;
            new_node->right = huff_tree->next;
            new_node->freq = new_node->left->freq + new_node->right->freq;

            while(tmp->next!= nullptr && tmp->next->freq<new_node->freq){
                tmp = tmp->next;
            }

            new_node->next = tmp->next;
            tmp->next = new_node;
            huff_tree = huff_tree->next->next;
        }
        return huff_tree;
    };

    void get_prefix_free_codes(symbol_data* &syms, size_t& eff_alph){
        //sort by code length
        std::sort(syms, syms+eff_alph, [](auto const& lhs, auto const& rhs){
            return lhs.len < rhs.len;
        });

        syms[0].code = 0;
        for(size_t i=1;i<eff_alph;i++){
            syms[i].code = (syms[i-1].code + 1)<< (syms[i].len - syms[i-1].len);
        }
    };

    void build_huff_rep(symbol_data* &syms, size_t& eff_alph){
        size_t max_symbol=0;
        for(size_t i=0;i<eff_alph;i++){
            if(syms[i].symbol>max_symbol){
                max_symbol = syms[i].symbol;
            }
        }
        assert(max_symbol!=0);

        symbols.width(sdsl::bits::hi(max_symbol)+1);
        symbols.resize(eff_alph+1);
        symbols[0] = 0;
        for(size_t i=1;i<eff_alph+1;i++){
            symbols[i] = syms[i-1].symbol;
        }

        max_depth = syms[eff_alph-1].len;
        heights = (uint64_t*) malloc(sizeof(uint64_t)*8*(max_depth+1));
        first = (uint64_t*) malloc(sizeof(uint64_t)*8*(max_depth+1));

        for(size_t i=0;i<max_depth+1;i++) heights[i] = 0;

        for(size_t i=eff_alph;i-->0;){
            heights[syms[i].len] = i+1;
            first[syms[i].len] = syms[i].code;
        }

        for(size_t i=max_depth+1;i-->1;){
            if(heights[i]==0){
                heights[i] = heights[i+1];
                first[i] = first[i+1]/2;
            }
        }
    }

    [[nodiscard]] inline size_t binary_search(size_t value, size_t left, size_t right) const {
        size_t mid, left_bound, right_bound;
        while(true){
            mid = left + ((right-left)>>1UL);
            left_bound = first[mid]<<(max_depth-mid);
            right_bound = mid<max_depth ? (first[mid+1]<<(max_depth-mid-1)) : 1UL<<max_depth;

            if(left_bound<=value && value<right_bound) break;

            if(value<left_bound){
                right = mid-1;
            }else{
                left = mid+1;
            }
        }
        return mid;
    }

    inline size_t decode_symbol(size_t& pos) const{
        auto block = bit_str.read(pos, pos+max_depth-1);
        size_t len = binary_search(block, 1, max_depth);
        block = block>>(max_depth-len);
        pos +=len;
        return symbols[heights[len] + block - first[len]];
    };

    /*void compute_sampled_pos(){
        sdsl::int_vector<> samp_p_tmp(1 + ((m_size - 1) / m_samp_window));
        size_t pos=0,smp_pos=0;
        for(size_t i=0;i<m_size;i++){
            if(i % m_samp_window == 0){
                samp_p_tmp[smp_pos++] = pos;
            }
            auto block = bit_str.read(pos, pos+max_depth-1);
            size_t len = binary_search(block, 1, max_depth);
            pos +=len;
        }
        sdsl::util::bit_compress(samp_p_tmp);
        sampled_pointers.swap(samp_p_tmp);
    }*/

public:

    huff_vector(): max_depth(0), m_size(0), m_samp_window(0), m_written_bits(0){};

    explicit huff_vector(arr_t &input_array): max_depth(0), m_size(0), m_written_bits(0){

        std::filesystem::path path(input_array.filename());
        sdsl::cache_config config(false, path.parent_path());
        size_t max_symbol=0;

        //represent the alphabet as a bit vector
        for (auto const &sym : input_array) if (sym > max_symbol) max_symbol = sym;
        sdsl::bit_vector alph_bv(max_symbol+1, false);
        sdsl::bit_vector::rank_1_type alph_bv_rs;
        for (auto const &sym : input_array) alph_bv[sym]=true;
        sdsl::util::init_support(alph_bv_rs, &alph_bv);
        size_t eff_alph = alph_bv_rs(alph_bv.size());
        //

        //allocate space for the huffman tree
        size_t huff_bytes = sizeof(h_leaf)*eff_alph + sizeof(h_int_node)*eff_alph;
        void * huff_memory = malloc(huff_bytes);
        memset((char *)huff_memory, 0, huff_bytes);
        //

        //sort the symbols according their frequencies
        size_t idx;
        auto huff_leaves = reinterpret_cast<h_leaf *>(huff_memory);
        for (auto const &sym : input_array) {
            idx = alph_bv_rs(sym);
            huff_leaves[idx].is_leaf = true;//indicate node is a leaf
            huff_leaves[idx].symbol = sym;
            huff_leaves[idx].freq++;
        }
        std::sort(huff_leaves, huff_leaves+ eff_alph, [](auto &lhs, auto &rhs) {
            return lhs.freq < rhs.freq;
        });
        //

        //compute the statistical entropy
        /*double h0=0;
        double n = input_array.size();
        for(size_t i=0;i<eff_alph;i++){
            double ns = huff_leaves[i].freq;
            h0+=(ns/n)*log2(n/ns);
        }
        std::cout<<std::fixed<<"    Emp. entropy: "<<h0<< " ("<<(((h0*n)/8)/1000000)<<" megabytes) "<<input_array.size()<<std::endl;
        std::cout<<std::fixed<<"    Worst case entropy: "<<log2(max_symbol)<< " ("<<(((log2(max_symbol)*n)/8)/1000000)<<" megabytes) "<<input_array.size()<<std::endl;
         */
        //

        //build the huffman codes
        h_int_node* huff_tree = build_tree(huff_memory, eff_alph);
        ivb code_data(sdsl::cache_file_name("code_len", config), std::ios::out, BUFFER_SIZE);
        compute_code_lengths(code_data, eff_alph, huff_tree, 0);
        free(huff_memory);

        auto syms = reinterpret_cast<symbol_data *>(malloc(eff_alph*sizeof(symbol_data)));
        for(size_t i=0, j=0; i < code_data.size(); i+=2, j++){
            syms[j].symbol = code_data[i];
            syms[j].len = code_data[i + 1];
        }
        code_data.close(true);
        get_prefix_free_codes(syms, eff_alph);
        //

        //build the succinct representation of the canonical huffman codes
        build_huff_rep(syms, eff_alph);

        //resort the symbols according their original values
        std::sort(syms, syms+eff_alph, [](auto &lhs, auto &rhs) {
            return lhs.symbol < rhs.symbol;
        });

        bit_str.stream_size = BUFFER_SIZE / sizeof(size_t);
        bit_str.stream = reinterpret_cast<size_t *>(malloc(BUFFER_SIZE));
        std::string str_file = sdsl::cache_file_name("stream", config);

        m_samp_window = 100 / smp_per;
        std::string samp_file = sdsl::cache_file_name("samples", config);
        ivb samp_buff(samp_file, std::ios::out, BUFFER_SIZE);

        //store the stream of codes
        {
            std::ofstream str_f(str_file, std::ios::binary);
            m_written_bits = 0;

            size_t i = 0, j, pos = 0, max_bits = BUFFER_SIZE * 8;
            for (auto const &sym : input_array) {

                if (pos % m_samp_window == 0){
                    samp_buff.push_back(m_written_bits);
                }

                auto data = syms[alph_bv_rs(sym)];
                j = i + data.len - 1;

                if (j >= max_bits) {
                    if(i%max_bits!=0){//there is a remainder
                        str_f.write((char *) bit_str.stream, (bit_str.stream_size-1) * sizeof(size_t));
                        bit_str.stream[0] = bit_str.stream[bit_str.stream_size - 1];
                    }else{ //buffer was completely used
                        str_f.write((char *) bit_str.stream, bit_str.stream_size * sizeof(size_t));
                    }

                    i = i % inv_bitstream<size_t>::word_bits;
                    j = i + data.len - 1;
                }

                bit_str.write(i, j, data.code);
                i = j + 1;
                pos++;
                m_written_bits +=data.len;
            }

            size_t bytes_last_chunk = INT_CEIL(i + max_depth, inv_bitstream<size_t>::word_bits) * sizeof(size_t);
            if (bytes_last_chunk > BUFFER_SIZE) {
                bit_str.stream = reinterpret_cast<size_t *>(realloc(bit_str.stream, bytes_last_chunk));
            }
            str_f.write((char *) bit_str.stream, bytes_last_chunk);
            str_f.close();
        }
        free(syms);

        //load the complete stream from disk
        bit_str.stream_size = INT_CEIL(m_written_bits+max_depth, inv_bitstream<size_t>::word_bits);
        bit_str.stream = reinterpret_cast<size_t *>(realloc(bit_str.stream, bit_str.stream_size*sizeof(size_t)));
        std::ifstream str_f(str_file, std::ios::in | std::ios::binary);
        str_f.read((char *)bit_str.stream, bit_str.stream_size*sizeof(size_t));
        if(remove(str_file.c_str())){
            std::cout<<"Error trying to remove file "<<str_file<<std::endl;
            exit(1);
        }
        //

        //load the sampled pointer from disk
        assert(samp_buff[samp_buff.size()-1]!=0);
        samp_buff.close();
        sdsl::load_from_file(sampled_pointers, samp_file);
        sdsl::util::bit_compress(sampled_pointers);
        if(remove(samp_file.c_str())){
            std::cout<<"Error trying to remove file "<<samp_file<<std::endl;
            exit(1);
        }
        //
        m_size = input_array.size();

        /*std::cout<<sdsl::size_in_bytes(symbols)/1000000.0<<std::endl;
        std::cout<<sdsl::size_in_bytes(bit_str)/1000000.0<<std::endl;
        std::cout<<sdsl::size_in_bytes(sampled_pointers)/1000000.0<<std::endl;
        std::cout<<sdsl::size_in_bytes(*this)/1000000.0<<std::endl;*/
    }

    huff_vector(huff_vector<arr_t, smp_per> && other) noexcept {
        first = std::exchange(other.first, nullptr);
        heights = std::exchange(other.heights, nullptr);
        symbols.swap(other.symbols);
        bit_str.swap(other.bit_str);
        sampled_pointers.swap(other.sampled_pointers);
        max_depth = other.max_depth;
        m_size = other.m_size;
        m_samp_window = other.m_samp_window;
        m_written_bits = other.m_written_bits;
    }

    huff_vector& operator=(huff_vector<arr_t, smp_per> && other) noexcept {
        first = std::exchange(other.first, nullptr);
        heights = std::exchange(other.heights, nullptr);
        symbols.swap(other.symbols);
        bit_str.swap(other.bit_str);
        sampled_pointers.swap(other.sampled_pointers);
        max_depth = other.max_depth;
        m_size = other.m_size;
        m_samp_window = other.m_samp_window;
        m_written_bits = other.m_written_bits;
        return *this;
    }

    huff_vector& operator=(const huff_vector<arr_t, smp_per> & other) noexcept {
        if(this!=&other) {
            first = other.first;
            heights = other.heights;
            symbols = other.symbols;
            bit_str = other.bit_str;
            sampled_pointers = other.sampled_pointers;
            max_depth = other.max_depth;
            m_size = other.m_size;
            m_samp_window = other.m_samp_window;
            m_written_bits = other.m_written_bits;
        }
        return *this;
    }

    inline void get_range(size_t start, size_t end, std::vector<size_t> &res) const {
        size_t pos = sampled_pointers[start / m_samp_window];
        for(size_t i= (start / m_samp_window) * m_samp_window; i < start; i++){
            decode_symbol(pos);
        }
        size_t idx=0;
        for(size_t i=start;i<=end;i++){
            res[idx++]= decode_symbol(pos);
        }
    }

    [[nodiscard]] inline size_t size() const{
        return m_size;
    }

    inline size_t operator[](size_t index) const{
        size_t pos = sampled_pointers[index / m_samp_window];
        for(size_t i= (index / m_samp_window) * m_samp_window; i < index; i++){
            auto block = bit_str.read(pos, pos+max_depth-1);
            size_t len = binary_search(block, 1, max_depth);
            pos +=len;
        }
        return decode_symbol(pos);
    }

    [[nodiscard]] inline size_t get_bit_pos(size_t index) const {
        size_t pos;
        if(index<size()){
            pos = sampled_pointers[index / m_samp_window];
            for(size_t i= (index / m_samp_window) * m_samp_window; i < index; i++){
                decode_symbol(pos);
            }
        }else{
            pos = m_written_bits;
        }
        return pos;
    }

    ~huff_vector(){
        if(heights!= nullptr){
            free(heights);
        }
        if(first!= nullptr){
            free(first);
        }
        if(bit_str.stream!=nullptr){
            free(bit_str.stream);
        }
    }

    const_iterator begin() const{
        return const_iterator(*this, 0);
    }

    const_iterator end() const{
        return const_iterator(*this, size()+1);
    }

    void load(std::istream &in){
        symbols.load(in);
        bit_str.load(in);
        sampled_pointers.load(in);
        sdsl::read_member(max_depth, in);
        sdsl::read_member(m_size, in);
        sdsl::read_member(m_samp_window, in);
        sdsl::read_member(m_written_bits, in);

        heights = (uint64_t*) malloc(sizeof(uint64_t)*8*(max_depth+1));
        first = (uint64_t*) malloc(sizeof(uint64_t)*8*(max_depth+1));
        for(size_t i=0;i<max_depth+1;i++){
            sdsl::read_member(first[i], in);
        }
        for(size_t i=0;i<max_depth+1;i++){
            sdsl::read_member(heights[i], in);
        }
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += symbols.serialize(out, child, "symbols");
        written_bytes += bit_str.serialize(out, child, "bit_stream");
        written_bytes += sampled_pointers.serialize(out, child, "sampled_pointers");
        written_bytes += sdsl::write_member(max_depth, out, child, "max_depth");
        written_bytes += sdsl::write_member(m_size, out, child, "m_size");
        written_bytes += sdsl::write_member(m_samp_window, out, child, "m_samp_window");
        written_bytes += sdsl::write_member(m_written_bits, out, child, "m_written_bits");
        for(size_t i=0;i<max_depth+1;i++){
            written_bytes += sdsl::write_member(first[i], out, child, "first_"+std::to_string(i));
        }
        for(size_t i=0;i<max_depth+1;i++){
            written_bytes += sdsl::write_member(heights[i], out, child, "max_depth_"+std::to_string(i));
        }
        return written_bytes;
    }
};

#endif //CDT_HUFF_VECTOR_HPP
