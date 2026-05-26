//
// Created by diego on 27-08-20.
//

#ifndef CDS_INT_ARRAY_H
#define CDS_INT_ARRAY_H

#include <iostream>
#include <limits>
#include <cassert>
#include <utility>
#include "memory_handler.hpp"
#include "bit_stream.h"

template<class word_t=uint64_t, uint8_t w_bits=std::numeric_limits<word_t>::digits>
struct int_array{

    static_assert(w_bits<=std::numeric_limits<word_t>::digits);
    typedef size_t size_type;
    typedef word_t value_type;
    typedef bit_stream<word_t, w_bits> stream_t;

    size_t m_cap = 0;
    size_t m_size = 0;
    uint8_t m_width = w_bits;
    stream_t bits;

    class proxy {

        int_array & m_arr;
        size_t m_idx;

    public:
        proxy(int_array & _arr, size_t _idx) : m_arr(_arr), m_idx(_idx) {
            assert(m_idx<m_arr.m_cap);
        }

        operator value_type() const {
            return m_arr.bits.read(m_idx * m_arr.m_width, (m_idx + 1) * m_arr.m_width - 1);
        }

        proxy& operator=(value_type value) {
            m_arr.write(m_idx, value);
            return *this;
        }

        proxy& operator=(const proxy& other) {
            m_arr.write(m_idx, static_cast<value_type>(other));
            return *this;
        }

        proxy& operator&=(size_t val) {
            size_t elm = m_arr.bits.read(m_idx * m_arr.m_width, (m_idx + 1) * m_arr.m_width - 1);
            m_arr.write(m_idx, elm & val);
            return *this;
        }

        proxy& operator++() {
            auto val = static_cast<value_type>(*this);
            m_arr.write(m_idx, val+1);
            return *this;
        }

        value_type operator++(int) {
            auto val = static_cast<value_type>(*this);
            ++(*this);
            return val;
        }

        proxy& operator--() {
            auto val = static_cast<value_type>(*this);
            m_arr.write(m_idx, val-1);
            return *this;
        }

        value_type operator--(int) {
            auto val = static_cast<value_type>(*this);
            --(*this);
            return val;
        }

        proxy& operator+=(const value_type x) {
            auto val = static_cast<value_type>(*this);
            m_arr.write(m_idx, val+x);
            return *this;
        }

        proxy& operator-=(const value_type x) {
            auto val = static_cast<value_type>(*this);
            m_arr.write(m_idx, val-x);
            return *this;
        }

        bool operator==(const proxy& x)const {
            return value_type(*this) == value_type(x);
        }

        bool operator<(const proxy& x)const {
            return value_type(*this) < value_type(x);
        }
    };

    //simple constructor (cap is the capacity of the vector) — disabled for w_bits==1
    template<uint8_t B = w_bits, std::enable_if_t<B != 1, int> = 0>
    int_array(size_t cap_, uint8_t width_) : m_cap(cap_), m_width(width_){
        assert(m_width <= w_bits);
        bits.stream_size = INT_CEIL(m_cap * width_, stream_t::word_bits);
        bits.stream = mem::allocate<word_t>(bits.stream_size);
        memset(bits.stream, 0, bits.stream_size*sizeof(word_t));
    }

    //bit_array constructor: capacity + default bit value (only for w_bits==1)
    template<uint8_t B = w_bits, std::enable_if_t<B == 1, int> = 0>
    int_array(size_t cap_, value_type def_val) : m_cap(cap_), m_width(1){
        assert(def_val <= 1);
        bits.stream_size = INT_CEIL(m_cap, stream_t::word_bits);
        bits.stream = mem::allocate<word_t>(bits.stream_size);
        m_size = m_cap;
        initialize(def_val, m_cap);
    }

    //initialize from allocated memory
    int_array(word_t *ptr, size_t size_, uint8_t width_): m_cap(size_), m_size(size_), m_width(width_){
        assert(m_width <= w_bits);
        bits.stream_size = INT_CEIL(m_size*m_width, stream_t::word_bits);
        bits.stream = ptr;
    }

    //initialize with the default value and size
    int_array(size_t size_, value_type def_val, uint8_t width_) : m_cap(size_), m_size(size_), m_width(width_){
        assert(m_width <= w_bits);
        bits.stream_size = INT_CEIL(m_size * width_, stream_t::word_bits);
        bits.stream = mem::allocate<word_t>(bits.stream_size);
        initialize(def_val, m_size);
    }

    //initializer-list constructor
    int_array(const std::initializer_list<size_t> list) {
        assert(list.size()>0);
        size_t max=0;
        for(auto const& sym : list) if(sym>max) max=sym;

        m_width = sym_width(max);
        assert(m_width<=w_bits);
        m_size = list.size();
        m_cap = m_size;

        bits.stream_size = n_words();
        bits.stream = mem::allocate<word_t>(bits.stream_size);
        size_t i=0;
        for(auto const& val : list){
            write(i++,val);
        }
        mask_tail();
    }

    ~int_array() {
        mem::deallocate(bits.stream);
    }

    //default constructor
    int_array()=default;

    //move constructor
    explicit int_array(int_array&& other) noexcept {
        bits.stream=nullptr;
        bits.stream_size=0;
        move(std::forward<int_array>(other));
    };

    //copy constructor
    int_array(const int_array& other) noexcept {
        copy(other);
    };

    void move(int_array&& other) noexcept {
        m_size = std::exchange(other.m_size, 0);
        m_cap = std::exchange(other.m_cap, 0);
        m_width = std::exchange(other.m_width, 0);
        bits.stream_size = std::exchange(other.bits.stream_size, 0);
        mem::deallocate(bits.stream);
        bits.stream = std::exchange(other.bits.stream, nullptr);
    }

    void copy(const int_array &other) {
        m_size = other.m_size;
        m_cap = other.m_cap;
        m_width = other.m_width;
        bits.stream_size = cap_words();
        mem::deallocate(bits.stream);
        bits.stream = mem::allocate<word_t>(bits.stream_size);
        memcpy(bits.stream, other.bits.stream, n_bytes());
    }


    void clear() {
        m_size=0;
    }

    void erase() {
        m_size=0;
        m_cap = 0;
        mem::deallocate(bits.stream);
        bits.stream = nullptr;
        bits.stream_size = 0;
    }

    //copy assignment operator
    int_array& operator=(const int_array & other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    //move assignment operator
    int_array& operator=(int_array && other) noexcept{
        if(this!=&other){
            move(std::forward<int_array>(other));
        }
        return *this;
    }

    bool operator==(const int_array &other) const {
        if(m_size!=other.m_size || m_width!=other.m_width) return false;
        size_t tot_bits = m_size * m_width;
        return bits.compare_chunk(other.bits.stream, 0, tot_bits);
    }

    [[nodiscard]] void* data() const{
        return bits.stream;
    }

    [[nodiscard]] word_t* stream() const{
        return bits.stream;
    }

    void set_data(word_t* new_data, size_t size) {
        bits.stream = new_data;
        bits.stream_size = size;
    }

    [[nodiscard]] size_t width() const{
        return m_width;
    }

    [[nodiscard]] size_t n_words() const{
        return m_size==0? 0: INT_CEIL(n_bits(), stream_t::word_bits);
    }

    [[nodiscard]] size_t cap_words() const{
        return m_cap==0? 0: INT_CEIL(cap_bits(), stream_t::word_bits);
    }

    [[nodiscard]] size_t n_bytes() const{//number of bytes used by the data (ceil)
        return INT_CEIL(n_bits(), 8);
    }

    [[nodiscard]] size_t cap_bytes() const{//number of bits used by the data (exact)
        return INT_CEIL(cap_bits(), 8);
    }

    [[nodiscard]] size_t n_bits() const{//number of bits used by the data (exact)
        return m_size * m_width;
    }

    [[nodiscard]] size_t cap_bits() const{//number of bits used by the data (exact)
        return m_cap * m_width;
    }

    void mask_tail() {
        size_t tot_bits = m_size * m_width;
        size_t n_cells = tot_bits / stream_t::word_bits;//floor

        if (tot_bits > n_cells * stream_t::word_bits) {
            bits.write(tot_bits, (n_cells + 1) * stream_t::word_bits - 1, 0);
        }
    }

    void push_back(value_type value) {
        if((m_size+1)>m_cap){
            m_cap = m_cap==0? 2: m_cap*2;
            reserve(m_cap);
        }
        size_t i = m_size * m_width;
        size_t j = (m_size+1) * m_width - 1;
        bits.write(i, j , value);
        m_size++;
    }

    void pop_back() {
        if(m_size>=1) m_size--;
    }

    void reserve(size_t new_cap) {
        //reserve memory for new_cap number of elements
        size_t new_stream_size = 0;
        if(new_cap>0){
            new_stream_size = INT_CEIL(new_cap * m_width, stream_t::word_bits);
        }

        if(new_stream_size>bits.stream_size){
            bits.stream = mem::reallocate<word_t>(bits.stream, new_stream_size);
            bits.stream_size = new_stream_size;
        }
        m_cap = new_cap;
    }

    void resize(size_t new_size) {
        //reserve memory for new_size number of elements
        size_t new_buffer_size = INT_CEIL(new_size * m_width, stream_t::word_bits);
        bits.stream = mem::reallocate<word_t>(bits.stream, new_buffer_size);
        bits.stream_size = new_buffer_size;
        m_size = new_size;
        m_cap = m_size;
    }

    [[nodiscard]] size_t size() const{
        return m_size;
    }

    [[nodiscard]] size_t capacity() const{
        return m_cap;
    }

    [[nodiscard]] bool empty() const{
        return m_size==0;
    }

    value_type back() const{
        assert(m_size>0);
        return bits.read((m_size-1)*m_width, m_size*m_width-1);
    }

    void write(size_t idx, value_type value){
        assert(idx<m_cap);
        if(idx>=m_size) m_size = idx+1;
        bits.write(idx * m_width, (idx + 1) * m_width - 1, value);
    }

    value_type read(const size_t idx) const {
        assert(idx<m_size);
        return bits.read(idx * m_width, (idx + 1) * m_width - 1);
    }

    value_type operator[](const size_t idx) const {
        assert(idx<m_size);
        return bits.read(idx * m_width, (idx + 1) * m_width - 1);
    }

    proxy operator[](const size_t idx) {
        return proxy(*this, idx);
    }

    size_type serialize(std::ostream &out) const{
        size_t written_bytes = bits.serialize(out);
        written_bytes +=serialize_elm(out, m_size);
        written_bytes +=serialize_elm(out, m_cap);
        written_bytes +=serialize_elm(out, m_width);
        return written_bytes;
    }

    void load(std::istream &in){
        bits.load(in);
        load_elm(in, m_size);
        load_elm(in, m_cap);
        load_elm(in, m_width);
    }

    void initialize(value_type val, size_t n_elems){
        assert(n_elems<=m_cap);
        size_t n_cells = INT_CEIL(n_elems * m_width, stream_t::word_bits);
        assert(n_cells <= bits.stream_size);

        if(val == 0){
            memset(bits.stream, 0, n_cells * sizeof(word_t));
            if (n_elems>m_size) m_size = n_elems;
            return;
        }

        // find GCD(m_width, word_bits) to compute the smallest repeating pattern
        size_t a = m_width, b = stream_t::word_bits;
        while(b){ a %= b; std::swap(a, b); }
        size_t pat_words = m_width / a;   // words in one pattern = LCM/word_bits
        size_t pat_elems = stream_t::word_bits / a;  // elements in one pattern = LCM/m_width

        // build pattern into a stack buffer (pat_words <= m_width <= w_bits)
        word_t pat_buf[w_bits] = {};
        stream_t tmp;
        tmp.stream = pat_buf;
        tmp.stream_size = pat_words;
        for(size_t i = 0; i < pat_elems; i++){
            tmp.write(i * m_width, (i + 1) * m_width - 1, val);
        }

        // fill stream by exponentially expanding the pattern
        size_t pat_bytes = pat_words * sizeof(word_t);
        size_t total_bytes = n_cells * sizeof(word_t);
        memcpy(reinterpret_cast<char*>(bits.stream), pat_buf, pat_bytes);
        size_t filled = pat_bytes;
        while(filled + filled <= total_bytes){
            memcpy(reinterpret_cast<char*>(bits.stream) + filled, bits.stream, filled);
            filled += filled;
        }
        if(filled < total_bytes){
            memcpy(reinterpret_cast<char*>(bits.stream) + filled, bits.stream, total_bytes - filled);
        }

        // mask tail bits beyond n_elems * m_width in the last word
        size_t rem = (n_elems * m_width) & (stream_t::word_bits - 1);
        if(rem > 0){
            bits.stream[n_cells - 1] &= stream_t::masks[rem];
        }
        if (n_elems>m_size) m_size = n_elems;
    }

    void set_to_value(value_type value) {
        initialize(value, m_size);
    }

    void swap(int_array& other) noexcept {
        std::swap(m_size, other.m_size);
        std::swap(m_cap, other.m_cap);
        std::swap(m_width, other.m_width);
        bits.swap(other.bits);
    }

    void set_width(const uint8_t new_width){
        assert(new_width<=w_bits);
        m_width = new_width;
    }
};

typedef int_array<uint64_t, 1> bit_array;
#endif //CDS_INT_ARRAY_H
