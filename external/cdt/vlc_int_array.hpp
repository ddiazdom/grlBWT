//
// Created by Diego Diaz on 9/3/19.
//

#ifndef MULTI_BWT_VLC_INT_VECTOR_HPP
#define MULTI_BWT_VLC_INT_VECTOR_HPP

#include <iostream>
#include <cstring>
#include <cassert>
#include "integer_encoders.hpp"
#include "inv_bitstream.hpp"

//!pushing new elements will invalidate an iterator instance
template<class vlc_arr_t>
class vlc_array_iterator{

    typedef vlc_array_iterator<vlc_arr_t>  self_type;
    typedef typename vlc_arr_t::size_type  value_type;
    typedef typename vlc_arr_t::coder_type coder_type;

private:
    const vlc_arr_t& array;
    size_t           index;
    value_type       value;
    size_t           len;
public:

    typedef value_type& reference;
    typedef value_type* pointer;
    typedef std::forward_iterator_tag iterator_category;

    vlc_array_iterator(const vlc_arr_t& array_, size_t index_) : array(array_), index(index_){
        if(index<array.m_written_bits){
            auto res = coder_type::decode(array.stream.stream, index);
            value = res.first;
            len = res.second-index;
        }else{
            index = array.m_written_bits;
        }
    };

    vlc_array_iterator(vlc_array_iterator&& other) noexcept : array(other.array),
                                                              index(other.index),
                                                              value(other.value),
                                                              len(other.len){};

    inline vlc_array_iterator& operator++(){
        size_t tmp_index = index+len;
        if(tmp_index<array.m_written_bits){
            auto res = coder_type::decode(array.stream.stream, tmp_index);
            value = res.first;
            index = tmp_index;
            len = res.second-index;
        }else{
           index = array.m_written_bits;
        }
        return *this;
    }

    inline vlc_array_iterator operator++(int){
        vlc_array_iterator old(*this);
        ++(*this);
        return old;
    }

    inline reference& operator*() {
        return value;
    }

    inline bool operator==(const self_type& rhs) const {
        return index == rhs.index;
    }

    inline bool operator!=(const self_type& rhs) const {
        return index != rhs.index;
    }

    vlc_array_iterator& operator=(self_type&& other) noexcept {
        move(std::forward<self_type>(other));
        return *this;
    }

    void move(self_type&& other) {
        std::swap(array, other.array);
        std::swap(index, other.index);
        std::swap(value, other.value);
    }
};

template<class coder_t>
class vlc_int_array {

public:
    typedef typename coder_t::size_type                size_type;
    typedef vlc_array_iterator<vlc_int_array<coder_t>> iterator_t;

private:
    friend                iterator_t;
    typedef coder_t       coder_type;
    inv_bitstream<size_t> stream{nullptr,0};
    size_type             last_val=0;
    size_t                m_size=0;
    size_t                m_written_bits=0;

    void move(vlc_int_array<coder_t>&& other){
        std::swap(last_val, other.last_val);
        std::swap(m_size, other.m_size);
        std::swap(m_written_bits, other.m_written_bits);
        stream.swap(other.stream);
    }

    void copy(const vlc_int_array<coder_t>& other){
        last_val = other.last_val;
        m_size = other.m_size;
        m_written_bits = other.m_written_bits;

        if(stream.stream!=other.stream.stream){
            stream.stream_size = other.stream.stream_size;
            size_t n_bytes = stream.stream_size*sizeof(size_t);
            if(stream.stream!=nullptr){
                free(stream.stream);
            }
            stream.stream = reinterpret_cast<size_t*>(malloc(n_bytes));
            memcpy(stream.stream, other.stream.stream, n_bytes);
        }
    }

public:

    //default constructor
    vlc_int_array() {
        stream.stream_size = 2;
        stream.stream = reinterpret_cast<size_t*>(malloc(sizeof(size_t)*2));
    };

    //copy constructor
    vlc_int_array(const vlc_int_array<coder_t>& other) {
        copy(other);
    };

    //move constructor
    vlc_int_array(vlc_int_array<coder_t>&& other) noexcept {
        move(std::forward<vlc_int_array<coder_t>>(other));
    };

    vlc_int_array<coder_t>& operator=(const vlc_int_array<coder_t>& other) {
        if(&other!=this){
            copy(other);
        }
        return *this;
    };

    vlc_int_array<coder_t>& operator=(vlc_int_array<coder_t>&& other) noexcept {
        if(&other!=this){
            move(std::forward<vlc_int_array<coder_t>>(other));
        }
        return *this;
    };

    ~vlc_int_array(){
        if(stream.stream!=nullptr){
            free(stream.stream);
        }
    };

    void swap(vlc_int_array<coder_t>& other){
        move(std::move(other));
    }

    iterator_t begin(){
        return iterator_t(*this, 0);
    };

    iterator_t end(){
        return iterator_t(*this, m_written_bits);
    };

    inline size_type back() const {
        return last_val;
    };

    inline size_t read(size_t idx) const{
        assert(idx<size());
        size_t i=0, pos=0;
        while(i<idx){
            auto res = coder_type::decode(stream.stream, pos);
            pos = res.second;
            i++;
        }
        auto res = coder_type::decode(stream.stream, pos);
        return res.first;
    }

    void push_back(size_type value) {
        auto code = coder_t::encode(value);
        if(INT_CEIL((m_written_bits+code.second), inv_bitstream<size_t>::word_bits)>=stream.stream_size){
            stream.stream_size *=1.5;
            stream.stream = reinterpret_cast<size_t*>(realloc(stream.stream, stream.stream_size*sizeof(size_t)));
        }

        stream.write(m_written_bits, m_written_bits+code.second-1, code.first);
        m_size++;
        m_written_bits+=code.second;
        last_val = value;
    };

    inline size_t size() const{
        return m_size;
    };

    void shrink_to_fit(){
        stream.stream_size = INT_CEIL(m_written_bits, inv_bitstream<size_t>::word_bits);
        stream.stream = reinterpret_cast<size_t*>(realloc(stream.stream, stream.stream_size*sizeof(size_t)));
    };

    size_t buffer_bytes() const{
        return stream.stream_size*sizeof(size_t);
    }
};


/*template<class coder_t, class allocator_t>
inline void vlc_int_array<coder_t, allocator_t>::push_back(size_type value) {

    //! this way of handling memory is not okay
    if(m_buffer_size==0 && m_buffer== nullptr){
        m_buffer_size = 20;
        m_buffer = m_allocator.allocate(m_buffer_size);
    }

    auto code = coder_t::encode(value);
    size_t n_words = 1 + (((written_bits+code.second)-1) / ULONG_WIDTH);
    if(n_words>=m_buffer_size){
        size_type * tmp = m_allocator.allocate(m_buffer_size*2);
        memcpy(tmp, m_buffer, m_buffer_size * 8);
        std::swap(tmp, m_buffer);
        m_allocator.deallocate(tmp, m_buffer_size);
        m_buffer_size*=2;
    }
    //!

    write_bits(code.first,
               written_bits,
               (written_bits+code.second-1));

    written_bits+=code.second;
    m_size++;
    last_val = value;
}

template<class coder_t, class allocator_t>
inline size_t vlc_int_array<coder_t, allocator_t>::size() const {
    return m_size;
}

template<class coder_t, class allocator_t>
void vlc_int_array<coder_t, allocator_t>::shrink_to_fit() {
    size_t n_words = 1 + ((written_bits-1) / ULONG_WIDTH);
    if(n_words!=m_buffer_size){
        size_type * tmp = m_allocator.allocate(n_words);
        memcpy(tmp, m_buffer, n_words * 8);
        std::swap(tmp, m_buffer);
        m_allocator.deallocate(tmp, m_buffer_size);
        m_buffer_size=n_words;
    }
}

template<class coder_t, class allocator_t>
void vlc_int_array<coder_t, allocator_t>::write_bits(size_type elm, size_t i, size_t j) {
    //! bits are read backward, i.e., from left to right as in a regular array
    size_t cell_i, cell_j;
    cell_i = i/ULONG_WIDTH;
    cell_j = j/ULONG_WIDTH;

    assert(cell_i<=cell_j);
    assert(cell_j<m_buffer_size);

    if(cell_i==cell_j){
        size_t cell_i_pos = ULONG_WIDTH - (j & (ULONG_WIDTH - 1UL)) - 1;
        m_buffer[cell_j] &= ~(magic_numbers::masks[j-i+1] << cell_i_pos);
        m_buffer[cell_j] |= elm << cell_i_pos;
    }else{
        size_t cell_i_pos = ULONG_WIDTH - (i & (ULONG_WIDTH - 1UL))-1; //equivalent to (i % ULONG_WIDTH)
        size_t cell_j_pos = ULONG_WIDTH - (j & (ULONG_WIDTH - 1UL))-1;

        m_buffer[cell_i] = (m_buffer[cell_i] & ~magic_numbers::masks[cell_i_pos+1] ) |
                           (elm >> (ULONG_WIDTH - cell_j_pos));
        m_buffer[cell_j] = (m_buffer[cell_j] & magic_numbers::masks[cell_j_pos] ) |
                           (elm << (cell_j_pos));
    }
}

template<class coder_t, class allocator_t>
inline std::pair<typename coder_t::size_type, size_t> vlc_int_array<coder_t, allocator_t>::read_bits(size_t bit_index) const {
    if(bit_index>=written_bits){
        return {0,written_bits+1};
    }else{
        return coder_t::decode(m_buffer, bit_index);
    }
}*/

#endif //MULTI_BWT_VLC_INT_VECTOR_HPP
