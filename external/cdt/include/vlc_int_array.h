//
// Created by Diego Diaz on 9/3/19.
//

#ifndef CDT_VLC_INT_VECTOR_HPP
#define CDT_VLC_INT_VECTOR_HPP

#include <iostream>
#include <cstring>
#include <cassert>
#include <vector>
#include "integer_encoders.hpp"
#include "inv_bitstream.hpp"
#include "macros.h"

template<class arr_t>
class vlc_int_array_iterator{

    typedef vlc_int_array_iterator<arr_t>  self_type;
    typedef size_t                         value_type;
    typedef typename arr_t::coder_type     coder_type;

private:
    arr_t*           array=nullptr;
    size_t           bit_pos{};
    value_type       value{};
    size_t           len{};
    decode_t         code{};

    void move(self_type&& other) {
        std::swap(array, other.array);
        std::swap(bit_pos, other.bit_pos);
        std::swap(value, other.value);
        std::swap(len, other.len);
        std::swap(code, other.code);
    }

public:

    typedef value_type& reference;
    typedef value_type* pointer;
    typedef std::forward_iterator_tag iterator_category;

    vlc_int_array_iterator()=default;

    vlc_int_array_iterator(arr_t* array_, size_t index_) : array(array_), bit_pos(index_){
        if(bit_pos < array->m_written_bits){
            coder_type::decode(array->stream.stream, bit_pos, code);
            value = code.value-2;
            len = code.next_pos - bit_pos;
        }else{
            bit_pos = array->m_written_bits;
        }
    };

    vlc_int_array_iterator(vlc_int_array_iterator&& other) noexcept {
        move(std::forward<vlc_int_array_iterator>(other));
    };

    void swap(self_type& other){
        move(std::forward<self_type>(other));
    }

    inline vlc_int_array_iterator& operator++(){
        size_t tmp_index = bit_pos + len;
        if(tmp_index<array->m_written_bits){
            coder_type::decode(array->stream.stream, tmp_index, code);
            value = code.value-2;
            bit_pos = tmp_index;
            len = code.next_pos - bit_pos;
        }else{
            bit_pos = array->m_written_bits;
        }
        return *this;
    }

    inline vlc_int_array_iterator operator++(int){
        vlc_int_array_iterator old(*this);
        ++(*this);
        return old;
    }

    inline reference& operator*() {
        return value;
    }

    inline bool operator==(const self_type& rhs) const {
        return bit_pos == rhs.bit_pos;
    }

    inline bool operator!=(const self_type& rhs) const {
        return bit_pos != rhs.bit_pos;
    }

    vlc_int_array_iterator& operator=(self_type&& other) noexcept {
        move(std::forward<self_type>(other));
        return *this;
    }

};

template<size_t samp_den=128, class coder_t=elias_delta>
class vlc_int_array {

public:
    typedef typename coder_t::size_type size_type;
    typedef vlc_int_array_iterator<vlc_int_array<samp_den, coder_t>> iterator_t;
    typedef coder_t coder_type; //type of encoder storing the variable-length elements
    friend  iterator_t;

private:
    inv_bitstream<size_t> stream{nullptr, 0};//binary stream where the codes are stored
    size_type last_elm = 0;
    size_t m_size = 0; //size of the vector.
    size_t m_written_bits = 0; //sum of the code lengths
    // Every sample[i] stores the bit position of V[i*sample_den] in the binary stream
    std::vector<size_t> samples;

    void move(vlc_int_array<samp_den, coder_t> &&other) {
        stream.swap(other.data);
        std::swap(m_size, other.m_size);
        std::swap(m_written_bits, other.m_written_bits);
        std::swap(last_elm, other.last_elm);
        std::swap(samp_den, other.samp_den);
        std::swap(samples, other.samples);
    }

    void copy(const vlc_int_array<samp_den, coder_t> &other) {
        last_elm = other.last_val;
        m_size = other.m_size;
        m_written_bits = other.m_written_bits;

        if (stream.stream != other.stream.stream) {
            stream.stream_size = other.stream.stream_size;
            size_t n_bytes = stream.stream_size * sizeof(size_t);
            if (stream.stream != nullptr) {
                free(stream.stream);
            }
            stream.stream = reinterpret_cast<size_t *>(malloc(n_bytes));
            memcpy(stream.stream, other.stream.stream, n_bytes);
        }
    }

public:

    //default constructor
    explicit vlc_int_array() {
        stream.stream_size = 2;
        stream.stream = reinterpret_cast<size_t *>(malloc(sizeof(size_t) * 2));
    };

    //copy constructor
    vlc_int_array(const vlc_int_array<samp_den, coder_t> &other) {
        copy(other);
    };

    //move constructor
    vlc_int_array(vlc_int_array<samp_den, coder_t> &&other) noexcept {
        move(std::forward<vlc_int_array<samp_den, coder_t>>(other));
    };

    vlc_int_array<samp_den, coder_t> &operator=(const vlc_int_array<samp_den, coder_t> &other) {
        if (&other != this) {
            copy(other);
        }
        return *this;
    };

    vlc_int_array<samp_den, coder_t> &operator=(vlc_int_array<samp_den, coder_t> &&other) noexcept {
        if (&other != this) {
            move(std::forward<vlc_int_array<samp_den, coder_t>>(other));
        }
        return *this;
    };

    ~vlc_int_array() {
        if (stream.stream != nullptr) {
            free(stream.stream);
        }
    };

    void swap(vlc_int_array<samp_den, coder_t> &other) {
        move(std::move(other));
    }

    inline iterator_t begin() {
        return iterator_t(this, 0);
    };

    inline iterator_t end() {
        return iterator_t(this, m_written_bits);
    };

    inline size_type back() const {
        return last_elm;
    };

    [[nodiscard]] inline size_t operator[](size_t idx) const {
        assert(idx < size());
        size_t i = (idx/samp_den)*samp_den, pos = samples[idx/samp_den];
        decode_t res{};
        while (i < idx) {
            coder_type::decode(stream.stream, pos, res);
            pos = res.next_pos;
            i++;
        }
        coder_type::decode(stream.stream, pos, res);
        return res.value-2;
    }

    iterator_t begin(size_t idx) const {
        assert(idx<m_size);
        size_t cell = idx / samp_den;
        size_t i=cell*samp_den, pos=samples[cell];
        decode_t res{};

        while(i<idx){
            coder_type::decode(stream.stream, pos, res);
            pos = res.next_pos;
            i++;
        }
        return iterator_t(const_cast<vlc_int_array *>(this), pos);
    }

    void push_back(size_type value) {
        auto code = coder_t::encode(value+2);

        if (INT_CEIL((m_written_bits + code.second), inv_bitstream<size_t>::word_bits) >= stream.stream_size) {
            stream.stream_size *= 1.5;
            stream.stream = reinterpret_cast<size_t *>(realloc(stream.stream, stream.stream_size * sizeof(size_t)));
        }
        stream.write(m_written_bits, m_written_bits + code.second - 1, code.first);

        if(m_size % samp_den ==0){
            samples.push_back(m_written_bits);
        }

        m_size++;
        m_written_bits += code.second;
        last_elm = value;
    };

    [[nodiscard]] inline size_t size() const {
        return m_size;
    };

    void shrink_to_fit() {
        stream.stream_size = INT_CEIL(m_written_bits, inv_bitstream<size_t>::word_bits);
        stream.stream = reinterpret_cast<size_t *>(realloc(stream.stream, stream.stream_size * sizeof(size_t)));
        samples.shrink_to_fit();
    };

    void load(std::istream &in) {
        stream.load(in);
        load_elm(in, m_size);
        load_elm(in, m_written_bits);
        load_elm(in, last_elm);
        load_plain_vector(in, samples);
        size_t s_sp;
        load_elm(in, s_sp);
        assert(s_sp==samp_den);
    }

    size_t serialize(std::ostream &out) const {
        size_t written_bytes = stream.serialize(out);
        written_bytes += serialize_elm(out, m_size);
        written_bytes += serialize_elm(out, m_written_bits);
        written_bytes += serialize_elm(out, last_elm);
        written_bytes += serialize_plain_vector(out, samples);
        written_bytes += serialize_elm(out, samp_den);
        return written_bytes;
    }
};
#endif //CDT_VLC_INT_VECTOR_HPP