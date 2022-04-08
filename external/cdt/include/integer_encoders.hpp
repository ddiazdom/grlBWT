//
// Created by diediaz on 04-09-19.
//
#ifndef MULTI_BWT_INTEGER_ENCODERS_HPP
#define MULTI_BWT_INTEGER_ENCODERS_HPP

#include <climits>
#include <iostream>
#include <limits>
#include <cassert>

const size_t masks[65]={0x0,
                        0x1,0x3, 0x7,0xF,
                        0x1F,0x3F, 0x7F,0xFF,
                        0x1FF,0x3FF, 0x7FF,0xFFF,
                        0x1FFF,0x3FFF, 0x7FFF,0xFFFF,
                        0x1FFFF,0x3FFFF, 0x7FFFF,0xFFFFF,
                        0x1FFFFF,0x3FFFFF, 0x7FFFFF,0xFFFFFF,
                        0x1FFFFFF,0x3FFFFFF, 0x7FFFFFF,0xFFFFFFF,
                        0x1FFFFFFF,0x3FFFFFFF, 0x7FFFFFFF,0xFFFFFFFF,
                        0x1FFFFFFFF,0x3FFFFFFFF, 0x7FFFFFFFF,0xFFFFFFFFF,
                        0x1FFFFFFFFF,0x3FFFFFFFFF, 0x7FFFFFFFFF,0xFFFFFFFFFF,
                        0x1FFFFFFFFFF,0x3FFFFFFFFFF, 0x7FFFFFFFFFF,0xFFFFFFFFFFF,
                        0x1FFFFFFFFFFF,0x3FFFFFFFFFFF, 0x7FFFFFFFFFFF,0xFFFFFFFFFFFF,
                        0x1FFFFFFFFFFFF,0x3FFFFFFFFFFFF, 0x7FFFFFFFFFFFF,0xFFFFFFFFFFFFF,
                        0x1FFFFFFFFFFFFF,0x3FFFFFFFFFFFFF, 0x7FFFFFFFFFFFFF,0xFFFFFFFFFFFFFF,
                        0x1FFFFFFFFFFFFFF,0x3FFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFF,
                        0x1FFFFFFFFFFFFFFF,0x3FFFFFFFFFFFFFFF, 0x7FFFFFFFFFFFFFFF,0xFFFFFFFFFFFFFFFF};

struct decode_t{
    size_t value;
    size_t next_pos;
};

class elias_gamma{
private:
public:
    typedef size_t size_type;
    static constexpr uint8_t word_bits  = std::numeric_limits<size_type>::digits;
    static constexpr uint8_t word_shift = __builtin_ctz(word_bits);

    //! value : integer to be encoded
    //! return : a pair (cod, n_bit). The resulting code is stored in
    //! the first 'n_bit' bits of 'cod'
    static std::pair<size_type, size_type> encode(size_type value);
    //! bit_stream : bit stream storing the value to be decoded
    //! bit_index : position where the code starts in the stream
    //! returns a pair (value, new_index) where value is the decoded
    //! integer and new_index is the position in the stream immediately
    //! after the code
    static void decode(const size_type *bit_stream, size_t bit_index, decode_t& res);
};

class elias_delta{
public:
    typedef size_t size_type;
    static constexpr uint8_t word_bits  = std::numeric_limits<size_type>::digits;
    static constexpr uint8_t word_shift = __builtin_ctz(word_bits);
    //! value : integer to be encoded
    //! return : a pair (cod, n_bit). The resulting code is stored in
    //! the first 'n_bit' bits of 'cod'
    static std::pair<size_type, size_type> encode(size_type value);
    //! bit_stream : bit stream storing the value to be decoded
    //! bit_index : position where the code starts in the stream
    //! returns a pair (value, new_index) where value is the decoded
    //! integer and new_index is the position in the stream immediately
    //! after the code
    static void decode(const size_type *bit_stream, size_t bit_index, decode_t& res);
};

/*
template<size_t l>
class rice{
    typedef unsigned int size_type;
    static std::pair<size_type, size_type> encode(size_t value);
    static std::pair<size_type, size_type> decode(const size_type *bit_stream,
                                                  size_t bit_index);
};

class vbyte{
    static std::pair<size_t, size_t> encode(size_t value);
};

class simple9{
    static std::pair<size_t, size_t> encode(size_t value);
};

class pfordelta{
    static std::pair<size_t, size_t> encode(size_t value);
};*/

inline std::pair<elias_gamma::size_type, elias_gamma::size_type> elias_gamma::encode(size_type value) {
    size_t len = word_bits - __builtin_clzl(value);
    return {value,(2*len) -1};
}

inline void elias_gamma::decode(const size_type *bit_stream, size_t bit_index, decode_t& res) {

    size_t code_len, n_zeroes;
    size_t cell_start = bit_index >> word_shift;
    size_t bit_pos = bit_index & (word_bits-1UL);

    size_type tmp = bit_stream[cell_start] << bit_pos;

    if(tmp==0){
        size_t n_zeroes_next = __builtin_clzl(bit_stream[cell_start + 1]);
        n_zeroes = (word_bits - bit_pos) + n_zeroes_next;
        code_len = (n_zeroes*2)+1;
        tmp = bit_stream[cell_start + 1] >> (word_bits - (n_zeroes_next + n_zeroes + 1U));

    }else{
        n_zeroes = __builtin_clzl(tmp);
        code_len = (n_zeroes*2)+1;
        size_t cell_end = (bit_index + code_len-1) >> word_shift;

        if(cell_start==cell_end){
            tmp >>= word_bits - code_len;
        }else{
            size_t left = (word_bits-(bit_pos+n_zeroes));
            size_t right = n_zeroes+1 - left;
            tmp = (bit_stream[cell_start] & masks[left]) << right;
            tmp |= bit_stream[cell_end]>>(word_bits-right);
        }
    }

    res.value = tmp;
    res.next_pos = bit_index + code_len;
}

inline std::pair<elias_delta::size_type , elias_delta::size_type > elias_delta::encode(size_type value) {
    elias_gamma::size_type len = word_bits - __builtin_clzl(value);
    auto g_code = elias_gamma::encode(len);
    assert((len+g_code.second-1)<=64);
    return {value + ((1UL<<(len-1))*(g_code.first-1)),len + g_code.second-1};
}

inline void elias_delta::decode(const size_type *bit_stream, size_t bit_index, decode_t& res) {

    elias_gamma::decode(bit_stream, bit_index, res);

    if(res.value==1){
        res.next_pos = bit_index+1;
        return;
    }

    size_t word_pos, bit_pos;
    size_type tmp;

    bit_index = res.next_pos;
    word_pos = bit_index >> word_shift;
    bit_pos = bit_index & (word_bits-1UL);

    if((bit_pos+res.value-1) > word_bits){

        size_t left = word_bits - bit_pos;
        size_t right = res.value-left-1;

        tmp = (bit_stream[word_pos] & masks[left])<<right;
        tmp |= bit_stream[word_pos+1] >> (word_bits - right);
        tmp |= 1UL << (res.value-1UL);

    }else{
        tmp = bit_stream[word_pos] >> (word_bits - (bit_pos + res.value - 1));
        tmp &= masks[res.value-1];
        tmp |= 1UL << (res.value-1UL);
    }

    res.next_pos = bit_index+res.value-1;
    res.value = tmp;
}
#endif //MULTI_BWT_INTEGER_ENCODERS_HPP
