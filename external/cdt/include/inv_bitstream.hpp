//
// Created by diediaz on 03-12-19.
//

#ifndef CDT_BIT_STREAM_HPP
#define CDT_BIT_STREAM_HPP
#include <iostream>
#include <limits>
#include "cdt_common.hpp"

template<class word_t>
struct inv_bitstream  {
    constexpr static uint8_t word_bits = std::numeric_limits<word_t>::digits;
    constexpr static uint8_t word_shift = __builtin_ctz(word_bits);
    const static size_t masks[65];
    using size_type = size_t;
    word_t *stream= nullptr;
    size_t stream_size{};
public:

    [[nodiscard]] inline size_t read(size_t i, size_t j) const {
        size_t cell_i = i >> word_shift;
        size_t cell_j = j >> word_shift;
        size_t len = j-i+1UL;

        if(cell_i == cell_j){
            size_t offset = word_bits - (j & (word_bits - 1UL)) -1UL;
            return (stream[cell_j] >> offset) & masks[len];
        }else{
            size_t l_off = word_bits - (i & (word_bits - 1UL));
            size_t r_off = len-l_off;
            return (((stream[cell_i] & masks[l_off]) << r_off) | (stream[cell_j] >> (word_bits - r_off)));
        }
    }

    inline void write(size_t i, size_t j, size_t value){

        size_t cell_i = i >> word_shift;
        size_t cell_j = j >> word_shift;
        size_t len = j-i+1UL;

        if(cell_i==cell_j){
            size_t offset = word_bits - (j & (word_bits - 1UL)) - 1UL;
            stream[cell_j] &= ~(masks[len] << offset);
            stream[cell_j] |= value << offset;
        }else{
            size_t l_off = word_bits - (i & (word_bits - 1UL));
            size_t r_off = len-l_off;
            stream[cell_i] = (stream[cell_i] & ~masks[l_off]) | (value >> r_off);
            stream[cell_j] = (stream[cell_j] & ~(masks[r_off]<<(word_bits-r_off))) | (value << (word_bits-r_off));
        }
    }

    size_type serialize(std::ostream &out) const{
        size_t written_bytes = serialize_elm(out, stream_size);
        out.write((char *)stream, sizeof(word_t)*stream_size);
        return written_bytes + (sizeof(word_t)*stream_size);
    }

    void load(std::istream &in){
        load_elm(in, stream_size);
        stream = reinterpret_cast<word_t *>(malloc(sizeof(word_t)*stream_size));
        in.read((char *)stream, sizeof(word_t)*stream_size);
    }

    inline inv_bitstream<word_t>& swap(inv_bitstream<word_t>& other) {
        std::swap(stream, other.stream);
        std::swap(stream_size, other.stream_size);
        return *this;
    }
};

template<class word_t>
const size_t inv_bitstream<word_t>::masks[65]={0x0,
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

#endif //CDT_BIT_STREAM_HPP
