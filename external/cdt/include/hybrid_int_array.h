//
// Created by Diaz, Diego on 14.3.2022.
//

#ifndef CDT_HYBRID_INT_ARRAY_H
#define CDT_HYBRID_INT_ARRAY_H

#include "hash_table.hpp"
#ifdef __ARM_NEON__
#include "sse2neon.h"
#else
#include <immintrin.h>
#endif

template<class arr_t>
class hyb_int_array_iterator{

    typedef hyb_int_array_iterator<arr_t> self_type;
    typedef size_t                        value_type;

private:
    const arr_t * array{};
    size_t        index{};
    value_type    value{};

    void move(self_type&& other) {
        std::swap(array, other.array);
        std::swap(index, other.index);
        std::swap(value, other.value);
    }

    void copy(const self_type& other) {
        array = other.array;
        index = other.index;
        value = other.value;
    }

public:

    typedef value_type&                reference;
    typedef value_type*                pointer;
    typedef std::forward_iterator_tag  iterator_category;

    hyb_int_array_iterator()=default;

    hyb_int_array_iterator(const arr_t * array_, size_t index_) : array(array_), index(index_){
        if(index < array->size()){
            value = (*array)[index];
        }
    };

    hyb_int_array_iterator(hyb_int_array_iterator&& other) noexcept {
        move(std::forward<hyb_int_array_iterator>(other));
    };

    hyb_int_array_iterator(const hyb_int_array_iterator& other) noexcept {
        copy(other);
    };

    void swap(hyb_int_array_iterator& other){
        move(std::forward<hyb_int_array_iterator>(other));
    }

    inline hyb_int_array_iterator& operator++(){
        if(index<(array->size()-1)){
            index++;
            value = (*array)[index];
        }else{
            index = array->size();
        }
        return *this;
    }

    inline hyb_int_array_iterator& operator--(){
        assert(index>0);
        index--;
        value = (*array)[index];
        return *this;
    }

    inline hyb_int_array_iterator operator++(int){
        hyb_int_array_iterator old(*this);
        ++(*this);
        return old;
    }

    inline size_t operator*() const {
        return value;
    }

    inline bool operator==(const self_type& other) const {
        return index == other.index;
    }

    inline bool operator!=(const self_type& other) const {
        return index != other.index;
    }

    hyb_int_array_iterator& operator=(self_type&& other) noexcept {
        move(std::forward<self_type>(other));
        return *this;
    }
};

template<typename size_type,
         typename hash_table_type=bit_hash_table<size_t, sizeof(size_t)*8, size_t, 32, true>>
class hybrid_int_array {

private:
    std::vector<size_type> m_vector;
    hash_table_type        m_ht;
    const size_t           max_val = std::numeric_limits<size_type>::max();

    void move(hybrid_int_array<size_type>&& other) noexcept {
        m_vector.swap(other.m_vector);
        m_ht.swap(other.m_ht);
    }

    inline void copy(const hybrid_int_array<size_type> &other) {
        m_vector = other.m_vector;
        m_ht = other.m_ht;
    }

    /*
    inline uint16_t sum_32_neon(const uint8_t a[32]) const {
        uint8x16_t acc = vdupq_n_u8(0); // clear accumulators
        uint64x2_t sum0 = _mm_sad_epu8( acc, vld1q_u8(a));
        uint64x2_t sum1 = _mm_sad_epu8( acc, vld1q_u8(a+16));
        sum0 = vaddq_u64(sum0, sum1);
        uint32x2_t vec64a = vget_low_u32(sum0);
        uint32x2_t vec64b = vget_high_u32(sum0);
        return vget_lane_u32(vadd_u32(vec64a, vec64b), 0); // extract lanes and
    }

    inline uint16_t sum_32_intel(const uint8_t a[32]) const {
        __m128i zero = _mm_xor_si128(zero, zero);
        __m128i sum0 = _mm_sad_epu8( zero, _mm_load_si128((__m128i*)a));
        __m128i sum1 = _mm_sad_epu8( zero, _mm_load_si128((__m128i*)(a+16)));
        __m128i sum2 = _mm_add_epi16(sum0, sum1);
        __m128i totalsum = _mm_add_epi16(sum2, _mm_shuffle_epi32(sum2, 2));
        return _mm_extract_epi16(totalsum, 0);
    }

    inline uint16_t sum_16_neon(const uint8_t a[16]) const {
        uint8x16_t acc = vdupq_n_u8(0); // clear accumulators
        uint64x2_t sum0 = _mm_sad_epu8( acc, vld1q_u8(a));
        uint32x2_t vec64a = vget_low_u32(sum0);
        uint32x2_t vec64b = vget_high_u32(sum0);
        return vget_lane_u32(vadd_u32(vec64a, vec64b), 0); // extract lanes and
    }

    inline uint16_t sum_16_intel(const uint8_t a[32]) const {
        __m128i zero = _mm_xor_si128(zero, zero);
        __m128i sum0 = _mm_sad_epu8( zero, _mm_load_si128((__m128i*)a));
        __m128i totalsum = _mm_add_epi16(sum0, _mm_shuffle_epi32(sum0, 2));
        return _mm_extract_epi16(totalsum, 0);
    }*/

public:

    const std::vector<size_type>& container = m_vector;

    typedef hyb_int_array_iterator<hybrid_int_array<size_type, hash_table_type>> iterator_t;

    hybrid_int_array()=default;

    hybrid_int_array(size_t n_elms, size_t def_val){
        assert(def_val<=max_val);
        m_vector = std::vector<size_type>(n_elms, def_val);
    }

    [[nodiscard]] inline size_t size() const{
        return m_vector.size();
    }

    [[nodiscard]] inline bool empty() const{
        return m_vector.empty();
    }

    [[nodiscard]] inline size_t back() const{
        assert(!m_vector.empty());
        return read(m_vector.size()-1);
    }

    [[nodiscard]] inline size_t read(size_t idx) const{
        assert(idx<m_vector.size());
        size_t val = m_vector[idx];

        if(val!=0) return val-1;

        auto res = m_ht.find(&idx, 64);
        assert(res.second);
        m_ht.get_value_from(res.first, val);
        return val-1;
    }

    inline size_t operator[](size_t idx) const{
        return read(idx);
    }

    void push_back(size_t val){
        if(val>=max_val){
            size_t idx = m_vector.size();
            m_ht.insert(&idx, 64, val+1);
            m_vector.push_back(0);
        }else{
            m_vector.push_back(val+1);
        }
    }

    inline void write(size_t idx, size_t val) {
        assert(idx<m_vector.size());
        if(val>=max_val){
            m_ht.insert(&idx, 64, val+1);
            m_vector[idx] = 0;
        }else{
           m_vector[idx] = val+1;
        }
    }

    iterator_t begin() {
        return iterator_t(this, 0);
    };

    iterator_t begin(size_t idx) const {
        return iterator_t(this, idx);
    };

    iterator_t end() {
        return iterator_t(this, m_vector.size());
    };

    //move constructor
    hybrid_int_array(hybrid_int_array&& other) noexcept {
        move(std::forward<hybrid_int_array<size_type>>(other));
    };

    //copy constructor
    hybrid_int_array(const hybrid_int_array& other) noexcept {
        copy(other);
    };

    //copy assignment operator
    hybrid_int_array& operator=(const hybrid_int_array& other){
        if(this!=&other){
            copy(other);
        }
        return *this;
    }

    //move assignment operator
    hybrid_int_array& operator=(hybrid_int_array && other) noexcept{
        if(this!=&other){
            move(std::forward<hybrid_int_array<size_type>>(other));
        }
        return *this;
    }

    size_type serialize(std::ostream &out) const{
        size_t written_bytes =serialize_plain_vector(out, m_vector);
        written_bytes += m_ht.serialize(out);
        return written_bytes;
    }

    void load(std::istream &in){
        load_plain_vector(in, m_vector);
        m_ht.load(in);
    }

    void shrink_to_fit(){
        m_vector.shrink_to_fit();
        m_ht.shrink_databuff();
    }

    //number of values exceeding the max_value threshold (i.e., number of values with overflow)
    [[nodiscard]] inline size_t overflow() const {
        return m_ht.size();
    }

    void swap(hybrid_int_array& other){
        move(std::forward<hybrid_int_array>(other));
    }

    [[nodiscard]] inline size_t vec_sum_red(size_t l, size_t r, const uint8_t* ptr) const {

        __m128i zero = _mm_setzero_si128();
        __m128i sum0 = _mm_sad_epu8( zero, _mm_load_si128(reinterpret_cast<const __m128i*>(ptr+l)));
        l+=16;
        while(l<=r){
            __m128i sum1 = _mm_sad_epu8( zero, _mm_load_si128(reinterpret_cast<const __m128i*>(ptr+l)));
            sum0 = _mm_add_epi32(sum0, sum1);
            l+=16;
        }

        __m128i totalsum = _mm_add_epi32(sum0, _mm_shuffle_epi32(sum0, 2));
        size_t res = _mm_cvtsi128_si32(totalsum);

        r++;
        while(r<l) res-=ptr[r++];

        return res;
    }

    [[nodiscard]] size_t test_scalar_sum_red(size_t l, size_t r, const uint8_t* ptr) const {
        size_t acc=0;
        for(size_t i=l;i<=r;i++){
            acc+=ptr[i];
        }
        return acc;
    }

    [[nodiscard]] size_t scl_sum_red(size_t l, size_t r) const {
        size_t acc=0;
        for(size_t i=l;i<=r;i++){
            acc+=m_vector[i];
            if(m_vector[i]==0){
                auto res = m_ht.find(&i, 64);
                size_t val=0;
                m_ht.get_value_from(res.first, val);
                acc+=val;
            }
        }
        return acc-(r-l+1);
    }
};

#endif //CDT_HYBRID_INT_ARRAY_H
