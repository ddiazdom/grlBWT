//
// Created by diegodiaz on 12-10-20.
//

#ifndef LPG_COMPRESSOR_SI_INT_ARRAY_H
#define LPG_COMPRESSOR_SI_INT_ARRAY_H

#include <vector>
#include "integer_encoders.hpp"
#include "inv_bitstream.hpp"
#include "vlc_int_array.h"

//forward iterator for the class si_int_array
template<class arr_t, class arr_it_t>
class si_int_array_iterator{

    typedef si_int_array_iterator<arr_t, arr_it_t> self_type;
    typedef size_t                                 value_type;

private:
    arr_t*     array= nullptr;
    size_t     value={};
    size_t     idx={};
    arr_it_t   m_iterator;

    void move(self_type&& other) {
        std::swap(array, other.array);
        std::swap(idx, other.idx);
        m_iterator.swap(other.m_iterator);
    }

    void copy(self_type&& other) {
        array = other.array;
        idx = other.index;
        m_iterator = other.m_iterator;
    }

public:

    typedef value_type& reference;
    typedef value_type* pointer;
    typedef std::forward_iterator_tag iterator_category;

   explicit si_int_array_iterator() = default;

   si_int_array_iterator(arr_t& array_, size_t idx_) : array(&array_), idx(idx_) {
       if(idx<array->size()){
           /*size_t cell = idx / array->m_samp_den;
           size_t i=cell*array->m_samp_den;
           value=array->samples[cell];
           m_iterator = array->m_vector.begin(i);
           while(i<idx){
               value += *m_iterator;
               i++;
               ++m_iterator;
           }
           value +=*m_iterator;*/
           value = (*array)[idx];
       }else{
           idx = array->size();
       }
   };

    si_int_array_iterator(si_int_array_iterator&& other)  noexcept {
        move(std::forward<si_int_array_iterator>(other));
    };

    si_int_array_iterator(si_int_array_iterator& other)  noexcept {
        copy(other);
    };

    si_int_array_iterator& operator=(const si_int_array_iterator& other) {
        if(&other!=this){
            copy(other);
        }
        return *this;
    }

    si_int_array_iterator& operator=(si_int_array_iterator&& other)  noexcept {
        if(&other!=this){
            move(std::forward<si_int_array_iterator>(other));
        }
        return *this;
    }

    inline si_int_array_iterator& operator++(){
        /*if(idx<array->size()){
            ++m_iterator;
            value +=*m_iterator;
            idx++;
        }else{
            idx = array->size();
        }*/
        if(idx<(array->size()-1)){
            value += array->m_vector[++idx];
        }else{
            idx = array->size();
        }
        return *this;
    }

    inline si_int_array_iterator operator++(int){
        si_int_array_iterator old(*this);
        ++(*this);
        return old;
    }

    inline size_t operator*() const {
        return value;
    }

    inline bool operator==(const self_type& other) const {
        return idx == other.idx;
    }

    inline bool operator!=(const self_type& other) const {
        return idx != other.idx;
    }
};

//vector V storing a list of strictly increasing integer values.
// The element V[i], with i>0, is encoded as a delta V[i]-V[i-1].
// coder_t: the underlying encoder to store the deltas.
template<size_t samp_den=128,
        class vector_type=vlc_int_array<samp_den>>
class si_int_array{

private:
    size_t              last_elm={};
    size_t m_samp_den = samp_den;
    static constexpr size_t width = __builtin_ctz(samp_den);
    vector_type         m_vector;
    std::vector<size_t> samples;


    void move(si_int_array&& other){
        m_vector.swap(other.m_vector);
        samples.swap(other.samples);
        std::swap(last_elm, other.last_elm);
        std::swap(m_samp_den, other.m_samp_den);
        assert(m_samp_den==samp_den);
    }

    void copy(si_int_array& other){
        samples = other.samples;
        m_vector = other.m_vector;
        last_elm = other.last_elm;
        m_samp_den = other.samp_den;
        assert(m_samp_den==samp_den);
    }


public:
    typedef si_int_array_iterator<si_int_array<samp_den, vector_type>, typename vector_type::iterator_t> iterator_t;
    typedef size_t                                                                                       size_type;
    friend                                                                                               iterator_t;

    explicit si_int_array()=default;

    si_int_array(si_int_array&& other) noexcept{
        move(std::forward<si_int_array>(other));
    }

    si_int_array(si_int_array& other){
        copy(other);
    }

    iterator_t begin(){
        return iterator_t(*this, 0);
    }

    iterator_t begin(size_t idx) {
        assert(idx<m_vector.size());
        return iterator_t(*this, idx);
    }

    iterator_t end() {
        return iterator_t(*this, m_vector.size());
    }

    si_int_array<samp_den, vector_type>& operator=(const si_int_array<samp_den, vector_type>& other) {
        if(other!=*this){
            copy(other);
        }
        return *this;
    };

    si_int_array<samp_den, vector_type>& operator=(si_int_array<samp_den, vector_type>&& other) noexcept {
        if(other!=*this){
            move(std::forward<si_int_array<samp_den, vector_type>>(other));
        }
        return *this;
    };

    inline si_int_array& push_back(size_t value){
        size_t elm;
        if(m_vector.size()==0){
            elm = value;
        }else{
            assert(value>last_elm);
            elm = value-last_elm;
        }

        if(m_vector.size() % m_samp_den == 0){
            samples.push_back(last_elm);
        }
        m_vector.push_back(elm);
        last_elm = value;
        return *this;
    }

    [[nodiscard]] inline size_t operator[](size_t idx) const {
        assert(idx<m_vector.size());

        //auto start = std::chrono::steady_clock::now();
        //m_vector.vec_sum_red(((idx>>width)<<width), idx);
        //auto end = std::chrono::steady_clock::now();
        //std::cout << "Elapsed time in nanoseconds: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns" << std::endl;
        //start = std::chrono::steady_clock::now();
        //m_vector.sum_red_scalar(cell*m_samp_den, idx);
        //end = std::chrono::steady_clock::now();
        //std::cout << "Elapsed time in nanoseconds: " << std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count() << " ns\n" << std::endl;

        return samples[(idx>>width)] + m_vector.scl_sum_red(((idx>>width)<<width), idx);
    }

    [[nodiscard]] inline size_t size() const{
        return m_vector.size();
    }

    [[nodiscard]] inline size_t empty() const{
        return m_vector.empty();
    }

    [[nodiscard]] inline size_t back() const{
        return last_elm;
    }

    inline void shrink_to_fit(){
        m_vector.shrink_to_fit();
        samples.shrink_to_fit();
    }

    void load(std::istream &in){
        load_elm(in, last_elm);
        load_elm(in, m_samp_den);
        m_vector.load(in);
        load_plain_vector(in, samples);
    }

    size_t serialize(std::ostream &out) const {
        size_t written_bytes = serialize_elm(out, last_elm);
        written_bytes += serialize_elm(out, m_samp_den);
        written_bytes += m_vector.serialize(out);
        written_bytes += serialize_plain_vector(out, samples);
        return written_bytes;
    }

    void swap(si_int_array<samp_den, vector_type>& other){
        move(std::move(other));
    }
};

#endif //LPG_COMPRESSOR_SI_INT_ARRAY_H
