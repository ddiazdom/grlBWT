//
// Created by diegodiaz on 12-10-20.
//

#ifndef LPG_COMPRESSOR_SI_INT_ARRAY_H
#define LPG_COMPRESSOR_SI_INT_ARRAY_H
#include "vlc_int_array.hpp"

template<class arr_t>
class si_int_array_iterator{
private:
    typedef si_int_array_iterator<arr_t> self_type;
    typename arr_t::iterator_t           iter;
    size_t                               idx;
    arr_t&                               arr;
    size_t                               value;

public:
    typedef size_t      value_type;
    typedef value_type& reference;
    typedef value_type* pointer;
    typedef std::forward_iterator_tag iterator_category;

    si_int_array_iterator(typename arr_t::iterator_t iter_, arr_t& arr_, size_t idx_, size_t value_) : iter(std::move(iter_)),
                                                                                                       idx(idx_),
                                                                                                       arr(arr_),
                                                                                                       value(value_){};

    inline self_type& operator++(){
        if(idx<arr.size()){
            ++iter;
            value+=*iter;
            idx++;
        }
        return *this;
    }

    inline bool operator==(const self_type& rhs) const {
        return idx == rhs.idx;
    }

    inline bool operator!=(const self_type& rhs) const {
        return idx != rhs.idx;
    }

    inline reference& operator*() {
        return value;
    }
};

class si_int_array{

private:
    typedef vlc_int_array<elias_delta> arr_type;
    size_t                             last_elm{};
    arr_type                           data;

    void move(si_int_array&& other){
        std::swap(last_elm, other.last_elm);
        data.swap(other.data);
    }

    void copy(si_int_array& other){
        last_elm = other.last_elm;
        data = other.data;
    }
public:
    typedef si_int_array_iterator<arr_type> iterator_t;
    typedef size_t                          size_type;
    friend                                  iterator_t;

    si_int_array()=default;

    si_int_array(si_int_array&& other) noexcept{
        move(std::forward<si_int_array>(other));
    }

    si_int_array(si_int_array& other){
        copy(other);
    }

    iterator_t begin(){
        return iterator_t(data.begin(), data, 0, data.read(0));
    }

    iterator_t end(){
        return iterator_t(data.end(), data, data.size(), 0);
    }

    inline si_int_array& push_back(size_t value){
        if(empty()){
            data.push_back(value);
        }else{
            assert(value>last_elm);
            data.push_back(value-last_elm);
        }
        last_elm = value;
        return *this;
    }

    inline size_t size() const{
        return data.size();
    }

    inline size_t empty() const{
        return data.size()==0;
    }

    inline size_t back() const{
        return last_elm;
    }

    inline void shrink_to_fit(){
        data.shrink_to_fit();
    }

    inline size_t buffer_bytes() const{
        return data.buffer_bytes();
    }
};

#endif //LPG_COMPRESSOR_SI_INT_ARRAY_H
