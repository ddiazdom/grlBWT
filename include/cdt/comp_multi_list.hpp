//
// Created by diego on 20-04-20.
//

#ifndef LPG_COMPRESSOR_COMP_MULTI_LIST_HPP
#define LPG_COMPRESSOR_COMP_MULTI_LIST_HPP

#include <sdsl/bit_vectors.hpp>

template<class arr_t>
class comp_multi_list {

public:
    typedef sdsl::bit_vector                            bit_vector_t;
    typedef sdsl::int_vector_buffer<>                   buffer_t;
    typedef sdsl::int_vector<>                          vector_t;
    typedef size_t                                      size_type;

    bit_vector_t                                        has_list_bv;
    bit_vector_t::rank_1_type                           has_list_bv_rs;

    bit_vector_t                                        list_limits;
    bit_vector_t::select_1_type                         list_limits_ss;
    arr_t                                               arr;

    comp_multi_list()=default;

    //elements are placed either backward or forward in every list
    comp_multi_list(buffer_t& _elm_buffer, bit_vector_t& universe, bool rev_placement){
        auto tmp_arr = build_lists(_elm_buffer, universe, rev_placement);
        arr.swap(tmp_arr);
        _elm_buffer.close(true);
    }

    //return the {start,len} for a list
    inline std::pair<size_t, size_t> get_list_locus(size_t idx) const{
        assert(has_list_bv[idx]);
        size_t list = has_list_bv_rs(idx) + 1;
        size_t start = list_limits_ss(list);
        size_t end = list_limits_ss(list + 1) - 1;
        return {start, end-start+1};
    }

    inline bool empty_list(const size_t idx) const{
        return !has_list_bv[idx];
    }

    inline size_t tot_lists() const{
        return has_list_bv_rs(has_list_bv.size());
    };

    //maps the idx to a list in the suffix pointers
    inline size_t idx2list(const size_t idx) const{
        assert(has_list_bv[idx]);
        return has_list_bv_rs(idx);
    };

    vector_t build_lists(buffer_t& _elm_buffer, bit_vector_t& _has_list_bv, bool rev_placement){

        //buffer is (idx, value)
        has_list_bv.swap(_has_list_bv);

        sdsl::util::init_support(has_list_bv_rs, &has_list_bv);

        vector_t list_counts(tot_lists(), 0);

        auto s_it = _elm_buffer.begin();
        auto s_end = _elm_buffer.end();
        while (s_it != s_end) {
            list_counts[idx2list(*s_it)]++;
            s_it += 2;
        }
        sdsl::util::bit_compress(list_counts);

        vector_t _tmp_array((_elm_buffer.size() / 2), 0);
        bit_vector_t list_limits_tmp(bit_vector_t(_tmp_array.size() + 1, false));
        list_limits.swap(list_limits_tmp);

        size_t pos = 0, c_pos = 0;
        while (pos < _tmp_array.size()) {
            list_limits[pos] = true;
            pos += list_counts[c_pos];
            list_counts[c_pos++] = 0;
        }
        list_limits[pos] = true;
        sdsl::util::init_support(list_limits_ss, &list_limits);

        auto s_it2 = _elm_buffer.begin();
        size_t list;
        if(rev_placement){
            size_t l_end;
            while (s_it2 != s_end) {
                list = idx2list(*s_it2);
                l_end = list_limits_ss(list + 2)-1;
                s_it2++;
                _tmp_array[l_end - list_counts[list]++] = *s_it2;
                s_it2++;
            }
        }else{
            size_t l_start;
            while (s_it2 != s_end) {
                list = idx2list(*s_it2);
                l_start = list_limits_ss(list + 1);
                s_it2++;
                _tmp_array[l_start + list_counts[list]++] = *s_it2;
                s_it2++;
            }
        }
        return _tmp_array;
    }

    void print_lists(){
        for(size_t i=0;i<has_list_bv.size();i++){
            if(has_list_bv[i]){
                std::cout<<(i+1)<<" : ";
                auto locus =  get_list_locus(i);
                for(size_t j=locus.first;j<(locus.first+locus.second);j++){
                    std::cout<<arr[j]<<" ";
                }
                std::cout<<" "<<std::endl;
            }
        }
    }

    comp_multi_list<arr_t>& operator=(comp_multi_list<arr_t>&& other) noexcept {
        has_list_bv.swap(other.has_list_bv);
        has_list_bv_rs.swap(other.has_list_bv_rs);
        list_limits.swap(other.list_limits);
        list_limits_ss.swap(other.list_limits_ss);
        arr.swap(other.arr);

        has_list_bv_rs.set_vector(&has_list_bv);
        list_limits_ss.set_vector(&list_limits);
    }

    void clear(){
        sdsl::memory_manager::clear(has_list_bv);
        sdsl::memory_manager::clear(list_limits);
        sdsl::memory_manager::clear(arr);
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += has_list_bv.serialize(out, child, "has_list_bv");
        written_bytes += has_list_bv_rs.serialize(out, child, "has_list_bv_rs");
        written_bytes += list_limits.serialize(out, child, "list_limits");
        written_bytes += list_limits_ss.serialize(out, child, "list_limits_ss");
        written_bytes += arr.serialize(out, child, "arr");
        return written_bytes;
    }

    void load(std::istream &in){
        has_list_bv.load(in);
        has_list_bv_rs.load(in);
        list_limits.load(in);
        list_limits_ss.load(in);
        arr.load(in);

        has_list_bv_rs.set_vector(&has_list_bv);
        list_limits_ss.set_vector(&list_limits);
    }
};

template<class arr_t>
class delta_multi_list: public comp_multi_list<arr_t>{
    typedef size_t                                      size_type;
    typedef sdsl::bit_vector                            bit_vector_t;
    typedef sdsl::int_vector_buffer<>                   buffer_t;

public:
    delta_multi_list(): comp_multi_list<arr_t>(){};

    delta_multi_list(buffer_t& _elm_buffer, bit_vector_t& universe, bool rev_placement): comp_multi_list<arr_t>(){
        auto tmp_arr = (*this).build_lists(_elm_buffer, universe, rev_placement);
        _elm_buffer.close(true);
        size_t tmp1, tmp2=0;

        for (size_t i = 0; i < tmp_arr.size(); i++) {
            if((*this).list_limits[i]){
                tmp2 = tmp_arr[i];
            }else{
                tmp1 = tmp_arr[i];
                tmp_arr[i] -=tmp2;
                tmp2 = tmp1;
            }
        }
        sdsl::util::assign((*this).arr, arr_t(tmp_arr));
    }

    delta_multi_list<arr_t>& operator=(delta_multi_list<arr_t>&& other) noexcept {
        (*this).has_list_bv.swap(other.has_list_bv);
        (*this).has_list_bv_rs.swap(other.has_list_bv_rs);
        (*this).list_limits.swap(other.list_limits);
        (*this).list_limits_ss.swap(other.list_limits_ss);
        (*this).arr.swap(other.arr);
        (*this).has_list_bv_rs.set_vector(&(*this).has_list_bv);
        (*this).list_limits_ss.set_vector(&(*this).list_limits);
        return *this;
    }
};


#endif //LPG_COMPRESSOR_COMP_MULTI_LIST_HPP
