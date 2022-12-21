//
// Created by Diaz, Diego on 20.12.2022.
//

#ifndef GRLBWT_SIMPLE_HASH_TABLE_H
#define GRLBWT_SIMPLE_HASH_TABLE_H

#include <cstring>
#include <iostream>

template<class key_t,
         class value_t,
         uint8_t key_bytes,
         uint8_t val_bytes>
struct bucket_type{
    static_assert(key_bytes<=8 && val_bytes<=8);

    uint8_t data[key_bytes+val_bytes+1]={0};
    static constexpr size_t key_mask = (1UL <<(key_bytes*8UL))-1UL;
    static constexpr size_t value_mask = (1UL <<(val_bytes*8UL))-1UL;

    [[nodiscard]] inline key_t get_key() const {
        size_t tmp = *(std::reinterpret_pointer_cast<size_t*>(data));
        return (tmp & key_mask);
    }

    [[nodiscard]] inline value_t get_value() const {
        size_t tmp = *(std::reinterpret_pointer_cast<size_t*>(data[key_bytes] ));
        return (tmp & value_mask);
    }

    inline void init_bucket(key_t key, value_t val) {
        data[key_bytes+val_bytes] = true;
    }

    inline void modify_value(value_t new_value) {
        size_t * tmp = std::reinterpret_pointer_cast<size_t*>(data[key_bytes]);
        *tmp &= ~value_mask;
        *tmp |= new_value;
    }

    inline bool unused(){
        return !data[key_bytes+val_bytes];
    }

    inline bool mark_as_unused(){
        return data[key_bytes+val_bytes] = false;
    }
};

template<class key_t, class value_t, uint8_t key_bytes=sizeof(key_t), uint8_t val_bytes=sizeof(value_t)>
class simple_hash_table{
public:
    typedef bucket_type<key_t, value_t, key_bytes, val_bytes> bucket_t;

private:
    bucket_t * table = nullptr;
    size_t n_buckets{};
    float load_factor{};
    float max_load_factor{};
    size_t n_elms=0;

    const std::function<size_t(key_t&)>& hasher;
    const std::function<size_t(key_t&)>& key_comp;
    const std::function<size_t(size_t&)>& probe;

    simple_hash_table(size_t size,
                      std::function<size_t(key_t&)>& hasher_,
                      std::function<bool(key_t&)>& key_comp_,
                      float max_load_factor_=0.6): hasher(hasher_),
                                               key_comp(key_comp_),
                                               max_load_factor(max_load_factor_){
        resize_table(size);
    }

    std::pair<bucket_t *, bool> insert(key_t key, value_t value){

        size_t bucket = hasher(key) & n_buckets;
        bool new_insert=false;

        if(table[bucket].unused()){
            table[bucket].init_bucket(key, value);
            new_insert = true;
        }else{

            bool found = false;
            while(!table[bucket].unsued()){
                if(key_comp(key, table[bucket].get_key())){
                    found = true;
                    break;
                }else{
                    bucket = probe(bucket);
                }
            }

            if(!found){
                table[bucket].init_bucket(key, value);
                new_insert = true;
            }
        }

        if(new_insert){
            n_elms++;
            load_factor = double(n_elms)/double(n_buckets);
            if(load_factor>max_load_factor){
                bucket = rehash(bucket);
            }
        }

        return {&table[bucket], new_insert};
    }

    bucket_t * find(key_t key){
        size_t bucket = hasher(key) & n_buckets;
        while(!table[bucket].empty()){
            if(key_comp(key, table[bucket])){
                return &table[bucket];
            }
            bucket = probe(bucket);
        }
        return nullptr;
    }

    size_t rehash(size_t bucket){
        resize_table(n_buckets<<1UL);
        return 0;
    }

    void resize_table(size_t new_size){
        if(table==nullptr){
            table = (bucket_t *)malloc(new_size*sizeof(bucket_t));
            memset((char *)table, 0, new_size*sizeof(bucket_t));
        }else{
            table = realloc(table, new_size*sizeof(bucket_t));
        }
        n_buckets = new_size;
    }

    void reserve(size_t new_buckets){
        for(size_t i=0;i<n_buckets;i++){
            table[i].mark_as_unused();
        }
        resize_table(new_buckets);
    }

    ~simple_hash_table(){
        free(table);
    }
};
#endif //GRLBWT_SIMPLE_HASH_TABLE_H
