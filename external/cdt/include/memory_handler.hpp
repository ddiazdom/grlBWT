//
// Created by Diego Diaz on 2/10/20.
//

#ifndef LMS_COMPRESSOR_MEMORY_HANDLER_HPP
#define LMS_COMPRESSOR_MEMORY_HANDLER_HPP

#endif //LMS_COMPRESSOR_MEMORY_HANDLER_HPP

class allocator{
public:
    template<class value_t>
    inline static value_t *allocate(size_t n_vals, bool init = true, value_t def_val = 0) {
        auto all_mem = (value_t *)malloc(sizeof(value_t) * n_vals);
        if(init){
            for(size_t i=0; i < n_vals;i++) all_mem[i] = def_val;
        }
        return all_mem;
    }

    template<class value_t>
    static value_t *
    reallocate(value_t *pointer, size_t old_size, size_t new_size, bool init = true, value_t def_val = 0) {
        auto re_all_mem = (value_t*)realloc(pointer, sizeof(value_t)*new_size);
        if(new_size>old_size && init){
            for(size_t i=old_size;i<new_size;i++){
                re_all_mem[i] = def_val;
            }
        }
        return re_all_mem;
    }

    template<class value_t>
    static void deallocate(value_t * pointer){
        if(pointer!= nullptr){
            free(pointer);
        }
    }
};
