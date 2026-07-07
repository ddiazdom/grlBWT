//
// Created by Diaz, Diego on 23.11.2021.
//

#ifndef GRLBWT_COMMON_H
#define GRLBWT_COMMON_H

#include <vector>
#include <string>
#include <sdsl/bit_vectors.hpp>
#include "hash_table.hpp"
#include "int_array.h"
#include "rank_support.h"
#include "file_streams.hpp"

#define BUFFER_SIZE 8388608 //8MB of buffer

typedef sdsl::bit_vector                        bv_t;
typedef sdsl::bit_vector::rank_1_type           bv_rs_t;
typedef int_array<size_t>                       vector_t;
typedef int_array<size_t>                       string_t;
typedef hash_table<size_t, 44>                  phrase_map_t;
typedef buffered_hash_table<size_t, 44>         buffered_map_t;
typedef typename phrase_map_t::buff_t           ht_buff_t;

// the phrases are stored in a bit-compressed hash table:
// this wrapper reinterprets the bits back as phrases
struct key_wrapper{
    size_t width;
    size_t d_bits;//bits used to describe the string
    const bitstream<ht_buff_t>& stream;

    //offset is the bit where the key description starts
    [[nodiscard]] inline size_t read(size_t offset, size_t idx) const {
        return stream.read(offset + d_bits + (idx * width),
                           offset + d_bits + ((idx + 1) * width - 1));
    }

    //offset is the bit where the key description starts
    [[nodiscard]] inline size_t size(size_t offset) const{
        return stream.read(offset, offset + d_bits - 1) / width;
    }

    //compare two phrases are different positions of the bit stream
    [[nodiscard]] inline bool compare(size_t a, size_t b) const{

        size_t a_bits = stream.read(a, a + d_bits - 1);
        size_t b_bits = stream.read(b, b + d_bits - 1);
        size_t min_bits = std::min<size_t>(a_bits, b_bits);

        size_t a_pos = a+d_bits+a_bits-1;
        size_t b_pos = b+d_bits+b_bits-1;
        size_t rm_diff = stream.inv_com_segments(a_pos, b_pos, min_bits);

        if(rm_diff < min_bits){
            a_pos = a+d_bits+(((a_bits - rm_diff-1) / width) * width);
            b_pos = b+d_bits+(((b_bits - rm_diff-1) / width) * width);
            size_t sym_a = stream.read(a_pos, a_pos+width-1);
            size_t sym_b = stream.read(b_pos, b_pos+width-1);
            return sym_a<sym_b;
        }else{
            return a_bits>b_bits;
        }
    }
};

// The global phrase table split into n_parts independent sub-tables by a hash of
// the phrase. A phrase always lands in the same part, so the parts are disjoint
// and can be filled in parallel with no locking; metasymbols are assigned from
// the sorted dictionary, so the split does not change the output. With n_parts==1
// this is just a single phrase_map_t. Consumers that iterate the table (the
// dictionary constructor and the metasymbol assignment) walk parts[0..n_parts)
// in order; lookups route to the owning part.
struct partitioned_map_t {
    std::vector<phrase_map_t> parts;
    size_t n_parts=1;

    partitioned_map_t()=default;
    partitioned_map_t(size_t n_parts_, float max_load_factor, size_t d_bits): n_parts(n_parts_) {
        parts.reserve(n_parts);
        for(size_t i=0;i<n_parts;i++) parts.emplace_back(max_load_factor, d_bits);
    }

    //the part a phrase belongs to (same function used to build and to look up)
    static inline size_t part_of(const void* key, size_t key_bits, size_t n_parts) {
        return XXH3_64bits(key, INT_CEIL(key_bits, 8)) % n_parts;
    }

    //pure size queries; nodiscard makes the compiler warn if the result is ignored
    [[nodiscard]] inline size_t size() const {
        size_t s=0; for(auto const& m: parts) s+=m.size(); return s;
    }
    [[nodiscard]] inline size_t description_bits() const { return parts[0].description_bits(); }

    inline bool key2value(const void* key, size_t key_bits, size_t& value) const {
        return parts[part_of(key, key_bits, n_parts)].key2value(key, key_bits, value);
    }

    //dump/reload/free the parts (one file per part, suffixed .<i>)
    void store_data_to_file(const std::string& prefix) {
        for(size_t i=0;i<parts.size();i++) parts[i].store_data_to_file(prefix+"."+std::to_string(i));
    }
    void load_data_from_file(const std::string& prefix) {
        for(size_t i=0;i<parts.size();i++) parts[i].load_data_from_file(prefix+"."+std::to_string(i));
    }
    void remove_data_files(const std::string& prefix) {
        for(size_t i=0;i<parts.size();i++) remove((prefix+"."+std::to_string(i)).c_str());
    }
    void destroy_table(){ for(auto& m: parts) m.destroy_table(); }
    void destroy_data(){ for(auto& m: parts) m.destroy_data(); }
    void shrink_databuff(){ for(auto& m: parts) m.shrink_databuff(); }
};
#endif //GRLBWT_COMMON_H
