//
// Created by Diaz, Diego on 3.2.2023.
//
#include <vector>
#include "hash_table.hpp"
#include "int_array.h"
#include "mxm3.h"

#include <chrono>
#include <random>

int main(int argc, char** argv) {
    std::ifstream ifs("data_ht_bug.txt");
    bit_hash_table<size_t> ht;
    long long duration_total = 0;
    long long duration_total2 = 0;

    std::random_device dev;
    std::mt19937_64 rng(dev());
    std::uniform_int_distribution<std::mt19937_64::result_type> dist_long(0,1099511627775ULL);

    std::random_device dev2;
    std::mt19937 rng2(dev2());
    std::uniform_int_distribution<std::mt19937::result_type> dist_short(2,100);

    std::vector<std::vector<size_t>> list(1000000, std::vector<size_t>());
    for(auto & str : list){
        size_t length = dist_short(rng2);
        for(size_t j=0;j<length;j++){
            str.push_back(dist_long(rng));
        }
    }

    int_array<size_t> wide_str(2, 40);
    for(auto & str : list){
        wide_str.clear();
        for(auto & sym : str){
            wide_str.push_back(sym);
        }
        wide_str.mask_tail();

        auto start = std::chrono::steady_clock::now();
        XXH3_64bits(wide_str.data(), wide_str.n_bits());
        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
        duration_total += duration;

        auto start2 = std::chrono::steady_clock::now();
        mx3::hash(reinterpret_cast<const uint8_t *>((char *) wide_str.data()), INT_CEIL(wide_str.n_bits(), 8), 0);
        auto end2 = std::chrono::steady_clock::now();
        auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>( end2 - start2 ).count();
        duration_total2 += duration2;
    }
    std::cout <<"insertion: "<<duration_total<<" miliseconds "<<std::endl;
    std::cout <<"insertion: "<<duration_total2<<" miliseconds "<<std::endl;

    /*duration_total=0;
    for(auto & str : list){
        wide_str.clear();
        for(auto & sym : str){
            wide_str.push_back(sym);
        }
        wide_str.mask_tail();
        auto start = std::chrono::steady_clock::now();
        auto res = ht.find(wide_str.data(), wide_str.n_bits());
        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
        duration_total += duration;
        assert(res.second);
    }
    std::cout <<"find: "<<duration_total<<" miliseconds "<<std::endl;*/

    /*while (ifs >> a >> b) {
        auto start = std::chrono::steady_clock::now();
        //std::cout<<i++<<" "<<a<<" -> "<<b<<" "<<ht.max_bucket_dist()<<" "<<ht.size()<<" "<<ht.tot_buckets()<<" "<<ht.load_factor()<<std::endl;
        auto res = ht.insert(&a, sizeof(a) * 8, b);
        auto end = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>( end - start ).count();
        duration_total += duration;
    }
    std::cout <<duration_total<<" miliseconds "<<std::endl;
    ht.ht_stats(10);
    std::cout<<ht.size()<<std::endl;*/
}

