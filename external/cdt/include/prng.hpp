//
// Created by Diego Diaz on 11/6/20.
//

#ifndef LPG_COMPRESSOR_PRNG_HPP
#define LPG_COMPRESSOR_PRNG_HPP
#include <iostream>

class prng {
private:
    uint64_t state[2] = {0xdeadbeefcafebabe, 0x8badf00dbaada555};
    uint64_t splitmix64_x;

private:
    inline static uint64_t xoroshiro128plus(uint64_t s[2]) {
        uint64_t s0 = s[0];
        uint64_t s1 = s[1];
        uint64_t result = s0 + s1;
        s1 ^= s0;
        s[0] = ((s0 << 55U) | (s0 >> 9U)) ^ s1 ^ (s1 << 14U);
        s[1] = (s1 << 36U) | (s1 >> 28U);
        return result;
    }

    inline uint64_t splitmix64() {
        uint64_t z = (splitmix64_x += UINT64_C(0x9E3779B97F4A7C15));
        z = (z ^ (z >> 30U)) * UINT64_C(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27U)) * UINT64_C(0x94D049BB133111EB);
        return z ^ (z >> 31U);
    }

public:
    explicit prng(size_t seed){
        splitmix64_x = seed;
        state[0] = splitmix64();
        state[1] = splitmix64();
    }
    inline uint64_t rand(){
        return xoroshiro128plus(state);
    };
};
#endif //LPG_COMPRESSOR_PRNG_HPP
