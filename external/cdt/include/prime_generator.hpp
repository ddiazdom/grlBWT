//
// Created by Diego Diaz on 11/6/20.
//

#ifndef LPG_COMPRESSOR_PRIME_GENERATOR_HPP
#define LPG_COMPRESSOR_PRIME_GENERATOR_HPP
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <random>
#include "prng.hpp"
#include "mod_operations.hpp"

class prime_generator {
private:
    typedef unsigned long long int size_type;

public:
    inline static size_type get_next_prime(size_t number){

        //k = number of iterations
        size_type k=4, prime_candidate=2;

        if(!(number & 1UL)) number++;
        prime_candidate = number+2;

        std::random_device rd;
        prng prng_inst(rd());

        while(!miller_rabin_prim_test(prime_candidate, k, prng_inst)){
            prime_candidate+=2;
        }

        assert(mod_operations::is_prime(prime_candidate));
        return prime_candidate;
    };

    inline static bool witness(size_t a, size_t n){
        size_t t,u,x, prev_x;

        size_t tmp,pos =0;
        tmp = n-1;
        while(!(tmp & (1U<<pos))){
            pos++;
        }
        u = tmp >> pos;
        t = pos;

        x = mod_operations::mod_pow(a,u,n);
        prev_x = x;

        for(size_t i=0;i<t;i++){
            x = mod_operations::mod_pow(x ,2,n);
            if(x==1 && prev_x!=1 && prev_x != tmp){
                return true;
            }
            prev_x = x;
        }
        return x != 1;
    };

    inline static bool miller_rabin_prim_test(size_t n, size_t s, prng& prng_inst){
        size_t a;
        size_t tmp = n-1;
        for(size_t j=0;j<s;j++){
            a = 1 + prng_inst.rand()%tmp;
            if(witness(a, n)){
                return false;
            }
        }
        return true;
    };

};
#endif //LPG_COMPRESSOR_PRIME_GENERATOR_HPP
