//
// Created by Diego Diaz on 11/6/20.
//

#ifndef LPG_COMPRESSOR_MOD_OPERATIONS_HPP
#define LPG_COMPRESSOR_MOD_OPERATIONS_HPP
#include<iostream>
#include<unordered_set>

class mod_operations{
public:

    // C function for extended Euclidean Algorithm
    static long long int gcdExtended(long long int a, long long int b, long long int *x, long long int *y){
        // Base Case
        if (a == 0){
            *x = 0, *y = 1;
            return b;
        }

        signed long long int x1, y1; // To store results of recursive call
        long long int gcd = gcdExtended(b%a, a, &x1, &y1);

        // Update x and y using results of recursive
        // call
        *x = y1 - (b/a) * x1;
        *y = x1;

        return gcd;
    }

    // Function to find modulo inverse of a
    static size_t modInverse(size_t a, size_t m){
        long long int x, y;
        long long int g = gcdExtended((long long int)a, (long long int)m, &x, &y);
        if (g != 1) {
            std::cout << "Inverse doesn't exist";
            exit(EXIT_FAILURE);
        } else {
            // m is added to handle negative x
            return (x%m + m) % m;
        }
    }

    // Utility function to store prime factors of a number
    static void findPrimefactors(std::unordered_set<int> &s, int n){
        // Print the number of 2s that divide n
        while (n%2 == 0){
            s.insert(2);
            n = n/2;
        }

        // n must be odd at this point. So we can skip
        // one element (Note i = i +2)
        for (int i = 3; i <= sqrt(n); i = i+2) {
            // While i divides n, print i and divide n
            while (n%i == 0) {
                s.insert(i);
                n = n/i;
            }
        }

        // This condition is to handle the case when
        // n is a prime number greater than 2
        if (n > 2) s.insert(n);
    }

    // returns true if n is prime
    inline static bool is_prime(size_t n){
        // Corner cases
        if (n <= 1)  return false;
        if (n <= 3)  return true;

        // This is checked so that we can skip
        // middle five numbers in below loop
        if (n%2 == 0 || n%3 == 0) return false;

        for (size_t i=5; i*i<=n; i=i+6){
            if (n%i == 0 || n%(i+2) == 0){
                return false;
            }
        }
        return true;
    }

    // Function to find smallest primitive root of n
    static int findSmallestPrimitiveRoot(int n){
        std::unordered_set<int> s;

        // Check if n is prime or not
        if (!is_prime(n))
            return -1;

        // Find value of Euler Totient function of n
        // Since n is a prime number, the value of Euler
        // Totient function is n-1 as there are n-1
        // relatively prime numbers.
        int phi = n-1;

        // Find prime factors of phi and store in a set
        findPrimefactors(s, phi);

        // Check for every number from 2 to phi
        for (int r=2; r<=phi; r++)
        {
            // Iterate through all prime factors of phi.
            // and check if we found a power with value 1
            bool flag = false;
            for (auto it = s.begin(); it != s.end(); it++)
            {

                // Check if r^((phi)/primefactors) mod n
                // is 1 or not
                if(mod_pow(r, phi/(*it), n) == 1){
                    flag = true;
                    break;
                }
            }

            // If there was no power with value 1.
            if (!flag)
                return r;
        }

        // If no primitive root found
        return -1;
    }

    inline static size_t mod_pow(uint64_t x, uint64_t y, size_t p){
        uint64_t res = 1;// Initialize result
        x = x % p;//Update x if it is more than or

        // equal to p
        while (y > 0){
            // If y is odd, multiply x with result
            if (y & 1U)
                res = (res*x) % p;

            // y must be even now
            y = y>>1U; // y = y/2
            x = (x*x)%p;
        }
        return res;
    }
};
#endif //LPG_COMPRESSOR_MOD_OPERATIONS_HPP
