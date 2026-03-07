//
// Created by Diaz, Diego on 5.3.2022.
//
#include "cdt_common.hpp"

uint8_t sym_width(unsigned long val){
    if(val==0) return 0;
    return (sizeof(unsigned long)*8) - __builtin_clzl(val);
}
