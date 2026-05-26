//
// Created by Diaz, Diego on 28.4.2023.
//

#ifndef CDS_VBYTE_H
#define CDS_VBYTE_H
#include<iostream>
#include<cassert>

inline uint8_t vbyte_len(size_t value){
    uint8_t n_bits = (sizeof(unsigned long)<<3) - __builtin_clzl(value | 1);
    n_bits+=((n_bits+7)>>3);
    return (n_bits+7)>>3;
}

struct vbyte {

    static constexpr uint8_t shift[8]={0, 7, 14, 21, 28, 35, 42, 49};
    static constexpr uint64_t mask[8]={0x7F, 0x3FFF, 0x1FFFFF, 0xFFFFFFF,
                                       0x7FFFFFFFF, 0x3FFFFFFFFFF, 0x1FFFFFFFFFFFF, 0xFFFFFFFFFFFFFF};
    //typedef sym_t sym_type;
    /*inline static off_t forward(const uint8_t * ptr, off_t n_pos) {
        off_t n_bytes=0;
        for(off_t i=0;i<n_pos;i++){
            while(ptr[n_bytes]<128) n_bytes++;
            n_bytes++;
        }
        return n_bytes;
    }

    inline static off_t backward(uint8_t * ptr, const uint8_t *& l_boundary, off_t n_pos) {

        if(ptr==l_boundary) return 0;
        uint8_t *tmp = ptr-1;
        while(tmp!=l_boundary && n_pos>0){
            tmp--;
            while(tmp!=l_boundary && *tmp<128){
                --tmp;
            }
            n_pos--;
        }
        tmp++;
        assert(tmp!=l_boundary);
        return ptr-tmp;
    }

    //rightmost byte of the rightmost valid value before ptr:
    // For instance, 0 0 0 1* 0 0* 0 1 (1 is the last byte of each symbol).
    // The function returns 3 with pos = 5 (assuming uint32_t symbols of sizeof(uint32_t)= 4 bytes)
    inline static off_t mov_to_prev_valid_byte(uint8_t *&ptr, __attribute__((unused)) off_t& pos) {
        uint8_t * tmp = ptr;
        while(*ptr<128) --ptr;
        return tmp-ptr;
    }*/

    //static void write(uint8_t * stream, uint64_t x, size_t vb_len){
    //    switch (vb_len) {
    //        case 1:
    //            stream[0] = x + 128;
    //            break;
    //        case 2:
    //            stream[0] = x & 127;
    //            x>>=7;
    //            stream[1] = x + 128;
    //            break;
    //        case 3:
    //            stream[0] = x & 127;
    //            x>>=7;
    //            stream[1] = (x & 127);
    //            x>>=7;
    //            stream[2] = x + 128;
    //            break;
    //        case 4:
    //            stream[0] = x & 127;
    //            x>>=7;
    //            stream[1] = (x & 127);
    //            x>>=7;
    //            stream[2] = (x & 127);
    //            x>>=7;
    //            stream[3] = x + 128;
    //            break;
    //        case 5:
    //            stream[0] = x & 127;
    //            x>>=7;
    //            stream[1] = (x & 127);
    //            x>>=7;
    //            stream[2] = (x & 127);
    //            x>>=7;
    //            stream[3] = (x & 127);
    //            x>>=7;
    //            stream[4] = x + 128;
    //            break;
    //        case 6:
    //            stream[0] = x & 127;
    //            x>>=7;
    //            stream[1] = (x & 127);
    //            x>>=7;
    //            stream[2] = (x & 127);
    //            x>>=7;
    //            stream[3] = (x & 127);
    //            x>>=7;
    //            stream[4] = (x & 127);
    //            x>>=7;
    //            stream[5] = x + 128;
    //            break;
    //        case 7:
    //            stream[0] = x & 127;
    //            x>>=7;
    //            stream[1] = (x & 127);
    //            x>>=7;
    //            stream[2] = (x & 127);
    //            x>>=7;
    //            stream[3] = (x & 127);
    //            x>>=7;
    //            stream[4] = (x & 127);
    //            x>>=7;
    //            stream[5] = (x & 127);
    //            x>>=7;
    //            stream[6] = x + 128;
    //            break;
    //        default:
    //            stream[0] = x & 127;
    //            x>>=7;
    //            stream[1] = (x & 127);
    //            x>>=7;
    //            stream[2] = (x & 127);
    //            x>>=7;
    //            stream[3] = (x & 127);
    //            x>>=7;
    //            stream[4] = (x & 127);
    //            x>>=7;
    //            stream[5] = (x & 127);
    //            x>>=7;
    //            stream[6] = (x & 127);
    //            x>>=7;
    //            stream[7] = x + 128;
    //            break;
    //    }
    //    //size_t len;
    //    //size_t comp_x;
    //    //size_t c=0, b=0;
    //    //while(x>=128){
    //    //    c += (x & 127)<<b;
    //    //    b+=8;
    //    //    x>>=7;
    //    //}
    //    //comp_x = c+((x+128)<<b);
    //    //len = (b+8)>>3UL;//in bytes
    //    //memcpy(stream, &c, vb_len);
    //}

    static size_t write(uint8_t *stream, uint64_t x) {

        stream[0] = x & 0x7F;
        x >>= 7;
        if (x == 0) {
            stream[0] |= 0x80;
            return 1;
        }

        stream[1] = x & 0x7F;
        x >>= 7;
        if (x == 0) {
            stream[1] |= 0x80;
            return 2;
        }

        stream[2] = x & 0x7F;
        x >>= 7;
        if (x == 0) {
            stream[2] |= 0x80;
            return 3;
        }

        stream[3] = x & 0x7F;
        x >>= 7;
        if (x == 0) {
            stream[3] |= 0x80;
            return 4;
        }

        stream[4] = x & 0x7F;
        x >>= 7;
        if (x == 0) {
            stream[4] |= 0x80;
            return 5;
        }

        stream[5] = x & 0x7F;
        x >>= 7;
        if (x == 0) {
            stream[5] |= 0x80;
            return 6;
        }

        stream[6] = x & 0x7F;
        x >>= 7;
        if (x == 0) {
            stream[6] |= 0x80;
            return 7;
        }

        stream[7] = x | 0x80;
        return 8;
    }

    static size_t read(const uint8_t * ptr, uint64_t& x) {
        /*uint64_t i=0;
        x=ptr[0];
        std::cout<<x<<std::endl;
        while(i<7 && ptr[i]<128){
            uint64_t tmp = ptr[++i];
            x |= tmp<<shift[i];
            std::cout<<x<<std::endl;
        }
        x &= mask[i];
        return i+1;*/

        /*const uint64_t block = *reinterpret_cast<const uint64_t *>(ptr);
        const unsigned i = (__builtin_ctzll(block & 0x8080808080808080ULL)) >> 3;

        x = ((block & 0x7F)      ) |
            ((block & 0x7F00)    ) >> 1 |
            ((block & 0x7F0000)  ) >> 2 |
            ((block & 0x7F000000)) >> 3 |
            ((block & 0x7F00000000ULL)) >> 4 |
            ((block & 0x7F0000000000ULL)) >> 5 |
            ((block & 0x7F000000000000ULL)) >> 6 |
            ((block & 0x7F00000000000000ULL)) >> 7;
        x &= mask[i];
        return i + 1;*/

        /*const uint64_t block = *reinterpret_cast<const uint64_t *>(ptr);
        const unsigned i = (__builtin_ctzll(block & 0x8080808080808080ULL)) >> 3;
        switch (i) {
            case 0:
                x = (block & 0x7F);
                break;
            case 1:
                x = ((block & 0x7F)      ) |
                    ((block & 0x7F00)    ) >> 1;
                break;
            case 2:
                x = ((block & 0x7F)      ) |
                    ((block & 0x7F00)    ) >> 1 |
                    ((block & 0x7F0000)  ) >> 2;
                break;
            case 3:
                x = ((block & 0x7F)      ) |
                    ((block & 0x7F00)    ) >> 1 |
                    ((block & 0x7F0000)  ) >> 2 |
                    ((block & 0x7F000000)) >> 3;
                break;
            case 4:
                x = ((block & 0x7F)      ) |
                    ((block & 0x7F00)    ) >> 1 |
                    ((block & 0x7F0000)  ) >> 2 |
                    ((block & 0x7F000000)) >> 3 |
                    ((block & 0x7F00000000ULL)) >> 4;
                break;
            case 5:
                x = ((block & 0x7F)      ) |
                    ((block & 0x7F00)    ) >> 1 |
                    ((block & 0x7F0000)  ) >> 2 |
                    ((block & 0x7F000000)) >> 3 |
                    ((block & 0x7F00000000ULL)) >> 4 |
                    ((block & 0x7F0000000000ULL)) >> 5;
                break;
            case 6:
                x = ((block & 0x7F)      ) |
                    ((block & 0x7F00)    ) >> 1 |
                    ((block & 0x7F0000)  ) >> 2 |
                    ((block & 0x7F000000)) >> 3 |
                    ((block & 0x7F00000000ULL)) >> 4 |
                    ((block & 0x7F0000000000ULL)) >> 5 |
                    ((block & 0x7F000000000000ULL)) >> 6;
                break;
            case 7:
                x = ((block & 0x7F)      ) |
                    ((block & 0x7F00)    ) >> 1 |
                    ((block & 0x7F0000)  ) >> 2 |
                    ((block & 0x7F000000)) >> 3 |
                    ((block & 0x7F00000000ULL)) >> 4 |
                    ((block & 0x7F0000000000ULL)) >> 5 |
                    ((block & 0x7F000000000000ULL)) >> 6 |
                    ((block & 0x7F00000000000000ULL)) >> 7;
                break;
            default:
                throw std::runtime_error("Invalid vbyte");
        }
        x &= mask[i];
        return i + 1;*/

        //manually unrolled
        uint64_t b = ptr[0];
        uint64_t tmp = b & 0x7f;
        if (b >= 128) {
            x = tmp;
            return 1;
        }

        b = ptr[1];
        tmp |= (b & 0x7f) << 7;
        if (b >= 128) {
            x = tmp;
            return 2;
        }

        b = ptr[2];
        tmp |= (b & 0x7f) << 14;
        if (b >= 128) {
            x = tmp;
            return 3;
        }

        b = ptr[3];
        tmp |= (b & 0x7f) << 21;
        if (b >= 128) {
            x = tmp;
            return 4;
        }

        b = ptr[4];
        tmp |= (b & 0x7f) << 28;
        if (b >= 128) {
            x = tmp;
            return 5;
        }

        b = ptr[5];
        tmp |= (b & 0x7f) << 35;
        if (b >= 128) {
            x = tmp;
            return 6;
        }

        b = ptr[6];
        tmp |= (b & 0x7f) << 42;
        if (b >= 128) {
            x = tmp;
            return 7;
        }

        b = ptr[7];
        tmp |= (b & 0x7f)  << 49;
        x = tmp;
        return 8;
    }

    //this assumes ptr points to the rightmost byte of a symbol in the stream
    //notice this function changes the memory address
    static off_t read_backwards(uint8_t *&ptr, const uint8_t *& l_boundary, uint64_t& sym) {
        assert(ptr!=l_boundary);
        assert(*ptr>=128);

        uint8_t * tmp = ptr;
        sym = *ptr-128;
        ptr--;
        while(ptr!=l_boundary && *ptr<128){
            sym = (sym<<7) + *ptr;
            --ptr;
        }
        ptr++;
        return tmp-ptr;
    }

    //it gets the rightmost symbol that differs from the symbol in ptr
    //it scans "prefix" number of bytes preceding ptr to look for the mismatching symbol
    static off_t mov_to_prev_diff_sym(uint8_t *& ptr, const uint8_t *l_boundary, uint64_t& sym) {
        //get the length in bytes of the current symbol
        uint8_t * tmp = ptr;
        while(tmp!=l_boundary && *tmp<128) tmp++;
        off_t bytes = tmp-ptr+1;
        sym = 0;
        off_t dist = ptr-l_boundary-1;//number of bytes we can visit

        off_t n_comp = dist/bytes;
        tmp = ptr;
        ptr-=bytes;
        while(n_comp>0 && memcmp(tmp, ptr, bytes)==0){
            ptr -= bytes;
            n_comp--;
        }

        ptr+=bytes-1;
        if(ptr!=l_boundary){
            //this extra loop is to cover a weird corner case.
            // E.g., 92 131, 131, 157.
            // In this case, the comparison of rightmost occurrence of 131 will finish in 92,
            // but the backward read has to start from the leftmost occurrence of 131
            while(*ptr<128) ptr++;
            assert(tmp-ptr>0);
            read_backwards(ptr, l_boundary, sym);
        }
        return tmp-ptr;
    }

    //get the leftmost symbol that differs from the symbol in ptr
    //it will stop checking once the scan reaches the boundary address
    //the function assumes ptr and boundary are addresses within the same memory area
    static off_t mov_to_next_diff_sym(uint8_t *& ptr, const uint8_t* r_boundary, uint64_t& sym) {

        //get the length in bytes of the current symbol
        uint8_t * tmp = ptr;
        while(*tmp<128) tmp++;
        off_t bytes = tmp-ptr+1;

        sym = 0;
        tmp = ptr;
        off_t dist = r_boundary-ptr-1;//number of bytes we can visit
        off_t n_comp = dist/bytes;
        ptr+=bytes;
        while(n_comp>0 && memcmp(tmp, ptr, bytes)==0){
            ptr+=bytes;
            n_comp--;
        }

        if(ptr!=r_boundary){
            read(ptr, sym);
        }
        return ptr-tmp;
    }
};

template<typename sym_t>
struct inv_vbyte {
    static void write(uint8_t * stream, sym_t x, uint8_t vb_len){
        size_t c;
        switch (vb_len) {
            case 1:
                c = (x & 127)+128;
                break;
            case 2:
                c = (x & 127)+128;
                x>>=7;
                c += (x & 127)<<8;
                break;
            case 3:
                c = (x & 127)+128;
                x>>=7;
                c += (x & 127)<<8;
                x>>=7;
                c += (x & 127)<<16;
                break;
            case 4:
                c = (x & 127)+128;
                x>>=7;
                c += (x & 127)<<8;
                x>>=7;
                c += (x & 127)<<16;
                x>>=7;
                c += (x & 127)<<24;
                break;
            default:
                exit(1);
        }
        memcpy(stream, &c, vb_len);
    }

    //this assumes ptr points to the rightmost byte of a vbyte-encoded symbol in the stream
    static off_t read(uint8_t *ptr,  sym_t& sym) {
        uint8_t * tmp = ptr;
        sym=0;
        while(*ptr<128){
            sym = (sym<<7) + *ptr;
            --ptr;
        }
        sym = (sym<<7) + (*ptr-128);
        return tmp-ptr+1;
    }
};
#endif //CDS_VBYTE_H
