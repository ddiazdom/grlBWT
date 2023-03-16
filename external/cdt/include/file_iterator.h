//
// Created by Diaz, Diego on 15.3.2023.
//

#ifndef CDT_FILE_ITERATOR_H
#define CDT_FILE_ITERATOR_H

#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <cassert>
#include <string>
#include <iostream>

template<typename sym_t>
class file_iterator{

    int fd=-1;
    static constexpr uint8_t w_bytes = sizeof(sym_t);
    sym_t *buffer=nullptr;
    off_t buff_size;
    std::string file;
    size_t b_pos=0;
    size_t curr_pos;
    size_t last_pos;

public:
    explicit file_iterator(std::string &_file, off_t _buff_size=1024*1024, off_t start=0, off_t end=0) : file(_file),
                                                                                                           buff_size(_buff_size),
                                                                                                           curr_pos(start){
        struct stat st{};
        stat(file.c_str(), &st);
        off_t fsz = st.st_size;
        off_t tot_syms = fsz / w_bytes;
        assert(tot_syms>0);

        if(start<tot_syms){

            if(end==0) end = tot_syms;
            assert(start<end && end<=tot_syms);

            fd = open(file.c_str(), O_RDONLY);
            assert(fd>0);

#ifdef __linux__
            ssize_t tot_bytes = (end-start)*w_bytes;
            posix_advise(fd, start*w_bytes, tot_bytes, POSIX_FADV_SEQUENTIAL);
#endif
            buff_size = std::min(buff_size, fsz);
            buffer = (sym_t *)malloc(buff_size);

            lseek(fd, start*w_bytes, SEEK_SET);
            size_t read_elm = read(fd, buffer, buff_size);
            assert(read_elm>0 && read_elm<=buff_size);
            last_pos = end-1;
        }else{
            curr_pos = tot_syms;
            last_pos = tot_syms;
        }
    }

    inline sym_t operator *() const {
        assert(curr_pos<=last_pos);
        return buffer[b_pos];
    }

    inline void operator++(){
        b_pos++;
        curr_pos++;
        if(b_pos==buff_size && curr_pos<last_pos){
            ssize_t read_elm = read(fd, buffer, buff_size);
            assert(read_elm>0 && read_elm<=buff_size);
            buff_size = read_elm;
            b_pos = 0;
        }
    }

    [[nodiscard]] inline bool consumed() const {
        return curr_pos>last_pos;
    }
};
#endif //CDT_FILE_ITERATOR_H
