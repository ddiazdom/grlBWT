//
// Created by Diaz, Diego on 14.3.2023.
//
#include "file_streams.hpp"

template<class sym_t>
int i_file_stream<sym_t>::fd_advice(int fd, off_t offset, off_t length, int advice) {
#ifdef __linux__
    return posix_fadvice(fd, offset, length, advice);
#else
    return -1;
#endif
}
