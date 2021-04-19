//
// Created by diegodiaz on 25-09-20.
//

#ifndef LPG_COMPRESSOR_SE_INT_QUEUE_H
#define LPG_COMPRESSOR_SE_INT_QUEUE_H
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <utility>

//semi-external integer queue
template<class value_t>
class se_int_queue{

private:
    std::fstream fs;
    size_t buff_bytes=0;
    value_t *buffer=nullptr;
    long long idx_head=0;
    long long idx_tail=-1;
    long long n_cells=0;
    size_t elm_on_disk=0;
    size_t g_start=0;
    std::string buffer_file;

    void copy(se_int_queue<value_t>& other){
        fs = other.fs;
        buff_bytes = other.buff_bytes;
        buffer = reinterpret_cast<value_t*>(malloc(buff_bytes));
        memcpy(buffer, other.buffer, buff_bytes);
        idx_head = other.idx_head;
        idx_tail = other.idx_tail;
        n_cells = other.n_cells;
        elm_on_disk = other.elm_on_disk;
        g_start = other.g_start;
        buffer_file = other.buffer_file;
    }

    void move(se_int_queue<value_t>&& other){
        std::swap(fs, other.fs);
        buff_bytes = other.buff_bytes;
        std::swap(buffer, other.buffer);
        idx_head = other.idx_head;
        idx_tail = other.idx_tail;
        n_cells = other.n_cells;
        elm_on_disk = other.elm_on_disk;
        g_start = other.g_start;
        buffer_file = other.buffer_file;
    }

public:

    se_int_queue(const se_int_queue<value_t>& other){
        copy(other);
    }

    se_int_queue(const se_int_queue<value_t>&& other) noexcept {
        move(std::forward<se_int_queue<value_t>>(other));
    }

    se_int_queue& operator=(const se_int_queue<value_t>& other){
        if(this != &other){
            copy(other);
        }
        return *this;
    }

    se_int_queue& operator=(se_int_queue<value_t>&& other) noexcept{
        if(this != &other){
            move(std::forward<se_int_queue<value_t>>(other));
        }
        return *this;
    }

    se_int_queue(std::string file, size_t buff_bytes_) : idx_head(0),
                                                         idx_tail(-1),
                                                         elm_on_disk(0),
                                                         buffer_file(std::move(file)) {
        assert(buff_bytes_ % sizeof(value_t)==0);
        n_cells = buff_bytes_ / sizeof(value_t);
        assert(n_cells % 2 ==0);
        fs.open(buffer_file, std::fstream::out);
        fs.close();
        fs.open(buffer_file, std::fstream::in | std::fstream::out| std::fstream::binary);
        assert(fs.good());
        buff_bytes = buff_bytes_;
        buffer = reinterpret_cast<value_t*>(malloc(buff_bytes));
    }

    void push_back(value_t elm){
        if((idx_tail+1) >= n_cells){
            long long n_elms = size()+1;
            if(n_elms>=n_cells){
                fs.seekg(g_start + (elm_on_disk*sizeof(value_t)) );
                fs.write((char *)&buffer[n_cells/2], (buff_bytes / 2));
                assert(fs.good());
                idx_tail = (n_cells / 2)-1;
                elm_on_disk += n_cells / 2;
            }else{
                memmove(buffer, &buffer[idx_head], (n_elms-1)*sizeof(value_t));
                idx_head=0;
                idx_tail = n_elms-2;
            }
        }
        buffer[++idx_tail] = elm;
    }

    value_t pop(){
        size_t val = buffer[idx_head];
        idx_head++;
        if(idx_head >= (n_cells/2) ){
            if(elm_on_disk>0){
                size_t n_elms = std::min<size_t>(elm_on_disk, (n_cells / 2));
                size_t read_bytes = n_elms*sizeof(value_t);
                fs.seekg(g_start);
                fs.read((char *)buffer, read_bytes);
                assert(fs.good());
                elm_on_disk -= n_elms;
                g_start += read_bytes;
                idx_head = 0;
            }else{
                long long n_elms = idx_tail-idx_head+1;
                memmove(buffer, &buffer[idx_head], n_elms*sizeof(value_t));
                idx_head=0;
                idx_tail = n_elms-1;
            }
        }else if(idx_head>idx_tail){
            idx_head = 0;
            idx_tail = -1;
        }
        return val;
    }

    inline value_t front() const {
        assert(idx_tail>=idx_head);
        return buffer[idx_head];
    }

    inline value_t tail() const {
        assert(idx_tail>=idx_head);
        return buffer[idx_tail];
    }

    inline size_t size() const {
        if(idx_tail<idx_head){
            return 0;
        }else{
            return idx_tail-idx_head+1+elm_on_disk;
        }
    }

    inline bool empty() const {
        return idx_tail<idx_head;
    }

    void close(bool remove_file){
        fs.close();
        if(remove_file){
            if(remove(buffer_file.c_str())){
                std::cout<<"Error trying to remove file "<<buffer_file;
            }
        }
    }

    ~se_int_queue(){
        if(buffer!= nullptr) free(buffer);
        close(false);
    }
};

#endif //LPG_COMPRESSOR_SE_INT_QUEUE_H
