//
// Created by diego on 17-02-20.
//

#ifndef LMS_COMPRESSOR_FILE_BUFFER_H
#define LMS_COMPRESSOR_FILE_BUFFER_H

#include <fstream>
#include <cassert>
#include <iostream>
#include <cstring>
#include <limits>
#include "bitstream.h"

template<class sym_t>
struct i_file_stream{

    std::ifstream text_i;
    size_t tot_cells{};
    size_t block_bg{};
    bitstream<sym_t> buffer;

    std::string file;
    static constexpr size_t w_bytes = sizeof(sym_t);

    i_file_stream(i_file_stream<sym_t> &&other) noexcept{
        move(std::forward<i_file_stream<sym_t>>(other));
    }

    i_file_stream(const std::string& i_file, size_t buff_size_){
        file =  i_file;
        text_i = std::ifstream(file, std::ifstream::binary);

        text_i.seekg (0, std::ifstream::end);
        tot_cells = text_i.tellg()/w_bytes;
        text_i.seekg (0, std::ifstream::beg);

        size_t buff_size = std::min<size_t>(buff_size_/w_bytes, tot_cells);
        buffer.stream_size = buff_size;
        buffer.stream = new sym_t[buff_size];

        block_bg=0;
        text_i.read((char*) buffer.stream, buff_size*w_bytes);
    }

    i_file_stream& operator=(i_file_stream<sym_t> &&other) noexcept{
        move(std::forward<i_file_stream<sym_t>>(other));
        return *this;
    }

    void move(i_file_stream<sym_t> &&other){
        std::swap(file, other.file);
        std::swap(text_i, other.text_i);
        std::swap(tot_cells, other.tot_cells);
        std::swap(block_bg, other.block_bg);
        buffer.swap(other.buffer);
    }

    inline size_t size() const{
        return tot_cells;
    }

    inline const std::string& filename() const{
        return file;
    }

    void close(bool rem=false){
        if(buffer.stream!= nullptr){
            delete [] buffer.stream;
            buffer.stream = nullptr;
        }
        if(text_i.is_open()) text_i.close();
        if(rem){
            if(remove(file.c_str())){
                std::cout<<"Error trying to remove file"<<file<<std::endl;
            }
        }
    }

    ~i_file_stream(){
        close();
    }

    inline size_t read(size_t i) {
        assert(i<tot_cells);
        if(i<block_bg || (block_bg+buffer.stream_size)<=i){
            block_bg = (i/buffer.stream_size)*buffer.stream_size;
            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
        }
        return buffer.stream[i-block_bg];
    }

    inline void read_chunk(void *dst, size_t i, size_t j) {
        assert(i<=j);

        size_t cell_i = i/(w_bytes*8);
        size_t cell_j = j/(w_bytes*8);

        assert(cell_j<tot_cells && (cell_j-cell_i+1)<=buffer.stream_size);

        if(cell_i<block_bg || (block_bg+buffer.stream_size)<=cell_j){
            block_bg = cell_i;
            text_i.seekg(cell_i*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
        }

        size_t start = i - block_bg*w_bytes*8;
        buffer.read_chunk(dst, start, start+(j-i));
    }
};

template<class sym_t>
struct o_file_stream{

    typedef long long int size_type;
    std::ofstream ofs;
    std::ifstream ifs;
    std::string file;

    size_type w_bytes=0;
    size_type block_bg=0;
    size_type rm_pos=-1;
    size_type lm_pos=std::numeric_limits<size_type>::max();
    size_type last_pos=-1;
    bitstream<sym_t> buffer;
    bool modified=false;

    o_file_stream()=default;

    o_file_stream(o_file_stream<sym_t> &&other) noexcept{
        move(std::forward<o_file_stream<sym_t>>(other));
    }

    o_file_stream(const std::string& o_file, size_t buff_size_, std::ios::openmode mode){
        file = o_file;

        ofs.open(file, mode | std::ios::out | std::ios::binary);
        assert(ofs.good());

        ifs.open(file,  std::ios::in | std::ios::binary);
        assert(ifs.good());

        w_bytes = sizeof(sym_t);

        size_t buff_size = buff_size_/w_bytes;
        buffer.stream_size = buff_size;
        buffer.stream = new sym_t[buff_size];

        ifs.seekg (0, std::ifstream::end);
        last_pos = (ifs.tellg()/w_bytes)-1;
        ifs.seekg (0, std::ifstream::beg);

        block_bg=0;
        read_block();
    }

    o_file_stream<sym_t>& operator=(o_file_stream<sym_t> &&other) noexcept{
        if(this!=&other){
            move(std::forward<o_file_stream<sym_t>>(other));
        }
        return *this;
    }

    void move(o_file_stream<sym_t>&& other){
        std::swap(ifs, other.ifs);
        std::swap(ofs, other.ofs);
        std::swap(file, other.file);
        std::swap(w_bytes, other.w_bytes);
        std::swap(block_bg, other.block_bg);
        std::swap(rm_pos, other.rm_pos);
        std::swap(lm_pos, other.lm_pos);
        std::swap(last_pos, other.last_pos);
        std::swap(modified, other.modified);
        buffer.swap(other.buffer);
    }

    void push_back(sym_t val){

        if((block_bg + (size_type)buffer.stream_size) <= (++last_pos)){
            //write data
            write_block();

            //reset buffer
            memset(buffer.stream, 0, buffer.stream_size*sizeof(sym_t));
            block_bg = (last_pos/buffer.stream_size)*buffer.stream_size;
            rm_pos=0;
            lm_pos=std::numeric_limits<size_type>::max();
        }

        size_type buff_pos = last_pos-block_bg;
        buffer.stream[buff_pos] = val;
        rm_pos = buff_pos;
        if(buff_pos<lm_pos) lm_pos = buff_pos;
        modified=true;
    }

    inline size_t size() const{
        return (last_pos+1);
    }

    const std::string& filename() const{
        return file;
    }

    size_t read(size_t i) {
        assert((size_type)i<=last_pos);
        if((size_type)i<block_bg || (block_bg+(size_type)buffer.stream_size)<=(size_type)i){

            if(modified){
                //write data
                write_block();
            }

            block_bg = (i/buffer.stream_size)*buffer.stream_size;
            read_block();
        }
        return buffer.stream[i-block_bg];
    }

    void write_block(){
        assert(lm_pos<=rm_pos);
        ofs.seekp((block_bg+lm_pos)*w_bytes);
        ofs.write((char*) &buffer.stream[lm_pos], (rm_pos-lm_pos+1)*w_bytes);
        //TODO caveat: it might produce desynchronization with the ifs in some platforms (OSx)
        // if I don't flush the ofs buffer. However, flusing the buffer is expensive.
        // Maybe I should define an if macro here.
        assert(ofs.good());
        modified=false;
    }

    void read_block(){
        assert(!modified);
        ifs.seekg(block_bg*w_bytes);
        ifs.read((char*) buffer.stream, buffer.stream_size*w_bytes);
        if(ifs.gcount()<(size_type)buffer.stream_size*w_bytes){
            ifs.clear();
        }
        assert(ifs.good());
        rm_pos=-1;
        lm_pos=std::numeric_limits<size_type>::max();
    }

    void flush(){
        if(modified){
            write_block();
            ofs.flush();
            //init the buffer
            block_bg = 0;
            read_block();
        }else{
            ofs.flush();
        }
    }

    void write(size_type i, sym_t val) {

        if(i<block_bg || (block_bg+(size_type)buffer.stream_size)<=i){
            //write data
            if(modified){
                write_block();
                modified = false;
            }

            //read the new segment
            block_bg = (i/buffer.stream_size)*buffer.stream_size;
            read_block();
        }
        size_type buff_pos = i-block_bg;

        if(i >= last_pos) last_pos = i;
        if(buff_pos>rm_pos){
            rm_pos = buff_pos;
        }

        if(buff_pos<lm_pos){
            lm_pos = buff_pos;
        }

        buffer.stream[buff_pos] = val;
        modified = true;
    }

    void close(bool rem=false){
        flush();
        if(buffer.stream!=nullptr){
            delete [] buffer.stream;
            buffer.stream=nullptr;
        }
        if(ofs.is_open()){
            ofs.close();
        }
        if(ifs.is_open()){
            ifs.close();
        }
        if(rem){
            if(remove(file.c_str())){
                std::cout<<"Error trying to remove file"<<file<<std::endl;
            }
        }
    }

    ~o_file_stream(){
        close();
    }
};
#endif //LMS_COMPRESSOR_FILE_BUFFER_H
