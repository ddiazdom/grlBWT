//
// Created by Diaz, Diego on 23.1.2022.
//

#ifndef GBWT_BWT_IO_H
#define GBWT_BWT_IO_H
#include <iostream>
#include <fstream>

class bwt_array{
    char * buffer;
    uint8_t sb;
    uint8_t fb;
    uint8_t bpr;
    size_t n_bytes;

public:

    bwt_array(uint8_t _sb, uint8_t _fb, size_t capacity) : sb(_sb),
                                                           fb(_fb),
                                                           bpr(sb+fb),
                                                           n_bytes(capacity*bpr) {
        buffer = (char *)malloc(n_bytes);
    }

    void write_run(size_t id, size_t sym, size_t freq){
    }

    void read_run(size_t id, size_t& sym, size_t& freq) const {
    }

    [[nodiscard]] inline size_t size() const {
        return n_bytes/bpr;
    }

    ~bwt_array(){
        free(buffer);
    }
};

class bwt_buff_reader{

    char * buffer =nullptr;
    char * sec_buffer;

    size_t bpr=0; //bytes per run
    size_t sb=0; // bytes for a run symbol
    size_t fb=0; // bytes for a run freq
    size_t block_bg=0; //start of the current block
    size_t tot_runs =0; //total number of runs in the file
    size_t buffer_size=0; //size of the buffer
    size_t offset =sizeof(size_t)*2;
    std::ifstream ifs;

public:

    explicit bwt_buff_reader(const std::string& input_file, size_t buff_size= 1024 * 1024){
        ifs = std::ifstream(input_file, std::ifstream::binary);
        buffer = (char *)malloc(buff_size);
        ifs.read((char *)&sb, sizeof(size_t));
        ifs.read((char *)&fb, sizeof(size_t));
        buffer_size = buff_size;
        bpr = sb + fb;
        sec_buffer = (char *) malloc(bpr);

        block_bg=0;
        ifs.seekg (0, std::ifstream::end);
        tot_runs = ifs.tellg();
        tot_runs = (tot_runs-offset)/bpr;
        ifs.seekg ((std::streamoff)offset);
        ifs.read(buffer, (std::streamsize) buffer_size);
    }

    inline void read(size_t i, size_t& sym, size_t& freq) {
        assert(i<tot_runs);
        sym = 0;
        freq = 0;
        size_t start = (i*bpr);
        size_t end = start + bpr-1;
        size_t buff_start = (start & (buffer_size-1));

        if(start<block_bg || (block_bg+buffer_size)<=end) {

            size_t b_start = start / buffer_size;
            size_t b_end = end / buffer_size;

            block_bg = (end / buffer_size) * buffer_size;
            ifs.seekg((std::streamoff)(offset+block_bg));

            if (b_start < b_end) {
                size_t left = buffer_size-buff_start;
                memcpy(sec_buffer, buffer + buff_start, left);
                ifs.read((char *) buffer, (std::streamsize) buffer_size);
                memcpy(sec_buffer+left, buffer, bpr-left);
                memcpy(&sym, sec_buffer, sb);
                memcpy(&freq, sec_buffer+sb, fb);
                return;
            }else{
                ifs.read((char *) buffer, (std::streamsize) buffer_size);
            }
        }
        memcpy(&sym, buffer+buff_start, sb);
        memcpy(&freq, buffer+buff_start+sb, fb);
    }

    inline size_t size() const {
       return tot_runs;
    }

    void close(){
        if(buffer!= nullptr){
            free(buffer);
            buffer = nullptr;
            free(sec_buffer);
            sec_buffer = nullptr;
        }
        if(ifs.is_open()) ifs.close();
    }
};

class bwt_buff_writer {

    char * buffer =nullptr;
    char * sec_buff =nullptr;
    size_t bpr=0; //bytes per run
    size_t sb=0; // bytes for a run symbol
    size_t fb=0; // bytes for a run freq
    size_t block_bg=0; //start of the current block
    size_t tot_runs =0; //total number of runs in the file
    size_t buffer_size=0; //size of the buffer
    std::ofstream ofs;

public:

    explicit bwt_buff_writer(const std::string& input_file, uint8_t _sb, uint8_t _fb, size_t buff_size= 1024 * 1024){

        ofs.open(input_file, std::ios::binary);
        assert(ofs.good());

        buffer_size = buff_size;
        buffer = (char *)malloc(buffer_size);
        sb = _sb;
        fb = _fb;
        bpr = sb+fb;
        sec_buff = (char *)malloc(bpr);
        ofs.write((char *)&sb, sizeof(size_t));
        ofs.write((char *)&fb, sizeof(size_t));

        block_bg=0;
    }

    inline void push_back(size_t& sym, size_t& freq) {

        size_t start = (tot_runs*bpr);
        size_t end = start + bpr-1;
        size_t buff_start = (start & (buffer_size - 1));

        if(start>=block_bg && (block_bg + buffer_size) <= end){
            size_t b_start = start/buffer_size;
            size_t b_end = end/buffer_size;
            if(b_start<b_end) {
                size_t left = buffer_size - buff_start;
                memcpy(sec_buff, &sym, sb);
                memcpy(sec_buff + sb, &freq, fb);
                memcpy(buffer+buff_start, sec_buff, left);
                ofs.write((char *) buffer, (std::streamsize) buffer_size);
                block_bg = (end/buffer_size)*buffer_size;
                memcpy(buffer, sec_buff+left, bpr-left);
                tot_runs++;
                return;
            }else{
                ofs.write((char *) buffer, (std::streamsize) buffer_size);
                block_bg = (end/buffer_size)*buffer_size;
            }
        }
        memcpy(buffer+buff_start, &sym, sb);
        memcpy(buffer+buff_start+sb, &freq, fb);
        tot_runs++;
    }

    void close(){
        if(buffer!= nullptr){
            size_t start = (tot_runs*bpr);
            size_t buff_start = (start & (buffer_size - 1));
            if(buff_start!=0){
                ofs.write((char *) buffer, (std::streamsize) buff_start);
            }else{
                ofs.write((char *) buffer, (std::streamsize) buffer_size);
            }
            free(buffer);
            buffer = nullptr;
        }
        if(sec_buff!= nullptr){
            free(sec_buff);
            sec_buff = nullptr;
        }
        if(ofs.is_open()) ofs.close();
    }
};

#endif //GBWT_BWT_IO_H
