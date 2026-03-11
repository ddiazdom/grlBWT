//
// Created by Diaz, Diego on 9.3.2026.
//

#ifndef CDT_RLBWT_IO_H
#define CDT_RLBWT_IO_H

#include <cassert>
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <cdt/vbyte.h>
#include <cdt/memory_handler.hpp>
#include <cstring>

class rle_vbyte_buff_reader {
    size_t buff_cap=0;
    size_t buff_len=0;
    size_t buff_pos=0;
    size_t file_size=0;
    size_t tot_read_bytes=0;
    size_t last_vbyte_pos=0;
    size_t read_runs=0;
    uint8_t *buffer=nullptr;
    std::ifstream ifs;

    void update_last_vbyte_pos() {
        if(buff_len == 0){
            throw std::runtime_error("Empty or malformed vbyte stream");
        }
        assert(buffer!=nullptr);
        last_vbyte_pos = buff_len-1;
        while(last_vbyte_pos>0 && buffer[last_vbyte_pos]<128) last_vbyte_pos--;
        if(buffer[last_vbyte_pos]<128)  throw std::runtime_error("Malformed vbyte stream");
    }

    void decode_elm(uint64_t& elm) {
        if(buff_pos>last_vbyte_pos) {
            const size_t rem_bytes = buff_len - buff_pos;
            memmove(buffer, buffer+buff_pos, rem_bytes);
            buff_pos = 0;
            ifs.read(reinterpret_cast<char *>(buffer+rem_bytes),
                     static_cast<std::streamsize>(buff_cap-rem_bytes));
            buff_len = static_cast<size_t>(ifs.gcount())+rem_bytes;
            update_last_vbyte_pos();
        }
        size_t read_vbyte = vbyte::read(buffer+buff_pos, elm);
        buff_pos+=read_vbyte;
        tot_read_bytes+=read_vbyte;
    }

public:
    explicit rle_vbyte_buff_reader(const std::string& input_file, size_t buff_bytes= 1024 * 1024) {
        assert(std::filesystem::exists(input_file));
        file_size = std::filesystem::file_size(input_file);
        ifs.open(input_file, std::ios::in | std::ios::binary);
        assert(ifs.good());
        buff_cap = std::max<size_t>(16, std::min(buff_bytes, file_size));
        buffer = mem::allocate<uint8_t>(buff_cap+8);//the +8 is to decode vbyte fast
        memset(buffer + buff_cap, 0x80, 8);//fill the tail with termination (0x80) vbyte symbols
        ifs.read(reinterpret_cast<char *>(buffer), static_cast<std::streamsize>(buff_cap));
        buff_len = static_cast<size_t>(ifs.gcount());
        update_last_vbyte_pos();
    }

    void next_run(uint64_t& sym, uint64_t&len) {
        decode_elm(sym);
        decode_elm(len);
        read_runs++;
    }

    bool has_next() const {
        return tot_read_bytes<file_size;
    }

    size_t n_read_runs() const {
        return read_runs;
    }

    rle_vbyte_buff_reader(const rle_vbyte_buff_reader&) = delete;
    rle_vbyte_buff_reader& operator=(const rle_vbyte_buff_reader&) = delete;

    ~rle_vbyte_buff_reader() noexcept {
        ifs.close();
        mem::deallocate(buffer);
    }
};

//reads a run-length-vbyte-encoded sequence and updated the runs in-place.
//this class works like an interator
//invariant: it can update the last accessed run, but the new values
//have always to be equal or smaller than the previous
class rle_vbyte_buff_updater {
    size_t buff_cap=0;//buffer capacity
    size_t buff_len=0;//valid bytes load in the buffer
    size_t buff_i_pos=0;//byte position of the read operation
    size_t buff_o_pos=0;//byte position of the write operation (update)
    size_t file_size=0;//original file size
    size_t f_read_pos=0;//bytes read from the file
    size_t f_write_pos=0;//bytes written in the file prefix (always less than f_read_pos)
    size_t tot_dec_vbytes=0;//total vbytes decoded
    size_t last_vbyte_pos=0;//the last byte of the last fully valid vbyte in the buffer
    size_t last_sym=0;//last decoded run symbol
    size_t last_len=0;//last decoded run length
    size_t read_runs=0;//number of reads decoded
    uint8_t *buffer=nullptr;//buffer
    std::string input_file;//original file
    bool closed=false;
    static constexpr size_t buff_offset=16;//minimum difference between buffers. We need this offset to avoid overlaps in edge cases

    std::fstream fs;

    void update_last_vbyte_pos() {
        if(buff_len == 0){
            throw std::runtime_error("Empty or malformed vbyte stream");
        }
        assert(buffer!=nullptr);
        last_vbyte_pos = buff_len-1;
        while(last_vbyte_pos>0 && buffer[last_vbyte_pos]<128) last_vbyte_pos--;
        if(buffer[last_vbyte_pos]<128)  throw std::runtime_error("Malformed vbyte stream");
    }

    void decode_elm(uint64_t& elm) {
        if(buff_i_pos>last_vbyte_pos) {
            refill_buffer(buff_len-buff_i_pos);
        }
        const size_t read_vbyte = vbyte::read(buffer+buff_i_pos, elm);
        buff_i_pos+=read_vbyte;
        tot_dec_vbytes+=read_vbyte;
    }

    void compress_previous() {
        if (read_runs==0) return;
        const size_t vlen_sym = vbyte_len(last_sym);
        const size_t vlen_len = vbyte_len(last_len);
        assert(vlen_sym + vlen_len <= buff_cap);
        vbyte::write(buffer+buff_o_pos, last_sym, vlen_sym);
        buff_o_pos+=vlen_sym;
        vbyte::write(buffer+buff_o_pos, last_len, vlen_len);
        buff_o_pos+=vlen_len;
        assert(buff_o_pos<=buff_i_pos);
    }

    void flush_modified_prefix() {
        if (buff_o_pos==0) return;
        fs.clear();
        fs.seekp(static_cast<std::streamoff>(f_write_pos), std::ios::beg);
        fs.write(reinterpret_cast<char *>(buffer),
                 static_cast<std::streamsize>(buff_o_pos));
        assert(fs.good());

        f_write_pos+=buff_o_pos;
        buff_o_pos=0;
    }

    void refill_buffer(const size_t rem_bytes) {
        flush_modified_prefix();
        memmove(buffer+buff_offset, buffer+buff_i_pos, rem_bytes);

        fs.clear();
        fs.seekg(static_cast<std::streamoff>(f_read_pos), std::ios::beg);
        fs.read(reinterpret_cast<char *>(buffer+buff_offset+rem_bytes), static_cast<std::streamsize>(buff_cap-buff_offset-rem_bytes));

        auto read_bytes = static_cast<size_t>(fs.gcount());
        f_read_pos+=read_bytes;

        buff_len = buff_offset+rem_bytes + read_bytes;
        buff_i_pos=buff_offset;
        update_last_vbyte_pos();
    }

public:

    explicit rle_vbyte_buff_updater(const std::string& input_file_, size_t buff_bytes= 1024 * 1024) {
        input_file = input_file_;
        assert(std::filesystem::exists(input_file));
        file_size = std::filesystem::file_size(input_file);
        fs.open(input_file, std::ios::in | std::ios::out | std::ios::binary);
        assert(fs.good());

        buff_cap = std::max<size_t>(16, std::min(buff_bytes, file_size))+buff_offset;
        buffer = mem::allocate<uint8_t>(buff_cap+8);//we add+8 bytes to decode vbytes fast
        memset(buffer + buff_cap, 0x80, 8);//fill the tail with termination vbyte symbols

        assert(buffer!=nullptr);
        refill_buffer(0);
    }

    rle_vbyte_buff_updater(const rle_vbyte_buff_updater&) = delete;

    rle_vbyte_buff_updater& operator=(const rle_vbyte_buff_updater&) = delete;

    void next_run(uint64_t& sym, uint64_t&len) {
        assert(!closed);
        compress_previous();
        decode_elm(sym);
        decode_elm(len);
        last_sym = sym;
        last_len = len;
        read_runs++;
    }

    void update_sym(const uint64_t sym) {
        assert(sym<=last_sym && read_runs>0);
        last_sym = sym;
    }

    void update_len(const uint64_t len) {
        assert(len<=last_len && read_runs>0);
        last_len = len;
    }

    bool has_next() const {
        return tot_dec_vbytes<file_size;
    }

    size_t n_read_runs() const {
        return read_runs;
    }

    void close() {
        if(closed) return;
        compress_previous();
        flush_modified_prefix();
        fs.flush();
        fs.close();
        assert(f_write_pos<=file_size);
        std::filesystem::resize_file(input_file, f_write_pos);
        closed=true;
    }

    ~rle_vbyte_buff_updater() noexcept {
        close();
        mem::deallocate(buffer);
    }
};

class rle_vbyte_buff_writer {

    size_t last_sym=0;
    size_t last_len=0;
    size_t n_runs=0;
    size_t buff_len=0;
    size_t buff_pos=0;
    size_t seq_size=0;
    uint8_t *buffer=nullptr;
    bool closed=false;
    std::ofstream ofs;

    void compress_previous() {
        if (n_runs==0) return;

        const size_t vlen_sym = vbyte_len(last_sym);
        const size_t vlen_len = vbyte_len(last_len);
        assert(vlen_sym + vlen_len <= buff_len);

        //store to file if we exceed buffer size
        if((buff_pos+vlen_sym+vlen_len)>=buff_len) {
            ofs.write(reinterpret_cast<char *>(buffer),
                      static_cast<std::streamsize>(buff_pos));
            buff_pos = 0;
        }

        //write vbytes and move the offset
        vbyte::write(buffer+buff_pos, last_sym, vlen_sym);
        buff_pos+=vlen_sym;
        vbyte::write(buffer+buff_pos, last_len, vlen_len);
        buff_pos+=vlen_len;
    }

public:

    void push_back(size_t sym_, size_t len_) {
        compress_previous();
        last_sym = sym_;
        last_len = len_;
        n_runs++;
        seq_size+=len_;
    }

    size_t sym() const {
        return last_sym;
    }

    size_t len() const {
        return last_len;
    }

    size_t size() const {
        return seq_size;
    }

    size_t tot_runs() const {
        return n_runs;
    }

    void update_sym(size_t sym_) {
        last_sym = sym_;
    }

    void increase_len(size_t len_) {
        last_len += len_;
        seq_size += len_;
    }

    void decrease_len(size_t len_) {
        assert(len_<last_len);
        last_len -= len_;
        seq_size -= len_;
    }

    explicit rle_vbyte_buff_writer(const std::string& input_file,
                               size_t buff_bytes= 1024 * 1024) {
        ofs.open(input_file,  std::ios::out | std::ios::binary);
        assert(ofs.good());
        buff_len = std::max<size_t>(16, buff_bytes);
        buffer = mem::allocate<uint8_t>(buff_len);
    }

    rle_vbyte_buff_writer(const rle_vbyte_buff_writer&) = delete;

    rle_vbyte_buff_reader& operator=(const rle_vbyte_buff_reader&) = delete;

    void close() {
        if (!closed) {
            compress_previous();
            if(buff_pos>0) {
                ofs.write(reinterpret_cast<char *>(buffer), static_cast<std::streamsize>(buff_pos));
                buff_pos=0;
                ofs.flush();
            }
            closed = true;
        }
    }

    ~rle_vbyte_buff_writer() noexcept {
        close();
        buff_len = 0;
        ofs.close();
        mem::deallocate(buffer);
    }
};
#endif //CDT_RLBWT_IO_H