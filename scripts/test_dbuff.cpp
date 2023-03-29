//
// Created by Diaz, Diego on 5.3.2023.
//
#include <iostream>
#include <vector>
#include <cassert>
#include "utils.h"
#include "../include/common.h"
#include "../include/parsing_strategies.h"
#include <thread>
#include <queue>
#include <mutex>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <condition_variable>

template<class sym_t>
struct text_chunk {
    typedef sym_t sym_type;
    size_t id{};
    off_t bytes{};
    size_t max_len{};
    sym_t * buffer = nullptr;

    std::vector<long> str_ptr;

    explicit text_chunk(off_t chunk_bytes): bytes(chunk_bytes),
                                            max_len(chunk_bytes/sizeof(sym_t)){
        buffer = (sym_t*) malloc(bytes);
    }

    text_chunk() =default;

    text_chunk(text_chunk&& other)  noexcept {
        std::swap(id, other.id);
        std::swap(bytes, other.bytes);
        std::swap(buffer, other.buffer);
        std::swap(str_ptr, other.str_ptr);
        std::swap(max_len, other.max_len);
    }

    ~text_chunk(){
        if(buffer!= nullptr){
            free(buffer);
        }
    }

    inline sym_type operator[](size_t idx) const {
        if(idx>=max_len){
            std::cout<<idx<<" "<<max_len<<" "<<bytes<<" "<<id<<std::endl;
        }
        assert(idx<max_len);
        return buffer[idx];
    }
};

template<typename T>
class unbounded_queue {

public:
    unbounded_queue() = default;
    ~unbounded_queue() {
        done();
    };

    void push(const T& item) {
        {
            std::unique_lock guard(m_queue_lock);
            m_queue.push(item);
        }
        m_condition.notify_one();
    }

    void push(T&& item) {
        {
            std::unique_lock guard(m_queue_lock);
            m_queue.push(std::move(item));
        }
        m_condition.notify_one();
    }

    bool pop(T& item) {
        std::unique_lock guard(m_queue_lock);
        m_condition.wait(guard, [&]() { return !m_queue.empty() || m_done; });
        if(m_done) return false;
        item = std::move(m_queue.front());
        m_queue.pop();
        return true;
    }

    std::size_t size() const {
        std::unique_lock guard(m_queue_lock);
        return m_queue.size();
    }

    bool empty() const {
        std::unique_lock guard(m_queue_lock);
        return m_queue.empty();
    }

    void done() {
        {
            std::unique_lock guard(m_queue_lock);
            m_done = true;
        }
        m_condition.notify_all();
    }

private:
    using queue_t = std::queue<T>;
    queue_t m_queue;
    mutable std::mutex m_queue_lock;
    std::condition_variable m_condition;
    bool m_done = false;
};

template<class text_chunk_t>
void read_chunk_from_file(int fd, text_chunk_t& chunk, off_t& rem_text_bytes){

    using sym_t = typename text_chunk_t::sym_type;

    off_t chunk_bytes = chunk.bytes<rem_text_bytes ? chunk.bytes : rem_text_bytes;
    chunk.bytes = chunk_bytes;

    off_t acc_bytes = 0;
    off_t read_bytes;
    off_t fd_buff_bytes = 8388608;// 8MB buffer
    sym_t sep_symbol = 10;

    off_t sym_bytes = sizeof(sym_t);
    off_t buff_pos=0;
    off_t end_pos;

    chunk.str_ptr.clear();
    chunk.str_ptr.push_back(0);

    while(chunk_bytes>0) {

        fd_buff_bytes = fd_buff_bytes<chunk_bytes ? fd_buff_bytes : chunk_bytes;
        read_bytes = read(fd, &chunk.buffer[buff_pos], fd_buff_bytes);
        assert(read_bytes>0);

        end_pos = buff_pos + (read_bytes/sym_bytes);

        for(off_t i=buff_pos;i<end_pos;i++){
            if(chunk.buffer[i]==sep_symbol){
                chunk.str_ptr.emplace_back(i+1);
            }
        }
        buff_pos = end_pos;

        chunk_bytes-=read_bytes;
        acc_bytes+=read_bytes;
    }

    if(chunk.str_ptr.empty()){
        //keep reading until compleating at least one string
    }

    assert(!chunk.str_ptr.empty());
    off_t offset = acc_bytes-(chunk.str_ptr.back()*sizeof(sym_t));

    std::cout<<chunk.id<<" "<<acc_bytes<<" "<<chunk.str_ptr.back()*sizeof(sym_t)<<" "<<offset<<" "<<chunk.str_ptr.back()-1<<std::endl;

    lseek(fd, offset*-1, SEEK_CUR);

    rem_text_bytes-=chunk.str_ptr.back()*sizeof(sym_t);
    std::cout<<"rem bytes: "<<rem_text_bytes<<std::endl;
    assert(chunk_bytes==0);
}

int main (int argc, char** argv){

    if(argc!=2){
        std::cout<<"usage: ./text_stats text"<<std::endl;
        std::cout<<"text_stats: a plain text in one-string-per-line format"<<std::endl;
        exit(0);
    }

    std::string input_file = std::string(argv[1]);

    std::cout<<" Running in parallel "<<std::endl;
    size_t active_chunks=20;
    off_t chunk_size = 1024*1024*80;

    unbounded_queue<size_t> in_strings;
    unbounded_queue<size_t> out_strings;
    std::mutex mutex;

    std::vector<text_chunk<uint8_t>> text_chunks;
    bool all_strings_submitted=false;

    auto file_reader = [&]() -> void{

        int fd = open(input_file.c_str(), O_RDONLY);

        struct stat st{};
        if(stat(input_file.c_str(), &st) != 0)  return;
        off_t rem_bytes = st.st_size;

#ifdef __linux__
        posix_fadvise(fd, 0, rem_bytes, POSIX_FADV_SEQUENTIAL);
#endif

        size_t chunk_id=0;
        text_chunks.resize(active_chunks);
        while(chunk_id<active_chunks && rem_bytes>0){
            chunk_size = std::min(chunk_size, rem_bytes);
            text_chunks[chunk_id].bytes = chunk_size;
            text_chunks[chunk_id].max_len = chunk_size/sizeof(uint8_t);
            text_chunks[chunk_id].buffer = (uint8_t *)malloc(chunk_size);
            text_chunks[chunk_id].id = chunk_id;
            read_chunk_from_file<text_chunk<uint8_t>>(fd, text_chunks[chunk_id], rem_bytes);
            in_strings.push(chunk_id);
            chunk_id++;
        }

        size_t buff_idx;
        while(rem_bytes>0){
            out_strings.pop(buff_idx);
            text_chunks[buff_idx].id = chunk_id++;
            read_chunk_from_file<text_chunk<uint8_t>>(fd, text_chunks[buff_idx], rem_bytes);
            in_strings.push(buff_idx);
        }
        while(!in_strings.empty());

        while(!out_strings.empty()){
            out_strings.pop(buff_idx);
        }

        in_strings.done();
        out_strings.done();

        close(fd);

        {
            std::lock_guard<std::mutex> guard(mutex);
            all_strings_submitted = true;
        }
    };


    auto string_consumer = [&](){
        phrase_map_t map;

        auto hasher=[&](string_t& phrase){
            phrase.mask_tail();
            map.increment_value(phrase.data(), phrase.n_bits(), 1);
        };

        size_t buff_id;
        bool res;

        while(true){
            res = in_strings.pop(buff_id);
            assert(text_chunks[buff_id].bytes>0);
            if(!res) break;
            auto init_str=[&](size_t idx)-> std::pair<long, long>{
                return {text_chunks[buff_id].str_ptr[idx], text_chunks[buff_id].str_ptr[idx+1]-1};
            };
            lms_parsing<text_chunk<uint8_t>, string_t, true>()(text_chunks[buff_id], 0, text_chunks[buff_id].str_ptr.size()-2, 90,  hasher, init_str);
            out_strings.push(buff_id);
        }

        /*{
            std::lock_guard<std::mutex> guard(mutex);
            for(size_t i=0;i<256;i++){
                sym_freqs[i]+=t_sym_freqs[i];
            }
        }*/
    };

    size_t n_threads=4;
    std::vector<std::thread> threads;
    threads.emplace_back(file_reader);
    for(size_t i=0;i<n_threads;i++){
        threads.emplace_back(string_consumer);
    }

    for(auto & thread : threads){
        thread.join();
    }

    /*size_t acc=0;
    for(size_t i=0;i<256;i++){
        if(sym_freqs[i]!=0){
            acc+=sym_freqs[i];
            std::cout<<(char)i<<" "<<sym_freqs[i]<<std::endl;
        }
    }
    std::cout<<"Total symbols "<<acc<<std::endl;*/
    std::vector<text_chunk<uint8_t>>().swap(text_chunks);

#ifdef __linux__
    posix_fadvise(fd, 0, st.st_size, POSIX_FADV_DONTNEED);
#endif

    /*uint8_t sep_sym;
    {
        std::ifstream ifs(input_file, std::ios::binary);
        ifs.seekg(-1, std::ios::end);
        ifs.read((char *)&sep_sym, 1);//we assume the last symbol is the sep symbol
        ifs.seekg(0, std::ios::beg);
    }

    dbuff_i_file_stream<uint8_t, 1024,4> db_if(input_file);
    auto start = high_resolution_clock::now();
    size_t sym_frq[256]={0};
    uint8_t sym;
    size_t cont=0, pos=0, str_len;
    size_t longest_string=0;
    size_t shortest_string=std::numeric_limits<size_t>::max();

    for(size_t i=0;i<db_if.size();i++){
        sym = db_if.value();
        sym_frq[sym]++;
        cont++;
        if(sym==sep_sym){
            str_len = cont - pos;
            if(str_len>longest_string) longest_string=str_len;
            if(str_len<shortest_string) shortest_string=str_len;
            pos = cont;
        }
        db_if.next();
    }

    std::vector<uint8_t> alphabet;
    std::vector<size_t> sym_freqs;
    size_t n_symbols=0;
    for(size_t i=0;i<256;i++){
        if(sym_frq[i]!=0){
            alphabet.push_back(i);
            sym_freqs.push_back(sym_frq[i]);
            n_symbols+=sym_frq[i];
        }
    }

    std::cout<<"Total symbols:         "<<n_symbols<<std::endl;
    std::cout<<"Alphabet:              "<<alphabet.size()<<std::endl;
    for(size_t i=0;i<alphabet.size();i++){
        std::cout<<"    symbol: "<<(int)alphabet[i]<<" freq: "<<sym_freqs[i]<<std::endl;
    }
    std::cout<<"Smallest symbol:       "<<(int)alphabet[0]<<std::endl;
    std::cout<<"Greatest symbol:       "<<(int)alphabet.back()<<std::endl;
    std::cout<<"Number of strings:     "<<sym_freqs[0]<<std::endl;
    std::cout<<"Longest string:        "<<longest_string<<std::endl;
    std::cout<<"Shortest string:       "<<shortest_string<<std::endl;
    std::cout<<"Average string length: "<<n_symbols/sym_freqs[0]<<std::endl;
    std::cout<<"Sep. symbol:           "<<(int)sep_sym<<std::endl;

    if(sep_sym!=alphabet[0]){
        std::cerr<<"Warning: the separator symbol is not the smallest symbol in the text"<<std::endl;
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    std::cout<<"Elapsed time: "<<duration.count()<<" microseconds"<<std::endl;
    return 0;*/
}