//
// Created by Diaz, Diego on 5.3.2023.
//
#include <iostream>
#include <vector>
#include <cassert>
#include "utils.h"
#include <thread>
#include <queue>
#include <mutex>
#include <fcntl.h>

using namespace std::chrono;
off_t acc=0;

template<class sym_t>
struct string_buffer {
    sym_t * string;
    size_t len;
    size_t id;
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

template<class sym_t>
void next_string_from_file(int fd, off_t buff_len, string_buffer<sym_t>& string){

    off_t len = string.len*sizeof(sym_t);
    auto *tmp = (char *)string.string;
    off_t read_bytes;

    while(len>0){
        buff_len = std::min(buff_len, len);
        read_bytes = read(fd, tmp, buff_len);
        len-=read_bytes;
        assert(read_bytes>0);
        tmp+=read_bytes;
    }
    assert(string.string[string.len-1]=='\n');
}

int main (int argc, char** argv){

    if(argc!=2){
        std::cout<<"usage: ./text_stats text"<<std::endl;
        std::cout<<"text_stats: a plain text in one-string-per-line format"<<std::endl;
        exit(0);
    }

    std::string input_file = std::string(argv[1]);

    std::cout<<" Computing collection stats "<<std::endl;
    str_collection str_col = collection_stats(input_file);
    str_col.str_ptrs.push_back((long)str_col.n_char);

    std::cout<<"File has "<<str_col.n_strings<<" strings "<<std::endl;
    std::cout<<" Running in parallel "<<std::endl;

    size_t active_strings=20;
    unbounded_queue<size_t> in_strings;
    unbounded_queue<size_t> out_strings;
    std::mutex print_mutex;
    auto string_buffs = (string_buffer<uint8_t> *)malloc(sizeof(string_buffer<uint8_t>)*active_strings);
    bool all_strings_submitted=false;
    size_t sym_freqs[256]={0};

    auto file_reader = [&](){

        int fd = open(input_file.c_str(), O_RDONLY);
        size_t proc_str=0;

        size_t str_id;
        for(size_t i=0;i<active_strings;i++){
            string_buffs[i].len = str_col.str_ptrs[i+1]-str_col.str_ptrs[i];
            string_buffs[i].string = (uint8_t*) malloc(sizeof(uint8_t)*str_col.longest_string);
            string_buffs[i].id = i;
            next_string_from_file<uint8_t>(fd, 1024*1024, string_buffs[i]);

            in_strings.push(i);
        }
        str_id = active_strings;

        while(str_id<str_col.n_strings){
            size_t buff_idx;
            out_strings.pop(buff_idx);

            proc_str++;
            string_buffs[buff_idx].len = str_col.str_ptrs[str_id+1]-str_col.str_ptrs[str_id];
            string_buffs[buff_idx].id = str_id++;
            next_string_from_file<uint8_t>(fd, 1024*1024, string_buffs[buff_idx]);

            in_strings.push(buff_idx);
        }

        while(proc_str<str_col.n_strings){
            size_t buff_idx;
            out_strings.pop(buff_idx);
            proc_str++;
        }

        in_strings.done();
        out_strings.done();

        close(fd);
        {
            std::lock_guard<std::mutex> guard(print_mutex);
            all_strings_submitted = true;
        }
    };

    auto string_consumer = [&](){
        size_t buff_id;
        size_t t_sym_freqs[256]={0};
        bool res;

        while(true){
            res = in_strings.pop(buff_id);
            if(!res) break;

            /*{
                std::lock_guard<std::mutex> guard(print_mutex);
                std::cout <<"Thread "<<std::this_thread::get_id()<<" is consuming string buff_id: "<<buff_id<<" -> str_id: "<< string_buffs[buff_id].id <<" len: "<<string_buffs[buff_id].len<<" "<<in_strings.size()<<std::endl;
            }*/
            for(size_t i=0;i<string_buffs[buff_id].len;i++){
                t_sym_freqs[string_buffs[buff_id].string[i]]++;
            }
            out_strings.push(buff_id);
        }

        {
            std::lock_guard<std::mutex> guard(print_mutex);
            for(size_t i=0;i<256;i++){
                sym_freqs[i]+=t_sym_freqs[i];
            }
        }
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

    for(size_t i=0;i<256;i++){
        if(sym_freqs[i]!=0){
            std::cout<<(char)i<<" "<<sym_freqs[i]<<std::endl;
        }
    }

    for(size_t i=0;i<active_strings;i++){
        free(string_buffs[i].string);
    }
    free(string_buffs);

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