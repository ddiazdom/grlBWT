//
// Created by Diaz, Diego on 4.3.2023.
//

#ifndef CDT_DBUFF_FILE_STREAM_H
#define CDT_DBUFF_FILE_STREAM_H
#include <bitstream.h>
#include <thread>
#include <fcntl.h>
#include <unistd.h>
#include <mutex>
#include <condition_variable>

//std::mutex g_display_mutex;

template<class sym_t,
         size_t n_pages,
         uint8_t n_buff=2>
struct dbuff_i_file_stream{

    typedef sym_t sym_type;
    typedef bitstream<sym_t> stream_t;
    //std::ifstream text_i;
    size_t tot_syms{};
    size_t text_pos{};
    size_t buff_pos{};
    size_t curr_buff={};
    size_t buff_size={};
    int fp{};

    bool worker_thread_alive=true;
    std::mutex m;
    std::condition_variable cv;

    std::string data;
    bool needs_disk_acc = false;

    bitstream<sym_t> buffers[n_buff];
    size_t inv_buff=0;
    std::thread dw;

    std::string file;
    static constexpr size_t w_bytes = sizeof(sym_t);
    //constexpr static uint8_t w_bits = std::numeric_limits<sym_t>::digits;
    //constexpr static uint8_t w_shift = __builtin_ctz(w_bits);

    dbuff_i_file_stream(dbuff_i_file_stream &&other) noexcept{
        move(std::forward<dbuff_i_file_stream>(other));
    }

    explicit dbuff_i_file_stream(const std::string& i_file){

        file =  i_file;
        {
            std::ifstream text_i = std::ifstream(file, std::ifstream::binary);
            assert(text_i.good());
            text_i.seekg(0, std::ifstream::end);
            size_t file_size = text_i.tellg();
            text_i.seekg(0, std::ifstream::beg);
            tot_syms = INT_CEIL(file_size, w_bytes);
            text_i.close();
        }
        auto page_size = (size_t)sysconf(_SC_PAGESIZE);;
        int file_flags = O_RDONLY;

#ifdef __linux__
        file_flags|=O_DIRECT;
#endif
        fp = open(file.c_str(), file_flags);

        buff_pos = 0;
        text_pos = 0;
        curr_buff = 0;
        inv_buff = 0;

        assert(((page_size*n_buff) % sizeof(sym_t))==0);
        buff_size =  (page_size*n_pages)/sizeof(sym_t);

        for(size_t i=0;i<n_buff;i++){
            buffers[i].stream_size = buff_size;
            if (posix_memalign((void **)&buffers[i].stream, page_size, buff_size*sizeof(sym_t))) {
                perror("posix_memalign failed");
            }
            read(fp, (char *) buffers[i].stream, buff_size*sizeof(sym_t));
        }

        dw = std::thread([&](){
            std::unique_lock<std::mutex> lk(m);
            lk.unlock();

            while(worker_thread_alive){
                lk.lock();
                cv.wait(lk, [&]{return needs_disk_acc;});
                read(fp, (char *)buffers[inv_buff].stream, buff_size*sizeof(sym_t));
                needs_disk_acc = false;
                lk.unlock();
                cv.notify_one();
            }
        });
        assert(dw.joinable());

        /*std::thread::id this_id = std::this_thread::get_id();
        g_display_mutex.lock();
        std::cout << "ID main thread: " << this_id <<std::endl;
        g_display_mutex.unlock();*/

        /*size_t block_bytes = INT_CEIL(buff_size_, w_bytes)*w_bytes;
        buff->stream_size = block_bytes/w_bytes;
        //buff->stream = new sym_t[buffer->stream_size];
        block_bits = block_bytes*8;
        block_bg=0;
        //text_i.read((char*) buffer->stream, std::streamsize(std::min<size_t>(file_size, block_bytes)));*/
    }

    ~dbuff_i_file_stream(){
        needs_disk_acc=true;
        cv.notify_one();
        worker_thread_alive = false;
        dw.join();
        for(size_t i=0;i<n_buff;i++) free(buffers[i].stream);
        close_stream();
    }

    inline sym_t next(){
        buff_pos++;
        text_pos++;
        if(buff_pos>=buff_size){

            // wait for the worker
            {
                std::unique_lock<std::mutex> lk(m);
                cv.wait(lk, [&]{return !needs_disk_acc;});
            }

            buff_pos= 0;
            inv_buff = curr_buff;
            curr_buff = (curr_buff + 1) % n_buff;

            // send data to the worker thread
            {
                std::lock_guard<std::mutex> lk(m);
                needs_disk_acc = true;
            }
            cv.notify_one();
            //assert(!needs_disk_access);
            //std::cout<<"curr buff is :"<<curr_buff<<" "<<needs_disk_access<<std::endl;
        }
        return buffers[curr_buff].stream[buff_pos];
    }

    inline sym_t value() const{
        return buffers[curr_buff].stream[buff_pos];
    }

    dbuff_i_file_stream& operator=(dbuff_i_file_stream &&other) noexcept{
        move(std::forward<dbuff_i_file_stream>(other));
        return *this;
    }

    /*void move(dbuff_i_file_stream &&other){
        std::swap(file, other.file);
        //std::swap(text_i, other.text_i);
        std::swap(tot_syms, other.tot_cells);
        //std::swap(block_bg, other.block_bg);
        //buffer.swap(other.buffer);
    }*/

    [[nodiscard]] inline size_t size() const{
        return tot_syms;
    }

    [[nodiscard]] inline const std::string& filename() const{
        return file;
    }

    void close_stream(bool rem=false){
        //if(text_i.is_open()) text_i.close();
        close(fp);
        if(rem){
            if(remove(file.c_str())){
                std::cout<<"Error trying to remove file"<<file<<std::endl;
            }
        }
    }

    /*
    inline void read_bit_chunk(uint8_t *dst, size_t i, size_t j) {
        assert(i<=j);

        size_t cell_i = i >> w_shift;
        size_t cell_j = j >> w_shift;

        / *assert(cell_j<tot_cells && (cell_j-cell_i+1)<=buffer.stream_size);

        if(cell_i<block_bg || (block_bg+buffer.stream_size)<=cell_i){
            block_bg = (cell_i/buffer.stream_size)*buffer.stream_size;
            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
        }

        size_t start = i % block_bits;
        if((cell_i/buffer.stream_size) == (cell_j/buffer.stream_size)){ //same block
            return buffer.read_chunk(dst, start, start+(j-i));
        }else{ //different blocks

            buffer.read_chunk(dst, start, block_bits-1);
            size_t bits_read = block_bits-start;
            block_bg+=buffer.stream_size;

            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }

            size_t bytes_read = INT_CEIL(bits_read, 8);
            size_t offset = std::min((8*bytes_read), (j-i+1))-bits_read;
            dst+=bytes_read-1;
            if(offset>0){
                uint8_t tail = buffer.read(0, offset-1);
                uint8_t left = bits_read % 8;
                *dst &= ~(bitstream<uint8_t>::masks[offset]<<left);//clean the area
                *dst |= tail << left;//add new data in the area
                bits_read+=offset;
            }

            if(bits_read<(j-i+1)){
                dst++;
                size_t end = j % block_bits;
                buffer.read_chunk(dst, offset, end);
                bits_read+=end-offset+1;
            }
            assert(bits_read==(j-i+1));
        }* /
    }
    inline size_t read_bits(size_t i, size_t j) {

        / *assert(i<=j && (j-i+1)<=64);
        size_t cell_i = i >> w_shift;
        size_t cell_j = j >> w_shift;

        assert(cell_j<tot_cells && (cell_j-cell_i+1)<=buffer.stream_size);

        if(cell_i<block_bg || (block_bg+buffer.stream_size)<=cell_i){
            block_bg = (cell_i/buffer.stream_size)*buffer.stream_size;
            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
        }

        size_t start = i % block_bits;
        if((cell_i/buffer.stream_size) == (cell_j/buffer.stream_size)){ //same block
            return buffer.read(start, start+(j-i));
        }else{ //different blocks
            size_t data = buffer.read(start, block_bits-1);
            size_t delta = block_bits-start;
            block_bg+=buffer.stream_size;

            text_i.seekg(block_bg*w_bytes);
            text_i.read((char*) buffer.stream, buffer.stream_size*w_bytes);
            if((unsigned int)text_i.gcount()<(buffer.stream_size*w_bytes)){
                text_i.clear();
            }
            return (buffer.read(0, j % block_bits)<<delta) | data;
        }* /
    }*/
};
#endif //CDT_DBUFF_FILE_STREAM_H
