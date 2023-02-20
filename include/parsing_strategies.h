//
// Created by Diaz, Diego on 15.12.2022.
//

#ifndef GRLBWT_PARSING_STRATEGIES_H
#define GRLBWT_PARSING_STRATEGIES_H

#include <thread>
#include "common.h"

struct parsing_info{
    size_t lms_phrases=0; //number of LMS phrases in the current parsing round
    size_t tot_phrases=0; //total number of phrases generated in the current parsing round
    size_t p_round=0; //parsing round
    size_t prev_alph=0; //size of the alphabet in the previous parsing round
    size_t max_sym_freq=0; // most repeated symbol in the input text of the round
    size_t longest_str=0; //longest string in the parsing
    size_t active_strings=0; //number of active strings
    std::vector<long> str_ptrs;
};

template<class stream_t,
        class string_t,
        bool first_round>
struct lms_parsing{

    typedef stream_t                       stream_type;
    typedef typename stream_type::sym_type sym_type;

    static void compute_breaks(stream_t& ifs,
                               size_t f_string, size_t l_string,
                               const std::function<void(size_t)>& process_phrase,
                               const std::function<std::pair<long, long>(size_t)>& init_str) {

        sym_type curr_sym, prev_sym;
        size_t end_ps, start_ps;
        uint8_t type, rep;

        for(size_t str=l_string+1;str-->f_string;) {

            auto range = init_str(str);

            if(range.first<=range.second) { //if this is not true, it means the string was fully compressed

                start_ps = range.first;
                end_ps = range.second;

                prev_sym = ifs.read(end_ps);
                if constexpr (std::is_same<sym_type, uint8_t>::value){
                    rep = 3U;
                }else{
                    rep = prev_sym & 1U;
                    prev_sym >>=1UL;
                }

                type = 0;

                for (size_t i = end_ps; i-- > start_ps;) {

                    curr_sym = ifs.read(i);

                    if constexpr (!std::is_same<sym_type, uint8_t>::value){
                        rep = (rep << 1UL) | (curr_sym & 1UL);
                        curr_sym >>=1UL;
                    }

                    if (curr_sym != prev_sym) {
                        type = (type<<1UL) | (curr_sym < prev_sym);
                        if ((type & 3U) == 2 && (rep & 3U)==3U) {//LMS suffix
                            //process the previous phrase
                            process_phrase(i+1);
                        }
                    } else {
                        type = (type<<1UL) | (type & 1UL);
                    }
                    prev_sym = curr_sym;
                }
            }
        }
    }

    inline void operator()(stream_t& ifs,
                           size_t f_string, size_t l_string, size_t max_symbol,
                           std::function<void(string_t&)>&& process_phrase,
                           std::function<std::pair<long, long>(size_t)>&& init_str) const {

        sym_type curr_sym, prev_sym;
        string_t phrase(2, sym_width(max_symbol));
        size_t end_ps, start_ps;
        uint8_t type, rep;

        for(size_t str=l_string+1;str-->f_string;) {

            auto range = init_str(str);

            if(range.first<=range.second) { //if this is not true, it means the string was fully compressed

                start_ps = range.first;
                end_ps = range.second;

                prev_sym = ifs.read(end_ps);
                if constexpr (first_round){
                    rep = 3U;
                }else{
                    rep = prev_sym & 1U;
                    prev_sym >>=1UL;
                }

                phrase.push_back(prev_sym);
                type = 0;

                for (size_t i = end_ps; i-- > start_ps;) {

                    curr_sym = ifs.read(i);

                    if constexpr (!first_round){
                        rep = (rep << 1UL) | (curr_sym & 1UL);
                        curr_sym >>=1UL;
                    }

                    if (curr_sym != prev_sym) {
                        type = (type<<1UL) | (curr_sym < prev_sym);
                        if ((type & 3U) == 2 && (rep & 3U)==3U) {//LMS suffix
                            //process the previous phrase
                            assert(!phrase.empty());
                            process_phrase(phrase);

                            //create the new phrase
                            phrase.clear();
                            phrase.push_back(prev_sym);
                        }
                    } else {
                        type = (type<<1UL) | (type & 1UL);
                    }

                    phrase.push_back(curr_sym);
                    prev_sym = curr_sym;
                }

                assert(!phrase.empty());
                process_phrase(phrase);
                phrase.clear();
            }
        }
    }
};

template<class parser_type,
         template<class, class> class hash_functor,
         template<class, class, class> class parse_functor>
struct mt_parse_strat_t {//multi thread strategy

    typedef typename parser_type::stream_type istream_t;

    struct thread_worker_data_t{

        size_t               start_str{};
        size_t               end_str{};
        istream_t            ifs;
        std::string          o_file;
        phrase_map_t&        map;
        buffered_map_t       inner_map;
        std::vector<long>&   str_ptr;
        size_t               max_symbol{};
        size_t               active_strings=0;
        size_t               rb{};

        thread_worker_data_t(size_t start_, size_t end_, std::string& i_file, std::string& o_file_,
                             phrase_map_t& map_, std::vector<long>& str_ptr_,
                             const size_t &hb_size, void *hb_addr,
                             size_t max_symbol_): start_str(start_),
                                                  end_str(end_),
                                                  ifs(i_file, BUFFER_SIZE),
                                                  o_file(o_file_),
                                                  map(map_),
                                                  inner_map(hb_size, o_file + "_phrases", 0.7, hb_addr, map.description_bits()),
                                                  str_ptr(str_ptr_),
                                                  max_symbol(max_symbol_),
                                                  rb(str_ptr[end_str+1]-1){};
        thread_worker_data_t()=default;

        inline std::pair<long, long> str2range(size_t str) const {
            assert(start_str<=str && str<=end_str);
            size_t bg = str_ptr[str];
            size_t end = str ==end_str ? rb : str_ptr[str+1]-1;
            assert(bg<=end);
            return {bg, end};
        }
    };

    std::string         i_file;
    std::string         o_file;
    phrase_map_t        map;
    parsing_info&       p_info;
    std::vector<long>&  str_ptr;
    std::vector<thread_worker_data_t> threads_data;
    char *              buff_addr;
    size_t              text_size{};

    mt_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_, size_t hbuff_size,
                     size_t n_threads) : i_file(i_file_),
                                         o_file(o_file_),
                                         map(0.8, sym_width(INT_CEIL(p_info_.longest_str*sym_width(p_info_.tot_phrases),8)*8)),
                                         p_info(p_info_),
                                         str_ptr(p_info.str_ptrs) {

        std::vector<std::pair<size_t, size_t>> thread_ranges;
        size_t str_per_thread = INT_CEIL(p_info.str_ptrs.size()-1, n_threads);
        n_threads = INT_CEIL((p_info.str_ptrs.size()-1), str_per_thread);

        for(size_t i=0;i<n_threads;i++){
            thread_ranges.emplace_back(str_per_thread*i, std::min(str_per_thread*(i+1)-1, p_info.str_ptrs.size()-2));
        }
        threads_data.reserve(thread_ranges.size());

        // each thread has a hast table with a buffer of at least 8MB
        hbuff_size = std::max<size_t>(hbuff_size, BUFFER_SIZE*n_threads);

        // each thread has a buffer with the same size, and with an integral number of size_t words
        size_t hb_bytes = INT_CEIL(INT_CEIL(hbuff_size, n_threads), sizeof(size_t)) * sizeof(size_t);
        hbuff_size = hb_bytes*n_threads;

        //how many size_t cells we can fit in the buffer
        //size_t buff_cells = hbuff_size/sizeof(size_t);
        //number of bytes per thread
        //size_t hb_bytes = (buff_cells / thread_ranges.size()) * sizeof(size_t);

        buff_addr = (char *)malloc(hbuff_size);
        size_t k = 0;
        for (auto &range: thread_ranges) {
            std::string tmp_o_file = o_file.substr(0, o_file.size() - 5);
            tmp_o_file.append("_range_" + std::to_string(range.first) + "_" + std::to_string(range.second));
            threads_data.emplace_back(range.first, range.second, i_file, tmp_o_file, map, p_info.str_ptrs, hb_bytes,
                                      buff_addr+(hb_bytes*k), p_info.tot_phrases);
            k++;
        }
        assert(!threads_data.empty());
        text_size = threads_data[0].ifs.size();
    };

    mt_parse_strat_t()=default;

    std::pair<size_t, size_t> get_phrases() {

        std::vector<std::thread> threads(threads_data.size());
        hash_functor<thread_worker_data_t, parser_type> hf;

        for (size_t i = 0; i < threads_data.size(); i++) {
            threads[i] = std::thread(hf, std::ref(threads_data[i]));
        }

        for (size_t i = 0; i < threads_data.size(); i++) {
            threads[i].join();
        }

        p_info.active_strings = 0;
        for (size_t i = 0; i < threads_data.size(); i++) {
            threads_data[i].inner_map.flush();
            p_info.active_strings +=threads_data[i].active_strings;
        }

        free(buff_addr);
#ifdef __linux__
        malloc_trim(0);
#endif

        //join the different phrase files
        //std::cout <<"      Merging thread data into one single hash table" << std::flush;
        //auto start = std::chrono::steady_clock::now();
        auto join_res = join_thread_phrases();
        //auto end = std::chrono::steady_clock::now();
        //report_time(start, end, 3);
        return join_res;
    }

    std::pair<size_t, size_t> join_thread_phrases() {

        size_t dic_bits=0, freq, max_freq=0;
        std::string file;

        if(p_info.p_round>0){
            map.resize_table(prev_power_of_two(p_info.lms_phrases));
        }

        for(auto const& thread : threads_data) {

            file = thread.inner_map.dump_file();
            size_t tot_bytes = std::filesystem::file_size(file);
            if(tot_bytes==0) continue;

            size_t d_bits = thread.inner_map.description_bits();
            size_t value_bits = thread.inner_map.value_bits();
            size_t longest_key = thread.inner_map.longest_key();//in bits

            size_t buffer_size = std::max<size_t>(INT_CEIL(longest_key, 8), std::min<size_t>(tot_bytes, BUFFER_SIZE));
            //size_t buffer_size = std::max<size_t>(INT_CEIL(longest_key, 8), std::min<size_t>(tot_bytes, 1024));
            buffer_size = next_power_of_two(buffer_size);
            i_file_stream<size_t> data_disk_buffer(file, buffer_size);

            //TODO testing
            //std::ifstream text_i(file, std::ios_base::binary);
            //bitstream<ht_buff_t> bits;
            //bits.stream = (ht_buff_t*) malloc(INT_CEIL(tot_bytes,sizeof(ht_buff_t))*sizeof(ht_buff_t));
            //bits.stream_size = INT_CEIL(tot_bytes, sizeof(ht_buff_t));
            //text_i.read((char *)bits.stream, std::streamsize(tot_bytes));
            //auto * key_test= (uint8_t*) malloc(INT_CEIL(longest_key, bitstream<ht_buff_t>::word_bits)*sizeof(ht_buff_t));
            //

            //there is a region of the file that does not contain data
            // * * * | * 0 0 0 0 0 0 0 <- tot_bits
            //             | <- if next_bit falls in this region, the loop is still valid, but there is no more data
            size_t tot_bits = (tot_bytes*8)-8;
            size_t key_bits;
            size_t next_bit = 0;
            auto * key= (uint8_t*) malloc(INT_CEIL(longest_key, bitstream<ht_buff_t>::word_bits)*sizeof(ht_buff_t));

            while(next_bit<tot_bits) {

                key_bits = data_disk_buffer.read_bits(next_bit, next_bit+d_bits-1);

                //TODO testing
                //size_t key_bits_test = bits.read(next_bit, next_bit+d_bits-1);
                //assert(key_bits_test==key_bits);
                //

                assert(key_bits>0 && key_bits<=longest_key);

                key[INT_CEIL(key_bits, 8)-1] = 0;
                next_bit+=d_bits;
                data_disk_buffer.read_bit_chunk(key, next_bit, next_bit+key_bits-1);

                //TODO testing
                //key_test[INT_CEIL(key_bits_test, 8)-1] = 0;
                //bits.read_chunk(key_test, next_bit, next_bit+key_bits-1);
                //assert(memcmp(key, key_test, INT_CEIL(key_bits, 8))==0);
                //

                next_bit+=key_bits;
                freq = data_disk_buffer.read_bits(next_bit, next_bit+value_bits-1);

                //TODO testing
                //size_t freq_test = bits.read(next_bit, next_bit+value_bits-1);
                //assert(freq_test==freq);
                //
                next_bit+=value_bits;

                auto res = map.increment_value(key, key_bits, freq);
                if(res==freq) dic_bits+=key_bits;
                freq = res;
                if(freq>max_freq) max_freq = freq;

                //TODO this work for the opt BWT
                /*if(!res.second){
                    flag = (freq & 3UL);
                    freq>>=1UL;

                    size_t val;
                    map.get_value_from(res.first, val);

                    flag |= (val & 3UL);
                    freq += val>>2UL;
                    assert(flag<=3);
                    if(freq>max_freq) max_freq = freq;

                    val = (freq << 2UL) | flag;
                    map.insert_value_at(res.first, val);
                }else{
                    freq>>=2UL;
                    if(freq>max_freq) max_freq = freq;
                    dic_bits+=key_bits;
                }*/
            }

            //TODO testing
            //text_i.close();
            //free(bits.stream);
            //free(key_test);
            //

            data_disk_buffer.close(true);
            free(key);
        }
        map.shrink_databuff();
        return {dic_bits/ sym_width(p_info.tot_phrases), max_freq};
    }

    template<class o_sym_type>
    size_t parse_text() {

        std::vector<std::thread> threads(threads_data.size());
        parse_functor<thread_worker_data_t, parser_type, o_file_stream<o_sym_type>> pf;

        for(size_t i=0;i<threads_data.size();i++){
            threads[i] = std::thread(pf, std::ref(threads_data[i]));
        }
        for(size_t i=0;i<threads_data.size();i++) {
            threads[i].join();
        }

        //invert the chunks
        o_file_stream<o_sym_type> of(o_file, BUFFER_SIZE, std::ios::out);
        for(auto const& thread: threads_data){
            i_file_stream<o_sym_type> inv_chunk(thread.o_file, BUFFER_SIZE);
            //std::cout<<"whut? ---> "<<thread.o_file<<" "<<inv_chunk.size()<<std::endl;
            for(size_t i = inv_chunk.size();i-->0;){
                of.push_back(inv_chunk.read(i));
            }
            inv_chunk.close(true);
        }
        of.close();
        size_t psize = of.size();


        //join the phrases in one single file
        //size_t psize = join_parse_chunks();

        //update string pointers
        long acc=0, prev, str_len=0;
        for(size_t i=0;i<threads_data.size();i++){
            prev = p_info.str_ptrs[threads_data[i].start_str];
            for(size_t j=threads_data[i].start_str; j<=threads_data[i].end_str;j++){
                acc += (prev-p_info.str_ptrs[j]);
                prev = p_info.str_ptrs[j];

                p_info.str_ptrs[j] = acc;
                if(j>0 && (p_info.str_ptrs[j]-p_info.str_ptrs[j-1])>str_len){
                    str_len = p_info.str_ptrs[j]-p_info.str_ptrs[j-1];
                }
            }
            acc +=prev+1;
        }

        p_info.str_ptrs.back() = (long)psize;
        if((p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2])>str_len){
            str_len = (p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2]);
        }
        p_info.longest_str = str_len;

        //TODO testing
        /*for(size_t ctr=0;ctr<(p_info.str_ptrs.size()-1);ctr++){
            if(p_info.str_ptrs[ctr]<=(p_info.str_ptrs[ctr+1]-1)){
                std::cout<<p_info.str_ptrs[ctr]<<" - "<<(p_info.str_ptrs[ctr+1]-1)<<" "<<p_info.str_ptrs[ctr+1]-p_info.str_ptrs[ctr]<<" "<<p_info.longest_str<<std::endl;
            }
        }*/
        //

        //TODO testing
        /*if(p_info.p_round>4){
            for(size_t ctr=0;ctr<(p_info.str_ptrs.size()-1);ctr++){
                if(p_info.str_ptrs[ctr]<=(p_info.str_ptrs[ctr+1]-1)){
                    std::cout<<p_info.str_ptrs[ctr]<<" - "<<(p_info.str_ptrs[ctr+1]-1)<<" | ";
                    for(long k=p_info.str_ptrs[ctr];k<=(p_info.str_ptrs[ctr+1]-1);k++){
                        std::cout<<(ifs.read(k)>>1UL)<<" ";
                    }
                    std::cout<<""<<std::endl;
                }
            }
            size_t ctr=1;
            for(long i=0;i<(long)ifs.size();i++){
                if(ifs.read(i)<tot_phrases && phrase_desc[ifs.read(i)]){
                    if((i+1)!=p_info.str_ptrs[ctr]){
                        while((i+1)>p_info.str_ptrs[ctr]){
                            if(p_info.str_ptrs[ctr]>(p_info.str_ptrs[ctr+1]-1)){
                                std::cout<<" xxx "<<p_info.str_ptrs[ctr]<<" "<<p_info.str_ptrs[ctr+1]-1<<std::endl;
                            }else{
                                std::cout<<" *** "<<p_info.str_ptrs[ctr]<<" "<<p_info.str_ptrs[ctr+1]-1<<" | ";
                                for(long k=p_info.str_ptrs[ctr];k<=(p_info.str_ptrs[ctr+1]-1);k++){
                                    std::cout<<ifs.read(k)<<" ";
                                }
                                std::cout<<""<<std::endl;
                            }
                            ctr++;
                        }
                        if((i+1)==p_info.str_ptrs[ctr]){
                            std::cout<<(i+1)<<" "<<p_info.str_ptrs[ctr]<<" - "<<(p_info.str_ptrs[ctr+1]-1)<<" | ";
                            for(long k=p_info.str_ptrs[ctr];k<=(p_info.str_ptrs[ctr+1]-1);k++){
                                std::cout<<ifs.read(k)<<" ";
                            }
                            ctr++;
                        }else{
                            std::cout<<"wrong!"<<std::endl;
                        }
                    }else{
                        std::cout<<(i+1)<<" "<<p_info.str_ptrs[ctr]<<" - "<<p_info.str_ptrs[ctr+1]-1<<" | "<<ifs.size()<<" "<<ctr<<std::endl;
                        for(long k=p_info.str_ptrs[ctr];k<=(p_info.str_ptrs[ctr+1]-1);k++){
                            std::cout<<ifs.read(k)<<" ";
                        }
                        std::cout<<""<<std::endl;
                        ctr++;
                    }
                }
            }
        }*/
        //
        return psize;
    }

    size_t join_parse_chunks() {

        //concatenate the files
        std::ofstream of(o_file, std::ofstream::binary);
        size_t buff_size = BUFFER_SIZE/sizeof(size_t);
        size_t len, rem, to_read, start, end, tot_syms=0;
        auto *buffer = new size_t[buff_size];
        std::string file;

        for(auto const& thread: threads_data){

            file = thread.o_file;
            std::ifstream tmp_i_file(file, std::ifstream::binary);

            tmp_i_file.seekg (0, std::ifstream::end);
            len = tmp_i_file.tellg()/sizeof(size_t);
            tmp_i_file.seekg (0, std::ifstream::beg);

            if(len>0){
                rem=len;
                to_read = std::min<size_t>(buff_size, len);

                while(true){
                    tmp_i_file.seekg( (rem - to_read) * sizeof(size_t));
                    tmp_i_file.read((char *)buffer, sizeof(size_t)*to_read);
                    assert(tmp_i_file.good());

                    //invert data
                    start =0;end=to_read-1;
                    while(start<end){
                        std::swap(buffer[start++], buffer[end--]);
                    }

                    of.write((char *)buffer, sizeof(size_t)*to_read);
                    assert(of.good());

                    rem -= tmp_i_file.gcount()/sizeof(size_t);
                    to_read = std::min<size_t>(buff_size, rem);
                    if(to_read == 0) break;
                }
                tmp_i_file.close();
                tot_syms+=len;
            }

            if(remove(file.c_str())){
                std::cout<<"Error trying to remove temporal file"<<std::endl;
                std::cout<<"Aborting"<<std::endl;
                exit(1);
            }
        }
        delete[] buffer;
        of.close();
        return tot_syms;
    }

    /*void remove_files(){
        //remove remaining files
        for(size_t i=0;i<threads_data.size();i++){
            std::string tmp_file =  threads_data[i].ofs.file;
            if(remove(tmp_file.c_str())){
                std::cout<<"Error trying to delete file "<<tmp_file<<std::endl;
            }
        }
    }*/
};

template<class parser_type,
         template<class, class> class hash_functor,
         template<class, class, class> class parse_functor>
struct st_parse_strat_t {//parse data for single thread

    typedef parser_type                    parser_t;
    typedef typename parser_t::stream_type istream_t;

    istream_t           ifs;
    std::string         o_file;
    std::string         tmp_o_file;

    phrase_map_t        map;
    phrase_map_t&       inner_map;
    parsing_info&       p_info;
    std::vector<long>&  str_ptr;
    size_t              text_size{};
    size_t              max_symbol{};
    size_t              start_str;
    size_t              end_str;
    size_t              active_strings=0;

    st_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_): ifs(i_file_, BUFFER_SIZE),
                                             o_file(o_file_),
                                             map(0.8, sym_width(INT_CEIL(p_info_.longest_str*sym_width(p_info_.tot_phrases),8)*8)),
                                             inner_map(map),
                                             p_info(p_info_),
                                             str_ptr(p_info.str_ptrs),
                                             text_size(ifs.size()),
                                             max_symbol(p_info.tot_phrases),
                                             start_str(0),
                                             end_str(str_ptr.size()-2){
        tmp_o_file = o_file.substr(0, o_file.size() - 4);
        tmp_o_file.append("_inv");

        if(p_info.p_round>1){
            //auto n_buckets = std::max<size_t>(size_t(double(ifs.size())*0.2*0.4), 4);
            //n_buckets = next_power_of_two(n_buckets);
            //map.resize_table(n_buckets);
            size_t n_buckets = prev_power_of_two(p_info.lms_phrases);
            map.resize_table(n_buckets);
        }
    }

    [[nodiscard]] inline std::pair<long, long> str2range(size_t str) const {
        assert(start_str<=str && str<=end_str);
        size_t bg = str_ptr[str];
        size_t end = str_ptr[str+1]-1;
        assert(bg<=end);
        return {bg, end};
    }

    std::pair<size_t, size_t> get_phrases() {

        hash_functor<st_parse_strat_t, parser_type>()(*this);

        map.shrink_databuff();
        key_wrapper key_w{sym_width(max_symbol), map.description_bits(), map.get_data()};
        size_t n_syms=0, max_freq=0, freq;
        p_info.active_strings = active_strings;

        for(auto const &ptr : map){
            n_syms += key_w.size(ptr);
            freq = 0;
            map.get_value_from(ptr, freq);

            //TODO testing
            if(freq==0){
                std::cout<<"here I have a bug at "<<ptr<<" "<<freq<<std::endl;
            }
            //

            assert(freq>0);
            if(freq>max_freq) max_freq = freq;
        }
        return {n_syms, max_freq};
    }

    template<class o_sym_type>
    size_t parse_text() {

        size_t parse_size = parse_functor<st_parse_strat_t, parser_type, o_file_stream<o_sym_type>>()(*this);

        //update string pointers
        long acc=0, str_len=0;
        size_t prev = p_info.str_ptrs[start_str];
        for(size_t j=start_str; j<=end_str; j++){
            acc += (prev-p_info.str_ptrs[j]);
            prev = p_info.str_ptrs[j];
            p_info.str_ptrs[j] = acc;
            if(j>0 && (p_info.str_ptrs[j]-p_info.str_ptrs[j-1])>str_len){
                str_len = p_info.str_ptrs[j]-p_info.str_ptrs[j-1];
            }
        }
        p_info.str_ptrs.back() = (long)parse_size;
        if((p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2])>str_len){
            str_len = (p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2]);
        }
        p_info.longest_str = str_len;

        i_file_stream<o_sym_type> inv_parse(o_file, BUFFER_SIZE);
        o_file_stream<o_sym_type> final_parse(tmp_o_file, BUFFER_SIZE, std::ios::out);
        for(size_t i=inv_parse.size();i-->0;){
            final_parse.push_back(inv_parse.read(i));
        }
        inv_parse.close(true);
        final_parse.close();
        rename(tmp_o_file.c_str(), o_file.c_str());

        return parse_size;
    }
};

typedef i_file_stream<uint8_t>                        uint8t_i_stream;
typedef i_file_stream<uint16_t>                       uint16t_i_stream;
typedef i_file_stream<uint32_t>                       uint32t_i_stream;
typedef i_file_stream<uint64_t>                       uint64t_i_stream;

typedef lms_parsing<uint8t_i_stream, string_t, true>    char_parser_t;
typedef lms_parsing<uint8t_i_stream, string_t, false>   uint8t_parser_t;
typedef lms_parsing<uint16t_i_stream, string_t, false>  uint16t_parser_t;
typedef lms_parsing<uint32t_i_stream, string_t, false>  uint32t_parser_t;
typedef lms_parsing<uint64t_i_stream, string_t, false>  uint64t_parser_t;

#endif //GRLBWT_PARSING_STRATEGIES_H
