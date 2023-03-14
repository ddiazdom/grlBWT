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
    size_t text_size=0;
    std::vector<long> str_ptrs;
};

template<class stream_t,
        class string_t,
        bool first_round,
        bool forward=true>
struct lms_parsing {

public:
    typedef stream_t                       stream_type;
    typedef typename stream_type::sym_type sym_type;

private:
    static void forward_parsing(stream_t& ifs,
               size_t f_string, size_t l_string, size_t max_symbol,
               std::function<void(string_t&)>&& process_phrase,
               std::function<std::pair<long, long>(size_t)>&& init_str){

        bool phrase_break;
        sym_type curr_sym, prev_sym;
        string_t phrase(2, sym_width(max_symbol));
        size_t prev_i, tmp_sym, end_ps, start_ps;
        uint8_t rep;

        for(size_t str=f_string;str<=l_string;str++) {
            auto range = init_str(str);
            if(range.first<=range.second) { //if this is not true, it means the string was fully compressed

                start_ps = range.first;
                end_ps = range.second;
                prev_sym = 0;

                size_t i = start_ps;
                if constexpr (first_round){
                    rep = 3UL;
                }else{
                    rep = 0;
                }

                while(i<=end_ps) {

                    curr_sym = ifs[i];
                    if constexpr (first_round){
                        tmp_sym = curr_sym;
                    }else{
                        rep = (rep << 1UL) | (curr_sym & 1UL);
                        tmp_sym = curr_sym >> 1UL;
                    }

                    phrase.push_back(tmp_sym);
                    prev_i=i;
                    while(i<end_ps &&  curr_sym==ifs[i+1]) i++;

                    if constexpr (first_round){
                        phrase_break = (i==end_ps) || ((rep & 3UL) == 3 && prev_sym>tmp_sym && tmp_sym<ifs[i+1]);
                    }else{
                        phrase_break = (i==end_ps) || ((rep & 3UL) == 3 && prev_sym>tmp_sym && tmp_sym<(ifs[i+1]>>1UL));
                    }

                    if(phrase_break){
                        //process the phrase
                        assert(!phrase.empty());
                        process_phrase(phrase);

                        //create the new phrase
                        phrase.clear();
                        phrase.push_back(tmp_sym);
                    }

                    for(size_t j=prev_i;j<i;j++){
                        phrase.push_back(tmp_sym);
                    }

                    prev_sym = tmp_sym;
                    i++;
                }
                phrase.clear();
            }
        }

    }

    static void reverse_parsing(stream_t& ifs,
               size_t f_string, size_t l_string, size_t max_symbol,
               std::function<void(string_t&)>&& process_phrase,
               std::function<std::pair<long, long>(size_t)>&& init_str){

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

                    curr_sym = ifs[i];

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

public:
    inline void operator()(stream_t& ifs,
               size_t f_string, size_t l_string, size_t max_symbol,
               std::function<void(string_t&)>&& process_phrase,
               std::function<std::pair<long, long>(size_t)>&& init_str) const {

        if constexpr (forward){
            forward_parsing(ifs, f_string, l_string, max_symbol,
                            std::forward<std::function<void(string_t&)>>(process_phrase),
                            std::forward<std::function<std::pair<long, long>(size_t)>>(init_str));
        }else{
            reverse_parsing(ifs, f_string, l_string, max_symbol,
                            std::forward<std::function<void(string_t&)>>(process_phrase),
                            std::forward<std::function<std::pair<long, long>(size_t)>>(init_str));
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
        size_t               n_phrases=0;
        size_t               offset=0;
        size_t               rb{};

        thread_worker_data_t(size_t start_, size_t end_, std::string& i_file, std::string& o_file_,
                             phrase_map_t& map_, std::vector<long>& str_ptr_,
                             const size_t &hb_size, void *hb_addr,
                             size_t max_symbol_): start_str(start_),
                                                  end_str(end_),
                                                  ifs(i_file,
                                                      BUFFER_SIZE
#ifdef __linux__
, POSIX_FADV_SEQUENTIAL
#endif
                                                      ),
                                                  o_file(o_file_),
                                                  map(map_),
                                                  inner_map(hb_size, o_file +"_"+std::to_string(start_)+"_"+std::to_string(end_)+"_phrases", 0.7, hb_addr, map.description_bits()),
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

    phrase_map_t        map;
    parsing_info&       p_info;
    std::vector<thread_worker_data_t> threads_data;
    char *              buff_addr=nullptr;

    mt_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_, size_t hbuff_size,
                     size_t n_threads) : map(0.8, sym_width(INT_CEIL(p_info_.longest_str*sym_width(p_info_.tot_phrases),8)*8)),
                                         p_info(p_info_){

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

        buff_addr = (char *)malloc(hbuff_size);
        size_t k = 0;
        for (auto &range: thread_ranges) {
            threads_data.emplace_back(range.first, range.second, i_file_, o_file_, map, p_info.str_ptrs, hb_bytes,
                                      buff_addr+(hb_bytes*k), p_info.tot_phrases);
            k++;
        }
        assert(!threads_data.empty());
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

        auto join_res = join_thread_phrases();
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
            i_file_stream<size_t> data_disk_buffer(file, buffer_size
#ifdef __linux__
                                                  ,POSIX_ADV_SEQUENTIAL
#endif
);

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

            data_disk_buffer.close_reader(true);
            free(key);
        }
        map.shrink_databuff();
        return {dic_bits/ sym_width(p_info.tot_phrases), max_freq};
    }

    template<class o_sym_type>
    size_t parse_text() {

        size_t offset=0;
        for(size_t i=0;i<threads_data.size();i++){
            threads_data[i].offset = offset;
            offset += threads_data[i].n_phrases;
        }

        std::vector<std::thread> threads(threads_data.size());
        parse_functor<thread_worker_data_t, parser_type, o_file_stream<o_sym_type>> pf;

        for(size_t i=0;i<threads_data.size();i++){
            threads[i] = std::thread(pf, std::ref(threads_data[i]));
        }
        for(size_t i=0;i<threads_data.size();i++) {
            threads[i].join();
        }

        size_t psize = offset;
        p_info.str_ptrs.back() = (long)psize;

        size_t str_len;
        p_info.longest_str = 0;
        for(size_t i=0;i<p_info.str_ptrs.size()-1;i++){
            str_len = size_t(p_info.str_ptrs[i+1]-p_info.str_ptrs[i]);
            if(p_info.longest_str < str_len) p_info.longest_str = str_len;
        }

        return psize;
    }
};

template<class parser_type,
         template<class, class> class hash_functor,
         template<class, class, class> class parse_functor>
struct st_parse_strat_t {//parse data for single thread

    typedef parser_type                    parser_t;
    typedef typename parser_t::stream_type istream_t;

    istream_t           ifs;
    std::string         o_file;

    phrase_map_t        map;
    phrase_map_t&       inner_map;
    parsing_info&       p_info;
    std::vector<long>&  str_ptr;
    size_t              max_symbol{};
    size_t              start_str;
    size_t              end_str;
    size_t              n_phrases=0;
    size_t              offset=0;
    size_t              active_strings=0;

    st_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_): ifs(i_file_,
                                                 BUFFER_SIZE
#ifdef __linux__
, POSIX_ADV_SEQUENTIAL
#endif
                                                 ),
                                             o_file(o_file_),
                                             map(0.8, sym_width(INT_CEIL(p_info_.longest_str*sym_width(p_info_.tot_phrases),8)*8)),
                                             inner_map(map),
                                             p_info(p_info_),
                                             str_ptr(p_info.str_ptrs),
                                             max_symbol(p_info.tot_phrases),
                                             start_str(0),
                                             end_str(str_ptr.size()-2){
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
            assert(freq>0);
            if(freq>max_freq) max_freq = freq;
        }
        return {n_syms, max_freq};
    }

    template<class o_sym_type>
    size_t parse_text() {

        offset = 0;
        size_t parse_size = parse_functor<st_parse_strat_t, parser_type, o_file_stream<o_sym_type>>()(*this);
        p_info.str_ptrs.back() = (long)parse_size;

        size_t str_len;
        p_info.longest_str = 0;
        for(size_t j=0; j<str_ptr.size()-1; j++){
            str_len = size_t(p_info.str_ptrs[j+1]-p_info.str_ptrs[j]);
            if(str_len>p_info.longest_str) p_info.longest_str = str_len;
        }

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
