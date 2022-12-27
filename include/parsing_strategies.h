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
    std::vector<long> str_ptrs;
};

template<typename parse_data_t,
         typename parser_t>
struct hash_functor{
    void operator()(parse_data_t& data) {

        auto hash_phrase = [&](string_t& phrase) -> void {
            phrase.mask_tail();
            auto res = data.inner_map.insert(phrase.data(), phrase.n_bits(), 1);

            if(!res.second){
                size_t val;
                data.inner_map.get_value_from(res.first, val);
                val++;
                data.inner_map.insert_value_at(res.first, val);
            }
        };

        auto init_str = [&](size_t str) -> std::pair<long, long>{
            return {data.str_ptr[str], data.str_ptr[str+1]-1};
        };

        parser_t()(data.ifs, data.start, data.end, data.max_symbol, hash_phrase, init_str);
        //pthread_exit(nullptr);
    };
};

template<typename parse_data_t,
        typename parser_t>
struct counting_functor {
    void operator()(parse_data_t& data) {

        std::cout<< "    Computing the phrases in the text" << std::flush;
        //std::vector<std::pair<size_t, size_t>> tails;

        o_file_stream<size_t> breaks_list("breaks.txt", BUFFER_SIZE, std::ios::out);

        string_t suffix(2, sym_width(data.max_symbol));
        phrase_map_t map(0, "", 0.8, nullptr, 8);

        auto record_break = [&](size_t id) -> void {
            //tails.emplace_back(data.ifs.read(id-1), data.ifs.read(id));
            breaks_list.push_back(id);
            suffix.write(0, data.ifs.read(id-1));
            suffix.write(1, data.ifs.read(id));
            auto res = map.insert(suffix.data(), suffix.n_bits(), id<<1UL);
            if(!res.second){
                size_t first_occ = 0;
                map.get_value_from(res.first, first_occ);
                map.insert_value_at(res.first, first_occ | 1UL);
            }
        };

        auto init_str = [&](size_t str) -> std::pair<long, long>{
            return {data.str_ptr[str], data.str_ptr[str+1]-1};
        };

        parser_t::compute_breaks(data.ifs, data.start, data.end, record_break, init_str);
        breaks_list.close();

        o_file_stream<size_t> new_breaks_list("breaks_tmp.txt", BUFFER_SIZE, std::ios::out);
        i_file_stream<size_t> old_breaks_list("breaks.txt", BUFFER_SIZE);

        size_t pos, n_uniq=0, j=0;
        for(auto const &ptr : map ){
            pos=0;
            map.get_value_from(ptr, pos);
            if(!(pos & 1UL)){
                pos >>=1UL;
                while(j<old_breaks_list.size() && old_breaks_list.read(j)!=pos){
                    new_breaks_list.push_back(old_breaks_list.read(j));
                    j++;
                }
                j++;
                n_uniq++;
            }
        }
        while(j<old_breaks_list.size()){
            new_breaks_list.push_back(old_breaks_list.read(j));
            j++;
        }
        assert(new_breaks_list.size()==(old_breaks_list.size()-n_uniq));
        old_breaks_list.close();
        std::cout<<"There are "<<n_uniq<<" unnecessary breaks out of "<<new_breaks_list.size()<<std::endl;
        exit(0);

        /*std::sort(tails.begin(), tails.end(), [&](auto a, auto b){
            if(a.second!=b.second){
                return a.second < b.second;
            }else{
                return a.first < b.first;
            }
        });

        size_t l_sym = tails[0].second;
        size_t r_sym = tails[1].first;
        size_t freq=1, pos=0;

        for(size_t i=1;i<tails.size();i++) {
            if(tails[i].second!=l_sym){
                if(freq==1) tails[pos++] = tails[i-1];
                freq=0;
                r_sym = tails[i].first;
                l_sym = tails[i].second;
            }else{
                if(tails[i].first!=r_sym){
                    if(freq<=10) tails[pos++] = tails[i-1];
                    r_sym = tails[i].first;
                    freq=0;
                }
            }
            freq++;
            //std::cout<<tails[i].first<<" "<<tails[i].second<<std::endl;
        }
        if(freq==10) tails[pos++] = tails.back();
        tails.resize(pos);
        std::cout<<"\nThere are "<<tails.size()<<" we can avoid "<<std::endl;*/
    };
};

template<typename parse_data_t,
        typename parser_t>
struct parse_functor{

    void operator()(parse_data_t& data) {
        size_t str_len;

        auto phrase2symbol = [&](string_t& phrase){
            phrase.mask_tail();
            auto res = data.map.find(phrase.data(), phrase.n_bits());
            assert(res.second);
            size_t sym = 0;
            data.map.get_value_from(res.first, sym);

            if(phrase.size()<str_len){ //when (sym & 1UL) is true, it means there are > 1 copies of a string in the input
                data.ofs.push_back(sym);
            }
        };

        auto init_str = [&](size_t str) -> std::pair<long, long>{

            assert(str>=data.start && str<=data.end);
            size_t start = data.str_ptr[str];
            size_t end = data.str_ptr[str+1]-1;

            str_len = end-start+1;

            if((str+1)<=data.end){
                data.str_ptr[str+1] = data.ofs.size()-1;
            }
            return {start, end};
        };

        parser_t()(data.ifs, data.start, data.end, data.max_symbol, phrase2symbol, init_str);

        data.str_ptr[data.start] = data.ofs.size()-1;

        data.ofs.close();
        data.ifs.close();
        //pthread_exit(nullptr);
    };
};

template<class parser_type, class ostream_t>
struct mt_parse_strat_t {//multi thread strategy

    typedef typename parser_type::stream_type istream_t;
    typedef typename parser_type::sym_type    sym_type;

    struct thread_worker_data_t{

        size_t               start{};
        size_t               end{};
        istream_t            ifs;
        ostream_t            ofs;
        phrase_map_t&        map;
        phrase_map_t         inner_map;
        std::vector<long>&   str_ptr;
        size_t               max_symbol{};

        thread_worker_data_t(size_t start_, size_t end_, std::string& i_file, std::string& o_file,
                             phrase_map_t& map_, std::vector<long>& str_ptr_,
                             const size_t &hb_size, void *hb_addr,
                             size_t max_symbol_): start(start_),
                                                  end(end_),
                                                  ifs(i_file, BUFFER_SIZE),
                                                  ofs(o_file, BUFFER_SIZE, std::ios::out),
                                                  map(map_),
                                                  inner_map(hb_size, o_file + "_phrases", 0.7, hb_addr, map.description_bits()),
                                                  str_ptr(str_ptr_),
                                                  max_symbol(max_symbol_){};
        thread_worker_data_t()=default;
    };

    std::string         i_file;
    std::string         o_file;
    phrase_map_t        map;
    parsing_info&       p_info;
    std::vector<long>&  str_ptr;
    std::vector<thread_worker_data_t> threads_data;
    void *buff_addr=nullptr;
    size_t              text_size{};

    mt_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_, const size_t &hbuff_size,
                     size_t n_threads) : i_file(i_file_),
                                         o_file(o_file_),
                                         map(0, "", 0.8, nullptr, sym_width(INT_CEIL(p_info_.longest_str*sym_width(p_info_.tot_phrases),8)*8)),
                                         p_info(p_info_),
                                         str_ptr(p_info.str_ptrs){

        std::vector<std::pair<size_t, size_t>> thread_ranges;
        size_t str_per_thread = INT_CEIL(p_info.str_ptrs.size()-1, n_threads);

        for(size_t i=0;i<n_threads;i++){
            thread_ranges.emplace_back(str_per_thread*i, std::min(str_per_thread*(i+1)-1, p_info.str_ptrs.size()-2));
        }
        threads_data.reserve(thread_ranges.size());

        //how many size_t cells we can fit in the buffer
        size_t buff_cells = hbuff_size/sizeof(size_t);

        //number of bytes per thread
        size_t hb_bytes = (buff_cells / thread_ranges.size()) * sizeof(size_t);

        buff_addr = malloc(hbuff_size);
        auto tmp_addr = reinterpret_cast<char *>(buff_addr);
        size_t k = 0;
        for (auto &range: thread_ranges) {
            std::string tmp_o_file = o_file.substr(0, o_file.size() - 5);
            tmp_o_file.append("_range_" + std::to_string(range.first) + "_" + std::to_string(range.second));
            threads_data.emplace_back(range.first, range.second, i_file, tmp_o_file, map, p_info.str_ptrs, hb_bytes,
                                      tmp_addr + (k * hb_bytes), p_info.tot_phrases);
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

        for (size_t i = 0; i < threads_data.size(); i++) {
            threads_data[i].inner_map.flush();
        }
        free(buff_addr);

#ifdef __linux__
        malloc_trim(0);
#endif

        //join the different phrase files
        std::cout << "      Merging thread data into one single hash table" << std::flush;
        auto start = std::chrono::steady_clock::now();
        auto join_res = join_thread_phrases();
        auto end = std::chrono::steady_clock::now();
        report_time(start, end, 3);
        return join_res;
    }

    std::pair<size_t, size_t> join_thread_phrases() {

        size_t dic_bits=0, freq, max_freq=0;
        std::string file;

        for(auto const& thread : threads_data){

            file = thread.inner_map.dump_file();
            std::ifstream text_i(file, std::ios_base::binary);

            text_i.seekg (0, std::ifstream::end);
            size_t tot_bytes = text_i.tellg();
            text_i.seekg (0, std::ifstream::beg);

            if(tot_bytes==0) continue;

            auto buffer = reinterpret_cast<char *>(malloc(tot_bytes));

            text_i.read(buffer, std::streamsize(tot_bytes));

            bitstream<ht_buff_t> bits;
            bits.stream = reinterpret_cast<ht_buff_t*>(buffer);

            size_t d_bits = map.description_bits();
            size_t next_bit = d_bits;

            //there is a region of the file that does not contain data
            // * * * | * 0 0 0 0 0 0 0 <- tot_bits
            //             | <- if next_bit falls in this region, the loop is still valid, but there is no more data
            size_t tot_bits = (tot_bytes*8)-8;
            size_t key_bits;
            size_t value_bits=map.value_bits();

            void* key=nullptr;
            size_t max_key_bits=0;

            while(next_bit<tot_bits){

                key_bits = bits.read(next_bit-d_bits, next_bit-1);
                assert(key_bits>0);

                size_t n_bytes = INT_CEIL(key_bits, bitstream<ht_buff_t>::word_bits)*sizeof(ht_buff_t);
                if(key_bits>max_key_bits){
                    if(key==nullptr){
                        key = malloc(n_bytes);
                    }else {
                        key = realloc(key, n_bytes);
                    }
                    max_key_bits = key_bits;
                }

                char *tmp = reinterpret_cast<char*>(key);
                tmp[INT_CEIL(key_bits, 8)-1] = 0;

                bits.read_chunk(key, next_bit, next_bit+key_bits-1);
                next_bit+=key_bits;
                freq = bits.read(next_bit, next_bit+value_bits-1);
                next_bit+=value_bits+d_bits;

                auto res = map.insert(key, key_bits, freq);
                if(!res.second){
                    size_t val;
                    map.get_value_from(res.first, val);
                    val+=freq;
                    if(val>max_freq) max_freq = val;
                    map.insert_value_at(res.first, val);
                }else{
                    if(freq>max_freq) max_freq = freq;
                    dic_bits+=key_bits;
                }
            }
            text_i.close();

            if(remove(file.c_str())){
                std::cout<<"Error trying to remove temporal file"<<std::endl;
                std::cout<<"Aborting"<<std::endl;
                exit(1);
            }
            free(key);
            free(buffer);
        }
        map.shrink_databuff();

        return {dic_bits/ sym_width(p_info.tot_phrases), max_freq};
    }

    size_t parse_text() {

        std::vector<std::thread> threads(threads_data.size());
        parse_functor<thread_worker_data_t, parser_type> pf;
        for(size_t i=0;i<threads_data.size();i++){
            threads[i] = std::thread(pf, std::ref(threads_data[i]));
        }
        for(size_t i=0;i<threads_data.size();i++) {
            threads[i].join();
        }

        //join the phrases in one single file
        size_t psize = join_parse_chunks();

        //update string pointers
        long acc=0, prev, str_len=0;
        for(size_t i=0;i<threads_data.size();i++){
            prev = p_info.str_ptrs[threads_data[i].start];
            for(size_t j=threads_data[i].start; j<=threads_data[i].end;j++){
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

            file = thread.ofs.file;
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

    void remove_files(){
        //remove remaining files
        for(size_t i=0;i<threads_data.size();i++){
            std::string tmp_file =  threads_data[i].ofs.file;
            if(remove(tmp_file.c_str())){
                std::cout<<"Error trying to delete file "<<tmp_file<<std::endl;
            }
        }
    }
};

template<class parser_type, class ostream_t>
struct st_parse_strat_t {//parse data for single thread

    typedef parser_type                    parser_t;
    typedef typename parser_t::stream_type istream_t;
    typedef typename parser_type::sym_type sym_type;

    istream_t           ifs;
    ostream_t           ofs;
    std::string         o_file;

    phrase_map_t        map;
    phrase_map_t&       inner_map;
    parsing_info&       p_info;
    std::vector<long>&  str_ptr;
    size_t              text_size{};
    size_t              max_symbol{};
    size_t              start;
    size_t              end;

    st_parse_strat_t(std::string &i_file_, std::string& o_file_,
                     parsing_info& p_info_): ifs(i_file_, BUFFER_SIZE),
                                             o_file(o_file_),
                                             map(0, "", 0.8, nullptr, sym_width(INT_CEIL(p_info_.longest_str*sym_width(p_info_.tot_phrases),8)*8)),
                                             inner_map(map),
                                             p_info(p_info_),
                                             str_ptr(p_info.str_ptrs),
                                             text_size(ifs.size()),
                                             max_symbol(p_info.tot_phrases),
                                             start(0),
                                             end(str_ptr.size()-2){

        std::string tmp_o_file = o_file.substr(0, o_file.size() - 5);
        tmp_o_file.append("_range_" + std::to_string(start) + "_" + std::to_string(end));
        ofs = ostream_t(tmp_o_file, BUFFER_SIZE, std::ios::out);
    }

    std::pair<size_t, size_t> get_phrases() {

        auto s = std::chrono::steady_clock::now();
        hash_functor<st_parse_strat_t, parser_type>()(*this);
        auto e = std::chrono::steady_clock::now();
        report_time(s, e, 3);

        map.shrink_databuff();
        key_wrapper key_w{sym_width(max_symbol), map.description_bits(), map.get_data()};
        size_t n_syms=0, max_freq=0, freq;
        for(auto const &ptr : map){
            n_syms += key_w.size(ptr);
            freq = 0;
            map.get_value_from(ptr, freq);
            assert(freq>0);
            if(freq>max_freq) max_freq = freq;
        }
        return {n_syms, max_freq};
    }

    size_t parse_text() {

        parse_functor<st_parse_strat_t, parser_type>()(*this);

        //update string pointers
        long acc=0, str_len=0;
        size_t prev = p_info.str_ptrs[start];
        for(size_t j=start; j<=end; j++){
            acc += (prev-p_info.str_ptrs[j]);
            prev = p_info.str_ptrs[j];
            p_info.str_ptrs[j] = acc;
            if(j>0 && (p_info.str_ptrs[j]-p_info.str_ptrs[j-1])>str_len){
                str_len = p_info.str_ptrs[j]-p_info.str_ptrs[j-1];
            }
        }
        p_info.str_ptrs.back() = (long)ofs.size();
        if((p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2])>str_len){
            str_len = (p_info.str_ptrs.back() - p_info.str_ptrs[p_info.str_ptrs.size()-2]);
        }
        p_info.longest_str = str_len;

        size_t tot_sym = reverse_text();

        return tot_sym;
    }

    size_t reverse_text(){

        std::ofstream of(o_file, std::ofstream::binary);
        size_t buff_size = BUFFER_SIZE/sizeof(size_t);

        size_t len, rem, to_read, p_start, p_end;
        auto *buffer = new size_t[buff_size];

        std::string file = ofs.file;
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
                p_start =0;
                p_end=to_read-1;
                while(p_start<p_end){
                    std::swap(buffer[p_start++], buffer[p_end--]);
                }

                of.write((char *)buffer, sizeof(size_t)*to_read);
                assert(of.good());

                rem -= tmp_i_file.gcount()/sizeof(size_t);
                to_read = std::min<size_t>(buff_size, rem);
                if(to_read == 0) break;
            }
            tmp_i_file.close();
        }

        if(remove(file.c_str())){
            std::cout<<"Error trying to remove temporal file"<<std::endl;
            std::cout<<"Aborting"<<std::endl;
            exit(1);
        }

        delete[] buffer;
        of.close();
        return len;
    }

    void remove_files(){
        //remove remaining files
        std::string tmp_file =  ofs.file;
        if(remove(tmp_file.c_str())){
            std::cout<<"Error trying to delete file "<<tmp_file<<std::endl;
        }
    }
};

typedef i_file_stream<uint8_t>                        byte_i_stream;
typedef i_file_stream<size_t>                         int_i_stream;
typedef o_file_stream<size_t>                         int_o_stream;

typedef lms_parsing<byte_i_stream, string_t>          byte_parser_t;
typedef lms_parsing<int_i_stream, string_t>           int_parser_t;
typedef mt_parse_strat_t<byte_parser_t, int_o_stream> mt_byte_parse_strategy;
typedef mt_parse_strat_t<int_parser_t, int_o_stream>  mt_int_parse_strategy;
typedef st_parse_strat_t<byte_parser_t, int_o_stream> st_byte_parse_strategy;
typedef st_parse_strat_t<int_parser_t, int_o_stream>  st_int_parse_strategy;
#endif //GRLBWT_PARSING_STRATEGIES_H
