//
// Created by Diego Diaz on 11/25/20.
//

#ifndef LPG_COMPRESSOR_PARALLEL_STRING_SORT_HPP
#define LPG_COMPRESSOR_PARALLEL_STRING_SORT_HPP
#include <pthread.h>

struct sorting_data{
    sdsl::bit_vector                 b_limits{};
    sdsl::bit_vector::select_1_type  b_limits_ss{};
    size_t                           n_buckets{};
    size_t                           qs_phrases{};
    size_t                           tot_phrases{};
    std::string                      list_file{};
};

template<typename comparator>
struct sorting_thread_data{
    const size_t                          start;
    const size_t                          end;

    const sdsl::bit_vector::select_1_type &bucket_limits_ss;
    sdsl::int_vector<>                    thread_labels{};
    const comparator&                     comp;

    sorting_thread_data(size_t start_, size_t end_,
                        const sdsl::bit_vector::select_1_type& bl_ss_,
                        sdsl::int_vector<>&& labels_,
                        const comparator& comp_): start(start_),
                                                  end(end_),
                                                  bucket_limits_ss(bl_ss_),
                                                  thread_labels(labels_),
                                                  comp(comp_){};

    void operator()(){
        size_t b_start, b_end, len, offset = bucket_limits_ss(start);
        size_t n_buckets=0, biggest=0, n_phrases=0;

        for(size_t j=start; j <= end; j++){
            b_start = bucket_limits_ss(j);
            b_end = bucket_limits_ss(j+1)-1;
            len = b_end-b_start+1;
            auto it = thread_labels.begin() + (b_start-offset);
            if(len>1){
                std::sort(it, it+len,
                          [&](size_t l, size_t r){
                              return comp(l,r);
                          });
                n_buckets++;
                if(len>biggest) biggest = len;
                n_phrases+=len;
            }
        }
    }
};

template<typename comparator>
void * sort_p(void *data){
    auto t_data = reinterpret_cast<sorting_thread_data<comparator>*>(data);
    (*t_data)();
    pthread_exit(nullptr);
}

template<typename accessor>
void light_sort(std::string& file, accessor& acc, size_t alph_size, sorting_data& sdata, sdsl::cache_config& config,
                std::string& bck_c_file){

    sdsl::bit_vector to_sort(alph_size, false);
    size_t max_freq;
    sdsl::int_vector_buffer<> list(file, std::ios::in);
    {
        sdsl::int_vector<> bucket_counts;
        if(bck_c_file.empty()){
            bucket_counts = sdsl::int_vector<>(alph_size + 1, 0, sdsl::bits::hi(list.size()) + 1);
            for (auto const &elm : list) {
                bucket_counts[acc(elm, 0)]++;
            }
        }else{
            sdsl::load_from_file(bucket_counts, bck_c_file);
        }

        max_freq = 0;
        size_t i = 0;
        sdata.tot_phrases=0;
        for (auto const &freq : bucket_counts) {
            if (freq > 0) {
                if (freq > max_freq) max_freq = freq;
                to_sort[i] = true;
                sdata.tot_phrases +=freq;
                sdata.n_buckets++;
                if (freq > 1) sdata.qs_phrases += freq;
            }
            i++;
        }
        sdata.list_file = sdsl::cache_file_name("presorted_list", config);
        sdata.b_limits = sdsl::bit_vector(sdata.tot_phrases + 1, false);

        size_t accu_freq = 0;
        for (auto const &freq : bucket_counts) {
            if (freq > 0) {
                sdata.b_limits[accu_freq] = true;
                accu_freq += freq;
            }
        }
        sdata.b_limits[accu_freq] = true;
    }

    sdsl::util::init_support(sdata.b_limits_ss, &sdata.b_limits);
    sdsl::bit_vector::rank_1_type to_sort_rs;
    sdsl::util::init_support(to_sort_rs, &to_sort);
    sdsl::int_vector<> bucket_pos(sdata.n_buckets+1, 0, sdsl::bits::hi(max_freq) + 1);
    sdsl::int_vector<> sorted_list(sdata.tot_phrases, 0, list.width());

    size_t rank, bucket;
    for(auto const& elm: list){
        rank = acc(elm,0);
        bucket = to_sort_rs(rank)+1;
        sorted_list[sdata.b_limits_ss(bucket)+bucket_pos[bucket]]=elm;
        bucket_pos[bucket]++;
    }
    list.close();
    sdsl::store_to_file(sorted_list, sdata.list_file);
}

template<typename comparator,
         typename accessor>
void parallel_str_sort(std::string& file, comparator comp,
                       accessor acce, size_t alph_size,
                       size_t n_threads, sdsl::cache_config& config,
                       std::string bck_counts= ""){

    sorting_data sdata;
    light_sort(file, acce, alph_size, sdata, config, bck_counts);

    std::vector<sorting_thread_data<comparator>> threads_data;
    if(sdata.qs_phrases>0){

        size_t len, acc=0, f_bucket=1, s, e;
        std::vector<std::pair<size_t, size_t>> t_ranges;
        size_t phrases_per_thread = INT_CEIL(sdata.qs_phrases, n_threads);
        size_t sum=0;

        for(size_t j=1;j<=sdata.n_buckets;j++){
            s = sdata.b_limits_ss(j);
            e = sdata.b_limits_ss(j+1)-1;
            len = e-s+1;
            if(len>1){
                acc+=len;
                sum+=len;
                if(acc>=phrases_per_thread){
                    t_ranges.emplace_back(f_bucket, j);
                    acc=0;
                    f_bucket=j+1;
                }
            }
        }

        if(f_bucket<=sdata.n_buckets){
            t_ranges.emplace_back(f_bucket, sdata.n_buckets);
        }

        sdsl::int_vector_buffer<> symbols_buff(sdata.list_file, std::ios::in);
        for(auto const& pair : t_ranges){
            s = sdata.b_limits_ss(pair.first);
            e = sdata.b_limits_ss(pair.second+1)-1;
            sdsl::int_vector<> thread_symbols(e - s + 1, 0, symbols_buff.width());
            for(size_t i=s,k=0;i<=e;i++,k++){
                //assert(symbols_buff[i] != 0);
                thread_symbols[k] = symbols_buff[i];
            }
            threads_data.emplace_back(pair.first, pair.second,
                                      sdata.b_limits_ss, std::move(thread_symbols), comp);
        }
        symbols_buff.close(true);

        std::vector<pthread_t> threads(threads_data.size());
        for (size_t j = 0; j < threads_data.size(); j++) {
            int ret = pthread_create(&threads[j],
                                     nullptr,
                                     &sort_p<comparator>,
                                     (void *) &threads_data[j]);
            if (ret != 0) {
                printf("Error: pthread_create() failed\n");
                exit(EXIT_FAILURE);
            }
        }
        for (size_t j = 0; j < threads_data.size(); j++) {
            pthread_join(threads[j], nullptr);
        }
    }else{
        sdsl::int_vector<> thread_label;
        sdsl::load_from_file(thread_label, sdata.list_file);
        threads_data.emplace_back(1, sdata.n_buckets, sdata.b_limits_ss,  std::move(thread_label), comp);
        if(remove(sdata.list_file.c_str())){
            std::cout<<"Error trying to remove file "<<sdata.list_file<<std::endl;
        }
    }

    sdsl::int_vector_buffer<> symbols_buff(file, std::ios::in);
    size_t i = 0;
    for(auto const& data : threads_data){
        for(auto && label : data.thread_labels){
            symbols_buff[i] = label;
            i++;
        }
    }
    symbols_buff.close();
}

#endif //LPG_COMPRESSOR_PARALLEL_STRING_SORT_HPP
