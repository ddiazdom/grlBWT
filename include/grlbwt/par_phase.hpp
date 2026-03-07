//
// Created by Diaz, Diego on 26.1.2023.
//

#ifndef GRLBWT_PAR_PHASE_H
#define GRLBWT_PAR_PHASE_H

#ifdef USE_MALLOC_COUNT
#include "malloc_count.h"
#endif

#include <cdt/utils.h>
#include "LMS_induction.hpp"
#include "grlbwt_common.h"
#include "par_strategies.h"

template<typename parse_data_t, typename parser_t>
struct hash_functor {

    void operator()(parse_data_t& data) {

        auto init_str = [&](size_t str) -> std::pair<long, long>{
            //size_t start = data.str_ptr[str];
            //size_t end = data.str_ptr[str+1]-1;
            auto range = data.str2range(str);
            data.active_strings += (range.first<=range.second);
            return range;//{start, end};
        };

        parser_t()(data.ifs, data.start_str, data.end_str, data.max_symbol,
                //hash_phrase,
                [&](string_t& phrase) -> void {
                    phrase.mask_tail();
                    data.inner_map.increment_value(phrase.data(), phrase.n_bits(), 1);
                },
                init_str);
        //pthread_exit(nullptr);
    };
};

template<class parse_data_t,
         class parser_t,
         class o_stream_type>
struct parse_functor{

    size_t operator()(parse_data_t& data) {

        o_stream_type ofs(data.o_file, BUFFER_SIZE, std::ios::out);

        auto phrase2symbol = [&](string_t& phrase) -> void {
            phrase.mask_tail();
            size_t sym = 0;
            auto res = data.map.key2value(phrase.data(), phrase.n_bits(), sym);
            assert(res);
            ofs.push_back(sym);
        };

        auto init_str = [&](size_t str) -> std::pair<long, long>{

            //assert(str>=data.start && str<=data.end);
            //size_t start = data.str_ptr[str];
            //size_t end = data.str_ptr[str+1]-1;
            //assert(start<=end);
            auto range = data.str2range(str);
            if((str+1)<=data.end_str){
                data.str_ptr[str+1] = ofs.size()-1;
            }
            return range; //{start, end};
        };

        parser_t()(data.ifs, data.start_str, data.end_str, data.max_symbol, phrase2symbol, init_str);

        data.str_ptr[data.start_str] = ofs.size()-1;

        data.ifs.close();
        ofs.close();

        //pthread_exit(nullptr);
        return ofs.size();
    };
};

struct dictionary {

    typedef size_t size_type;
    size_t alphabet{};      //alphabet of the dictionary
    size_t prev_alphabet{};  //size of the previous alphabet
    size_t n_phrases{};     //number of LMS phrases in the dictionary
    size_t t_size{};        //size of the text from which the dictionary was generated
    size_t max_sym_freq{};  //maximum symbol frequency in the parse from which the dictionary was created
    size_t end_str_dummy{}; //symbol to mark the end of a string
    size_t bwt_dummy{};     //dummy symbol indicating that the induction is from the BWT i+1
    size_t hocc_dummy{};    //dummy symbol indicating that the induction is from the hocc array
    size_t metasym_dummy{}; //dummy metasymbol for the compressed suffixes
    vector_t dict;          //list of phrases in the dictionary
    vector_t freqs;         //frequency of every dictionary phrase in the original text
    bv_t d_lim;
    bv_t phrases_has_hocc; //mark the phrases with hidden occurrences
    bv_t *desc_bv = nullptr;

    dictionary() = default;

    dictionary(phrase_map_t &mp_map, size_t dict_syms,
               size_t max_freq, bv_t &is_suffix_bv, size_t _t_size, size_t _p_alph_size,
               size_t _max_sym_freq) : alphabet(is_suffix_bv.size()),
                                       prev_alphabet(_p_alph_size),
                                       n_phrases(mp_map.size()),
                                       t_size(_t_size),
                                       max_sym_freq(_max_sym_freq),
                                       end_str_dummy(alphabet),
                                       bwt_dummy(alphabet + 1),
                                       hocc_dummy(alphabet + 2),
                                       dict(dict_syms, 0, sym_width( alphabet)),//the +3 corresponds to the extra symbols I use for the BWT induction
                                       freqs(n_phrases, 0, sym_width(max_freq)),
                                       d_lim(dict_syms, false),
                                       desc_bv(&is_suffix_bv) {

        key_wrapper key_w{sym_width(alphabet), mp_map.description_bits(), mp_map.get_data()};
        size_t j = 0, k = 0, freq;

        //TODO testing
        //std::vector<std::vector<size_t>> plain_dict;
        //

        for (auto const &ptr: mp_map) {
            //TODO
            //std::vector<size_t> phrase;
            //

            for (size_t i = key_w.size(ptr); i-- > 0;) {
                dict[j] = key_w.read(ptr, i);
                //phrase.push_back(dict[j]);
                d_lim[j++] = false;
            }
            //plain_dict.emplace_back(phrase);
            d_lim[j - 1] = true;

            freq = 0;
            mp_map.get_value_from(ptr, freq);
            assert(freq <= max_freq);
            freqs[k++] = freq;
        }

        //TODO
        //std::cout<<"storing the dictionary in plain format, delete this"<<std::endl;
        //std::sort(plain_dict.begin(), plain_dict.end(), [](auto a, auto b){
        //    for(size_t i=0;i<std::min(a.size(), b.size()); i++){
        //        if(a[i]!=b[i]) return a[i]<b[i];
        //    }
        //    return a.size()<b.size();
        //});

        //std::ifstream ifs("serialized_dict_"+std::to_string(t_size), std::ios::in | std::ios::binary);
        //std::vector<std::vector<size_t>> correct_dict;
        //for(size_t i=0;i<plain_dict.size();i++){
        //    std::vector<size_t> correct_phrase;
        //    load_plain_vector(ifs, correct_phrase);
        //    if(correct_phrase.size()!=plain_dict[i].size()){
        //        std::cout<<i<<" -> length corr "<<correct_phrase.size()<<" length malo "<<plain_dict[i].size()<<std::endl;
        //        std::cout<<"bueno: "<<std::endl;
        //        for(size_t u=0;u<correct_phrase.size();u++){
        //            std::cout<<correct_phrase[u]<<" ";
        //        }
        //        std::cout<<"malo: "<<std::endl;
        //        for(size_t u=0;u<plain_dict[i].size();u++){
        //            std::cout<<plain_dict[i][u]<<" ";
        //        }
        //    }
        //    assert(correct_phrase.size()==plain_dict[i].size());
        //    for(size_t u=0;u<correct_phrase.size();u++){
        //        assert(correct_phrase[u]==plain_dict[i][u]);
        //    }
        //}
        //ifs.close();
        //for(auto const &vector : plain_dict ){
        //    serialize_plain_vector(ofs, vector);
        //}
        //
        assert(j == dict_syms);
    }

    [[nodiscard]] inline bool is_suffix(size_t sym) const {
        return (*desc_bv)[sym];
    };

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v = nullptr, std::string name = "") const {
        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_type written_bytes = sdsl::write_member(alphabet, out, child, "alphabet");
        written_bytes += sdsl::write_member(prev_alphabet, out, child, "p_alpha_size");
        written_bytes += sdsl::write_member(n_phrases, out, child, "n_phrases");
        written_bytes += sdsl::write_member(t_size, out, child, "t_size");
        written_bytes += sdsl::write_member(max_sym_freq, out, child, "max_sym_freq");
        written_bytes += sdsl::write_member(end_str_dummy, out, child, "end_str_dummy");
        written_bytes += sdsl::write_member(bwt_dummy, out, child, "bwt_dummy");
        written_bytes += sdsl::write_member(hocc_dummy, out, child, "hocc_dummy");
        written_bytes += sdsl::write_member(metasym_dummy, out, child, "metasym_dummy");
        dict.serialize(out);
        phrases_has_hocc.serialize(out, child);
        return written_bytes;
    }

    void load(std::istream &in) {
        sdsl::read_member(alphabet, in);
        sdsl::read_member(prev_alphabet, in);
        sdsl::read_member(n_phrases, in);
        sdsl::read_member(t_size, in);
        sdsl::read_member(max_sym_freq, in);
        sdsl::read_member(end_str_dummy, in);
        sdsl::read_member(bwt_dummy, in);
        sdsl::read_member(hocc_dummy, in);
        sdsl::read_member(metasym_dummy, in);
        dict.load(in);
        phrases_has_hocc.load(in);
    }
};

template<class sa_type>
void produce_grammar(dictionary& dict, sa_type& s_sa,  phrase_map_t& new_phrases_ht, bv_t& phr_marks, parsing_info& p_info, tmp_workspace& ws) {

    string_t phrase(2, sym_width(dict.alphabet));

    dict.alphabet+=3;
    dict.metasym_dummy= dict.alphabet+s_sa.size()+1;

    size_t alphabet = dict.alphabet;
    size_t meta_sym_dummy = dict.metasym_dummy;

    size_t pos, new_size=0, l_sym, r_sym;
    vector_t new_dict((s_sa.size()+1)*2, 0, sym_width(meta_sym_dummy));

    //TODO testing
    //std::cout<<"\ndist stats"<<std::endl;
    //new_phrases_ht.ht_stats(10);
    //

    for(size_t u=0;u<s_sa.size();u++) {

        pos = s_sa[u];
        assert(pos<dict.dict.size() && (!dict.d_lim[pos] || dict.is_suffix(dict.dict[pos])));

        if(dict.d_lim[pos]){
            assert(dict.is_suffix(dict.dict[pos]));
            new_dict.write(new_size++, meta_sym_dummy);
            new_dict.write(new_size++, dict.dict[pos]);
        }else{
            pos++;
            while(!phr_marks[pos] && !dict.d_lim[pos]) pos++;

            l_sym = dict.dict.read(pos-1);
            assert(l_sym<alphabet);

            if(phr_marks[pos]){

                //TODO
                /*if(p_info.p_round==11){
                    std::cout<<"\n"<<pos<<std::endl;
                }*/
                //

                phrase.clear();
                do{
                    phrase.push_back(dict.dict[pos]);
                }while(!dict.d_lim[pos++]);
                phrase.mask_tail();

                //TODO
                /*if(p_info.p_round==11){
                    for(size_t i=0;i<phrase.size();i++){
                        std::cout<<phrase[i]<<" ";
                    }
                    std::cout<<""<<std::endl;
                }*/
                //

                auto res = new_phrases_ht.find(phrase.data(), phrase.n_bits());
                assert(res.second);

                r_sym = 0;
                new_phrases_ht.get_value_from(res.first, r_sym);
                r_sym+=dict.alphabet;

                new_dict.write(new_size++, l_sym);
                new_dict.write(new_size++, r_sym);
            }else{
                r_sym = dict.dict.read(pos);
                new_dict.write(new_size++, meta_sym_dummy);
                new_dict.write(new_size++, dict.is_suffix(r_sym) ? r_sym : l_sym);
            }
        }
    }

    dict.dict.swap(new_dict);
    dict.n_phrases = s_sa.size();

    sdsl::util::clear(dict.d_lim);
    std::string dict_file = ws.get_file("dict_lev_"+std::to_string(p_info.p_round));
    sdsl::store_to_file(dict, dict_file);
}

template<class sa_type>
size_t produce_pre_bwt(dictionary &dict, sa_type &sa, phrase_map_t &new_phrases_ht,
                       bv_t &phr_marks, parsing_info &p_info, tmp_workspace &ws) {

    string_t phrase(2, sym_width(dict.alphabet));

    size_t n_breaks=0, n_phrases=0, acc_freq=0, null_sym=std::numeric_limits<size_t>::max();
    size_t pos, pl_sym=null_sym, bg_range=0, rank=0, l_sym;
    bool is_full_phrase, is_valid_suffix;

    bv_rs_t d_lim_rs(&dict.d_lim);

    std::string pre_bwt_file = ws.get_file("pre_bwt_lev_"+std::to_string(p_info.p_round));
    size_t sb = INT_CEIL(sym_width(dict.hocc_dummy), 8);
    size_t fb = INT_CEIL(sym_width(dict.t_size),8);
    bwt_buff_writer pre_bwt(pre_bwt_file, std::ios::out, sb, fb);

    vector_t ranks(dict.n_phrases, 0, sym_width(dict.dict.size())+1);
    dict.phrases_has_hocc.resize(dict.dict.size());
    sdsl::util::set_to_value(dict.phrases_has_hocc, false);

    size_t u = 0 ;
    while(u<sa.size()) {

        pos = (sa[u]>>1UL) - 1;
        is_valid_suffix = !dict.d_lim[pos] || dict.is_suffix(dict.dict[pos]);

        if(!is_valid_suffix) {
            u++;
            while(sa[u] & 1UL) u++;
            bg_range = u;
            continue;
        }

        is_full_phrase = (pos == 0 || dict.d_lim[pos - 1]);
        size_t freq = dict.freqs[d_lim_rs(pos)];
        if(is_full_phrase){
            ranks[d_lim_rs(pos)] = rank<<1UL | (freq>1);
            l_sym = dict.bwt_dummy;
        }else{
            l_sym =  dict.dict[pos-1];
        }

        n_breaks +=pl_sym!=l_sym;
        n_phrases += is_full_phrase;
        acc_freq += freq;

        if((u+1)==sa.size() || !(sa[u+1] & 1UL)) {

            if(n_breaks>1 || n_phrases==1) {

                l_sym = dict.bwt_dummy;
                if((u-bg_range+1) > 1) {

                    size_t samp_pos = pos;
                    phrase.clear();
                    do{
                        phrase.push_back(dict.dict[samp_pos]);
                    } while(!dict.d_lim[samp_pos++]);

                    phrase.mask_tail();
                    auto res = new_phrases_ht.insert(phrase.data(), phrase.n_bits(), rank);
                    assert(res.second);
                    dict.phrases_has_hocc[rank] = true;

                    for (size_t j = bg_range; j <= u; j++) {
                        phr_marks[(sa[j] >> 1UL) - 1] = true;
                    }
                    l_sym = dict.hocc_dummy;
                }
                sa[rank] = pos;
                rank++;
            }

            if (pre_bwt.size() > 1 && pre_bwt.last_sym() == l_sym) {
                pre_bwt.inc_freq_last(acc_freq);
            } else {
                pre_bwt.push_back(l_sym, acc_freq);
            }

            bg_range = u+1;
            l_sym = null_sym;
            n_breaks = 0;
            n_phrases = 0;
            acc_freq = 0;
        }
        pl_sym = l_sym;
        u++;
    }

    pre_bwt.close();
    sa.resize(rank);
    dict.phrases_has_hocc.resize(rank);

    sdsl::util::clear(d_lim_rs);
    dict.freqs.erase();
    store_to_file(ws.get_file("phr_ranks"), ranks);
    ranks.erase();
    //new_phrases_ht.shrink_databuff();

#ifdef __linux__
    malloc_trim(0);
#endif
    return rank;
}

template<class vector_type, class sa_type>
size_t process_dictionary_int(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

    sa_type sa;
    if constexpr (std::is_same<sa_type, vector_t>::value) {
        size_t width = sym_width(dict.dict.size())+2;
        sa.set_width(width);
        sa.resize(dict.dict.size());
        sa.initialize(0, sa.size());
    }else{
        sa.resize(dict.dict.size());
        memset(sa.data(), 0, sa.size()*sizeof(typename sa_type::value_type));
    }

    std::cout <<"    Sorting the dictionary and constructing the preliminary BWT" << std::flush;
    auto start = std::chrono::steady_clock::now();
    suffix_induction<dictionary, sa_type, vector_type>(dict, sa);

    phrase_map_t new_phrases_ht;
    bv_t phr_marks(dict.dict.size(), false);
    size_t n_phrases = produce_pre_bwt<sa_type>(dict, sa, new_phrases_ht, phr_marks, p_info, ws);
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 2);

    //malloc_count_print_status();
    //malloc_count_reset_peak();

    std::cout << "    Compressing the dictionary" << std::flush;
    start = std::chrono::steady_clock::now();
    produce_grammar<sa_type>(dict, sa, new_phrases_ht, phr_marks, p_info, ws);
    end = std::chrono::steady_clock::now();
    report_time(start, end, 35);

#ifdef __linux__
    malloc_trim(0);
#endif
    return n_phrases;
}

inline size_t process_dictionary(dictionary &dict, parsing_info &p_info, tmp_workspace &ws) {

    uint8_t width = sym_width(dict.dict.size()) + 1;
    size_t n_phrases;

    if (width <= 8) {
        using uint8_vector_t = std::vector<uint8_t, mallocator<uint8_t>>;
        n_phrases = process_dictionary_int<uint8_vector_t, uint8_vector_t>(dict, p_info, ws);
    } else if (width <= 16) {
        using uint16_vector_t = std::vector<uint16_t, mallocator<uint16_t>>;
        n_phrases = process_dictionary_int<uint16_vector_t, uint16_vector_t>(dict, p_info, ws);
    } else if (width <= 32) {
        using uint32_vector_t = std::vector<uint32_t, mallocator<uint32_t>>;
        n_phrases = process_dictionary_int<uint32_vector_t, uint32_vector_t>(dict, p_info, ws);
    } else {
        using uint64_vector_t = std::vector<uint64_t, mallocator<uint64_t>>;
        n_phrases = process_dictionary_int<uint64_vector_t, vector_t>(dict, p_info, ws);
    }
    return n_phrases;
}


template<class parse_strategy_t>
size_t par_round(parse_strategy_t &p_strategy, parsing_info &p_info, bv_t &phrase_desc, tmp_workspace &ws) {

#ifdef __linux__
    malloc_trim(0);
#endif
    std::cout << "    Computing the dictionary of LMS phrases" << std::flush;
    auto start = std::chrono::steady_clock::now();
    auto res = p_strategy.get_phrases();
    auto end = std::chrono::steady_clock::now();
    report_time(start, end, 22);


    store_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);
    std::vector<long>().swap(p_info.str_ptrs);

#ifdef __linux__
    malloc_trim(0);
#endif

    phrase_map_t &map = p_strategy.map;
    size_t psize;//<- for the iter stats
    assert(map.size() > 0);

    size_t dict_sym = res.first;
    size_t max_freq = res.second;

    //save a copy of the hash table into a file
    std::string ht_file = ws.get_file("ht_data");
    map.store_data_to_file(ht_file);

    //temporal unload of the hash table (not the data)
    map.destroy_table();
    size_t tot_phrases;
    {
        //create a dictionary from where the ids will be computed
        std::cout << "    Compacting the dictionary" << std::flush;
        start = std::chrono::steady_clock::now();
        dictionary dict(map, dict_sym, max_freq, phrase_desc,
                        p_strategy.text_size, p_info.prev_alph, p_info.max_sym_freq);
        end = std::chrono::steady_clock::now();
        map.destroy_data();
        report_time(start, end, 36);

        //process the dictionary
        tot_phrases = process_dictionary(dict, p_info, ws);
        p_info.prev_alph = dict.alphabet;
    }

    //reload the hash table
    map.load_data_from_file(ht_file);
    ws.remove_file("ht_data");

    {
        std::cout << "    Assigning metasymbols to the LMS phrases" << std::flush;
        start = std::chrono::steady_clock::now();
        bv_t new_phrase_desc(tot_phrases, false);
        key_wrapper key_w{sym_width(p_info.tot_phrases), map.description_bits(), map.get_data()};

        size_t j = 0;
        vector_t ranks;
        std::string ranks_file = ws.get_file("phr_ranks");
        load_from_file(ranks_file, ranks);
        for (auto const &ptr: map) {
            size_t val = ranks[j++];
            map.insert_value_at(ptr, val);

            //the first bit marks if the phrase is repeated or not.
            // We need to shift it to get the real id
            new_phrase_desc[(val >> 1UL)] = phrase_desc[key_w.read(ptr, 0)];
        }
        ws.remove_file("phr_ranks");
        std::string suffix_file = ws.get_file("suffix_file");
        sdsl::store_to_file(new_phrase_desc, suffix_file);
        end = std::chrono::steady_clock::now();
        report_time(start, end, 21);
    }

    std::cout << "    Creating the parse of the text" << std::flush;
    start = std::chrono::steady_clock::now();
    load_pl_vector(ws.get_file("str_ptr"), p_info.str_ptrs);

    size_t bps = sym_width(tot_phrases)+1;
    if(bps<=8){
        psize = p_strategy.template parse_text<uint8_t>();
    }else if(bps<=16){
        psize = p_strategy.template parse_text<uint16_t>();
    } else if(bps<=32){
        psize = p_strategy.template parse_text<uint32_t>();
    }else{
        psize = p_strategy.template parse_text<uint64_t>();
    }
    assert(psize>=map.size());//the parse can't be smaller than the number of phrases

    end = std::chrono::steady_clock::now();
    report_time(start, end, 31);

    {
        //keep track of the phrases that have to be rephrased
        std::string suffix_file = ws.get_file("suffix_file");
        bv_t new_phrase_desc;
        sdsl::load_from_file(new_phrase_desc, suffix_file);
        phrase_desc.swap(new_phrase_desc);
    }

    p_info.max_sym_freq = max_freq;
    p_info.lms_phrases = map.size();
    p_info.tot_phrases = tot_phrases;
    p_info.p_round++;

    std::cout << "    Stats:" << std::endl;
    std::cout << "      Parsing phrases:                  " << p_info.lms_phrases << std::endl;
    std::cout << "      Number of symbols in the phrases: " << dict_sym << std::endl;
    std::cout << "      Number of unsolved BWT blocks:    " << p_info.tot_phrases << std::endl;
    std::cout << "      Parse size:                       " << psize << std::endl;

    map.destroy_data();
    map.destroy_table();

#ifdef __linux__
    malloc_trim(0);
#endif
    return (p_info.str_ptrs.size()-1) == psize ? 0 : p_info.tot_phrases;
}

template<class par_type>
size_t par_phase_int(std::string& i_file, std::string& o_file, parsing_info& p_info,
                     size_t hbuff_size, size_t n_threads, bv_t& symbol_desc, tmp_workspace& ws){
    size_t alph_size;
    if (n_threads > 1) {
        using par_strat_t = mt_parse_strat_t<par_type, hash_functor, parse_functor>;
        par_strat_t p_strat(i_file, o_file, p_info, hbuff_size, n_threads);
        alph_size = par_round<par_strat_t>(p_strat, p_info, symbol_desc, ws);
    } else {
        using par_strat_t = st_parse_strat_t<par_type, hash_functor, parse_functor>;
        par_strat_t p_strat(i_file, o_file, p_info);
        alph_size = par_round<par_strat_t>(p_strat, p_info, symbol_desc, ws);
    }
    return alph_size;
}

template<class sym_type>
size_t par_phase(std::string &i_file, size_t n_threads, float hbuff_frac, tmp_workspace &ws) {

    std::cout<<"Reading the file"<<std::endl;
    str_collection str_coll = collection_stats<sym_type>(i_file);
    std::cout<<"Stats: "<<std::endl;
    std::cout<<"  Smallest symbol               : "<<str_coll.min_sym<<std::endl;
    std::cout<<"  Greatest symbol               : "<<str_coll.max_sym<<std::endl;
    std::cout<<"  Number of symbols in the file : "<<str_coll.n_syms<<std::endl;
    std::cout<<"  Number of strings             : "<<str_coll.n_strings<<std::endl;

    auto hbuff_size = std::max<size_t>(64 * n_threads, size_t(ceil(float(str_coll.n_syms) * hbuff_frac)));

    std::cout << "Parsing the text:    " << std::endl;
    if(n_threads>1){
        std::cout<<"  Running with up to "<<n_threads<<" working threads "<<std::endl;
        std::cout<<"  Using "<<report_space((off_t)hbuff_size)<<" for the thread hash tables ("<< report_space(off_t(hbuff_size/n_threads))<<" each)"<<std::endl;
    }

    std::string output_file = ws.get_file("tmp_output");
    std::string tmp_i_file = ws.get_file("tmp_input");

    // mark which symbols represent string boundaries
    bv_t symbol_desc(str_coll.max_sym + 1, false);
    symbol_desc[str_coll.min_sym] = true;

    parsing_info p_info;
    p_info.max_sym_freq = str_coll.max_sym_freq;
    p_info.tot_phrases = str_coll.max_sym + 1;
    p_info.str_ptrs.swap(str_coll.str_ptrs);
    p_info.str_ptrs.push_back((long) str_coll.n_syms);
    p_info.str_ptrs.shrink_to_fit();
    p_info.longest_str = str_coll.longest_string;
    p_info.active_strings = str_coll.n_strings;

    size_t iter = 1;
    size_t n_syms;
    using f_parser_t = lms_parsing<i_file_stream<sym_type>, string_t, true>;

    std::cout << "  Parsing round " << iter++ << std::endl;
    auto start = std::chrono::steady_clock::now();
    n_syms = par_phase_int<f_parser_t>(i_file, tmp_i_file, p_info, hbuff_size, n_threads, symbol_desc, ws);
    auto end = std::chrono::steady_clock::now();

    report_time(start, end, 4);
#ifdef USE_MALLOC_COUNT
    malloc_count_print_status();
    malloc_count_reset_peak();
#endif

    while (n_syms > 0) {
        start = std::chrono::steady_clock::now();
        std::cout << "  Parsing round " << iter++ << std::endl;
        size_t bps = sym_width(n_syms)+1;
        if(bps<=8){
            n_syms = par_phase_int<uint8t_parser_t>(tmp_i_file, output_file, p_info,
                                                    hbuff_size, n_threads, symbol_desc, ws);
        }else if(bps<=16){
            n_syms = par_phase_int<uint16t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, symbol_desc, ws);
        } else if(bps<=32){
            n_syms = par_phase_int<uint32t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, symbol_desc, ws);
        } else{
            n_syms = par_phase_int<uint64t_parser_t>(tmp_i_file, output_file, p_info,
                                                     hbuff_size, n_threads, symbol_desc, ws);
        }
        end = std::chrono::steady_clock::now();
        report_time(start, end, 4);

        remove(tmp_i_file.c_str());
        rename(output_file.c_str(), tmp_i_file.c_str());

#ifdef USE_MALLOC_COUNT
        malloc_count_print_status();
        malloc_count_reset_peak();
#endif
    }

    sdsl::util::clear(symbol_desc);
    ws.remove_file("suffix_file");
    return iter - 2;
}
#endif //GRLBWT_PAR_PHASE_H
