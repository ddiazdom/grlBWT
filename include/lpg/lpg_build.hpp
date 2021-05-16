//
// Created by diego on 06-03-20.
//

#ifndef LG_COMPRESSOR_LMS_ALGO_H
#define LG_COMPRESSOR_LMS_ALGO_H

#include <unordered_map>
#include <vector>

#include "cdt/file_streams.hpp"
#include "cdt/int_array.h"
#include "cdt/hash_table.hpp"

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support_v.hpp>

#define L_TYPE false
#define S_TYPE true
#define BUFFER_SIZE 8388608 //8MB of buffer

class lpg_build {

    typedef sdsl::bit_vector                             bv_t;
    typedef sdsl::int_vector_buffer<1>                   bvb_t;
    typedef sdsl::int_vector_buffer<>                    ivb_t;
    typedef int_array<size_t>                            string_t;
    typedef bit_hash_table<bool,1>                       string_map_t;
    typedef bit_hash_table<size_t,44>                    phrase_map_t;
    typedef typename string_map_t::buff_t                buff_t;

public:
    typedef std::vector<std::pair<uint8_t, size_t>>      alpha_t;

    struct plain_grammar_t{
        uint8_t                            sigma{}; // terminal's alphabet
        size_t                             r{}; //r: number of rules
        size_t                             c{}; //c: length of the right-hand of the start symbol
        size_t                             g{}; //g: sum of the rules' right-hand sides
        size_t                             max_tsym{}; //maximal terminal symbol
        std::unordered_map<size_t,uint8_t> sym_map; //map terminal symbols to their original byte symobols of the text
        std::string                        rules_file; // rules are concatenated in this array
        std::string                        rules_lim_file; // bit vector marking the last symbol of every right-hand
        std::vector<size_t>                rules_per_level; // number of rules generated in every parsing round
        std::string                        is_rl_file; //bit vector that marks which rules are run-length compressed
        std::string                        lvl_breaks_file; //file with the LMS breaks per level
        uint8_t                            lms_rounds{}; // rounds of LMS parsing

        plain_grammar_t() = default;
        plain_grammar_t(std::string& r_file_,
                        std::string& r_lim_file_,
                        std::string& rl_rules_file_,
                        std::string& lvl_breaks_file_): rules_file(r_file_),
                                                        rules_lim_file(r_lim_file_),
                                                        is_rl_file(rl_rules_file_),
                                                        lvl_breaks_file(lvl_breaks_file_){}

        void save_to_file(std::string& output_file);
        void load_from_file(std::string &g_file);
    };

    // the phrases are stored in a bit compressed hash table:
    // this wrapper reinterprets the bits back as phrases
    struct key_wrapper{
        size_t width;
        size_t d_bits;//bits used to describe the string
        const bitstream<buff_t>& stream;

        //offset is the bit where the key description starts
        inline size_t read(size_t offset, size_t idx) const {
            return stream.read(offset + d_bits + (idx * width),
                               offset + d_bits + ((idx + 1) * width - 1));
        }

        //offset is the bit where the key description starts
        inline size_t size(size_t offset) const{
            return stream.read(offset, offset + d_bits - 1) / width;
        }

        //compare two phrases are different positions of the bit stream
        inline bool compare(size_t a, size_t b) const{

            size_t a_bits = stream.read(a, a + d_bits - 1);
            size_t b_bits = stream.read(b, b + d_bits - 1);
            size_t min_bits = std::min<size_t>(a_bits, b_bits);

            size_t a_pos = a+d_bits+a_bits-1;
            size_t b_pos = b+d_bits+b_bits-1;
            size_t rm_diff = stream.inv_com_segments(a_pos, b_pos, min_bits);

            if(rm_diff < min_bits){
                a_pos = a+d_bits+(((a_bits - rm_diff-1) / width) * width);
                b_pos = b+d_bits+(((b_bits - rm_diff-1) / width) * width);
                size_t sym_a = stream.read(a_pos, a_pos+width-1);
                size_t sym_b = stream.read(b_pos, b_pos+width-1);
                return sym_a<sym_b;
            }else{
                return a_bits>b_bits;
            }
        }
    };

    /***
     *
     * @param i_file : input text file
     * @param p_gram_file : file where the plain grammar will be stored
     * @param n_threads : number of working threads
     * @param config : temporal files handler
     * @param hbuff_size : buffer size for the hashing step
     * @param sep_symbol : string delimiter in the input text
     */
    static void compute_LPG(std::string &i_file, std::string &p_gram_file, size_t n_threads, sdsl::cache_config &config,
                            size_t hbuff_size, alpha_t &alphabet);

    /***
     * check if the grammar is correct
     * @param g_file : file with the plain grammar
     * @param uncomp_file : original input text
     */
    static void check_plain_grammar(plain_grammar_t& p_gram, std::string& uncomp_file);
private:

    template<class sym_type>
    struct lms_info {

        i_file_stream<sym_type>    ifs;
        o_file_stream<size_t>      ofs;
        const sdsl::int_vector<2>& phrase_desc;

        phrase_map_t&              m_map;
        size_t                     start;
        size_t                     end;
        const uint8_t              sym_width;
        string_map_t               thread_map;

        lms_info(std::string &i_file_, std::string &o_file_, phrase_map_t &m_map_,
                 size_t start_, size_t end_,
                 const size_t &alph,
                 const size_t &hb_size, void *hb_addr,
                 const sdsl::int_vector<2> &phrase_desc_) : ifs(i_file_, BUFFER_SIZE),
                                                            ofs(o_file_, BUFFER_SIZE, std::ios::out),
                                                            phrase_desc(phrase_desc_),
                                                            m_map(m_map_),
                                                            start(start_),
                                                            end(end_),
                                                            sym_width(sdsl::bits::hi(alph)+1),
                                                            thread_map(hb_size, o_file_+"_phrases", 0.8, hb_addr) {

            //TODO for the moment the input string has to have a sep_symbol appended at the end
            //TODO assertion : sep_symbols cannot be consecutive
        };

        inline void hash_phrase(string_t& phrase) {
            phrase.mask_tail();
            auto res = thread_map.insert(phrase.data(), phrase.n_bits(), false);
        };

        inline void store_phrase(string_t& phrase){
            phrase.mask_tail();
            auto res = m_map.find(phrase.data(), phrase.n_bits());
            if(res.second){
                ofs.push_back(res.first.value()>>1UL);
            }else{
                assert(phrase.size()==1 && is_suffix(phrase[0]));
                ofs.push_back(phrase[0]);
            }
        };

        inline bool is_suffix(sym_type symbol) const{
            return phrase_desc[symbol] & 2;
        }
    };

    /***
     * Find the previous S* suffix before S[idx]
     * @tparam sym_t : type of symbol in S
     * @param idx : the scan starts from this position
     * @param ifs : input file stream for S
     * @param phrase_desc : input array with the class of each symbol S
     * @return
     */
    template<class sym_t>
    static long long prev_lms_sym(long long idx, i_file_stream<sym_t>& ifs, sdsl::int_vector<2>& phrase_desc) {

        bool type, prev_type;
        size_t sym, prev_sym = ifs.read(idx);

        //get the suffix type of the symbol at idx
        if ((phrase_desc[prev_sym] & 2UL)) {
            return idx;
        } else {
            long long pos = idx;
            sym = ifs.read(++pos);
            while (!(phrase_desc[sym] & 2UL) && prev_sym == sym) {
                sym = ifs.read(++pos);
            }
            prev_type = prev_sym < sym;
        }

        //reach the next LMS-symbol
        while(true) {
            idx--;
            if(idx<0) return idx;

            sym = ifs.read(idx);
            type = sym==prev_sym? prev_type : sym < prev_sym;

            if(phrase_desc[sym] & 2U) return idx;
            if(type==L_TYPE && prev_type==S_TYPE) return idx+1;

            prev_sym = sym;
            prev_type = type;
        }
    }

    template<class sym_type>
    static std::vector<std::pair<size_t, size_t>>
    compute_thread_ranges(size_t n_threads, std::string& i_file, sdsl::int_vector<2>& phrase_desct);

    template<class sym_type>
    static size_t
    compute_LPG_int(std::string &i_file, std::string &o_file, size_t n_threads, size_t hbuff_size,
                    plain_grammar_t &p_gram, ivb_t &rules, bvb_t &rules_lim,
                    sdsl::int_vector<2> &phrase_desc, sdsl::cache_config &config);
    static void
    assign_ids(phrase_map_t &mp_map, size_t max_sym, key_wrapper &key_w, ivb_t &r, bvb_t &r_lim,
               size_t n_threads, sdsl::cache_config &config);

    static void join_parse_chunks(const std::string &output_file,
                                  std::vector<std::string> &chunk_files);
    static void join_thread_phrases(phrase_map_t& mp_map, std::vector<std::string> &chunk_files);

    template<class sym_t>
    static void * hash_phrases(void * data);
    template<class sym_t>
    static void * record_phrases(void *data);

    //mark the nonterminals that can be removed from the grammar
    static bv_t mark_nonterminals(plain_grammar_t& p_gram);
    static void simplify_grammar(lpg_build::plain_grammar_t &p_gram, bv_t &rem_nts, bv_t::rank_1_type &rem_nts_rs);
    static void decomp(size_t nt, sdsl::int_vector<> &rules, bv_t::select_1_type &rlim_ss, bv_t &rem_nt,
                       bv_t::rank_1_type &rem_nt_rs,
                       ivb_t &dec);
    //these functions are to build the index
    static void create_lvl_breaks(plain_grammar_t &p_gram, bv_t &rem_nts, bv_t::rank_1_type &rem_nts_rs);
    static std::vector<uint8_t>
    rec_dc(size_t nt, uint8_t lev, sdsl::int_vector<> &rules, bv_t &rem_nts, bv_t::select_1_type &r_lim_ss);
    static void rec_dc_int(size_t nt, uint8_t lev, size_t& pos, bool rm, sdsl::int_vector<>& rules,
                           bv_t& rem_nts, bv_t::select_1_type &r_lim_ss,
                           std::vector<uint8_t> &breaks);
    static void colex_nt_sort(plain_grammar_t &p_gram);
    static void run_length_compress(plain_grammar_t& p_gram, sdsl::cache_config& config);
};
#endif //LG_COMPRESSOR_LMS_ALGO_H
