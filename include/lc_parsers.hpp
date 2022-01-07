//
// Created by Diaz, Diego on 24.11.2021.
//

#ifndef LPG_COMPRESSOR_LC_PARSERS_HPP
#define LPG_COMPRESSOR_LC_PARSERS_HPP

#define L_TYPE false
#define S_TYPE true

template<class stream_t,
         class string_t>
struct lms_parsing{

    typedef stream_t                       stream_type;
    typedef typename stream_type::sym_type sym_type;

    const sdsl::int_vector<2>& phrase_desc;

    explicit lms_parsing(const sdsl::int_vector<2>& pr_desc): phrase_desc(pr_desc){};

    inline bool is_suffix(sym_type symbol) const{
        return phrase_desc[symbol] & 2;
    }

    /***
    * Find the greatest position j such that ifs[j]<=ifs[idx] and ifs[j] is S*
    * @param idx : the scan starts from this position
    * @param ifs : input file stream
    * @return
    */
    long long prev_break(long long idx, stream_t& ifs) {
        bool type, prev_type;
        size_t sym, prev_sym = ifs.read(idx);

        //get the suffix type of the symbol at idx
        if (is_suffix(prev_sym)) {
            return idx;
        } else {
            long long pos = idx;
            sym = ifs.read(++pos);
            while (!is_suffix(sym) && prev_sym == sym) {
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

            if(is_suffix(sym)) return idx;
            if(type==L_TYPE && prev_type==S_TYPE) return idx+1;

            prev_sym = sym;
            prev_type = type;
        }
    }

    void operator()(stream_t& ifs,
                    size_t start, size_t end,
                    std::function<void(string_t&, bool)> task) {


        bool s_type, prev_s_type = S_TYPE;
        sym_type curr_sym, prev_sym;

        string_t curr_lms(2, sdsl::bits::hi(phrase_desc.size())+1);
        prev_sym = ifs.read(end);
        curr_lms.push_back(prev_sym);

        /*if(end>=447492){
            for(size_t i=447488;i<=447500;i++){
                std::cout<<ifs.read(i)<<" ";
            }
            std::cout<<" "<<std::endl;
        }*/

        for (size_t i = end; i-- > start;) {

            curr_sym = ifs.read(i);

            //                                     L_TYPE   S_TYPE*
            //                                        ---- ----
            //this is a junction between two strings = ...$ $...
            if (is_suffix(curr_sym)) {
                bool full_str = curr_lms.size()==1 && is_suffix(curr_lms[0]);
                if (!curr_lms.empty()) {
                    task(curr_lms, full_str);
                }
                curr_lms.clear();
                s_type = S_TYPE;
            } else {
                if (curr_sym < prev_sym) {//S_TYPE type
                    s_type = S_TYPE;
                } else if (curr_sym == prev_sym) {
                    s_type = prev_s_type;
                } else {//L_TYPE type
                    s_type = L_TYPE;

                    if (prev_s_type == S_TYPE) {//Left-most suffix
                        curr_lms.pop_back();
                        if (!curr_lms.empty()) {
                            task(curr_lms, false);
                            curr_lms.clear();
                        }
                        curr_lms.push_back(prev_sym);
                    }
                }
            }
            curr_lms.push_back(curr_sym);
            prev_sym = curr_sym;
            prev_s_type = s_type;
        }
        assert(curr_lms[0]!=1);
        bool full_str = curr_lms.size()==1 &&
                        is_suffix(curr_lms[0]) &&
                        (start == 0 || is_suffix(ifs.read(start - 1)));
        if(!curr_lms.empty()) task(curr_lms, full_str);
    }
};

#endif //LPG_COMPRESSOR_LC_PARSERS_HPP
