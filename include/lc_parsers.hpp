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

    std::vector<size_t>& str_ptr;
    const size_t dummy_sym;
    const size_t dummy_bd;//dummy symbol that indicate a masked suffix of a string

    explicit lms_parsing(size_t d_sym, std::vector<size_t>& str_ptr_): str_ptr(str_ptr_),
                                                                       dummy_sym(d_sym),
                                                                       dummy_bd(d_sym+2){};

    /***
    inline bool is_suffix(sym_type symbol) const{
        return phrase_desc[symbol];
    }
    * Find the greatest position j such that ifs[j]<=ifs[idx] and ifs[j] is S*
    * @param idx : the scan starts from this position
    * @param ifs : input file stream
    * @return
    *
    long long prev_break(long long idx, stream_t& ifs) const {

        if(idx<=0) return 0;
        size_t sym=ifs.read(idx), prev_sym=sym;
        size_t pos = idx;
        while(pos<ifs.size() && !is_suffix(sym) && prev_sym == sym) {
            prev_sym = ifs.read(pos++);
        }
        bool type, prev_type = sym < prev_sym;
        prev_sym = sym;

        while(true) {
            idx--;
            if(idx<0) return idx;

            sym = ifs.read(idx);
            type = sym==prev_sym? prev_type : sym < prev_sym;

            if(is_suffix(sym) || (type==L_TYPE && prev_type==S_TYPE)){
                return idx+1;
            }
            prev_sym = sym;
            prev_type = type;
        }
    }*/

    void operator()(stream_t& ifs,
                    size_t start, size_t end,
                    std::function<void(string_t&)> task) const {

        bool s_type, prev_s_type;
        sym_type curr_sym, prev_sym;
        string_t curr_lms(2, sdsl::bits::hi(dummy_sym)+1);

        for(size_t str=end;str>start;str--){

            size_t start_ps = str_ptr[str-1]+1;
            size_t end_ps = str_ptr[str];

            prev_sym = ifs.read(end);
            prev_s_type = S_TYPE;
            curr_lms.push_back(prev_sym);

            for (size_t i = end_ps; i-- > start_ps;) {

                curr_sym = ifs.read(i);
                if(curr_sym>=dummy_sym){//we reach a masked region x_1 # x_2

                    //if curr_sym==dummy_sym, then the left symbol (i-1) is S_type
                    //if curr_sym==(dummy_sym+1), then the left symbol (i-1) is L_type
                    s_type= (curr_sym==dummy_sym);

                    if(i>start){
                        curr_lms.push_back(dummy_sym);//skip the dummy
                        curr_sym = ifs.read(--i);

                        //the break for the next phrase occurs in the left end of the masked region (x_1).
                        // We can merge the active phrase with the next phrase in one single super phrase
                        if(i>start && s_type==S_TYPE && ifs.read(i-1)>curr_sym){
                            curr_lms.push_back(curr_sym);
                            curr_sym = ifs.read(--i);
                            s_type = L_TYPE;
                        }
                    }else{
                        break;
                    }
                } else {
                    if (curr_sym < prev_sym) {//S_TYPE type
                        s_type = S_TYPE;
                    } else if (curr_sym == prev_sym) {
                        s_type = prev_s_type;
                    } else {//L_TYPE type
                        s_type = L_TYPE;
                        if (prev_s_type == S_TYPE) {//LMS suffix

                            //the break for the next phrase occurs in the right end of the masked region (x_2)
                            // We can merge the active phrase with the next phrase in one single super phrase
                            if(i==start || (i>start && ifs.read(i-1)!=dummy_sym)){
                                if (curr_lms.size()>1) {
                                    task(curr_lms);
                                }
                                curr_lms.clear();
                                curr_lms.push_back(prev_sym);
                            }
                        }
                    }
                }
                curr_lms.push_back(curr_sym);
                prev_sym = curr_sym;
                prev_s_type = s_type;
            }

            if(!curr_lms.empty()){
                task(curr_lms);
            }
            curr_lms.clear();
        }
    }

    /*
    std::vector<std::pair<size_t, size_t>> partition_text(size_t n_chunks,
                                                          std::string& i_file) const {
        std::vector<std::pair<size_t, size_t>> thread_ranges;

        stream_type is(i_file, BUFFER_SIZE);
        size_t n_chars = is.tot_cells;
        assert(n_chars>0);
        size_t sym_per_chunk = INT_CEIL(n_chars, n_chunks);
        size_t start, end;
        size_t eff_threads = INT_CEIL(n_chars, sym_per_chunk);

        for(size_t i=0;i<eff_threads;i++){
            start = (i * sym_per_chunk);
            end = std::min<size_t>(((i + 1) * sym_per_chunk), n_chars-1);

            start = start==0? 0 : size_t(prev_break(start, is));
            if(end!=n_chars-1){
                long long tmp_end = prev_break(end, is);
                end = tmp_end<0?  0 : size_t(tmp_end);
                if(is_suffix(is.read(end-1))) end--;
            }else{
                assert(is_suffix(is.read(end)));
            }

            if(start<=end){
                thread_ranges.emplace_back(start, end);
            }
        }
        is.close();
        return thread_ranges;
    }*/
};

#endif //LPG_COMPRESSOR_LC_PARSERS_HPP
