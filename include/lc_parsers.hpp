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

    const bv_t& phrase_desc;

    explicit lms_parsing(const bv_t& pr_desc): phrase_desc(pr_desc){};

    inline bool is_suffix(sym_type symbol) const{
        return phrase_desc[symbol];
    }

    /***
    * Find the greatest position j such that ifs[j]<=ifs[idx] and ifs[j] is S*
    * @param idx : the scan starts from this position
    * @param ifs : input file stream
    * @return
    */
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
    }

    void operator()(stream_t& ifs,
                    size_t start, size_t end,
                    std::function<void(string_t&, bool)> task) const {


        bool s_type, prev_s_type = S_TYPE;
        sym_type curr_sym, prev_sym;

        string_t curr_lms(2, sdsl::bits::hi(phrase_desc.size())+1);
        prev_sym = ifs.read(end);
        curr_lms.push_back(prev_sym);

        for (size_t i = end; i-- > start;) {

            curr_sym = ifs.read(i);

            //                                     L_TYPE   S_TYPE*
            //                                        ---- ----
            //this is a junction between two strings = ...$ $...
            if (is_suffix(curr_sym)) {
                bool full_str = curr_lms.size()==1 && is_suffix(curr_lms[0]);
                if (!curr_lms.empty()) {
                    //TODO testing
                    /*if(!first){
                        std::cout<<"first -> ";
                        for(size_t j=curr_lms.size();j-->0;){
                            std::cout<<curr_lms[j]<<" ";
                        }
                        std::cout<<" "<<i<<" "<<end<<std::endl;
                        first = true;
                    }*/
                    //
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

                    if (prev_s_type == S_TYPE) {//Left-most S suffix
                        if (curr_lms.size()>1) {
                            //TODO testing
                            /*if(!first){
                                std::cout<<"ffirst -> ";
                                for(size_t j=curr_lms.size();j-->0;){
                                    std::cout<<curr_lms[j]<<" ";
                                }
                                std::cout<<" "<<std::endl;
                                first = true;
                            }*/
                            //
                            //std::cout<<i<<std::endl;
                            task(curr_lms, false);
                        }
                        curr_lms.clear();
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
        if(!curr_lms.empty()){
            //TODO testing
            /*for(size_t j=curr_lms.size();j-->0;){
                std::cout<<curr_lms[j]<<" ";
            }
            std::cout<<" "<<std::endl;*/
            //
            task(curr_lms, full_str);
        }
    }

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

            //std::cout<<"range: "<<start<<" "<<end<<std::endl;
            if(start<end){
                thread_ranges.emplace_back(start, end);
            }
        }
        is.close();
        return thread_ranges;
    }
};

#endif //LPG_COMPRESSOR_LC_PARSERS_HPP
