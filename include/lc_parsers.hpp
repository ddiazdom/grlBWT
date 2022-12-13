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
    const size_t dummy_sym;

    explicit lms_parsing(size_t d_sym): dummy_sym(d_sym){};

    void operator()(stream_t& ifs,
                    size_t f_string, size_t l_string,
                    std::function<void(string_t&)> process_phrase,
                    const std::function<std::pair<long, long>(size_t)>& init_str) const {

        bool s_type, prev_s_type;
        sym_type curr_sym, prev_sym;
        string_t phrase(2, sdsl::bits::hi(dummy_sym)+1);
        size_t end_ps, start_ps;

        for(size_t str=l_string+1;str-->f_string;) {

            auto range = init_str(str);

            if(range.first<=range.second) { //if this is not true, it means the string was fully compressed

                start_ps = range.first;
                end_ps = range.second;
                prev_sym = ifs.read(end_ps);
                prev_s_type = L_TYPE;
                phrase.push_back(prev_sym);

                for (size_t i = end_ps; i-- > start_ps;) {

                    curr_sym = ifs.read(i);

                    if(curr_sym>=dummy_sym){//we reach a masked region x_1 # x_2
                        //if curr_sym==dummy_sym, then the left symbol (i-1) is S_type
                        //if curr_sym==(dummy_sym+1), then the left symbol (i-1) is L_type
                        s_type= (curr_sym==dummy_sym);

                        if(i>start_ps){
                            phrase.push_back(dummy_sym);//skip the dummy
                            curr_sym = ifs.read(--i);

                            //the break for the next phrase occurs in the left end of the masked region (x_1).
                            // We can merge the active phrase with the next phrase in one single super phrase
                            if(i>start_ps && s_type==S_TYPE && ifs.read(i-1)>curr_sym){
                                phrase.push_back(curr_sym);
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
                                if(i==start_ps || (i>start_ps && ifs.read(i-1)!=dummy_sym)){
                                    assert(!phrase.empty());
                                    process_phrase(phrase);
                                    phrase.clear();
                                    phrase.push_back(prev_sym);
                                }
                            }
                        }
                    }

                    phrase.push_back(curr_sym);
                    prev_sym = curr_sym;
                    prev_s_type = s_type;
                }

                assert(!phrase.empty());
                process_phrase(phrase);
                phrase.clear();
            }
        }
    }
};

#endif //LPG_COMPRESSOR_LC_PARSERS_HPP
