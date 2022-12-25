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

    void operator()(stream_t& ifs,
                    size_t f_string, size_t l_string, size_t max_symbol,
                    const std::function<void(string_t&)>& process_phrase,
                    const std::function<std::pair<long, long>(size_t)>& init_str) const {

        bool s_type, prev_s_type, is_rep, prev_is_rep;
        sym_type curr_sym, prev_sym;
        string_t phrase(2, sdsl::bits::hi(max_symbol)+1);
        size_t end_ps, start_ps;

        //TODO this is for testing
        //size_t n_rep_sym, uniq_phrases=0;//, del_sym=0;
        //std::vector<bool> test_seq;
        //size_t run;
        //

        for(size_t str=l_string+1;str-->f_string;) {

            auto range = init_str(str);

            if(range.first<=range.second) { //if this is not true, it means the string was fully compressed

                start_ps = range.first;
                end_ps = range.second;

                prev_sym = ifs.read(end_ps);
                if constexpr (std::is_same<sym_type, uint8_t>::value){
                    prev_is_rep = true;
                }else{
                    prev_is_rep = prev_sym & 1UL;
                    prev_sym >>=1UL;
                }
                //n_rep_sym=prev_is_rep;
                //run=1;

                prev_s_type = L_TYPE;
                phrase.push_back(prev_sym);
                //test_seq.push_back(prev_is_rep);

                for (size_t i = end_ps; i-- > start_ps;) {

                    curr_sym = ifs.read(i);

                    if constexpr (std::is_same<sym_type, uint8_t>::value){
                        is_rep = true;
                    }else{
                        is_rep = curr_sym & 1UL;
                        curr_sym >>=1UL;
                    }

                    if (curr_sym < prev_sym) {//S_TYPE type
                        s_type = S_TYPE;
                    } else if (curr_sym == prev_sym) {
                        s_type = prev_s_type;
                    } else {//L_TYPE type
                        s_type = L_TYPE;
                        if (prev_s_type == S_TYPE && prev_is_rep && is_rep) {//LMS suffix

                            //process the previous phrase
                            assert(!phrase.empty());
                            process_phrase(phrase);

                            //TODO testing
                            /*for(size_t j=test_seq.size();j-->0;){
                                std::cout<<test_seq[j]<<" ";
                            }
                            std::cout<<" | "<<n_rep_sym<<std::endl;*/
                            //uniq_phrases+= (n_rep_sym<phrase.size());
                            //

                            //create the new phrase
                            phrase.clear();
                            //test_seq.clear();
                            phrase.push_back(prev_sym);
                            //test_seq.push_back(prev_is_rep);
                            //n_rep_sym = prev_is_rep;
                        }
                    }

                    /*if(is_rep!=prev_is_rep){
                        if(!prev_is_rep && run>1){
                            std::cout<<run<<std::endl;
                        }
                        run=0;
                    }*/

                    //run++;
                    phrase.push_back(curr_sym);
                    //test_seq.push_back(is_rep);
                    //n_rep_sym += is_rep;

                    prev_sym = curr_sym;
                    prev_s_type = s_type;
                    prev_is_rep = is_rep;
                }

                assert(!phrase.empty());
                process_phrase(phrase);

                /*if(!prev_is_rep && run>1){
                    std::cout<<run<<std::endl;
                }*/

                //TODO testing
                /*for(size_t j=test_seq.size();j-->0;){
                    std::cout<<test_seq[j]<<" ";
                }
                std::cout<<" | "<<n_rep_sym<<" fin string "<<std::endl;*/
                //
                //uniq_phrases+= (n_rep_sym<phrase.size());

                phrase.clear();
                //test_seq.clear();
            }
        }

        //std::cout<<"\n  There are "<<uniq_phrases<<" phrases that are not necessary to hash";
    }
};
#endif //LPG_COMPRESSOR_LC_PARSERS_HPP
