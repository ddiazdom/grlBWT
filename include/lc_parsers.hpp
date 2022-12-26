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

        sym_type curr_sym, prev_sym;
        string_t phrase(2, sdsl::bits::hi(max_symbol)+1);
        size_t end_ps, start_ps;
        uint8_t type, rep;

        for(size_t str=l_string+1;str-->f_string;) {

            auto range = init_str(str);

            if(range.first<=range.second) { //if this is not true, it means the string was fully compressed

                start_ps = range.first;
                end_ps = range.second;

                prev_sym = ifs.read(end_ps);
                if constexpr (std::is_same<sym_type, uint8_t>::value){
                    rep = 3U;
                }else{
                    rep = prev_sym & 1U;
                    prev_sym >>=1UL;
                }

                phrase.push_back(prev_sym);
                type = 0;

                for (size_t i = end_ps; i-- > start_ps;) {

                    curr_sym = ifs.read(i);

                    if constexpr (!std::is_same<sym_type, uint8_t>::value){
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
};
#endif //LPG_COMPRESSOR_LC_PARSERS_HPP
