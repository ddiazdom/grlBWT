//
// Created by Diaz, Diego on 23.1.2022.
//
#include "bwt_io.h"
#include <vector>

int main(int argc, char** argv) {

    std::string tmp_file = "some_test";
    for(size_t j=0;j<700;j++) {

        uint8_t sb= 1 + rand() % 7;
        uint8_t fb= 1 + rand() % 7;
        size_t size = rand() % 1000000;

        bwt_buff_writer bwt_out(tmp_file, std::ios::out, sb, fb, 64);
        std::vector<std::pair<size_t, size_t>> tmp_vector;
        size_t max_sym = (1UL << (sb*8));
        size_t max_freq = (1UL << (fb*8));
        for(size_t i=0;i<size;i++){

            size_t sym = rand() % max_sym;
            size_t freq = rand() % max_freq;
            bwt_out.push_back(sym, freq);
            tmp_vector.emplace_back(sym, freq);

            if(i==1000){//checking random reading
                size_t s;
                size_t f;
                std::cout<<"pasee por acaa"<<std::endl;
                for(size_t k=0;k<200;k++){
                    size_t pos = rand() % i;
                    s = bwt_out.read_sym(pos);
                    f = bwt_out.read_freq(pos);
                    assert(s==tmp_vector[pos].first);
                    assert(f==tmp_vector[pos].second);
                    pos = rand() % i;
                    bwt_out.read_run(pos, s, f);
                    assert(s==tmp_vector[pos].first);
                    assert(f==tmp_vector[pos].second);
                }
            }/* else if(i==2000){//checking modification of randon positions
                size_t f;
                for(size_t k=0;k<200;k++){
                    size_t pos = rand() % i;
                    f = bwt_out.read_freq(pos);
                    size_t new_val;
                    if(f!=0){
                        new_val = rand() % f;
                    }else{
                        new_val = 0;
                    }
                    bwt_out.dec_freq(pos, new_val);
                    tmp_vector[pos].second-=new_val;
                    assert(tmp_vector[pos].second == bwt_out.read_freq(pos));

                    size_t pos2 = rand() % i;
                    f = bwt_out.read_freq(pos2);
                    size_t limit = max_freq-f;
                    if(f==(max_freq-1)){
                        new_val = 0;
                    }else{
                        new_val = rand() % limit;
                        assert((new_val+f)<max_freq);
                    }
                    bwt_out.inc_freq(pos2, new_val);
                    tmp_vector[pos2].second+=new_val;
                    assert(tmp_vector[pos2].second == bwt_out.read_freq(pos2));
                }
            } else if(i==3000){//checking writing of random positions
                for(size_t k=0;k<200;k++){
                    size_t pos = rand() % i;
                    sym = rand() % max_sym;
                    freq = rand() % max_freq;
                    tmp_vector[pos].first = sym;
                    tmp_vector[pos].second = freq;
                    bwt_out.write_sym(pos, sym);
                    bwt_out.write_freq(pos, freq);
                    assert(tmp_vector[pos].first==bwt_out.read_sym(pos));
                    assert(tmp_vector[pos].second==bwt_out.read_freq(pos));
                }
            }*/
        }
        bwt_out.close();

        std::cout<<"termine aca"<<std::endl;
        bwt_buff_reader bwt_in("some_test", 64);
        assert(bwt_in.size()==size);
        size_t sym=0, freq=0;
        for(size_t i=0;i<bwt_in.size();i++){
            size_t pos = rand() % bwt_in.size();
            bwt_in.read_run(pos, sym, freq);
            if(tmp_vector[pos].first!=sym ||
               tmp_vector[pos].second!=freq){
                std::cout<<i<<" "<<sym<<" "<<freq<<std::endl;
                std::cout<<tmp_vector[i].first<<" "<<tmp_vector[i].second<<std::endl;
            }
            assert(sym==tmp_vector[pos].first);
            assert(freq==tmp_vector[pos].second);

            pos = rand() % bwt_in.size();
            size_t tmp_sym = bwt_in.read_sym(pos);
            assert(tmp_sym==tmp_vector[pos].first);
        }
        bwt_in.close();

        if(remove(tmp_file.c_str())){
            std::cout<<"Error trying to remove file"<<std::endl;
        }
    }
}
