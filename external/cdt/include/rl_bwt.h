//
// Created by Diaz, Diego on 2.3.2022.
//

#ifndef CDT_RL_FMI_H
#define CDT_RL_FMI_H
#include <sdsl/wavelet_trees.hpp>
#include "succ_pred.h"
#include "si_int_array.h"
#include "utils.h"

//Implementation of the RLBWT (Makinen and Navarro 2005; Makinen et al. 2010; Gagie et al. 2020)
template<class bwt_type, class bv_type>
class rl_bwt {

    sdsl::wt_huff<> m_heads; //run heads
    std::vector<size_t> m_r_C; // Every m_r_C[i] stores the number of runs in the BWT whose symbol is smaller than i
    std::vector<size_t> m_C; // Every m_C[i] stores the number of symbols in the BWT that are smaller than i
    std::vector<uint8_t> m_sym_map; //Maps a byte symbol to its compressed alphabet [1,\sigma]
    std::vector<uint8_t> m_inv_sym_map; //Maps a compressed symbol in [1,\sigma] to its original byte alphabet

    // A sparse bit vector run_lengths=B[1,n] in which B[i] = 1 if BWT[i] is the first symbol in a run.
    // B is encoded as an integer vector V = x_1, x_2,..., x_z, where every x_i is a position such that B[x_i]=1.
    // As V is strictly increasing, evey V[i] is encoded as V[i]-V[i-1] for i>1.
    bv_type run_lengths;

    // A vector R[1,r]=sorted_run_length storing the cumulative sums of the run lengths (stably sorted
    // according their symbols). Let R[i] be the value associated with the jth BWT run for symbol c \in [1,\sigma].
    // R[i] stores the sum of the lengths of the BWT runs whose symbols are lexicographically smaller than c plus
    // the sum of the lengths of the j-1 runs for symbol c that precede the BWT run of R[i].
    bv_type m_sorted_run_lengths;

    //data structure for successor and predecessor queries over B. The operation sp.predecessor(i), receives as input
    // an index position i \in [1,n], and returns a pair (rank, i'), where i' is the rightmost position i'<=i with
    // B[i']=1. The value for rank is the rank value for the bit in B[i'].
    succ_pred<bv_type> sp;

    void move(rl_bwt<bwt_type, bv_type>&& other){
        m_heads.swap(other.m_heads);
        m_r_C.swap(other.m_r_C);
        m_C.swap(other.m_C);
        m_sym_map.swap(other.m_sym_map);
        m_inv_sym_map.swap(other.m_inv_sym_map);
        run_lengths.swap(other.run_lengths);
        m_sorted_run_lengths.swap(other.m_sorted_run_lengths);
        sp.swap(other.sp);
        sp.set_vector(&run_lengths);
    }

    void copy(const rl_bwt<bwt_type, bv_type>& other){
        m_heads = other.m_heads;
        m_r_C = other.m_r_C;
        m_C = other.m_C;
        m_sym_map = other.m_sym_map;
        m_inv_sym_map = other.m_inv_sym_map;
        run_lengths = other.run_lengths;
        m_sorted_run_lengths= other.m_sorted_run_lengths;
        sp = other.sp;
        sp.set_vector(&run_lengths);
    }

public:

    const std::vector<size_t>& r_C          = m_r_C;
    const std::vector<size_t>& C            = m_C;
    const std::vector<uint8_t>& sym_map     = m_sym_map;
    const std::vector<uint8_t>& inv_sym_map = m_inv_sym_map;
    const bv_type& sorted_run_lengths       = m_sorted_run_lengths;

    rl_bwt()=default;

    rl_bwt(std::string& rl_bwt_file, tmp_workspace& ws) : m_inv_sym_map(255, 0) {
        std::cout<<"Creating a RLBWT data structure from the input file "<<rl_bwt_file<<std::endl;
        bwt_type rl_bwt_io(rl_bwt_file);
        size_t run_len, sym, k=0;
        std::string tmp_file = ws.get_file("plain_rl_bwt");
        uint8_t buffer[8192]={0};
        size_t sym_freqs[255]={0};
        std::ofstream ofs(tmp_file, std::ios::binary);

        size_t len_dist[65535]={0};
        std::cout<<"  Parsing the BWT runs"<<std::endl;
        {
            size_t rl_acc = 0;
            for (size_t i = 0; i < rl_bwt_io.size(); i++) {
                rl_bwt_io.read_run(i, sym, run_len);
                buffer[k++] = (uint8_t) sym;
                if (k == 8192) {
                    ofs.write((char *)buffer, 8192);
                    k = 0;
                }
                sym_freqs[sym]++;
                run_lengths.push_back(rl_acc);
                len_dist[run_len]++;
                rl_acc += run_len;
            }
            //run_lengths.push_back(rl_acc);
            run_lengths.shrink_to_fit();
            succ_pred tmp_sp(run_lengths);
            sp.swap(tmp_sp);
            ofs.write((char *)buffer, (std::streamsize)k);
            ofs.close();

            //TODO testing
            double acc_frac=0, frac;
            for(size_t i=0;i<65535;i++){
                if(len_dist[i]!=0){
                    frac = double(len_dist[i])/double(run_lengths.size());
                    acc_frac +=frac;
                }
                if(i==255 || i==65534){
                    std::cout<<"    "<<std::fixed<<std::setprecision(2)<<(acc_frac*100)<<"% of the BWT runs have a length shorter than "<<i<<std::endl;
                }
            }
            //
        }

        std::cout<<"  Constructing the wavelet tree with the run heads"<<std::endl;
        sdsl::construct(m_heads, tmp_file, 1);
        ws.remove_file("plain_rl_bwt");

        std::vector<size_t> tmp_r_C;
        size_t max_sym, c_sym=0;
        for(size_t i=0;i<255;i++){
            if(sym_freqs[i]>0){
                tmp_r_C.push_back(sym_freqs[i]);
                m_inv_sym_map[i] = c_sym++;
                m_sym_map.push_back(i);
                max_sym=i;
            }
        }
        m_sym_map.shrink_to_fit();
        m_inv_sym_map.resize(max_sym+1);
        m_inv_sym_map.shrink_to_fit();

        size_t acc=0, tmp;
        m_r_C.resize(tmp_r_C.size()+1);
        for(size_t i=0;i<tmp_r_C.size();i++){
            tmp = tmp_r_C[i];
            m_r_C[i] = acc;
            tmp_r_C[i] = acc;
            acc +=tmp;
        }
        m_r_C[tmp_r_C.size()] = acc;
        m_C = std::vector<size_t>(m_r_C.size(), 0);

        std::cout<<"  Stably sorting the run lengths"<<std::endl;
        std::vector<size_t> tmp_sr_lengths(run_lengths.size());
        for(size_t i=0;i<rl_bwt_io.size();i++){
            rl_bwt_io.read_run(i, sym, run_len);
            c_sym = m_inv_sym_map[sym];
            assert(tmp_r_C[c_sym]<tmp_sr_lengths.size());
            tmp_sr_lengths[tmp_r_C[c_sym]++] = run_len;
            m_C[c_sym]+=run_len;
        }

        acc=0;
        for(unsigned long & i : m_C){
            tmp = i;
            i = acc;
            acc +=tmp;
        }
        m_C[m_C.size()-1] = acc;

        acc=0;
        for(unsigned long tmp_sr_length : tmp_sr_lengths){
            tmp = tmp_sr_length;
            m_sorted_run_lengths.push_back(acc);
            acc+=tmp;
        }
        m_sorted_run_lengths.push_back(acc);
        m_sorted_run_lengths.shrink_to_fit();
    }

    [[nodiscard]] std::pair<uint8_t, size_t> lf(size_t idx) const{
        //the expensive operation is pred(idx) as I have to decode 64 delta codes in the worst case.
        auto pred_res = sp.pred(idx);// pred_res = {rank, pred_res} (rank is 1-based)
        assert(pred_res.first>=1);
        auto inv_select = m_heads.inverse_select(pred_res.first-1);
        size_t c_sym = m_inv_sym_map[inv_select.second];
        size_t val = m_sorted_run_lengths[m_r_C[c_sym] + inv_select.first];
        return {inv_select.second, val + (idx-pred_res.second)};
    }

    [[nodiscard]] uint8_t alphabet() const {
        return m_heads.sigma;
    }

    [[nodiscard]] uint8_t max_byte_symbol() const {
        return m_sym_map.back();
    }

    [[nodiscard]] uint8_t min_byte_symbol() const {
        return m_sym_map[0];
    }

    void interval_symbols(size_t i, size_t j, uint64_t& k, std::vector<uint8_t>& symbols,
                          std::vector<uint64_t>& ranks_i,
                          std::vector<uint64_t>& ranks_j) const{

        assert(i<j);
        size_t c_sym, i_sym, j_sym;
        auto res_i = sp.pred(i);
        auto res_j = sp.pred(j);

        if(res_i.first==res_j.first){//input range falls within the same BWT run
            k=1;
            auto inv_res = m_heads.inverse_select(res_i.first-1);
            c_sym = m_inv_sym_map[inv_res.second];
            symbols[0] = inv_res.second;
            ranks_i[0] =  m_sorted_run_lengths[m_r_C[c_sym] + inv_res.first]-m_C[c_sym]+(i-res_i.second);
            ranks_j[0] = ranks_i[0] + (j-i);
        }else{

            m_heads.interval_symbols(res_i.first-1, res_j.first-(res_j.second==j), k, symbols, ranks_i, ranks_j);

            i_sym = m_heads[res_i.first-1];
            j_sym = m_heads[res_j.first-1];
            for(size_t u=0;u<k;u++){
                c_sym = m_inv_sym_map[symbols[u]];

                ranks_i[u] = m_sorted_run_lengths[m_r_C[c_sym] + ranks_i[u]]-m_C[c_sym];
                ranks_j[u] = m_sorted_run_lengths[m_r_C[c_sym] + ranks_j[u]]-m_C[c_sym];

                if(symbols[u]==i_sym){
                    ranks_i[u]+=(i-res_i.second);
                }

                if(symbols[u]==j_sym && res_j.second!=j && j<C.back()){
                    ranks_j[u] -= run_lengths[res_j.first]-j;
                }
            }
        }
    }

    size_t serialize(std::ostream &out) const{
        size_t written_bytes=0;
        written_bytes += m_heads.serialize(out, nullptr, "rl_bwt");
        written_bytes += serialize_plain_vector(out, m_r_C);
        written_bytes += serialize_plain_vector(out, m_C);
        written_bytes += serialize_plain_vector(out, m_sym_map);
        written_bytes += serialize_plain_vector(out, m_inv_sym_map);
        written_bytes += run_lengths.serialize(out);
        written_bytes += m_sorted_run_lengths.serialize(out);
        written_bytes += sp.serialize(out);
        return written_bytes;
    }

    void load(std::istream &in){
        m_heads.load(in);
        load_plain_vector(in, m_r_C);
        load_plain_vector(in, m_C);
        load_plain_vector(in, m_sym_map);
        load_plain_vector(in, m_inv_sym_map);
        run_lengths.load(in);
        m_sorted_run_lengths.load(in);
        sp.load(in);
        sp.set_vector(&run_lengths);
    }

    void swap(rl_bwt<bwt_type, bv_type>& other){
        move(std::move(other));
    }
};


#endif //CDT_RL_FMI_H
