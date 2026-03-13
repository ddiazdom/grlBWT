//
// Created by Diaz, Diego on 10.3.2026.
//

#include<iostream>
#include <grlbwt/bwt_io.h>
#include <cdt/rle_vbyte_streams.h>
#include <cdt/file_streams.hpp>
#include <cdt/logger.h>

std::pair<size_t , size_t> speed_old(const std::string& old_rle_file) {
    TRACE_SCOPE();
    bwt_buff_reader bwt_reader(old_rle_file);
    size_t sym=0, len=0;
    for(size_t i=0;i<bwt_reader.size();i++) {
        bwt_reader.read_run(i, sym, len);
    }
    return std::make_pair(sym, len);
}

std::pair<uint64_t, uint64_t> speed_new(const std::string& new_rle_file) {
    TRACE_SCOPE();
    rle_vbyte_buff_reader bwt_reader(new_rle_file);
    uint64_t sym=0, len=0;
    while(bwt_reader.has_next()) {
        bwt_reader.next_run(sym, len);
    }
    return std::make_pair(sym, len);
}

void test_speed(const std::string& old_rle, const std::string& new_rle) {
    Logger::level = log_level::TRACE;
    auto res = speed_old(old_rle);
    std::cout<<"Old RLE: "<<res.first<<" "<<res.second<<std::endl;
    auto res2 = speed_new(new_rle);
    std::cout<<"New RLE: "<<res2.first<<" "<<res2.second<<std::endl;
}

void check_correctness_read_write(const std::string& input_file, const std::string& output_file) {

    bwt_buff_reader bwt_reader(input_file);
    std::cout<<"Parsing to the new RLE format"<<std::endl;
    rle_vbyte_buff_writer new_bwt_writer(output_file);
    size_t sym=0, len=0;
    for(size_t i=0;i<bwt_reader.size();i++) {
        bwt_reader.read_run(i, sym, len);
        new_bwt_writer.push_back(sym, len);
    }
    new_bwt_writer.close();

    std::cout<<"Checking they produce the same output"<<std::endl;
    rle_vbyte_buff_reader new_rlbwt_reader(output_file);
    uint64_t sym2, len2, i=0;
    while(new_rlbwt_reader.has_next()){
        new_rlbwt_reader.next_run(sym2, len2);
        bwt_reader.read_run(i++, sym, len);
        assert(sym==sym2);
        assert(len==len2);
    }
    assert(new_rlbwt_reader.n_read_runs()==bwt_reader.size());
    std::cout<<"Done, everything seems to be correct"<<std::endl;
}

void check_correctness_updater(const std::string& input_file, const std::string& output_file) {

    std::cout<<"Checking the updater works correctly"<<std::endl;

    {
        bwt_buff_reader bwt_reader(input_file);
        bwt_buff_writer bwt_writer(input_file+"_updated", std::ios::out, 5, 5);
        rle_vbyte_buff_updater rlbwt_updater(output_file);
        size_t sym=0, len=0;
        uint64_t sym2, len2;
        for(size_t i=0;i<bwt_reader.size();i++) {
            bwt_reader.read_run(i, sym, len);
            rlbwt_updater.next_run(sym2, len2);
            //std::cout<<i<<std::endl;
            assert(sym==sym2);
            if(sym>1) {
                sym=1;
                rlbwt_updater.update_sym(sym);
            }

            assert(len==len2);
            if (len>1) {
                len=1;
                rlbwt_updater.update_len(1);
            }
            bwt_writer.push_back(sym, len);
        }
        bwt_reader.close();
        bwt_writer.close();
        rlbwt_updater.close();
    }

    {
        bwt_buff_reader bwt_reader(input_file+"_updated");
        rle_vbyte_buff_reader rlbwt_reader(output_file);
        size_t sym=0, len=0;
        uint64_t sym2, len2;
        for(size_t i=0;i<bwt_reader.size();i++) {
            bwt_reader.read_run(i, sym, len);
            assert(rlbwt_reader.has_next());
            rlbwt_reader.next_run(sym2, len2);
            assert(sym==sym2);
            assert(len==len2);
        }
        assert(!rlbwt_reader.has_next());
    }
    std::cout<<"Done, no errors were found"<<std::endl;
    std::filesystem::remove(input_file+"_updated");
}

void check_correctness_updater_large_random_values(size_t n) {
    uint64_t lower_bound = 1;
    uint64_t upper_bound = 0xFFFFFFFFFFFFFF;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution distr(lower_bound, upper_bound);

    //create a sequence of large strings
    bwt_buff_writer bwt_writer("seq_old_enc", std::ios::out, 7, 7);
    rle_vbyte_buff_writer rlbwt_wt("seq_new_enc", 512);
    for(int i=0;i<n;i++) {
        uint64_t sym = distr(gen);
        uint64_t len = distr(gen);
        bwt_writer.push_back(sym, len);
        rlbwt_wt.push_back(sym, len);
    }
    rlbwt_wt.close();
    bwt_writer.close();
    test_speed("seq_old_enc", "seq_new_enc");
    //

    //update the sequence
    bwt_buff_reader bwt_rd("seq_old_enc");
    bwt_buff_writer bwt_wt_up("seq_old_enc_updated", std::ios::out, 7, 7);
    rle_vbyte_buff_updater rlbwt_updater("seq_new_enc");
    size_t sym=0, len=0;
    uint64_t sym2, len2;
    for(size_t i=0;i<bwt_rd.size();i++) {
        bwt_rd.read_run(i, sym, len);
        rlbwt_updater.next_run(sym2, len2);
        assert(sym==sym2);
        assert(len==len2);

        sym = distr(gen) % sym;
        len = distr(gen) % len;

        rlbwt_updater.update_sym(sym);
        rlbwt_updater.update_len(len);
        bwt_wt_up.push_back(sym, len);
    }
    bwt_rd.close();
    bwt_wt_up.close();
    rlbwt_updater.close();
    std::filesystem::remove("seq_old_enc");
    assert(!rlbwt_updater.has_next());
    //

    //compare the updated sequences
    bwt_buff_reader bwt_rd2("seq_old_enc_updated");
    rle_vbyte_buff_reader rlbwt_rd2("seq_new_enc");
    size_t sym3=0, len3=0;
    uint64_t sym4, len4;
    for(size_t i=0;i<bwt_rd2.size();i++) {
        bwt_rd2.read_run(i, sym3, len3);
        rlbwt_rd2.next_run(sym4, len4);
        assert(sym3==sym4);
        assert(len3==len4);
    }
    assert(!rlbwt_rd2.has_next());
}

int main(int argc, char** argv) {

    if (argc != 3) {
        std::cout << "usage: ./grlbwt2rle file.rlbwt output_prefix\n\n"
                     "file.rlbwt: the BCR BWT file generated by grlBWT\n"
                     "output_prefix: output prefix for the RLE files" << std::endl;
        exit(0);
    }

    const auto input_file = std::string(argv[1]);
    const auto new_file = std::string(argv[2])+"_new.rlbwt";

    //check_correctness_read_write(input_file, new_file);
    //test_speed(input_file, new_file);
    //check_correctness_updater(input_file, new_file);
    check_correctness_updater_large_random_values(100000000);
}