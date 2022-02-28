//
// Created by Diaz, Diego on 28.10.2021.
//

#ifndef LHTIGS_FASTQ_WRITER_H
#define LHTIGS_FASTQ_WRITER_H
#include <fstream>
#include <utility>
#include <zlib.h>
#include "common.h"
#include "dna_string.h"

#define insert_sym() \
sym_freqs[(uint8_t)out_buffer[k]]++; \
k++; \
if(k==buff_size){ \
    ofs.write(out_buffer, buff_size); \
    k=0; \
} \

template<bool is_gzipped>
struct fasta_parser {
    std::string input_file;
    std::string output_file;
    gzFile in_zfile=nullptr;
    std::ifstream ifs;
    std::ofstream ofs;
public:
    fasta_parser(std::string  ifile, std::string ofile): input_file(std::move(ifile)),
                                                         output_file(std::move(ofile)){};
    str_collection operator()();
};

template<bool is_gzipped>
class fastq_parser {
    std::string input_file;
    std::string output_file;
    gzFile in_zfile=nullptr;
    std::ifstream ifs;
    std::ofstream ofs;
public:
    fastq_parser(std::string  ifile, std::string ofile): input_file(std::move(ifile)),
                                                         output_file(std::move(ofile)){};
    str_collection operator()();
};


//transform an input FASTA/Q file in one-string-per-line format
str_collection fastx2plain(const std::string& input_file, std::string const& output_file);

template<bool is_gzipped>
str_collection fasta_parser<is_gzipped>::operator()(){

    if constexpr(is_gzipped){
        in_zfile = gzopen(input_file.c_str(), "rb");
    }else{
        ifs.open(input_file.c_str(), std::ios::binary | std::ios::in);
        assert(ifs.good());
    }

    std::streamsize k=0;
    std::streamsize buff_size = 8192;
    std::streamsize read_bytes;
    char * in_buffer, *out_buffer;
    size_t sym_freqs[255]={0};
    str_collection str_coll;

    in_buffer = (char *)calloc(buff_size, 1);
    out_buffer = (char *)calloc(buff_size, 1);

    ofs.open(output_file, std::ios::out | std::ios::binary);

    bool header=false;
    bool first = true;

    while (true) {

        if constexpr(is_gzipped) {
            read_bytes = (std::streamsize) gzread(in_zfile, in_buffer, buff_size);
        }else{
            ifs.read((char *)in_buffer, buff_size);
            read_bytes = ifs.gcount();
        }

        if (read_bytes > 0) {
            for(std::streamsize j=0;j<read_bytes;j++){

                if(in_buffer[j]=='>'){
                    if(!first){
                        out_buffer[k] = '\n';
                        insert_sym()
                    }
                    first=false;
                    header = true;
                    str_coll.n_strings++;
                }else if (in_buffer[j]=='\n'){
                    header = false;
                }else if(!header){
                    out_buffer[k] = (char )toupper(in_buffer[j]);
                    insert_sym()
                }
            }
        } else {
            break;
        }
    }

    out_buffer[k] = '\n';
    sym_freqs[(uint8_t)out_buffer[k]]++;
    k++;
    ofs.write(out_buffer, k);

    if constexpr(is_gzipped){
        gzclose(in_zfile);
    }else{
        ifs.close();
    }
    ofs.close();
    free(in_buffer);
    free(out_buffer);

    for(uint8_t i=0;i<255;i++){
        if(sym_freqs[i]!=0){
            str_coll.alphabet.push_back(i);
            str_coll.sym_freqs.push_back(sym_freqs[i]);
            str_coll.n_char+=sym_freqs[i];
        }
    }
    str_coll.min_sym = str_coll.alphabet[0];
    str_coll.max_sym = str_coll.alphabet.back();

    return str_coll;
}

template<bool is_gzipped>
str_collection fastq_parser<is_gzipped>::operator()() {

    if constexpr(is_gzipped){
        in_zfile = gzopen(input_file.c_str(), "rb");
    }else{
        ifs.open(input_file.c_str(), std::ios::binary | std::ios::in);
        assert(ifs.good());
    }

    std::streamsize read_bytes;
    std::streamsize k=0;
    std::streamsize buff_size = 8192;
    char *in_buffer, *out_buffer;
    in_buffer = (char *)calloc(buff_size, 1);
    out_buffer = (char *)calloc(buff_size, 1);
    ofs.open(output_file, std::ios::out | std::ios::binary);

    str_collection str_coll;
    size_t sym_freqs[255]={0};

    size_t line=0;
    while (true) {

        if constexpr(is_gzipped) {
            read_bytes = (std::streamsize) gzread(in_zfile, in_buffer, buff_size);
        }else{
            ifs.read((char *)in_buffer, buff_size);
            read_bytes = ifs.gcount();
        }

        if (read_bytes > 0) {
            for(std::streamsize i=0;i<read_bytes;i++){
                if(line==1){
                    out_buffer[k] = (char)in_buffer[i];
                    insert_sym()
                }
                if(in_buffer[i]=='\n'){
                    line = (line+1) & 3UL;
                    if(line==1) str_coll.n_strings++;
                }
            }
        } else {
            break;
        }
    }
    ofs.write(out_buffer, k);
    ofs.close();
    free(in_buffer);
    free(out_buffer);

    if constexpr(is_gzipped){
        gzclose(in_zfile);
    }else{
        ifs.close();
    }

    for(uint8_t i=0;i<255;i++){
        if(sym_freqs[i]!=0){
            str_coll.alphabet.push_back(i);
            str_coll.sym_freqs.push_back(sym_freqs[i]);
            str_coll.n_char+=sym_freqs[i];
        }
    }

    str_coll.min_sym = str_coll.alphabet[0];
    str_coll.max_sym = str_coll.alphabet.back();

    return str_coll;
}

#endif //LHTIGS_FASTQ_WRITER_H
