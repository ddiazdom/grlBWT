//
// Created by Diaz, Diego on 28.10.2021.
//

#include "fastx_handler.h"
#include "utils.h"

str_collection fastx2plain(const std::string& input_file, std::string const &output_file){

    bool is_fasta=false, is_fastq=false, is_gz=check_gzip(input_file);
    char f_sym;
    str_collection str_coll;

    //read the first symbol in the file to check if it is in FASTA or FASTQ format
    if(is_gz){
        gzFile zfile = gzopen(input_file.c_str(), "rb");
        gzread(zfile, &f_sym, 1);
        gzclose(zfile);
    }else{
        std::ifstream ifs(input_file, std::ios::binary);
        ifs.read(&f_sym, 1);
        ifs.close();
    }

    if(f_sym=='>'){
        is_fasta = true;
    }else if(f_sym=='@'){
        is_fastq = true;
    }

    if(is_fasta){
        if(is_gz){
            fasta_parser<true> fp(input_file, output_file);
            str_coll = fp();
        }else{
            fasta_parser<false> fp(input_file, output_file);
            str_coll = fp();
        }
    }else if(is_fastq){
        if(is_gz){
            fastq_parser<true> fp(input_file, output_file);
            str_coll = fp();
        }else{
            fastq_parser<false> fp(input_file, output_file);
            str_coll = fp();
        }
    }else{
        std::cout<<"Error trying to infer the input"<<std::endl;
    }
    return str_coll;
}
