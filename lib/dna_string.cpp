//
// Created by Diaz, Diego on 15.10.2021.
//
#include "../include/dna_string.h"

uint8_t dna_string::rev_comp[170] ={0,0,0,0,0,0,0,0,0,0,10,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,84,0,71,0,0,0,67,0,0,0,0,0,0,78,0,
                                    0,0,0,0,65,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,0,0,168,169,0,0,142,143,0,0,0,0,
                                    0,0,134,135,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                                    0,0,0,0,0,0,0,0,130,131};

std::ostream& operator<<(std::ostream& os, const dna_string& seq) {
    size_t i=0;
    for(auto const& sym : seq){
        if(sym<=84){
            os << char(sym);
        }else{
            os << char((sym >>1U));
        }
        if((sym & 1UL)){
            os << "->[" <<seq.get_run(i++)<<"]";
        }
    }
    return os;
}

void get_rev_comp(const std::string& plain_file, str_collection& str_coll){
    std::ifstream ifs;
    std::ofstream ofs;
    std::streampos buff_size = 8192, read_bytes;
    char *in_buffer, *out_buffer;
    in_buffer = (char *)calloc(buff_size, 1);
    out_buffer = (char *)calloc(buff_size, 1);

    ifs.open(plain_file, std::ios::binary | std::ios::in);
    ofs.open(plain_file, std::ios::binary | std::ios::out | std::ios::app);
    ifs.seekg(0, std::ios::end);
    ofs.seekp(0, std::ios::end);
    std::streampos n_bytes = ifs.tellg() - (std::streampos)1;
    while(true){
        read_bytes = std::min<std::streampos>(n_bytes, 8182);
        ifs.seekg(n_bytes - (std::streampos)read_bytes);
        ifs.read((char *)in_buffer, read_bytes);

        for(size_t j=read_bytes, i=0;j-->0;i++){
            out_buffer[i] = (char)dna_string::rev_comp[(uint8_t)in_buffer[j]];
        }
        ofs.write(out_buffer, read_bytes);
        n_bytes -= read_bytes;
        if(n_bytes==0) break;
    }
    out_buffer[0] = '\n';
    ofs.write(out_buffer, 1);
    ofs.close();
    ifs.close();
    free(in_buffer);
    free(out_buffer);

    str_coll.n_strings*=2;
    str_coll.n_char*=2;
    size_t inv_alph[255]={0};
    for(size_t i=0;i<str_coll.alphabet.size();i++){
        inv_alph[str_coll.alphabet[i]] = i;
    }

    std::vector<size_t> new_freqs(str_coll.alphabet.size(), 0);

    for(size_t i=0;i<str_coll.alphabet.size();i++){
        new_freqs[i] = str_coll.sym_freqs[i] + str_coll.sym_freqs[inv_alph[dna_string::rev_comp[str_coll.alphabet[i]]]];
    }
    std::swap(new_freqs, str_coll.sym_freqs);
}

