//
// Created by Diaz, Diego on 27.10.2021.
//
#include "utils.h"
#include <filesystem>
#ifdef __APPLE__
#include <unistd.h>
#endif
#include <random>
#include <zlib.h>

bool is_fastx(const std::string& input_file){

    bool is_gz=check_gzip(input_file);
    char f_sym;

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
    return f_sym=='>' || f_sym=='@';
}

bool file_exists(const std::filesystem::path& p, std::filesystem::file_status const& s){
    if(std::filesystem::status_known(s) ? std::filesystem::exists(s) : std::filesystem::exists(p)){
        return true;
    }else{
        return false;
    }
}

std::string random_string(size_t length){
    const std::string characters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
    std::random_device random_device;
    std::mt19937 generator(random_device());
    std::uniform_int_distribution<size_t> distribution(0, characters.size() - 1);
    std::string random_string;
    random_string.resize(length);
    for(std::size_t i = 0; i < length; i++) {
        random_string[i] = characters[distribution(generator)];
    }
    return random_string;
}

bool check_gzip(const std::string& file){
   std::ifstream ifs(file, std::ios::binary);
   char f_bytes[2]={0,0};
   auto ext = std::filesystem::path(file).extension();
   ifs.read(f_bytes, 2);
   ifs.close();
   return ext==".gz" && f_bytes[0]=='\x1F' && f_bytes[1]=='\x8B';
}

std::string append_to_path(std::string& orig_path, std::string const& suffix){
    return std::filesystem::path(orig_path) / suffix;
}

std::string remove_ext(std::string& input_string, std::vector<std::string>& exts){
    std::string res = input_string;
    for(auto const& ext: exts){
        if(ends_with(res, ext)){
            res.erase(res.size()-ext.size());
            break;
        }
    }
    return res;
}

std::string create_temp_folder(std::string& base_folder, std::string const & prefix){
    std::string tmp_path = base_folder + "/"+prefix+".XXXXXX";
    char temp[200] = {0};
    tmp_path.copy(temp, tmp_path.size() + 1);
    temp[tmp_path.size() + 1] = '\0';
    auto res = mkdtemp(temp);
    if (res == nullptr) {
        std::cout << "Error trying to create a temporal folder" << std::endl;
    }
    return {temp};
}

//a helper function to modify the output
bool ends_with(std::string const & value, std::string const & ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

void report_mem_peak(){

}

str_collection collection_stats(std::string& input_file){
    str_collection str_coll;
    size_t sym_frq[256] = {0};
    std::ifstream ifs(input_file, std::ios::binary);
    char buffer[8192]={0};
    size_t read_bytes;
    std::streampos buff_size=8192;

    while(true){
        ifs.read((char *)&buffer, buff_size);
        read_bytes = ifs.gcount();
        if(read_bytes>0){
            for(size_t i=0;i<read_bytes;i++){
                sym_frq[(uint8_t)buffer[i]]++;
            }
        }else{
            break;
        }
    }

    for(size_t i=0;i<255;i++){
        if(sym_frq[i]!=0){
            str_coll.alphabet.push_back(i);
            str_coll.sym_freqs.push_back(sym_frq[i]);
            str_coll.n_char+=sym_frq[i];
        }
    }

    str_coll.n_strings = sym_frq[10];
    str_coll.min_sym = str_coll.alphabet[0];
    str_coll.max_sym = str_coll.alphabet.back();
    ifs.close();
    return str_coll;
}
