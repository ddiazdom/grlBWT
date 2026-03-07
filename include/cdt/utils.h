//
// Created by Diaz, Diego on 27.10.2021.
//

#ifndef CDT_UTILS_H

#define CDT_UTILS_H
#include <chrono>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <cassert>
#include <sstream>
#include "utils.h"
#include <random>
#include <sys/resource.h>

#ifdef __APPLE__
#include <unistd.h>
#endif

struct str_collection {
    std::vector<long> str_ptrs;
    size_t n_strings=0;
    size_t longest_string=0;
    size_t max_sym=0;
    size_t min_sym=std::numeric_limits<size_t>::max();
    size_t n_syms=0;
    size_t max_sym_freq=0;
};

inline bool check_gzip(const std::string& file){
    std::ifstream ifs(file, std::ios::binary);
    char f_bytes[2]={0,0};
    auto ext = std::filesystem::path(file).extension();
    ifs.read(f_bytes, 2);
    ifs.close();
    return ext==".gz" && f_bytes[0]=='\x1F' && f_bytes[1]=='\x8B';
}

inline bool is_fastx(const std::string& input_file){
    bool is_gz=check_gzip(input_file);
    char f_sym;
    //read the first symbol in the file to check if it is in FASTA or FASTQ format
    if(is_gz){
        //gzFile zfile = gzopen(input_file.c_str(), "rb");
        //gzread(zfile, &f_sym, 1);
        //gzclose(zfile);
        std::cout<<"gzipped files are not supported"<<std::endl;
        exit(0);
    }else{
        std::ifstream ifs(input_file, std::ios::binary);
        ifs.read(&f_sym, 1);
        ifs.close();
    }
    return f_sym=='>' || f_sym=='@';
}

inline bool folder_exists(const std::filesystem::path& p) {
    return std::filesystem::is_directory(p);
}

inline std::string random_string(size_t length){
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

inline std::string append_to_path(std::string& orig_path, std::string const& suffix){
    return std::filesystem::path(orig_path) / suffix;
}

//a helper function to modify the output
inline bool ends_with(std::string const & value, std::string const & ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

inline std::string remove_ext(std::string& input_string, std::vector<std::string>& exts){
    std::string res = input_string;
    for(auto const& ext: exts){
        if(ends_with(res, ext)){
            res.erase(res.size()-ext.size());
            break;
        }
    }
    return res;
}

inline std::string create_temp_folder(std::string& base_folder, std::string const & prefix){

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

inline void report_mem_peak(){
    rusage usage{};
    getrusage(RUSAGE_SELF, &usage);
    std::cout<<"\n peak : "<<usage.ru_maxrss<<" megabytes "<<std::endl;
}

struct tmp_workspace{

    std::string tmp_folder;
    std::string ext;
    bool remove_all;

    explicit tmp_workspace(std::string const& base_folder=std::filesystem::temp_directory_path(),
                           bool rem_all=true,
                           std::string const& prefix="tmp") : remove_all(rem_all) {

        std::string tmp_path = std::filesystem::canonical(std::filesystem::path(base_folder)) / std::string(prefix+".XXXXXX");
        char temp[200] = {0};
        tmp_path.copy(temp, tmp_path.size() + 1);
        temp[tmp_path.size() + 1] = '\0';
        auto res = mkdtemp(temp);
        if (res == nullptr) {
            std::cout << "Error trying to create a temporal folder" << std::endl;
        }
        tmp_folder = std::string(temp);
        ext = random_string(3);
    }

    explicit tmp_workspace(std::string folder,
                           std::string ext_,
                           bool rem_all) : tmp_folder(folder),
                                           ext(ext_),
                                           remove_all(rem_all){
        assert(folder_exists(tmp_folder));
    }

    [[nodiscard]] std::string get_file(std::string const& prefix) const {
        return std::filesystem::path(tmp_folder) / std::string(prefix+"_"+ext);
    }

    void remove_file(std::string const& prefix) const {
        std::filesystem::path file =  std::filesystem::path(tmp_folder) / std::string(prefix+"_"+ext);
        bool res = remove(file);
        if(!res){
            std::cout<<"Error trying to remove "<<file<<std::endl;
            exit(1);
        }
    }

    ~tmp_workspace(){
        if(remove_all){
            std::filesystem::remove_all(tmp_folder);
        }
    }

    [[nodiscard]] std::string folder() const{
        return tmp_folder;
    }
};

template<class sym_type>
str_collection collection_stats(std::string& input_file){

    //We assume that the smallest symbol in the collection is the separator symbol.
    // This should be the last character in the file
    str_collection str_coll;

    size_t sym_frq[256] = {0};

    std::ifstream ifs(input_file, std::ios::binary);

    sym_type sep_sym;
    ifs.seekg(-1*sizeof(sym_type), std::ios::end);

    ifs.read((char *)&sep_sym, 1*sizeof(sym_type));
    str_coll.n_syms = ifs.tellg()/sizeof(sym_type);
    assert(str_coll.n_syms>0);
    str_coll.max_sym_freq = str_coll.n_syms;
    ifs.seekg(0, std::ios::beg);
    size_t pos = 0, cont=0;

    std::streampos buff_size = 8388608;
    std::streampos n_elms = buff_size/sizeof(sym_type);
    auto buffer = (sym_type *) calloc(n_elms, sizeof(sym_type));

    size_t read_bytes, str_len, read_syms;
    sym_type sym=0, min_sym=std::numeric_limits<sym_type>::max(), max_sym=0;

    while(true){
        ifs.read((char *)buffer, buff_size);

        read_bytes = ifs.gcount();

        if(read_bytes>0){

            read_syms = read_bytes/sizeof(sym_type);
            for(size_t i=0;i<read_syms;i++){
                sym = buffer[i];

                if constexpr (std::is_same<sym_type, uint8_t>::value){
                    sym_frq[sym]++;
                }else{
                    if(sym>max_sym) max_sym = sym;
                    if(sym<min_sym) min_sym = sym;
                }

                cont++;
                if(sym==sep_sym){
                    str_coll.str_ptrs.push_back((long)pos);
                    str_len = cont - pos;
                    if(str_len>str_coll.longest_string) str_coll.longest_string=str_len;
                    pos = cont;
                }
            }
        }else{
            break;
        }
    }

    str_coll.str_ptrs.shrink_to_fit();

    if constexpr (std::is_same<sym_type, uint8_t>::value) {
        size_t i=0;
        while(sym_frq[i] == 0) i++;
        min_sym = i;

        i=255;
        while(sym_frq[i]==0) i--;
        max_sym = i;

        size_t max_sym_freq = 0;
        for(unsigned long frq : sym_frq){
            if(frq>max_sym_freq) max_sym_freq = frq;
        }
        str_coll.max_sym_freq = max_sym_freq;
    }

    if(sep_sym!=min_sym || sym!=sep_sym){
        std::cout<<"Error: the file is ill formed"<<std::endl;
        exit(1);
    }

    str_coll.min_sym = min_sym;
    str_coll.max_sym = max_sym;
    str_coll.n_strings = str_coll.str_ptrs.size();

    ifs.close();
    free(buffer);
    return str_coll;
}

template<class time_t>
void report_time(time_t start, time_t end, size_t padding){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);

    for(size_t i=0;i<padding;i++) std::cout<<" ";
    if(h.count()>0){
        std::cout<<"Elapsed time (hh:mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<h.count()<<":"<<std::setfill('0')<<std::setw(2)<<m.count()<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(m.count()>0){
        std::cout<<"Elapsed time (mm:ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(m.count())<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(s.count()>0){
        std::cout<<"Elapsed time (ss.ms): "<<std::setfill('0')<<std::setw(2)<<size_t(s.count())<<"."<<ms.count()<<std::endl;
    }else{
        std::cout<<"Elapsed time (ms): "<<ms.count()<<std::endl;
    }
}

//from https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

inline std::string report_space(off_t bytes){
    if(bytes<1000){
        return std::to_string(bytes)+" bytes";
    }

    if(bytes<1000000){
        float b = float(bytes)/1000;
        return to_string_with_precision(b, 2)+" KBs";
    }

    if(bytes < 1000000000L){
        float b = float(bytes)/1000000;
        return to_string_with_precision(b, 2)+" MBs";
    }

    if(bytes < 1000000000000L){
        double b = double(bytes)/1000000000L;
        return to_string_with_precision(b, 2)+" GBs";
    }

    double b = double(bytes)/1000000000000L;
    return to_string_with_precision(b, 2)+" TBs";
}
#endif //CDT_UTILS_H
