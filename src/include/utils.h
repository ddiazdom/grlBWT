//
// Created by Diaz, Diego on 27.10.2021.
//

#ifndef GRLBWT_UTILS_H
#define GRLBWT_UTILS_H
#include <chrono>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <cassert>
#include <functional>

#ifdef __APPLE__
#include <unistd.h>
#endif

#define INT_CEIL(a,b) (a>0? 1+(a-1)/b : 0)

struct str_collection {
    std::vector<long> str_ptrs;
    size_t n_strings=0;
    size_t longest_string=0;
    size_t max_sym=0;
    size_t min_sym=std::numeric_limits<size_t>::max();
    size_t n_syms=0;
    size_t max_sym_freq=0;

    [[nodiscard]] std::string to_string(size_t padding) const {
        std::string pad(padding,' ');
        std::string stat = "Stats: \n"+
                           pad+"  Smallest symbol               : "+std::to_string(min_sym)+"\n"+
                           pad+"  Greatest symbol               : "+std::to_string(max_sym)+"\n"+
                           pad+"  Number of symbols in the file : "+std::to_string(n_syms)+"\n"+
                           pad+"  Number of strings             : "+std::to_string(n_strings)+"\n"+
                           pad+"  Max sym freq.                 : "+std::to_string(max_sym_freq);
        return stat;
    }
};

//check if the file is gzipped
bool check_gzip(const std::string& file);

//check if file exists
bool file_exists(const std::filesystem::path& p, std::filesystem::file_status const& s = std::filesystem::file_status{});

//create temporal folder
std::string create_temp_folder(std::string& base_folder, std::string const & prefix);

//a helper function to modify the output
bool ends_with(std::string const & value, std::string const & ending);

//remove extension
std::string remove_ext(std::string& input_string, std::vector<std::string>& exts);

//produce a random string
std::string random_string(size_t length);

//check if it is in FASTx
bool is_fastx(const std::string& file);


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
        assert(file_exists(tmp_folder));
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
str_collection collection_stats(std::string& input_file);

template<class time_t>
void report_time(time_t start, time_t end, size_t padding){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);

    for(size_t i=0;i<padding;i++) std::cout<<" ";
    if(h.count()>0){
        std::cout<<"Elapsed time (hh:mm:ss.ms): "<<h.count()<<":"<<m.count()<<":"<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(m.count()>0){
        std::cout<<"Elapsed time (mm:ss.ms): "<<size_t(m.count())<<":"<<s.count()<<"."<<ms.count()<<std::endl;
    }else if(s.count()>0){
        std::cout<<"Elapsed time (ss.ms): "<<size_t(s.count())<<"."<<ms.count()<<std::endl;
    }else{
        std::cout<<"Elapsed time (ms): "<<ms.count()<<std::endl;
    }
}

template<class time_t>
std::string time2str(time_t start, time_t end){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);

    std::stringstream ss;
    if(h.count()>0){
        ss <<"(hh:mm:ss.ms) "<<std::setfill('0')<<std::setw(2)<<h.count()<<":"<<std::setfill('0')<<std::setw(2)<<m.count()<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count();
    }else if(m.count()>0){
        ss<<"(mm:ss.ms) "<<std::setfill('0')<<std::setw(2)<<size_t(m.count())<<":"<<std::setfill('0')<<std::setw(2)<<s.count()<<"."<<ms.count();
    }else if(s.count()>0){
        ss<<"(ss.ms) "<<std::setfill('0')<<std::setw(2)<<size_t(s.count())<<"."<<ms.count();
    }else{
        ss<<"(ms) "<<std::setfill('0')<<std::setw(2)<<ms.count();
    }
    return ss.str();
}

//from https://stackoverflow.com/questions/16605967/set-precision-of-stdto-string-when-converting-floating-point-values
template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return std::move(out).str();
}

std::string report_space(off_t bytes);

// Define log levels
enum LogLevel {
    LOG_NONE = 0,
    LOG_ERROR = 1,
    LOG_WARNING = 2,
    LOG_INFO = 3,
};

class logger {
public:

    explicit logger(LogLevel lvl_): log_lvl(lvl_){};

    // Logging methods
    /*void log(const std::string& message) {
        std::lock_guard<std::mutex> guard(logMutex); // Ensure thread safety
        std::cout <<pad<< message << std::endl;
    }*/

    void error(const std::string& message) {
        if (log_lvl>=LOG_ERROR){
            std::cout <<pad<<"ERROR: "<< message << std::endl;
        }
    }

    void warning(const std::string& message) {
        if (log_lvl>=LOG_WARNING){
            std::cout <<pad<<"WARNING: "<< message << std::endl;
        }
    }

    inline void info(const std::string& message) const {
        if(log_lvl>=LOG_INFO){
            std::cout <<pad<< message << std::endl;
        }
    }

    inline void info_start(const std::string& message) const {
        if(log_lvl>=LOG_INFO){
            std::cout <<pad<< message << std::flush;
        }
    }
    inline void info_app(const std::string& message) const {
        if(log_lvl>=LOG_INFO){
            std::cout << message << std::flush;
        }
    }
    inline void info_end(const std::string& message) const {
        if(log_lvl>=LOG_INFO){
            std::cout << message << std::endl;
        }
    }


    // Measure and log the execution time of a block of code
    void measure_and_log_time(LogLevel level, const std::function<void()>& codeBlock, const std::string& description = "") {
        if (log_lvl>=level) {
            std::cout <<pad<< description << std::flush;
            auto start = std::chrono::high_resolution_clock::now();
            codeBlock();
            auto end = std::chrono::high_resolution_clock::now();
            //std::lock_guard<std::mutex> guard(logMutex); // Ensure thread safety
            report_time(start, end, pad.length());
        } else {
            codeBlock();
        }
    }

    [[nodiscard]] size_t get_padding() const {
        return pad.size();
    }

    void inc_pad(){
        pad.push_back(' ');
        pad.push_back(' ');
    }

    void dec_pad(){
        pad.pop_back();
        pad.pop_back();
    }

private:
    std::string pad;
    LogLevel log_lvl;
};

uint8_t sym_width(unsigned long val);

size_t next_power_of_two(unsigned long val);

size_t prev_power_of_two(unsigned long val);

bool is_power_of_two(unsigned long val);

template<class data_type>
void load_from_file(std::string const& file, data_type& dt){
    std::ifstream ifs(file, std::ios::binary);
    dt.load(ifs);
    ifs.close();
}

template<class data_type>
void store_to_file(std::string const& file, data_type& dt){
    std::ofstream ofs(file, std::ios::binary);
    dt.serialize(ofs);
    ofs.close();
}

template<class vector_t>
size_t serialize_plain_vector(std::ostream& ofs, vector_t& vector){
    size_t n = vector.size();
    ofs.write((char *)&n, sizeof(n));
    ofs.write((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*n));
    return sizeof(n)+ sizeof(typename vector_t::value_type)*n;
}

template<class size_type>
size_t serialize_raw_vector(std::ostream& ofs, size_type * vector, size_t len){
    ofs.write((char *)&len, sizeof(len));
    ofs.write((char *)vector, (std::streamsize)(sizeof(size_type)*len));
    return sizeof(len)+ sizeof(size_type)*len;
}

template<class val_type>
size_t serialize_elm(std::ostream& ofs, val_type value){
    ofs.write((char *)&value, sizeof(val_type));
    return sizeof(val_type);
}

template<class vector_t>
void load_plain_vector(std::istream& ifs, vector_t& vector){
    size_t n=0;
    ifs.read((char *)&n, sizeof(n));
    vector.resize(n);
    ifs.read((char *)vector.data(), (std::streamsize)(sizeof(typename vector_t::value_type)*n));
}

template<class size_type, class len_type>
void load_raw_vector(std::istream& ifs, size_type*& vector, len_type& len){
    ifs.read((char *)&len, sizeof(len_type));
    if(vector== nullptr){
        vector = (size_type *) malloc(sizeof(size_type)*len);
    }else{
        vector = (size_type *) realloc(vector, sizeof(size_type)*len);
    }
    ifs.read((char *)vector, (std::streamsize)(sizeof(size_type)*len));
}

template<class val_type>
void load_elm(std::istream& ifs, val_type& value){
    ifs.read((char *)&value, sizeof(val_type));
}

template<class vector_type>
void store_pl_vector(std::string const& file, vector_type& vector){
    std::ofstream ofs(file, std::ios::binary);
    serialize_plain_vector<vector_type>(ofs, vector);
    ofs.close();
}

template<class vector_type>
void load_pl_vector(std::string const& file, vector_type& vector){
    std::ifstream ifs(file, std::ios::binary);
    load_plain_vector<vector_type>(ifs, vector);
    ifs.close();
}

void report_mem_peak();
#endif //GRLBWT_UTILS_H
