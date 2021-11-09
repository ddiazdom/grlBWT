//
// Created by ale on 11-05-21.
//

#ifndef LPG_COMPRESSOR_UTILS_HPP
#define LPG_COMPRESSOR_UTILS_HPP

#include <chrono>
#include <iostream>
#include <vector>

//create temporal folder
std::string create_temp_folder(std::string& base_folder, std::string const & prefix);

//a helper function to modify the output
bool ends_with(std::string const & value, std::string const & ending);

//remove extension
std::string remove_ext(std::string& input_string,
                       std::vector<std::string>& exts);
//produce an output name
std::string get_output_name(std::string& output_p,
                            std::string& alternative,
                            std::vector<std::string>& in_ext,//extension to remove
                            std::vector<std::string>& out_ext);

template<class time_t>
void report_time(time_t start, time_t end){
    auto dur = end - start;
    auto h = std::chrono::duration_cast<std::chrono::hours>(dur);
    auto m = std::chrono::duration_cast<std::chrono::minutes>(dur -= h);
    auto s = std::chrono::duration_cast<std::chrono::seconds>(dur -= m);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(dur -= s);
    if(h.count()>0){
        std::cerr<<"Elapsed time (h:mm:ss.ms): "<<h.count()<<":"<<m.count()<<":"<<s.count()<<"."<<ms.count()<<std::endl;
    }else{
        std::cerr<<"Elapsed time (mm:ss.ms): "<<m.count()<<":"<<s.count()<<"."<<ms.count()<<std::endl;
    }
}


#endif //LPG_COMPRESSOR_UTILS_HPP
