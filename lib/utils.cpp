//
// Created by Diaz, Diego on 9.11.2021.
//

#include "utils.hpp"
#include <filesystem>
#ifdef __APPLE__
#include <unistd.h>
#endif

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

//produce an output name
std::string get_output_name(std::string& output_p,
                            std::string& alternative,
                            std::vector<std::string>& in_ext,//extension to remove
                            std::vector<std::string>& out_ext){ //extension to add

    std::string output;
    if(output_p.empty()){
        std::filesystem::path fs(alternative);
        output = fs.filename();
        for(auto const& ext : in_ext){
            if(ends_with(output, ext)){
                output.erase(output.size()-ext.size());
            }
        }
    }else{
        output = output_p;
    }

    bool has_extension = false;
    for(auto const& ext : out_ext){
        if(ends_with(output, ext)){
            has_extension = true;
        }
    }
    if(!has_extension) output+=out_ext[0];

    return output;
}
