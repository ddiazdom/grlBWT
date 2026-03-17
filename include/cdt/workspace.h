//
// Created by Diaz, Diego on 9.3.2026.
//

#ifndef GRLBWT_TMP_FOLDER_H
#define GRLBWT_TMP_FOLDER_H
#include <string>
#include <filesystem>
#include <iostream>

struct workspace{

    static std::filesystem::path system_tmp();
    std::string get_file_name(std::string const& prefix) const;
    void remove_file(std::string const& prefix) const;

    std::string workspace_folder() const;

    static void init(std::string const& base_path= system_tmp(),
                     std::string const& prefix="tmp",
                     bool rem_all=true);

    static workspace& get_instance();
    ~workspace();


private:
    std::filesystem::path workspace_path;
    std::string ext;
    bool remove_all;
    static std::unique_ptr<workspace> singleton;

    static std::string random_string(size_t len);

    void create_folder(std::string const& base_path, std::string const& prefix) {
        //iterate until finding a temporary folder that does not exist
        for (size_t i=0;i<100;i++){
            std::filesystem::path dir = base_path + "/" + prefix + "_" + random_string(6);
            if(!std::filesystem::exists(dir)) {
                workspace_path = dir;
                return;
            }
        }
        throw std::runtime_error("Failed to create temporary workspace");
    }

    explicit workspace(std::string const& base_path= system_tmp(),
                       std::string const& prefix="tmp",
                       const bool rem_all=true) : remove_all(rem_all) {

        //check bash_path is a valid directory
        assert(std::filesystem::is_directory(base_path));
        create_folder(base_path, prefix);
        std::filesystem::create_directories(workspace_path);
        ext = random_string(3);
    }

};

std::unique_ptr<workspace> workspace::singleton = nullptr;

inline std::filesystem::path workspace::system_tmp() {
    return std::filesystem::temp_directory_path();
}

inline workspace::~workspace(){
    if(remove_all){
        std::filesystem::remove_all(workspace_path);
    }
}

inline workspace& workspace::get_instance() {
    if (singleton == nullptr) {
        singleton = std::unique_ptr<workspace>(new workspace());
    }
    return *singleton;
}

inline void workspace::init(std::string const& base_path, std::string const& prefix, bool rem_all) {
    if(singleton != nullptr) {
        throw std::runtime_error("workspace already initialized");
    }
    singleton = std::unique_ptr<workspace>(new workspace(base_path, prefix, rem_all));
}

inline std::string workspace::random_string(const size_t len) {

    static const char alphabet[] =
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789";

    std::mt19937 rng(std::random_device{}());
    std::uniform_int_distribution<> dist(0, sizeof(alphabet)-2);
    std::string s;
    for (int i = 0; i < len; ++i) {
        s += alphabet[dist(rng)];
    }
    return s;
}

inline std::string workspace::get_file_name(std::string const& prefix) const {
    return workspace_path/(prefix+"_"+ext);
}

inline void workspace::remove_file(std::string const& prefix) const {
    std::filesystem::path file =  workspace_path/(prefix+"_"+ext);
    if(!remove(file)){
        std::cout<<"Error trying to remove "<<file<<std::endl;
        exit(1);
    }
}

inline std::string workspace::workspace_folder() const {
    return workspace_path.string();
}

#define SYSTEM_TMP workspace::system_tmp()
#define INIT_TMP_FOLDER(base, prefix, rem_all) workspace::init(base, prefix, rem_all)
#define TMP_WORKSPACE workspace::get_instance().workspace_folder()
#define TMP_FILE_NAME(prefix) workspace::get_instance().get_file_name(prefix)
#define TMP_REMOVE_FILE(prefix) workspace::get_instance().remove_file(prefix)

#endif //GRLBWT_TMP_FOLDER_H