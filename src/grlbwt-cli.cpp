#include "CLI11.hpp"
#include <grlbwt/grlbwt.hpp>
#include <version.h>//do not delete it (built dynamically to print the program version)

std::string version_string() {
    std::ostringstream out;
    out << PROJECT_NAME << " " << PROJECT_VERSION << "\n";
    out << "commit: " << GIT_COMMIT << "\n";
    out << "built: " << BUILD_DATE << " " << BUILD_TIME << "\n";
    out << "build: " << BUILD_TYPE << "\n";
    return out.str();
}

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads{};
    float hbuff_frac=0.5;
    uint8_t alph_bytes=1;
    log_level log_lvl=log_level::INFO;
};

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option *) const override { return ""; }
};

static void parse_app(CLI::App& app, arguments& args){

	const auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(25);
    app.formatter(fmt);

    app.add_option("TEXT",
                      args.input_file,
                      //"Input file in one-string-per-line or FASTA/Q format (automatically detected)"
                      "Input file in one-string-per-line format")->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output-file",
                      args.output_file,
                      "Output file")->type_name("STR");
    app.add_option("-a,--alphabet", args.alph_bytes, "Number of bytes for the alphabet (def. 1)")->
            check(CLI::IsMember({1, 2, 4, 8}))->default_val(1)->type_name("UINT");
    app.add_option("-t,--threads",
                      args.n_threads,
                      "Maximum number of working threads")->default_val(1)->type_name("UINT");
    app.add_option("-f,--hbuff",
                      args.hbuff_frac,
                      "Parsing will use INPUT_SIZE*f bytes. O=no limit (def. 0.5)")->check(CLI::Range(0.0,1.0))->default_val(0.15)->type_name("FLOAT");
    app.add_option("-T,--tmp",
                      args.tmp_dir,
                      "Temporary folder (def. /tmp/grl.bwt.xxxx)")->
            check(CLI::ExistingDirectory)->default_val(std::filesystem::temp_directory_path())->type_name("DIR");
    app.add_option("-l,--log-level",
        args.log_lvl,
        "Verbosity level (WARN=1, INFO=2, DEBUG=3, TRACE=4, def. 2)")->check(CLI::Range(0, 4))->type_name("UINT");
    app.set_version_flag("-v,--version",
        version_string(),
        "Print the software version and exit");
}

template<class sym_type>
void build_bwt_int(arguments& args){
    LOG_INFO("Alphabet type:    "+std::to_string(args.alph_bytes));
    LOG_INFO("BWT type:         BCR");
    grl_bwt_algo<sym_type>(args.input_file, args.output_file, args.n_threads,
                           args.hbuff_frac, args.tmp_dir, args.log_lvl);
    LOG_INFO("The resulting BCR BWT was stored in "+args.output_file);
}

int main(int argc, char** argv) {

    arguments args;

    CLI::App app("Repetition-aware BWT construction");
    parse_app(app, args);
    CLI11_PARSE(app, argc, argv);

    LOG_INFO("Input file:       "+args.input_file);
    if(args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
    args.output_file = std::filesystem::path(args.output_file).replace_extension(".rl_bwt");

    if(args.alph_bytes==1){
        build_bwt_int<uint8_t>(args);
    }else if(args.alph_bytes==2){
        build_bwt_int<uint16_t>(args);
    }else if(args.alph_bytes==4){
        build_bwt_int<uint32_t>(args);
    } else if(args.alph_bytes==8){
        build_bwt_int<uint64_t>(args);
    }
    return 0;
}
