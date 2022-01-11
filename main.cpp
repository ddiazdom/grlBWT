#include <thread>

#include "external/CLI11.hpp"
#include "utils.hpp"
#include "gbwt.hpp"

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads{};
    float hbuff_frac=0.5;
    bool ver=false;
};

class MyFormatter : public CLI::Formatter {
public:
    MyFormatter() : Formatter() {}
    std::string make_option_opts(const CLI::Option *) const override { return ""; }
};

static void parse_app(CLI::App& app, struct arguments& args){
    
	auto fmt = std::make_shared<MyFormatter>();

    fmt->column_width(23);
    app.formatter(fmt);

    app.add_option("TEXT",
                      args.input_file,
                      "Input text file")->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output-file",
                      args.output_file,
                      "Output file")->type_name("");
    app.add_option("-t,--threads",
                      args.n_threads,
                      "Maximum number of threads")->default_val(1);
    app.add_option("-f,--hbuff",
                      args.hbuff_frac,
                      "Hashing step will use at most INPUT_SIZE*f bytes. O means no limit (def. 0.5)")->
            check(CLI::Range(0.0,1.0))->default_val(0.15);
    app.add_option("-T,--tmp",
                      args.tmp_dir,
                      "Temporal folder (def. /tmp/gbwt_tmp.xxxx)")->
            check(CLI::ExistingDirectory)->default_val("/tmp");
    app.add_flag("-V,--version",
                 args.ver, "Print the software version and exit");

    app.footer("Report bugs to <diego.diaz@helsinki.fi>");
}

int main(int argc, char** argv) {

    arguments args;

    CLI::App app("Grammar-based construction of the BWT");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    std::cout << "Input file:        "<<args.input_file<<std::endl;
    std::cout << "Temporal folder:   "<<args.tmp_dir<<std::endl;
    std::string tmp_folder = create_temp_folder(args.tmp_dir, "gbwt_tmp");

    if(args.output_file.empty()){
        args.output_file = std::filesystem::path(args.input_file).filename();
        args.output_file += ".gram";
    }
    g_bwt_algo(args.input_file, args.output_file, tmp_folder, args.n_threads, args.hbuff_frac);

    return 0;
}
