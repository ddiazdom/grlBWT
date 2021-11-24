#include <thread>

#include "external/CLI11.hpp"
#include "utils.hpp"
#include "grammar_build.hpp"

struct arguments{
    std::string input_file;
    std::string output_file;
    std::string patter_list_file;
    std::vector<std::string> patterns;

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

    CLI::App *gram = app.add_subcommand("gram", "Create the LPG grammar");
    app.set_help_all_flag("--help-all", "Expand all help");

    app.add_flag("-V,--version",
                 args.ver, "Print the software version and exit");

    gram->add_option("TEXT",
                      args.input_file,
                      "Input text file")->check(CLI::ExistingFile)->required();
    gram->add_option("-o,--output-file",
                      args.output_file,
                      "Output file")->type_name("");
    gram->add_option("-t,--threads",
                      args.n_threads,
                      "Maximum number of threads")->default_val(1);
    gram->add_option("-f,--hbuff",
                      args.hbuff_frac,
                      "Hashing step will use at most INPUT_SIZE*f bytes. O means no limit (def. 0.5)")->
            check(CLI::Range(0.0,1.0))->default_val(0.5);
    gram->add_option("-T,--tmp",
                      args.tmp_dir,
                      "Temporal folder (def. /tmp/lpg_gram.xxxx)")->
            check(CLI::ExistingDirectory)->default_val("/tmp");
    app.require_subcommand(1);
    app.footer("Report bugs to <diego.diaz@helsinki.fi>");
}

int main(int argc, char** argv) {

	arguments args;

    CLI::App app("Grammar-based compression");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(app.got_subcommand("gram")) {

        std::cout << "Computing the grammar for the self-index" << std::endl;
        std::string tmp_folder = create_temp_folder(args.tmp_dir, "lpg_gram");

        if(args.output_file.empty()){
            args.output_file = std::filesystem::path(args.input_file).filename();
            args.output_file += ".gram";
        }
        build_gram(args.input_file, args.output_file, tmp_folder, args.n_threads, args.hbuff_frac);
    }else{

    }
    return 0;
}
