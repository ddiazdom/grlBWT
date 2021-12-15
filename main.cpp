#include <thread>

#include "external/CLI11.hpp"
#include "utils.hpp"
#include "grammar_build.hpp"

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads{};
    size_t b_buff=16;
    uint8_t comp_lvl=1;
    float hbuff_frac=0.5;
    bool ver=false;
    bool keep=false;
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

    CLI::App *gram = app.add_subcommand("gram", "Create a locally consistent grammar");
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
                      "Temporal folder (def. /tmp/lc_gram.xxxx)")->
            check(CLI::ExistingDirectory)->default_val("/tmp");
    gram->add_option("-L,--level-compression", args.comp_lvl, "Level of compression")->check(CLI::Range(1,2));

    CLI::App *dc = app.add_subcommand("decomp", "Decompress a locally consistent grammar to a file");
    dc->add_option("GRAM",
                     args.input_file,
                     "Input locally consistent grammar")->check(CLI::ExistingFile)->required();
    dc->add_option("-T,--tmp",
                     args.tmp_dir,
                     "Temporal folder (def. /tmp/lc_gram.xxxx)")->
                     check(CLI::ExistingDirectory)->default_val("/tmp");
    dc->add_option("-o,--output-file",
                     args.output_file,
                     "Output file")->type_name("");
    dc->add_option("-t,--threads",
                     args.n_threads,
                     "Number of threads")->default_val(1);
    dc->add_flag("-k,--keep",
                 args.keep,
                 "Keep the input grammar");
    dc->add_option("-B,--file-buffer",
                 args.b_buff,
                 "Size in MiB for the file buffer (def. 16 MiB)")->default_val(16);

    app.require_subcommand(1);
    app.footer("Report bugs to <diego.diaz@helsinki.fi>");
}

template<class vector_t>
void decompress_text(arguments& args){
    std::cout << "Decompressing the locally consistent grammar" << std::endl;
    std::string tmp_folder = create_temp_folder(args.tmp_dir, "lc_gram");

    if(args.output_file.empty()){
        args.output_file = std::filesystem::path(args.input_file).filename();
        args.output_file.resize(args.output_file.size()-5); //remove the ".gram" suffix
    }

    grammar<vector_t> gram;
    sdsl::load_from_file(gram, args.input_file);

    std::cout<<"Grammar size:                 "<<gram.gram_size()<<std::endl;
    std::cout<<"Number of terminals symbols:  "<<gram.ter()<<std::endl;
    std::cout<<"Number of nonterminals rules: "<<gram.nter()<<std::endl;
    std::cout<<"Number of run-length rules:   "<<gram.rl_rules()<<std::endl;
    std::cout<<"Number of SP rules:           "<<gram.sp_rules()<<std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    gram.se_decomp_str(0, gram.strings()-1,
                       args.output_file, tmp_folder,
                       args.n_threads, args.b_buff*1024*1024);
    auto end = std::chrono::high_resolution_clock::now();
    report_time(start, end);
}

int main(int argc, char** argv) {

    arguments args;

    CLI::App app("Grammar-based compression");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(app.got_subcommand("gram")) {

        std::cout << "Input file:        "<<args.input_file<<std::endl;
        std::cout << "Compression level: "<<(int)args.comp_lvl<<std::endl;
        std::cout << "Computing the grammar: "<<std::endl;
        std::string tmp_folder = create_temp_folder(args.tmp_dir, "lc_gram");

        if(args.output_file.empty()){
            args.output_file = std::filesystem::path(args.input_file).filename();
            args.output_file += ".gram";
        }

        build_gram(args.input_file, args.output_file, args.comp_lvl, tmp_folder, args.n_threads, args.hbuff_frac);

    }else if(app.got_subcommand("decomp")){

        std::ifstream ifs(args.input_file, std::ios::binary);
        uint8_t comp_level;
        ifs.read((char *)&comp_level, 1);
        ifs.close();

        if(comp_level==1){
            decompress_text<sdsl::int_vector<>>(args);
        }else if(comp_level==2){
            decompress_text<huff_vector<>>(args);
        }else{
            exit(1);
        }
    }
    return 0;
}
