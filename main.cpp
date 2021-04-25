#include <chrono>
#include <thread>

#include "third-party/CLI11.hpp"
#include "lpg/lpg.hpp"

struct arguments{
    std::string input_file;
    std::string output_file;
    std::string tmp_dir;
    size_t n_threads{};
    bool comp=false;
    bool decomp=false;
    float hbuff_frac=0.5;

    bool keep=false;
    bool verbose=false;
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
    int max_thread = (int)std::thread::hardware_concurrency();

    app.add_option("FILE",
                   args.input_file,
                   "Input file")->
                     check(CLI::ExistingFile)->required();

    auto *comp = app.add_flag("-c,--compress",
                 args.comp, "Compress the input file");
    auto *decom = app.add_flag("-d,--decompress",
                               args.decomp, "Decompress the input file")->excludes(comp);

    app.add_option("-o,--output-file",
                   args.output_file,
                   "Output file")->type_name("");

    app.add_flag("-k,--keep", args.keep, "Keep input file");

    app.add_option("-t,--threads",
                   args.n_threads,
                   "Maximum number of threads (def. [1-n_cores]=1)")->
                   check(CLI::Range(1, max_thread))->
                   default_val(1);

    app.add_option("-f,--hashing-buffer",
                   args.hbuff_frac,
                   "Hashing step will use at most INPUT_SIZE*f bytes. O means no limit (def. 0.5)")->
                   check(CLI::Range(0.0,1.0))->
                   default_val(0.5)->
                   excludes(decom);

    app.add_option("-T,--tmp",
                   args.tmp_dir,
                   "Temporal folder (def. /tmp/lpg.xxxx)")->
                     check(CLI::ExistingDirectory)->
                     default_val("/tmp");

    app.add_flag("-v,--verbose",
                 args.verbose, "Print extra information");
    app.add_flag("-V,--version",
                 args.ver, "Print the software version and exit");

    app.footer("By default, lpg will compress FILE if -c,-d or -b are not set\n\nReport bugs to <diediaz@dcc.uchile.cl>");
}

int main(int argc, char** argv) {

    arguments args;

    CLI::App app("A string compressor based on LMS induction");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(!args.comp && !args.decomp) args.comp = true;

    if(args.comp) {
        lpg<huff_vector<>> g(args.input_file, args.tmp_dir, args.n_threads, args.hbuff_frac);

        if(args.output_file.empty()){
            args.output_file = args.input_file;
        }

        args.output_file = args.output_file+".lg";

        sdsl::store_to_file(g, args.output_file);

        if(!args.keep){
            if(remove(args.input_file.c_str())){
                std::cout<<"There was an error trying to remove input file"<<std::endl;
            }
        }
    }else if(args.decomp){

        //TODO I have to include some mechanism to check if the file is ok
        std::cout<<"Loading the grammar"<<std::endl;
        lpg<huff_vector<>> g;
        sdsl::load_from_file(g, args.input_file);
        std::cout<<"  Terminals:                "<<(size_t)g.sigma<<std::endl;
        std::cout<<"  Nonterminals:             "<<g.tree.int_nodes()<<std::endl;
        std::cout<<"  Size of the comp. string: "<<g.tree.n_children(lpg<huff_vector<>>::root)<<std::endl;
        std::cout<<"  Grammar size:             "<<g.tree.nodes()<<std::endl;

        if(args.output_file.empty()){
            args.output_file = args.input_file.substr(0, args.input_file.size()-3);
        }

        g.decompress_text(args.tmp_dir, args.output_file, args.n_threads);

        if(!args.keep){
            if(remove(args.input_file.c_str())){
                std::cout<<"There was an error trying to remove input file"<<std::endl;
            }
        }
    }
    GET_MEM_PEAK();
    return 0;
}