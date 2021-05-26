#include <chrono>
#include <thread>

#include "third-party/CLI11.hpp"
#include "lpg/lpg_index.hpp"

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

    CLI::App *index = app.add_subcommand("index", "Create an LPG self-index");
    CLI::App *search = app.add_subcommand("search", "Search for a pattern in the index");
    app.set_help_all_flag("--help-all", "Expand all help");

    app.add_flag("-V,--version",
                 args.ver, "Print the software version and exit");

    index->add_option("TEXT",
                      args.input_file,
                      "Input text file")->check(CLI::ExistingFile)->required();
    index->add_option("-o,--output-file",
                      args.output_file,
                      "Output file")->type_name("");
    index->add_option("-t,--threads",
                      args.n_threads,
                      "Maximum number of threads")->default_val(1);
    index->add_option("-f,--hbuff",
                      args.hbuff_frac,
                      "Hashing step will use at most INPUT_SIZE*f bytes. O means no limit (def. 0.5)")->
            check(CLI::Range(0.0,1.0))->
            default_val(0.5);
    index->add_option("-T,--tmp",
                      args.tmp_dir,
                      "Temporal folder (def. /tmp/lpg_index.xxxx)")->
            check(CLI::ExistingDirectory)->default_val("/tmp");


    search->add_option("INDEX",
                       args.input_file,
                       "Input LPG index file")->check(CLI::ExistingFile)->required();

    search->add_option("-p,--patterns",
                                     args.patterns,
                                     "Pattern to search for in the index")->type_name("")->required(true);

    search->add_option("-F,--pattern-list",
                                  args.patter_list_file,
                    "File with a pattern list")->type_name("");
    search->add_option("-t,--threads",
                       args.n_threads,
                       "Maximum number of threads")->default_val(1);
    search->add_option("-o,--output-file",
                   args.output_file,
                   "Output file")->type_name("");
    app.require_subcommand(1);

    app.footer("By default, lpg_index will compress FILE if -c,-d or -b are not set\n\nReport bugs to <diediaz@dcc.uchile.cl>");
}

int main(int argc, char** argv) {

    arguments args;

    CLI::App app("A grammar-based self-index");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);

    if(app.got_subcommand("index")) {
        lpg_index g(args.input_file, args.tmp_dir, args.n_threads, args.hbuff_frac);

        if(args.output_file.empty()){
            args.output_file = args.input_file;
        }
        args.output_file = args.output_file+".idx";

        std::cout<<"Saving the self-index to file "<<args.output_file<<std::endl;
        sdsl::store_to_file(g, args.output_file);

    }else if(app.got_subcommand("search")){
        std::cout<<"Searching for patterns in the self-index"<<std::endl;
        lpg_index g;
        sdsl::load_from_file(g, args.input_file);

        std::cout<<"Index size:"<<sdsl::size_in_bytes(g)<<std::endl;
        std::cout<<"Index name:"<<args.input_file<<std::endl;
        std::set<std::string> patterns_set;
        std::vector<std::string> patterns;
        if (!args.patter_list_file.empty()){
            std::fstream in(args.patter_list_file,std::ios::in);
            if(in.good()){
                std::string ss;
                while(in >> ss){
                    patterns_set.insert(ss);
                    patterns.push_back(ss);
                }
#ifdef DEBUG_INFO
                std::cout<<"Patterns to search["<<patterns_set.size()<<"]"<<std::endl;
#endif
#ifdef DEBUG_PRINT
                for (const auto &s : patterns_set) std::cout<<s<<std::endl;
#endif

            }
        }
#ifdef CHECK_OCC
        std::string file; file.resize(args.input_file.size() - 3);
        std::copy(args.input_file.begin(),args.input_file.end()-3,file.begin());
        std::cout<<"file:"<<file<<std::endl;
#endif
//        std::cout<<"descomprimiendo gramatica"<<std::endl;
//        g.uncompress_grammar(args.input_file + ".ucmp");
//        if(!utils::compareFiles(args.input_file + ".ucmp",file))
//            std::cout<<"ERROR GRAMMAR DOES NOT REPRESENT THE TEXT\n";

        /*for(auto const& pattern : patterns){
            std::cout<<pattern<<std::endl;
            auto cuts = g.compute_pattern_cut(pattern);
            for(unsigned long j : cuts.first){
                std::cout<<j<<std::endl;
            }
        }*/
        std::cout<<"Searching for the patterns "<<std::endl;
        g.search(patterns
#ifdef CHECK_OCC
                 ,file
#endif
                 );
        //g.search(args.patterns);
        //g.search(args.patter_list_file);
    }
    return 0;
}