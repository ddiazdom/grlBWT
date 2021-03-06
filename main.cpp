#include <thread>

#include "external/CLI11.hpp"
#include "utils.h"
#include "grl_bwt.hpp"
#include "fastx_handler.h"

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads{};
    float hbuff_frac=0.5;
    bool ver=false;
    bool rev_comp=false;
    std::string version= "0.0.1";
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
                      "Input file in one-string-per-line or FASTA/Q format"
                      " (automatically detected)")->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output-file",
                      args.output_file,
                      "Output file")->type_name("");
    app.add_option("-t,--threads",
                      args.n_threads,
                      "Maximum number of working threads")->default_val(1);
    app.add_option("-f,--hbuff",
                      args.hbuff_frac,
                      "Hashing step will use at most INPUT_SIZE*f bytes. O means no limit (def. 0.5)")->
            check(CLI::Range(0.0,1.0))->default_val(0.15);
    app.add_option("-T,--tmp",
                      args.tmp_dir,
                      "Temporary folder (def. /tmp/grl.bwt.xxxx)")->
            check(CLI::ExistingDirectory)->default_val("/tmp");
    app.add_flag("-v,--version",
                 args.ver, "Print the software version and exit");
    app.add_flag("-R,--rev-comp", args.rev_comp, "Also consider the DNA reverse complements of the strings in TEXT");
    app.footer("Report bugs to <diego.diaz@helsinki.fi>");
}

int main(int argc, char** argv) {

    arguments args;

    CLI::App app("Repetition-aware BWT construction");
    parse_app(app, args);

    CLI11_PARSE(app, argc, argv);
    if(args.ver){
        std::cout<<args.version<<std::endl;
        exit(0);
    }

    tmp_workspace tmp_ws(args.tmp_dir, true, "grl.bwt");

    std::cout << "Input file:       "<<args.input_file<<std::endl;
    std::cout << "Temporary folder: "<<tmp_ws.folder()<<std::endl;

    if(args.output_file.empty()) args.output_file = args.input_file;
    args.output_file = std::filesystem::path(args.output_file).replace_extension(".rl_bwt");

    std::string input_collection = args.input_file;
    str_collection str_coll;

    if(is_fastx(args.input_file)){
        input_collection = tmp_ws.get_file("plain_input");
        str_coll = fastx2plain(args.input_file, input_collection);
        if(args.rev_comp) get_rev_comp(input_collection, str_coll);
    }else {
        input_collection = tmp_ws.get_file("plain_input");
        std::filesystem::copy(args.input_file, input_collection);
        str_coll = collection_stats(input_collection);
        if(args.rev_comp) get_rev_comp(input_collection, str_coll);
    }

    std::cout<<"Alphabet : "<<str_coll.alphabet.size()<<std::endl;
    std::cout<<std::setw(10)<<"Symbol"<<std::setw(12)<<"Frequency"<<std::endl;
    for(size_t i=0;i<str_coll.alphabet.size();i++){
        if(str_coll.alphabet[i]==10){
            std::cout<<std::setw(10)<<"'\\n'";
        }else{
            std::cout<<std::setw(8)<<"'"<<str_coll.alphabet[i]<<"'";
        }
        std::cout<<std::setw(10)<<str_coll.sym_freqs[i]<<std::endl;
    }
    std::cout<<"Number of strings :    "<<str_coll.n_strings<<std::endl;
    std::cout<<"Number of characters : "<<str_coll.n_char<<std::endl;

    grl_bwt_algo(input_collection, args.output_file, tmp_ws, args.n_threads,
                 str_coll, args.hbuff_frac);
    return 0;
}
