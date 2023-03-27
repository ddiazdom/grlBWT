#include <thread>

#include "CLI11.hpp"
#include "utils.h"
#include "grl_bwt.hpp"
#include "fastx_handler.h"

struct arguments{
    std::string input_file;
    std::string output_file;

    std::string tmp_dir;
    size_t n_threads{};
    uint8_t b_f_r=1;
    float hbuff_frac=0.5;
    bool ver=false;
    bool rev_comp=false;
    bool opt_bwt=false;
    uint8_t alph_bytes=1;
    size_t sep_symbol=0;
    std::string version= "v1.0.1 alpha";
};

struct CellWidthValidator : public CLI::Validator {
    CellWidthValidator() {
        name_ = "ValidCellWidth";
        func_ = [](const std::string &str) {
            bool valid = (str=="1") || (str=="2") || (str=="4") || (str=="8");
            if(!valid)
                return std::string(str+" is not a valid number of bytes for an native integer type");
            else
                return std::string();
        };
    }
};
const static CellWidthValidator ValidCellWidth;

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
                      //"Input file in one-string-per-line or FASTA/Q format (automatically detected)"
                      "Input file in one-string-per-line format")->check(CLI::ExistingFile)->required();
    app.add_option("-o,--output-file",
                      args.output_file,
                      "Output file")->type_name("");
    app.add_option("-a,--alphabet", args.alph_bytes, "Number of bytes for the alphabet (def. 1)")->
            check(CLI::Range(1, 8))->default_val(1)->check(ValidCellWidth);
    app.add_option("-t,--threads",
                      args.n_threads,
                      "Maximum number of working threads")->default_val(1);
    app.add_option("-f,--hbuff",
                      args.hbuff_frac,
                      "Hashing step will use at most INPUT_SIZE*f bytes. O means no limit (def. 0.5)")->
            check(CLI::Range(0.0,1.0))->default_val(0.15);
    app.add_option("-b,--run-len-bytes",
                   args.b_f_r,
                   "Max. number of bytes to encode the run lengths in the recursive BWTs (def. 1)")->
            check(CLI::Range(0,5))->default_val(1);
    app.add_option("-T,--tmp",
                      args.tmp_dir,
                      "Temporary folder (def. /tmp/grl.bwt.xxxx)")->
            check(CLI::ExistingDirectory)->default_val("/tmp");
    app.add_flag("-v,--version",
                 args.ver, "Print the software version and exit");
    //app.add_flag("-R,--rev-comp", args.rev_comp, "Also consider the DNA reverse complements of the strings in TEXT");
    //app.add_flag("-m,--min-bwt", args.opt_bwt, "Produce the optimal BCR BWT instead of the regular one");

    app.footer("Report bugs to <diego.diaz@helsinki.fi>");
}

template<class sym_type>
void run_int(std::string input_collection, arguments& args){

    tmp_workspace tmp_ws(args.tmp_dir, true, "grl.bwt");
    std::cout<< "Temporary folder: "<<tmp_ws.folder()<<std::endl;

    if(args.opt_bwt){
        //grl_bwt_algo<true>(input_collection, args.output_file, tmp_ws, args.n_threads, str_coll, args.hbuff_frac, args.b_f_r);
        //TODO this option is broken
        std::cout<<"BWT type:         BCR optimal"<<std::endl;
        std::cout<<"This option is broken"<<std::endl;
        exit(0);
    }else{
        std::cout<<"BWT type:         BCR exact"<<std::endl;
        grl_bwt_algo<sym_type, false>(input_collection, args.output_file, tmp_ws, args.n_threads, args.hbuff_frac, args.b_f_r);
    }
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

    std::cout << "Input file:       "<<args.input_file<<std::endl;
    if(args.output_file.empty()) args.output_file = std::filesystem::path(args.input_file).filename();
    args.output_file = std::filesystem::path(args.output_file).replace_extension(".rl_bwt");

    std::string input_collection = args.input_file;
    str_collection str_coll;

    //TODO I turned FASTA/Q support off for the moment
    /*if(is_fastx(args.input_file)){
        //input_collection = tmp_ws.get_file("plain_input");
        //std::cout<<"The input is in FASTA/Q format"<<std::endl;
        //std::cout<<"Creating a temporary file in plain format: "<<input_collection<<std::endl;
        //str_coll = fastx2plain_format(args.input_file, input_collection, args.rev_comp, '\n');
        std::cout<<"this option is not implemented yet"<<std::endl;
        exit(0);
    }else if(args.rev_comp) {
        //input_collection = tmp_ws.get_file("plain_input");
        //std::cout<<"The input is in plain format, but the BCR BWT computation requires the DNA reverse complements"<<std::endl;
        //std::cout<<"Creating a temporary file in plain format: "<<input_collection<<std::endl;
        //std::filesystem::copy(args.input_file, input_collection);
        //TODO get reverse complement in plain format
        std::cout<<"this option is not implemented yet"<<std::endl;
        exit(0);
    }else {
        str_coll = collection_stats<uint8_t>(input_collection);
    }*/

    if(args.alph_bytes>1){
        std::cout<<"Alphabet type:    integer"<<std::endl;
    }else{
        std::cout<<"Alphabet type:    byte"<<std::endl;
    }

    if(args.alph_bytes==1){
        run_int<uint8_t>(input_collection, args);
    }else if(args.alph_bytes==2){
        run_int<uint16_t>(input_collection, args);
    }else if(args.alph_bytes==4){
        run_int<uint32_t>(input_collection, args);
    } else if(args.alph_bytes==8){
        run_int<uint64_t>(input_collection, args);
    }
    return 0;
}
