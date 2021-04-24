//
// Created by diediaz on 11-12-19.
//

#ifndef LMS_GRAMMAR_REP_HPP
#define LMS_GRAMMAR_REP_HPP

#include "cdt/louds_tree.hpp"
#include "cdt/se_int_queue.h"

#include "lpg_build.hpp"
#include "lpg_iterator.hpp"
#include "lpg/suffpair_algo.hpp"

#include <pthread.h>

template<class label_arr_type>
class lpg{

private:
    typedef typename lpg_build::plain_grammar_t         plain_grammar_t;
    louds_tree                                          m_tree;
    label_arr_type                                      m_l_labels;
    sdsl::int_vector<8>                                 symbols_map;
    uint8_t                                             m_sigma{};
    size_t                                              m_max_degree{}; //only for decompression purposes
    friend struct grammar_iterator<lpg>;

    struct decomp_data{
        size_t start;
        size_t end;
        std::string tmp_file;
        lpg &grammar;
        decomp_data(size_t start_, size_t end_, std::string& tmp_file_, lpg& grammar_):
        start(start_), end(end_), tmp_file(tmp_file_), grammar(grammar_){};
    };


    static void build_grammar(std::string &i_file, std::string &g_file, uint8_t sep_symbol,
                              sdsl::cache_config &config, size_t n_threads, size_t hbuff_size){

        std::cout<<"Computing the LPG grammar"<<std::endl;
        auto start = std::chrono::high_resolution_clock::now();
        lpg_build::compute_LPG(i_file, g_file, n_threads, config, hbuff_size, sep_symbol);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;
        std::cout<<"  Elap. time (secs): "<<elapsed.count()<<std::endl;

        std::cout<<"Suffix-pairing nonterminals"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        suffpair(g_file, config, n_threads, hbuff_size);
        end = std::chrono::high_resolution_clock::now();
        elapsed = end - start;
        std::cout<<"  Elap. time (secs): "<<elapsed.count()<<std::endl;
    }

    void build_grammar_tree(plain_grammar_t &p_gram, sdsl::cache_config& config){

        sdsl::int_vector_buffer<1> tree_shape(sdsl::cache_file_name("gt_shape", config),
                                              std::ios::out);
        //TODO check the width of the g_pointers
        sdsl::int_vector_buffer<64> l_labs(sdsl::cache_file_name("l_labels", config),
                                              std::ios::out);
        m_max_degree=0;

        size_t tree_pos=0;
        tree_shape[tree_pos++] = true;
        tree_shape[tree_pos++] = false;

        {
            std::cout <<"  Visiting the grammar rules in level-order" << std::endl;
            bool int_node;
            size_t int_rank=1;
            sdsl::int_vector<> g_locus(p_gram.r, 0, sdsl::bits::hi(p_gram.r) + 1);
            std::string queue_file=sdsl::cache_file_name("queue_file", config);
            se_int_queue<size_t> g_queue(queue_file, BUFFER_SIZE);

            bv_t lmsg_as_sp_bv;
            bool lmsg_as_sp;
            sdsl::load_from_file(lmsg_as_sp_bv, p_gram.lms_as_sp_file);
            g_queue.push_back((p_gram.r - 1) << 1UL);

            //load the succinct (uncompressed) grammar
            sdsl::int_vector<> r;
            sdsl::bit_vector r_lim;
            sdsl::load_from_file(r, p_gram.rules_file);
            sdsl::load_from_file(r_lim, p_gram.rules_lim_file);
            sdsl::bit_vector::select_1_type r_lim_ss;
            sdsl::util::init_support(r_lim_ss, &r_lim);


            size_t start, gpos, len, tail;
            while (!g_queue.empty()) {

                auto tmp_sym = g_queue.front();
                g_queue.pop();

                lmsg_as_sp = (tmp_sym & 1UL);
                tmp_sym>>=1UL;

                if (tmp_sym < m_sigma) {//terminal
                    l_labs.push_back(tmp_sym);
                } else {
                    if (g_locus[tmp_sym] == 0 && !lmsg_as_sp) {//first visit of a nonterminal rule
                        int_node = true;
                        g_locus[tmp_sym] = int_rank;
                    } else {//it is not the first visit of a nonterminal rule
                        if(g_locus[tmp_sym] == 0){
                            //the locus is still undetermined
                            l_labs.push_back(p_gram.r + tmp_sym);
                        }else{
                            l_labs.push_back(m_sigma + g_locus[tmp_sym]);
                        }
                        int_node = false;
                    }

                    if (int_node) {
                        start = r_lim_ss(tmp_sym)+1;

                        gpos = start;
                        while(!r_lim[gpos]){
                            g_queue.push_back(r[gpos]<<1UL);
                            tree_shape[tree_pos++] = true;
                            gpos++;
                        }
                        len = gpos-start+1;

                        //check if the last symbol is an LMSg nonterminal that occurs in a SP context
                        tail = r[gpos]<<1UL;
                        if(lmsg_as_sp_bv[tmp_sym]) tail |=1UL;
                        g_queue.push_back(tail);
                        tree_shape[tree_pos++] = true;
                        //

                        if (tmp_sym != p_gram.r && len > max_degree) {
                            m_max_degree = len;
                        }
                        int_rank++;
                    }
                }
                tree_shape[tree_pos++] = false;
            }

            //update missing loci
            for(auto && g_pointer : l_labs){
                if(g_pointer > p_gram.r){
                    g_pointer = m_sigma + g_locus[g_pointer - p_gram.r];
                }
            }
            tree_shape.close();
            sdsl::util::clear(g_locus);
            sdsl::util::clear(r);
            sdsl::util::clear(r_lim);
            sdsl::util::clear(r_lim_ss);
            g_queue.close(true);
        }

        std::cout<<"  Compressing leaf labels"<<std::endl;
        m_l_labels = label_arr_type(l_labs, config);
        l_labs.close();

        std::cout<<"  Building LOUDS rep. of the tree shape"<<std::endl;
        sdsl::bit_vector bv;
        sdsl::load_from_file(bv, sdsl::cache_file_name("gt_shape", config));
        m_tree = louds_tree(std::move(bv));

        sdsl::register_cache_file("gt_shape", config);
        sdsl::register_cache_file("l_labels", config);
    };

    //static void check_uncomp_grammar(plain_grammar_t &r_data, std::string uncom_seq);

    void check_grammar(std::string &uncom_seq){
        std::cout<<"Checking the grammar produces the exact input string"<<std::endl;

        std::vector<uint32_t > tmp_decomp;
        std::vector<uint32_t > prev_dc_step;

        std::ifstream uc_vec(uncom_seq, std::ifstream::binary);

        uc_vec.seekg(0, std::ifstream::end);
        size_t len = uc_vec.tellg();
        uc_vec.seekg(0, std::ifstream::beg);

        char *buffer = new char[len];
        uc_vec.read(buffer, len);
        uc_vec.close();

        size_t n_children = m_tree.n_children(root);
        size_t pos=0, tmp_node;
        std::vector<uint8_t> buff;
        for(size_t i=1;i<=n_children;i++){
            tmp_node = m_tree.child(2, i);
            decompress_node(buff, tmp_node);
            for(unsigned char j : buff){
                assert(j==buffer[pos]);
                pos++;
            }
            buff.clear();
        }
        uc_vec.close();

        delete [] buffer;
        std::cout<<"Everything ok!"<<std::endl;
    }

    void swap(lpg&& other){
        std::swap(m_sigma, other.m_sigma);
        std::swap(m_max_degree, other.m_max_degree);
        std::swap(symbols_map, other.symbols_map);
        std::swap(m_tree, other.m_grammar_tree);
        std::swap(m_l_labels, other.m_l_labels);
    }

public:
    typedef size_t                        size_type;
    static constexpr size_t root          = 2;
    const louds_tree&       tree          = m_tree;
    const label_arr_type&   l_labels      = m_l_labels;
    const uint8_t&          sigma         = m_sigma;
    const size_t&           max_degree    = m_max_degree;

    lpg(std::string &input_file, std::string &tmp_folder,
        uint8_t sep_symbol, size_t n_threads, float hbuff_frac){

        std::cout<<"Input file: "<<input_file<<std::endl;

        //create a temporary folder
        std::string tmp_path = tmp_folder+"/lpg.XXXXXX";
        char temp[200]={0};
        tmp_path.copy(temp, tmp_path.size() + 1);
        temp[tmp_path.size()+1] = '\0';
        auto res = mkdtemp(temp);
        if(res==nullptr){
            std::cout<<"Error trying to create a temporal folder"<<std::endl;
        }
        std::cout<<"Temporal folder: "<<std::string(temp)<<std::endl;
        //

        sdsl::cache_config config(false, temp);

        std::string g_file = sdsl::cache_file_name("g_file", config);

        //build the uncompressed version of the grammar
        std::ifstream is (input_file, std::ifstream::binary);
        is.seekg (0, std::ifstream::end);
        auto n_chars = is.tellg();
        is.close();

        //maximum amount of RAM allowed to spend in parallel for the hashing step
        auto hbuff_size = size_t(std::ceil(float(n_chars)*hbuff_frac));

        std::cout<<"Input file contains "<<n_chars<<" symbols"<<std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        build_grammar(input_file, g_file, sep_symbol, config, n_threads, hbuff_size);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        //lpg_build::check_plain_grammar(g_file, input_file);

        //plain representation of the grammar
        plain_grammar_t plain_gram;
        plain_gram.load_from_file(g_file);

        std::cout <<"Resulting grammar stats: "<<std::endl;
        std::cout <<"  Terminals:                " << (size_t)plain_gram.sigma << std::endl;
        std::cout << "  Nonterminals:             " << plain_gram.r - plain_gram.sigma << std::endl;
        std::cout << "  Size of the comp. string: " << plain_gram.c << std::endl;
        std::cout << "  Grammar size:             " << plain_gram.g - plain_gram.sigma << std::endl;
        std::cout <<"  Elap. time (secs):        "<< elapsed.count()<<std::endl;

        m_sigma = plain_gram.sigma;
        symbols_map.resize(plain_gram.symbols_map.size());
        for(size_t i=0; i < plain_gram.symbols_map.size(); i++){
            symbols_map[i] = plain_gram.symbols_map[i];
        }

        std::cout<<"Building the grammar tree"<<std::endl;
        start = std::chrono::high_resolution_clock::now();
        build_grammar_tree(plain_gram, config);
        end = std::chrono::high_resolution_clock::now();

        elapsed = end - start;
        std::cout<<"  Grammar tree stats:"<<std::endl;
        std::cout <<"    Tot. nodes:         "<<m_tree.nodes() << std::endl;
        std::cout <<"    Int. nodes:         "<<m_tree.int_nodes() << std::endl;
        std::cout <<"    Leaves:             "<<m_tree.nodes() - m_tree.int_nodes() << std::endl;
        std::cout <<"    Space usage (MB):   "<<sdsl::size_in_mega_bytes(*this)<<std::endl;
        std::cout <<"    Compression ratio:  "<<double(n_chars)/sdsl::size_in_bytes(*this)<<std::endl;
        std::cout <<"  Elap. time (secs):  "<<elapsed.count()<<std::endl;
    }

    lpg(): m_sigma(0), m_max_degree(0){};

    lpg(const lpg &other){
        m_sigma = other.m_sigma;
        m_max_degree = other.m_max_degree;
        symbols_map = other.symbols_map;
        m_tree = other.m_grammar_tree;
        m_l_labels = other.m_l_labels;
    }

    lpg(lpg &&other) noexcept {
        swap(std::forward<lpg>(other));
    }

    lpg& operator=(lpg&& other) noexcept {
        swap(std::forward<lpg>(other));
        return *this;
    }

    grammar_iterator<lpg> g_begin() const {
        return grammar_iterator<lpg>{*this, root};
    }

    grammar_iterator<lpg> g_end() const {
        return grammar_iterator<lpg>{*this, m_tree.shape.size()};
    }

    //get an level-order iterator that start from the input node
    grammar_iterator<lpg> node(size_t node) const{
        assert(node>=root && node<m_tree.shape.size());
        return grammar_iterator<lpg>{*this, node};
    }

    void subtree_shape_int(size_t node_id, std::string& shape) const{
        if(m_tree.is_leaf(node_id)){
            size_t p_val = m_l_labels[m_tree.leaf_rank(node_id) - 1];
            if(p_val<=m_sigma){
                shape.push_back('(');
                shape.push_back(')');
                return;
            }else{
                node_id = m_tree.int_select(p_val - m_sigma);
            }
        }

        size_t n_children= m_tree.n_children(node_id);

        shape.push_back('(');
        for(size_t i=1;i<=n_children;i++){
            subtree_shape_int(m_tree.child(node_id, i), shape);
        }
        shape.push_back(')');
    }

    std::string parse_subtree_shape(size_t node_id) const{
        std::string shape;
        subtree_shape_int(node_id, shape);
        return shape;
    }

    std::string decompress_node(size_t node_id) const {
        std::string res;
        decompress_node<std::string>(res, node_id);
        return res;
    }

    template<class buff_t>
    void decompress_node(buff_t& buffer, size_t node_id) const {
        assert(node_id != 0);
        if(m_tree.is_leaf(node_id)){
            size_t lab = m_l_labels[m_tree.leaf_rank(node_id) - 1];
            if(lab<sigma){
                buffer.push_back(symbols_map[lab]);
                return;
            }else{
                node_id = m_tree.int_select(lab-m_sigma);
            }
        }

        size_t n_children = m_tree.n_children(node_id);
        for(size_t r=1;r<=n_children;r++){
            decompress_node(buffer, m_tree.child(node_id, r));
        }
    }

    static void * decompress_thread(void * data){

        size_t start = ((decomp_data*)data)->start;
        size_t end = ((decomp_data*)data)->end;
        lpg &grammar = ((decomp_data*)data)->grammar;

        size_t tmp_node;

        o_file_stream<uint8_t> of_s(((decomp_data*)data)->tmp_file, BUFFER_SIZE, std::ios::out);

        for(size_t i=start;i<=end;i++){
            tmp_node = grammar.m_tree.child(root, i);
            grammar.decompress_node<o_file_stream<uint8_t>>(of_s, tmp_node);
        }
        of_s.close();

        pthread_exit(nullptr);
    }

    void decompress_text(std::string& tmp_folder, std::string& o_file, size_t n_threads){
        std::cout<<"Decompressing the text"<<std::endl;

        if(tmp_folder.size()>1 && tmp_folder.back()=='/'){
            tmp_folder.pop_back();
        }
        //create a temporary folder
        std::string tmp_path = tmp_folder+"/lpg.XXXXXX";

        char *temp = (char *) malloc(tmp_path.size()+1);
        memcpy(temp, tmp_path.c_str(), tmp_path.size());
        temp[tmp_path.size()] = 0;

        if(mkdtemp(temp)== nullptr){
            std::cout<<"Error trying to set temporal folder : "<<errno<<std::endl;
            exit(1);
        }

        std::cout<<"Temporal folder: "<<std::string(temp)<<std::endl;
        //

        auto start_p = std::chrono::high_resolution_clock::now();

        size_t n_symbols = m_tree.n_children(2);
        size_t start, end;
        assert(n_symbols!=0 && n_threads>=1);

        size_t sym_per_thread = 1 + ((n_symbols - 1) / n_threads);

        std::vector<decomp_data> d_data;

        for(size_t i=0;i<n_threads;i++){

            start = (i * sym_per_thread) + 1;
            end = std::min<size_t>((i+1) * sym_per_thread, n_symbols);

            std::stringstream ss;
            ss<<temp<<"/chunk_"<<start<<"_"<<end<<".txt";
            std::string tmp_file = ss.str();

            d_data.emplace_back(start, end, tmp_file, *this);
        }

        std::vector<pthread_t> threads(n_threads);
        for(size_t i=0;i<n_threads;i++){

            int ret =  pthread_create(&threads[i],
                                      nullptr,
                                      &lpg::decompress_thread,
                                      (void*)&d_data[i]);
            if(ret != 0) {
                printf("Error: pthread_create() failed\n");
                exit(EXIT_FAILURE);
            }
        }

        for(size_t i=0;i<n_threads;i++) {
            pthread_join(threads[i], nullptr);
        }

        //concatenate elements into one file!!
        std::ofstream o_f(o_file, std::ios::out | std::ofstream::binary);
        assert(o_f.good());

        size_t len, rem, to_read;
        char *buffer = new char[BUFFER_SIZE];

        for(auto & chunk : d_data){

            std::ifstream i_f(chunk.tmp_file, std::ifstream::binary);

            i_f.seekg (0, std::ifstream::end);
            len = i_f.tellg();
            i_f.seekg (0, std::ifstream::beg);

            rem = len;
            to_read = std::min<size_t>(BUFFER_SIZE, len);

            while(true){

                i_f.read(buffer, to_read);
                assert(i_f.good());

                o_f.write(buffer, to_read);
                assert(o_f.good());
                //TODO if file not ok, then remove it

                rem -= i_f.gcount();
                to_read = std::min<size_t>(BUFFER_SIZE, rem);
                if(to_read == 0) break;
            }
            i_f.close();

            if(remove(chunk.tmp_file.c_str())){
                std::cout<<"Error trying to remove temporal file"<<std::endl;
                std::cout<<"Aborting"<<std::endl;
                exit(1);
            }
        }

        delete[] buffer;
        o_f.close();

        //delete temporary files
        std::cout<<"Deleting temporal data"<<std::endl;
        if(rmdir(temp)!=0){
            std::cout<<"temporary files couldn't be deleted properly"<<std::endl;
        }
        std::cout<<"Decompression completed"<<std::endl;
        auto end_p = std::chrono::high_resolution_clock::now();

        std::chrono::duration<double> elapsed = end_p - start_p;
        std::cout<<"Output file: "<<o_file<<std::endl;
        std::cout<<"Elap. time: "<<elapsed.count()<<" seconds"<<std::endl;
        free(temp);
    };

    //node label
    inline size_t label(size_t node_id) const {
        if(m_tree.is_leaf(node_id)){
            return m_l_labels[m_tree.leaf_rank(node_id) - 1];
        }else{
            return m_tree.int_rank(node_id)+m_sigma;
        }
    }

    //rth child in the parse tree
    inline size_t pt_child(size_t node_id, size_t r) const{
        assert(node_id>=2);
        if(m_tree.is_leaf(node_id)){
            size_t label = m_l_labels[m_tree.leaf_rank(node_id) - 1];
            assert(label>m_sigma);
            node_id = m_tree.int_select(label-m_sigma);
        }
        return m_tree.child(node_id, r);
    }

    //rightmost child in the parse tree
    inline size_t pt_rm_child(size_t node_id) const{
        assert(node_id>=2);
        if(m_tree.is_leaf(node_id)){
            size_t label = m_l_labels[m_tree.leaf_rank(node_id) - 1];
            assert(label>m_sigma);
            node_id = m_tree.int_select(label-m_sigma);
        }
        return m_tree.child(node_id, m_tree.n_children(node_id));
    }

    //leftmost child in the parse tree
    inline size_t pt_lm_child(size_t node_id) const{
        assert(node_id>=2);
        if(m_tree.is_leaf(node_id)){
            size_t label = m_l_labels[m_tree.leaf_rank(node_id) - 1];
            assert(label>m_sigma);
            node_id = m_tree.int_select(label-m_sigma);
        }
        return m_tree.child(node_id, 1);
    }

    void load(std::istream &in){
        m_tree.load(in);
        m_l_labels.load(in);
        symbols_map.load(in);
        sdsl::read_member(m_max_degree, in);
        sdsl::read_member(m_sigma, in);
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const{
        sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += m_tree.serialize(out, child, "lg_tree");
        written_bytes += m_l_labels.serialize(out, child, "lg_pointers");
        written_bytes += symbols_map.serialize(out, child, "symbols_map");
        written_bytes += sdsl::write_member(m_max_degree, out, child, "max_degree");
        written_bytes += sdsl::write_member(m_sigma, out, child, "sigma");
        return written_bytes;
    }
};

template <>
void lpg<sdsl::int_vector<>>::load(std::istream &in){

    m_tree.load(in);

    sdsl::cache_config config;
    std::string lab_file = sdsl::cache_file_name("labels", config);
    size_t lab_width = sdsl::bits::hi(m_tree.int_nodes())+1;
    size_t lab_space = m_tree.leaves()*lab_width;
    size_t buff_size = std::min<size_t>(BUFFER_SIZE, INT_CEIL(lab_space,8));
    sdsl::int_vector_buffer<> tmp(lab_file, std::ios::out, buff_size, lab_width);

    {
        huff_vector<> lab_data;
        lab_data.load(in);
        auto it = lab_data.begin();
        auto it_end = lab_data.end();
        while (it != it_end) {
            tmp.push_back(*it);
            ++it;
        }
        tmp.close();
    }

    sdsl::load_from_file(m_l_labels, lab_file);

    if(remove(lab_file.c_str())){
        std::cout<<"Error trying to remove file "<<lab_file<<std::endl;
    }

    symbols_map.load(in);
    sdsl::read_member(m_max_degree, in);
    sdsl::read_member(m_sigma, in);
}
#endif //LMS_GRAMMAR_REP_HPP
