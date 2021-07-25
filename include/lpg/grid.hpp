//
// Created by diego on 5/5/21.
//

#ifndef LPG_COMPRESSOR_GRID_H
#define LPG_COMPRESSOR_GRID_H


#include <sdsl/rrr_vector.hpp>
#include <sdsl/wt_int.hpp>
#include <sdsl/construct.hpp>
#include "macros.hpp"
#include "utils.hpp"

struct grid_point{

    size_t row{};
    size_t col{};
    size_t label{};
    uint32_t level{};

    grid_point() = default;

    grid_point(const size_t& r,const size_t& c,const size_t& l,const size_t& lv):row(r),col(c),label(l),level(lv){}

    void save(std::fstream& out) const {
        out.write((const char *)&row,sizeof (row));
        out.write((const char *)&col,sizeof (col));
        out.write((const char *)&label,sizeof (label));
        out.write((const char *)&level,sizeof (level));
    }
    void load(std::fstream& in){
        in.read(( char *)&row,sizeof (row));
        in.read(( char *)&col,sizeof (col));
        in.read(( char *)&label,sizeof (label));
        in.read(( char *)&level,sizeof (level));
    }

};

struct grid_query{
    size_t row1{};
    size_t col1{};
    size_t row2{};
    size_t col2{};
};
class grid {

public:


    typedef grid_point                                     point;
    typedef grid_query                                     query;
    typedef size_t                                     size_type;
    typedef sdsl::rrr_vector<>                             wt_bv;
    //typedef sdsl::wt_int<
    //        wt_bv,
    //        wt_bv::rank_1_type,
    //        wt_bv::select_1_type,
    //        wt_bv::select_0_type
    //        >                                               wt_s;
    typedef sdsl::wt_int<> wt_s;
    typedef sdsl::rrr_vector<>                              bv_x;
    typedef sdsl::int_vector<>                                vi;

protected:
    wt_s sb;
    vi labels;

//    bv_x xa;
    bv_x xb;

//    bv_x::rank_1_type xb_rank1;
    bv_x::select_1_type xb_sel1;
//    bv_x::select_0_type xb_sel0;
//    bv_x::rank_1_type xa_rank1;

public:

    grid() = default;
    grid( const grid& _g ):sb(_g.sb), labels(_g.labels),xb(_g.xb) {
        compute_rank_select_st();
    }
    grid(const std::vector<point>& _points) {
        build(_points);
    }

    virtual ~grid() = default;
    void build(const std::vector<point>& _points) {
        std::vector<point> level_points(_points.size());
        size_type n_cols = 0,n_rows = 0, n_points = 0;
        for (const auto & _point : _points){
            level_points[n_points] = _point;
            n_cols = (n_cols < _point.col)? _point.col:n_cols;
            n_rows = (n_rows < _point.row)? _point.row:n_rows;
            ++n_points;
        }
#ifdef DEBUG_INFO
        std::cout<<"GRID:n_points:"<<n_points<<std::endl;
        std::cout<<"GRID:n_cols:"<<n_cols<<std::endl;
        std::cout<<"GRID:n_rows:"<<n_rows<<std::endl;
#endif
        sort_points(level_points);

#ifdef DEBUG_INFO
        std::fstream fout("points",std::ios::out|std::ios::binary);
        fout.write((const char *)&n_points,sizeof (n_points));
        fout.write((const char *)&n_cols,sizeof (n_cols));
        fout.write((const char *)&n_rows,sizeof (n_rows));
        for (size_t i = 0; i < n_points ; ++i) {
            _points[i].save(fout);
        }
        std::cout<<"points saved....."<<std::endl;

        std::cout<<"sort_points"<<n_rows<<std::endl;
#endif
        build_bitvectors(level_points,n_cols,n_rows,n_points);

#ifdef DEBUG_INFO
        std::cout<<"build_bitvectors"<<n_rows<<std::endl;
#endif
        build_wt_and_labels(level_points,n_cols,n_rows,n_points);
#ifdef DEBUG_INFO
        std::cout<<"build_wt_and_labels"<<n_rows<<std::endl;
#endif
        compute_rank_select_st();
#ifdef DEBUG_INFO
        std::cout<<"compute_rank_select_st"<<n_rows<<std::endl;
#endif
    }
    void build(const std::vector<point>& _points,uint32_t level) {

        std::vector<point> level_points;
        size_type n_cols = 0,n_rows = 0, n_points = 0;
        for (const auto & _point : _points)
            if(_point.level == level){
                level_points.push_back(_point);
                /**
                 * compute max and min of col and row values...
                 * */
                n_cols = (n_cols < _point.col)? _point.col:n_cols;
                n_rows = (n_rows < _point.row)? _point.row:n_rows;
                ++n_points;
            }

        sort_points(level_points);
        build_bitvectors(level_points,n_cols,n_rows,n_points);
        build_wt_and_labels(level_points,n_cols,n_rows,n_points);
        compute_rank_select_st();
    }

    void breakdown_space() const {
        std::cout<<"wt_s:"<< sdsl::size_in_bytes(sb)<<std::endl;
        std::cout<<"vi:"<< sdsl::size_in_bytes(labels)<<std::endl;
        std::cout<<"bv_x:"<< sdsl::size_in_bytes(xb)<<std::endl;
    }

protected:

    void sort_points(std::vector<point>& _points){
        /*
         * Sort _points by rows (rules) then by cols(suffix)
         * */
        sort(_points.begin(), _points.end(), [](const point &a, const point &b) -> bool
        {
            if (a.row < b.row) return true;
            if (a.row > b.row) return false;
            return (a.col) < (b.col);
        });
    }

    size_type map(const size_type  & row) const{
        assert(row > 0);
        return xb_sel1(row)-row+1;
    }

    void build_bitvectors(const std::vector<point>& _points, const size_type& n_cols,const size_type& n_rows,const size_type& n_points){

        std::vector<size_type> card_rows(n_rows, 0);
//        std::vector<size_type> card_cols(n_cols, 0);

        /*
        * Computing the cardinal of every column and every row
        * */
        for (size_type i = 0; i < n_points; ++i) {
            if(_points[i].row != 0){
                card_rows.at(_points[i].row - 1)++;
            }else{
                std::cout<<"ZERO WARNING!!!!"<<std::endl;
                std::cout<<"row:"<<_points[i].row<<std::endl;
                std::cout<<"col:"<<_points[i].col<<std::endl;
                std::cout<<"label:"<<_points[i].label<<std::endl;
            }
//            card_cols[_points[i].col - 1]++;
        }

        /**
       * Building bit_vectors XA and XB
       *
       * */

        auto build_bv = [](const std::vector<size_type>& V,const size_type &n_points,const size_type &cardinality){
            sdsl::bit_vector X(n_points + cardinality + 1, 0);
            X[0] = 1;
            /**
             * Put a 0 in XB(XA) for each element in the ith row(col)
             * then add 1
             * */
            size_t p = 1;
            for (uint j = 0; j < cardinality; ++j) {
                X[p + V[j]] = true;
                p += V[j] + 1;
            }
            return X;
        };

        xb = bv_x(build_bv(card_rows,n_points,n_rows));
#ifdef DEBUG_PRINT
        std::cout<<"GRID:XB"<<std::endl;
        for (int i = 0; i < xb.size(); ++i) {
            std::cout<<xb[i];
        }
        std::cout<<std::endl;
#endif
//        xa = bv_x(build_bv(card_rows,n_points,n_rows));
    }

    void build_wt_and_labels(const std::vector<point>& _points, const size_type& n_cols,const size_type& n_rows,const size_type& n_points){

        /**
         * Build a wavelet_tree on SB( index of the columns not empty ) and plain representation for SL(labels)
         * */
        std::ofstream sb_file("sb_file", std::ios::binary);
        sdsl::int_vector<> _sl(n_points,0);
        sdsl::int_vector<> _sb(n_points,0);

        size_type j = 0;
        for (size_type  i = 0; i < n_points; ++i) {
            _sl[j] = _points[i].label;
            _sb[j] = _points[i].col;
            ++j;
        }

        sdsl::util::bit_compress(_sl);
        sdsl::util::bit_compress(_sb);

        labels = vi(_sl);
        sdsl::serialize(_sb,sb_file);
        sb_file.close();
        //todo refactor wt construction...
        std::string id_sb = sdsl::util::basename("sb_file") + "_";
        sdsl::cache_config file_conf_sb(false,"./",id_sb);
        sdsl::construct(sb,"sb_file",file_conf_sb,0);
#ifdef DEBUG_PRINT
        std::cout<<"GRID:SB"<<std::endl;
        for (int i = 0; i < sb.size(); ++i) {
            std::cout<<sb[i]<<" ";
        }
        std::cout<<std::endl;
        std::cout<<"GRID:LABELS"<<std::endl;
        for (int i = 0; i < labels.size(); ++i) {
            std::cout<<labels[i]<<" ";
        }
        std::cout<<std::endl;
#endif
    }

    void compute_rank_select_st(){
        xb_sel1 = bv_x::select_1_type(&xb);
    }


public:

    void load(std::istream &in) {

        sdsl::load(labels,in);
        sdsl::load(sb,in);
        sdsl::load(xb,in);
        sdsl::load(xb_sel1,in);

        compute_rank_select_st();
    }

    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {

//        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;

        written_bytes += sdsl::serialize(labels,out);
        written_bytes += sdsl::serialize(sb ,out);
        written_bytes += sdsl::serialize(xb,out);
        written_bytes += sdsl::serialize(xb_sel1,out);

        return written_bytes;
    }

    size_type first_label_col(const size_type  & col) const{
        return labels[sb.select(1,col)];
    }
    size_type size_cols()const{ return sb.size();}

    void search_2d(const query& q,std::vector<size_type>& R) const{
        size_t p1,p2;
        p1 = map(q.row1);
        p2 = map(q.row2+1)-1;
        if(p1 > p2) return;
        auto res = sb.range_search_2d2(p1,p2,q.col1,q.col2);
        R.resize(res.first,0);
        for ( size_type i = 0; i < R.size(); i++ ){
            R[i] = res.second[i].second;
        }
    }
};


class grid_t {

public:
    typedef grid                                            _grid;
    typedef sdsl::int_vector<>                                 vi;

    typedef size_t                                      size_type;
    typedef grid_point                                     point;
    typedef grid_query                                     query;


protected:
    std::vector<_grid> grid_levels{};
    vi  level_columns_map;


public:

    grid_t () = default;

    grid_t (const grid_t & _g ){
        for (uint32_t i = 0; i < _g.grid_levels.size() ; ++i) {
            grid_levels[i] = _g.grid_levels[i];
        }
        level_columns_map = _g.level_columns_map;
    }

    grid_t (const std::vector<point>& _points,const uint32_t &_l) {
        grid_levels.resize(_l);
        for (uint32_t i = 0; i < _l  ; ++i) {
            grid_levels[i].build(_points,i+1);
//            std::cout<<"grid-level-"<<i+1<<std::endl;
        }
        sdsl::int_vector<> V(_points.size(),0);
        for (size_type i = 0; i < _points.size() ; ++i) {
            V[i] = _points[i].level;
        }
        sdsl::util::bit_compress(V);
        level_columns_map = vi(V);

    }

    void breakdown_space() const {

        std::cout<<"level_columns_map,"<<sdsl::size_in_bytes(level_columns_map)<<std::endl;
        uint i = 0;
        for (const auto &item : grid_levels) {
            std::cout<<"grid_levels["<<i+1<<"],"<<sdsl::size_in_bytes(grid_levels[i])<<std::endl;
            ++i;
        }

    }


    virtual ~ grid_t() {}

    void load(std::istream &in) {

        sdsl::load(level_columns_map,in);
        size_type levels;
        sdsl::load(levels,in);
        grid_levels.resize(levels);
        for (uint32_t i = 0; i < levels; ++i) {
            grid_levels[i].load(in);
        }
    }



    size_type serialize(std::ostream &out, sdsl::structure_tree_node *v, std::string name) const {
//        sdsl::structure_tree_node *child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
        size_t written_bytes = 0;
        written_bytes += sdsl::serialize(level_columns_map,out);
        std::cout<<"written_bytes += sdsl::serialize(level_columns_map,out);"<<std::endl;
        size_type levels = grid_levels.size();
        written_bytes += sdsl::serialize(levels,out);
        for (uint32_t i = 0; i < levels; ++i) {
            written_bytes += sdsl::serialize(grid_levels[i],out);
            std::cout<<"written_bytes += sdsl::serialize(grid_levels["<<i<<"],out);"<<std::endl;
        }
        return written_bytes;
    }


    size_type get_preorder_node_from_suffix(const size_type & sfx, const uint32_t& level) const{
        // if the sfx is in level it should be a point
        if(level_columns_map[sfx] != level)
            return 0;
        return grid_levels[level-1].first_label_col(sfx);
    }

    /**
     * O(G) this is linear on G but can reduce the number of extraction rules expanded to compares
     */
    void map_suffixes_levels(const uint32_t& level,std::vector<size_type>&V)const{
        size_type len = grid_levels[level-1].size_cols();
        V.resize(len,0);
        size_type j = 0;
        for (size_type i = 0; i < level_columns_map.size(); ++i) {
            if(level_columns_map[i] == level){
                V[j] = i+1;
                ++j;
            }
        }
    }

    void search (const grid_query& q,const uint32_t & level, std::vector<size_type>&results) const {
        grid_levels[level-1].search_2d(q,results);
    }
    uint32_t get_levels()const{ return grid_levels.size();}
};

#endif //LPG_COMPRESSOR_GRID_H
