//
// Created by niqiu on 18-10-20.
//

#ifndef ABC_LAST_APPROX_H
#define ABC_LAST_APPROX_H

#include "simulation.h"
#include "util.h"
#include "algorithm"
#include <sstream>
#include <vector>
#include <random>

enum pattern_type{
    ALL_ONE, ALL_ZERO, PATTERN, COMPLTD_PATTERN, UNDECIDED
};

const string pattern_name[] = {
        "ALL_ONE", "ALL_ZERO", "PATTERN", "COMPLTD_PATTERN", "UNDECIDED"
};

int modify_constant0(Abc_Ntk_t* mffc);






class pattern{
public:
    pattern(int l, int r_size);
    ~pattern();

    double local_error(int ** dc, double** wt, int row_head, int row_tail, int* candidate);
    int* merge(int ** dc, double** wt, int row_head, int row_tail, int* candidates[], int candidate_size);
    int* p_cracker_helper(int ** dc, double** wt, int row_head, int row_tail);
    void pattern_cracker2(int** dc, double** wt);
    int** optimal_apparant_pt(int** dc, double** wt);
    double error_for_candidate(int** dc, double** wt, int* pt);

    bool percolate_right(int** dc, double** wt);
    bool percolate_up(int** dc, double** wt);
    double pt_entry_error(int position, int** dc, double** wt);
    double rowtype_entry_error(int position, int** dc, double** wt);
    double calculate_error(int ** dc, double** wt);

    void decide_rowtype(int** dc, double** wt);

    void print_info(bool detail);
    vector<int> decoder();

    double get_error(){return error;}
    int get_row_size(){return row_size;}
    int get_val(int i){return val[i];}
    void set_val(int i, int val_){this->val[i] = val_;}
    int get_length(){return length;}
    void copy_val(pattern* pt);
    pattern_type* get_pattern_type(){return row;};
private:
    int length;
    int row_size;
    int* val;
    double error;
    pattern_type* row;
};


class candidate{
public:
    pattern* pt;
    int* fSet;
    int* bSet;
    int bSetSize;
    int PiSize;
    int bi_index;
    double accu_error;
    candidate(pattern* pt_, int* fSet_, int* bSet_, int PiSize, double accu_error_, int bi_index_, int bSize=4);
    candidate();
    ~candidate();
    void print_info(bool isDetail);
    void consolidate();


};

Mffc* local_approx(Abc_Ntk_t* pNtk, Mffc* mffc, int PiSize, vector<vector<int>>& sim_result_ptr, int beamSize = 3, int min_level = 0);
int* choose_four(const int* arr, int size, int index);
int modify_mffc(Abc_Ntk_t* mffc, pattern* pt, int* bset, int* fset, int fSetSize, bool initial=true, int bSize = 4);
void write_back(Mffc* mffc, Abc_Ntk_t* original_ntk);
void amend_freeset(int** dec_chart, double** wt,
                   int fSetSize, int bSetSize, int aug_ID, int**& aug_dc, double**& aug_wt);


void print_array(int* arr, int sz);
void print_array(double* arr, int sz);

void write_back_input(Abc_Ntk_t* origin_Ntk, int origin_input, int replace_input,
                      int output_id, vector<int>& array);
#endif //ABC_LAST_APPROX_H
