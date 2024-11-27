//
// Created by frank on 18-10-17.
//

#ifndef ABC_LAST_SIMULATION_H
#define ABC_LAST_SIMULATION_H
#include "abc_api.h"
#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <iomanip>
#include <set>
#include <queue>
#include <algorithm>
#include "util.h"

using namespace std;

extern int FFC_count;
extern int MFFC_count;
extern set<int> area;

struct weight_compare{
    bool operator()(const int a, const int b){
        return false;
    }
};

struct firstcome{
    bool operator()(const int a, const int b) const {
        if(a < b) return true;
        else return false;
    }
};

class decomposed_chart{
public:
    int** decompose_chart;
    double ** weight;
    int* array_row;
    int* array_column;
    int row;
    int column;
    int row_size;
    int column_size;

    void print_decomposed_chart();
    void print_weight();
    void clear_weight();
    void change_set(int* bound_set, int* free_set, int new_bound_size = 0, int new_free_size = 0);
    void change_set1(int* bound_set, int* free_set);
    decomposed_chart(int** truth_table, int row, int column, int* array_row, int* array_column){
        this->decompose_chart=truth_table;
        this->column = column;
        this->row = row;
        this->array_row = new int[row];
        this->array_column = new int[column];
        for(int i = 0; i < column; i ++) {
            this->array_column[i] = array_column[i];
        }
        for(int i = 0; i < row; i ++) {
            this->array_row[i] = array_row[i];
        }
        this->row_size=pow(2, row);
        this->column_size=pow(2, column);
        weight=new double*[column_size];
        for(int i=0; i<column_size; i++){
            weight[i]=new double[row_size];
        }
        for(int i=0; i<column_size; i++){
            for(int j=0; j<row_size; j++){
                weight[i][j]=0;
            }
        }
    }
    ~decomposed_chart(){
        for(int i=0; i<column; i++){
            delete[] decompose_chart[i];
        }
        delete[] decompose_chart;
        delete[] array_row;
        delete[] array_column;
        for(int i=0; i<column; i++){
            delete[] weight[i];
        }
        delete[] weight;
    }
};

int inverse(int num);

void binary(int* array, int size, int num);

void simulate_node(Abc_Ntk_t* LUT_network, Abc_Obj_t* pnode);

int simulation(Abc_Ntk_t* LUT_network, int* bound_set, int* free_set,\
        int bound_size, int free_size, int input_row, int input_column, int min_level=0);

decomposed_chart* create_decomposed_chart(Abc_Ntk_t* LUT_network, int* bound_set,\
        int* free_set, int num, int min_level=0, int bound_size = 4);

int binary_inverse(int* array, int size);

void simulate_weight(decomposed_chart* decomposedChart, \
        Abc_Ntk_t* MFFC_network, vector<Abc_Obj_t*> input);

void create_weight(decomposed_chart* decomposedChart, \
        Abc_Ntk_t* MFFC_network, Abc_Ntk_t* LUT_network, vector<Abc_Obj_t*> input);

//void simulate_whole(Abc_Ntk_t* LUT_network);
//void Create_Weight(decomposed_chart* decomposedChart, Mffc MFFC_network, int* bound_set, int* free_set);

vector<vector<int>> simulate_whole(Abc_Ntk_t* LUT_network);
void Create_Weight(decomposed_chart* decomposedChart, Mffc* MFFC_network, int* bound_set, int* free_set,  vector<vector<int>>& array);

int** create_dc(Abc_Ntk_t* LUT_network, int* bound_set, int* free_set, int num, int min_level = 0, int bound_size = 4);

Mffc* create_ffc_new(Mffc* origin_MFFC, Abc_Ntk_t* origin_ntk, int upperbound);

void find_new_FFC(Abc_Ntk_t* MFFC, set<int, firstcome>& root_array, set<int, firstcome>& ID_of_FFC,
                  set<int, firstcome> old_root_array, set<int, firstcome>& possible_FFC, set<int, firstcome>& visited,set<int, firstcome>& real_root, int upperbound);

Abc_Ntk_t * Abc_NtkCreateffc( Abc_Ntk_t * pNtk, Abc_Obj_t * pNode, char * pNodeName, set<int, firstcome>& ID_of_FFC,  set<int, firstcome>& root_array );

void  create_ffc_set_new(Mffc* origin_FFC, Abc_Ntk_t* origin_ntk, int upperbound, int lowerbound, set<int, firstcome>& changed, vector<vector<Mffc*>>& ffc_set);

double constant_0(Mffc* testcase, vector<vector<int>>& sim_result_ptr, int bound_size = 4);

double constant_zero(Abc_Ntk_t* origin_Ntk, int upperbound, int lowerbound, vector<vector<int>>& sim_result, set<int, firstcome>& changed, set<int>& critical_path, double& upperbound_error);

double not_zero(Abc_Ntk_t* origin_Ntk, int upperbound, int lowerbound, double upperbound_error, vector<vector<int>>& sim_result, set<int, firstcome>& changed, set<int>& critical_path);

int check_level(Mffc* new_mffc);

void find_critical_path(Abc_Ntk_t* origin_Ntk, set<int>& critical, int level = 0);

void binary_new(int* array, int size, int num);

bool check_path(set<int>& critical_path, int* array, int size);

int reduce_input(Mffc* mffc, vector<vector<int>>& sim_result_ptr, int* pattern, int* boundset, double& error);

void amend_freeset(int** dec_chart, double** wt, int* fSet, int fSetSize, int* bSet, int bSetSize, int aug_ID, int**& aug_dc, double**& aug_wt);

//Mffc* local_approx(Abc_Ntk_t* pNtk, Mffc* mffc, int PiSize, vector<vector<int>>& sim_result_ptr, set<int>& critical_path, int min_level = 0);

void not_zero_sasimi(Abc_Ntk_t* origin_Ntk, vector<vector<int>>& sim_result, double& upperbound_error);

void replace_two_signal(Abc_Ntk_t* origin_Ntk, int signal_id_1, int signal_id_2);

double not_zero_input_replace(Abc_Ntk_t* origin_Ntk, vector<vector<int>>& sim_result, double& upper_bound);

#endif //ABC_LAST_SIMULATION_H
