//
// Created by niqiu on 19-11-5.
//

#ifndef ABC_LAST_ERROR_ESTIMATION_H
#define ABC_LAST_ERROR_ESTIMATION_H

#include <vector>
#include "abc_api.h"
#include "util.h"



class Pmatrix{
private:
    int M;
    int N;
    int O;
    int offset;
    //translate Abc_ObjId and the ordering in the matrix
    const std::vector<std::vector<int>>* sim;
    Abc_Ntk_t* pNtk;
    int misc;

public:
    int*** CPM;
    int*** BD;
    Pmatrix(int M, int N, int O);
    ~Pmatrix();


    void init(const std::vector<std::vector<int>>* sim_result, Abc_Ntk_t* ntk);
    //Effect: perform initialization and calculate BD matrix
    //Modifies: BD, offset, sim, pNtk;

    void calculateCPM();

    double error_propagation(int id_chg, const int* local_error) const;


    void test();

};

int* simulate_local_error(std::vector<std::vector<int>>& sim_result, Mffc* ffc, int root_id);


#endif //ABC_LAST_ERROR_ESTIMATION_H
