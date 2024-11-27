//
// Created by niqiu on 18-9-24.
//

#ifndef ABC_LAST_UTIL_H
#define ABC_LAST_UTIL_H

#include "abc_api.h"
#include <iostream>
#include <cstdio>
#include <vector>
#include <map>
#include <cmath>
#include <set>

using namespace std;


class Mffc{
public:
    Abc_Ntk_t* ref; //An exact copy of this Mffc
    vector<Abc_Obj_t*> Pi_origin; //Store the Pi for the MFFc in the original network
    Abc_Obj_t* root; // Root of this Mffc
    vector<int> ID_origin;
    Mffc* duplicate();
    Mffc():ref(nullptr),root(nullptr), error(0), reduced_num(0){}
    ~Mffc();
    double error;
    int reduced_num;
    bool hasConstant = false;
    int level_reduce=0;
    void RI_write_back(int reduce_ID, int* pattern);
};


//void Abc_NodeMffcConeSuppPrint( Abc_Obj_t * pNode );
void NodeMffcConeSupp(Abc_Obj_t* pNode, Vec_Ptr_t* vCone, Vec_Ptr_t* vSupp, bool print);

Abc_Ntk_t* ReadBlif(char* filename);

vector<Mffc*> getMffcNtk(Abc_Ntk_t* pNtk, int minPiSize = 1, int maxPiSize = 999);
//EFFECT: Get all the sub-Mffc-network with size larger than ConeSize

void printMffcSetInfo(vector<Mffc*> pNtkSet, bool isNode);

map<int ,int> getMffcSetNodeInfo(vector<Mffc*> pNtkSet);

map<int ,int> getMffcSetPiInfo(vector<Mffc*> pNtkSet);
//EFFECT: ix->first represents the size of the cone. ix->second represents
//Home many times the cone of this size appear in pNtkSet

string dec2bin(unsigned n, int num_of_digit);

Abc_Ntk_t* abc_map(Abc_Ntk_t* in);

void clean_up(Abc_Ntk_t* pNtk, Mffc* lt, set<int>& deleted_nodes);

bool check_multi_fanout(Mffc* lt);


void map_LUT(char* old_filename, char* new_filename);


#endif //ABC_LAST_UTIL_H

