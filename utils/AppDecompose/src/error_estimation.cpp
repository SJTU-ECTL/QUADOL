//
// Created by niqiu on 19-11-5.
//

#include <iostream>
#include "error_estimation.h"
#include "simulation.h"
#include "set"


// void Pmatrix::resize(int M, int N, int O) {
//     O = 1;
//     for(int i = 0; i < this->M; i++){
//         for(int j = 0; j < this->N; j++){
//             delete[] CPM[i][j];
//             //delete[] BD[i][j];
//         }
//         delete[] CPM[i];
//         //delete[] BD[i];
//     }
//     delete[] CPM;
//     delete[] BD;
//     this->M = M;
//     this->N = N + 10;
//     this->O = 1;
//     this->Index = new int[N + 10];
//     CPM = new int**[M];
//     BD = new int**[M];
//     for(int i = 0; i < M; i++){
//         CPM[i] = new int*[N + 10];
//         BD[i] = new int*[N + 10];
//         for(int j = 0; j < N + 10; j++){
//             CPM[i][j] = new int[O];
//             this->Index[j] = 0;
//             //BD[i][j] = new int[O];
//             for(int k = 0; k < 1; k++){
//                 CPM[i][j][k] = -1;
//                 //BD[i][j][k] = -1;
//             }
//         }
//     }
// }

Pmatrix::Pmatrix(int M, int N, int O) {
    O = 1;
    this->M = M;
    this->N = N + 1;
    this->O = O;
    CPM = new int**[M];
    BD = new int**[M];
    for(int i = 0; i < M; i++){
        CPM[i] = new int*[N + 1];
        BD[i] = new int*[N + 1];
        for(int j = 0; j < N + 1; j++){
            //CPM[i][j] = new int[O];
            CPM[i][j] = new int[O];
            //BD[i][j] = new int[O];
            for(int k = 0; k < O; k++){
                CPM[i][j][k] = -1;
                //BD[i][j][k] = -1;
            }
        }
    }
}


Pmatrix::~Pmatrix() {
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            delete[] CPM[i][j];
            //delete[] BD[i][j];
        }
        delete[] CPM[i];
        //delete[] BD[i];
    }
    delete[] CPM;
    delete[] BD;
}

void Pmatrix::init(const std::vector<std::vector<int>>* sim_result, Abc_Ntk_t* ntk){
    pNtk = ntk;
    sim = sim_result;

    Abc_Obj_t* pObj;
    int j;
    int min = 99999;
    Abc_NtkForEachNode(ntk, pObj, j){
            int id = Abc_ObjId(pObj);
            //cout << id << " " << Abc_ObjName(pObj) << endl;
            if(id < min) min = id;
        }
    this->offset = min;


    //int cnt = 0;
    for(int i = 0; i < M; i++){
        Abc_NtkForEachNode(ntk, pObj, j){
            int nodeid = Abc_ObjId(pObj);
            pObj->dTemp = (*sim_result)[nodeid][i];
        }
        Abc_NtkForEachNode(ntk, pObj, j){

                int numFanout = Abc_ObjFanoutNum(pObj);
                int nodeid = Abc_ObjId(pObj);
                string name_in = Abc_ObjName(pObj);

                BD[i][nodeid-offset] = new int[numFanout];


                int k;
                Abc_Obj_t* fanout;
                Abc_ObjForEachFanout(pObj, fanout, k){
                    string name_out = Abc_ObjName(fanout);
                    if(name_in == name_out){
                        continue;
                    }


                    int true_val = (*sim_result)[nodeid][i];

                    pObj->dTemp = 0;
                    simulate_node(ntk, fanout);
                    int result0 = (int)fanout->dTemp;

                    pObj->dTemp = 1;
                    simulate_node(ntk, fanout);
                    int result1 = (int)fanout->dTemp;

                    int result = (result0 != result1);
                    this->BD[i][nodeid-offset][k] = result;

                    pObj->dTemp = true_val;

                }
                //simulate_node(ntk, pObj);


            }

    }


}

void Pmatrix::calculateCPM() {
    struct cmp_level_t{
        bool operator()(Abc_Obj_t* in0, Abc_Obj_t* in1) const {
            return in0->Level > in1->Level;
        }
    }cmp_level;

    std::vector<Abc_Obj_t*> nodes_set;
    int j;
    Abc_Obj_t* pNode;
    Abc_NtkForEachNode(this->pNtk, pNode, j){
            nodes_set.emplace_back(pNode);
        }

    std::sort(nodes_set.begin(), nodes_set.end(), cmp_level);
    //int count2 = 0;
    //cout << nodes_set.size() << " " << N << endl;
    for(int i = 0; i < M; i++) {
        Abc_Obj_t* po;
        int k;

        Abc_NtkForEachPo(pNtk, po, k){
            Abc_Obj_t* node = Abc_ObjFanin0(po);
            int id = Abc_ObjId(node);
            CPM[i][id-offset][0] = 1;
        }



        for (auto &entry: nodes_set) {

            int id = Abc_ObjId(entry);
            string name_node = Abc_ObjName(entry);
            Abc_Obj_t* pFanout;
            int kk;


            Abc_NtkForEachPo(pNtk, po, kk) {
                int result = 0;
                if(CPM[i][id-offset][0] != -1) continue;
                Abc_ObjForEachFanout(entry, pFanout, k) {
                    string name_out = Abc_ObjName(pFanout);
                    int fanoutid = Abc_ObjId(pFanout);
                    if (name_node == name_out) continue;

                    int bd = BD[i][id-offset][k];
                    int cpm = CPM[i][fanoutid-offset][0];
                    assert(cpm != -1);
                    assert(bd != -1);
                    if(bd == 1 && cpm == 1){
                        result = 1;
                        break;
                    }

                }
                assert(result == 0 || result == 1);
                //assert(CPM[i][id-offset][kk] == -1);
                CPM[i][id-offset][0] = result;
                if(result == 1) break;


            }

            delete[] BD[i][id-offset];

        }

        delete[] BD[i];
    }

    //cout << "Num of 1: " << misc <<  " out of " << count2 << endl;
}


double Pmatrix::error_propagation(int id_chg, const int *local_error) const {
    int nPo = Abc_NtkPoNum(pNtk);
    auto aem = new double[nPo]{0};
    double er = 0;

    //int cnt = 0;

    for(int i = 0; i < this->M; i++){
        Abc_Obj_t* po;
        int k;
        int wrong = 0;
        if(local_error[i] == 1 &&
           CPM[i][id_chg-offset][0] == 1) er++;
        /*
        Abc_NtkForEachPo(this->pNtk, po, k){
            if(local_error[i] == 1 &&
               CPM[i][id_chg-offset][k] == 1){
                //er++;
                aem[k]++;
                wrong++;
            }
            //cout << CPM[i][id_chg-offset][k] << " " << k << endl;
        }
        //cout << "******************" << endl;
        if(wrong > 0) er++;
        */
    }

    //cout << "Num of 1 in CPM: " << cnt << endl;



    er /= M;
    delete[] aem;
    return er;
}


void set_level(Mffc* ffc);

/*------------------------------------------------*/
/* used to simulate one FFC basing on the given   */
/* input pattern information                      */
/*------------------------------------------------*/
int* simulate_local_error(vector<vector<int>>& sim_result, Mffc* ffc, int root_id){
    set_level(ffc);
    int count = 100000;
    int* result = new int[count]{0};
    //*result = 0;

    //int pout_id = Abc_ObjId(ffc->root);
    int pout_id = root_id;
    int level_max = ffc->ref->LevelMax;
    Abc_Obj_t* pin;
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    //base on the level of the network, I use a 2d Vector to store all the nodes with different level
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int j;
    Abc_NtkForEachNode(ffc->ref, pNode, j) {
            int level = pNode->Level - 1;
            if (level >= 0) {
                Vec_PtrPush(Nodes[level], pNode);
            }
        }
    for(int i=0; i<count; i++) {
        for (int k = 0; k < ffc->ref->vPis->nSize; k++) {
            pin = (Abc_Obj_t *) (ffc->ref->vPis->pArray[k]);
            pin->dTemp = sim_result[ffc->Pi_origin[k]->Id][i];
        }
        for (int n = 1; n < level_max; n++) {
            //simulate all the nodes base on their levels, we should simulate the low-level nodes first,
            // because the high-level nodes may take the low-level nodes as the Fanins
            int size = Nodes[n]->nSize;
            for (int l = 0; l < size; l++) {
                pin = (Abc_Obj_t *) Nodes[n]->pArray[l];
                simulate_node(ffc->ref, pin);
            }
        }
        pin = (Abc_Obj_t *) Nodes[level_max-1]->pArray[0];
        if(pin->dTemp != sim_result[pout_id][i]){
            result[i] = 1;
        }
    }
    for(int k=0;k<level_max;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
    return result;
}


void set_level(Mffc* ffc){
    Abc_Obj_t* pin;
    Abc_Obj_t* pNode;
    int max_level = 0;
    int j;
    Abc_NtkForEachNode(ffc->ref, pNode, j) {
            pNode->Level = 0;
        }
    int size = j;
    int count = 0;
    for (int k = 0; k < ffc->ref->vPis->nSize; k++) {
        pin = (Abc_Obj_t *) (ffc->ref->vPis->pArray[k]);
        pin->Level = 1;
        count ++;
    }
    j = 0;
    int count_old = 0;
    while(count != count_old) {
        count_old = count;
        Abc_NtkForEachNode(ffc->ref, pNode, j) {
                bool ready = true;
                int current_level = 0;
                if (pNode->Level == 0) {
                    for(int i = 0; i < pNode->vFanins.nSize; i++){
                        auto id = pNode->vFanins.pArray[i];
                        auto temp = Abc_NtkObj(ffc->ref, id);
                        if(temp->Level == 0){
                            ready = false;
                            break;
                        }
                        if(current_level < temp->Level + 1) current_level = temp->Level+1;
                    }
                    if(ready){
                        pNode->Level = current_level;
                        count ++;
                        if(current_level > max_level) max_level = current_level;
                    }
                }
            }
    }
    ffc->ref->LevelMax = max_level;
}

void Pmatrix::test() {
    int cnt = 0;
    for(int i = 0; i < M; i++){
        for(int j = 0; j < N; j++){
            for(int k = 0; k < O; k++){
                if(CPM[i][j][0] == 1)
                    cnt++;
            }
        }
        cout << cnt << endl;
    }

    int sz = 20;
    // auto arr = new arrBool(sz);

    // arr->setVal(12, 1);
    // //cout << (*arr)[12] << endl;
    // arr->setVal(0, 1);
    // //cout << (*arr)[0] << endl;
    // //cout << (*arr)[2] << endl;
    // arr->setVal(1, 1);
    // arr->setVal(3, 1);
    // arr->setVal(5, 1);
    // arr->setVal(0, 0);
    // for(int i = 0 ; i < sz; i++)
    //     cout << (*arr)[i] << " ";
    // cout << endl;

    // cout << cnt << " out of " << M*N*O << endl;

}


// int arrBool::operator[] (unsigned index) {
//     int d = sizeof(unsigned) * 8;
//     int aIndex = index/d;
//     unsigned phase = index - aIndex * d;
//     unsigned* num = this->val + aIndex;
//     unsigned temp = 1;
//     temp <<= phase;
//     int result = ((temp | (*num)) > 0);
//     return result;
// }

// void arrBool::setVal(unsigned index, unsigned val_) {
//     int d = sizeof(unsigned) * 8;
//     int aIndex = index/d;
//     unsigned phase = index - aIndex * d;
//     unsigned* num = this->val + aIndex;
//     // 0 0 -> 0, 0 1 -> 1, 1 1 -> 1, 1 0 -> 0
//     unsigned temp = val_;
//     if(val_ == 0){
//         temp <<= phase;
//         *num &= temp;
//     }
//     else{
//         temp <<= phase;
//         *num |= temp;
//     }



// }

// arrBool::arrBool(int size) {
//     this->size = size;
//     int d = sizeof(unsigned) * 8;
//     int aSize = size/d + 1;
//     this->val = new unsigned[aSize]{0};
// }

// arrBool::~arrBool() {
//     delete[] val;
// }
