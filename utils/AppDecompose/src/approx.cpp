//
// Created by niqiu on 18-10-20.
//

#include "approx.h"
#define PRINT print_array(bSet, 4);print_array(fSet, PiSize-4);

int modify_constant0(Abc_Ntk_t* mffc){
    //Modify the mffc using bound set and free set, return the ID of the bound node
    Abc_Obj_t* newpi;
    int l;
    cout<<"the ID of PI is as following:"<<endl;
    Abc_NtkForEachPi(mffc, newpi, l){
        cout<<newpi->Id<<endl;
    }
    int result;
    Abc_Obj_t *pNode, *qNode, *vicNode;
    int j, max = 0;

    Abc_NtkForEachNode(mffc, pNode, j) {
            if (Abc_ObjId(pNode) > max) max = Abc_ObjId(pNode);
        }
    Abc_NtkForEachNode(mffc, pNode, j) {
            if (Abc_ObjId(pNode) == max) continue;
            Abc_NtkDeleteObj(pNode);
        }


    pNode = Abc_NtkCreateNodeConst0(mffc);

    //Abc_NtkForEachPi(mffc, qNode, j){
    //    Abc_ObjAddFanin(pNode, qNode);
    //}
    //pNode = Abc_NtkCreateNodeConst1(mffc);

    string temp = std::to_string(Abc_ObjId(pNode));
    Nm_ManStoreIdName(mffc->pManName, Abc_ObjId(pNode), ABC_OBJ_NODE, (char *) "nn", (char *) temp.c_str());
    result = Abc_ObjId(pNode);
    //Asign name
    //cout << "New node has name " << Abc_ObjName(pNode) << endl;

    qNode = Abc_NtkPo(mffc, 0);
    Abc_ObjRemoveFanins(qNode);
    Abc_ObjAddFanin(qNode, pNode);
    //Connect to Po


    //vicNode = Abc_NtkObj(mffc, Abc_NtkPo(mffc, 0)->vFanins.pArray[0]);
    //cout << "Victim has node name: " << Abc_ObjName(vicNode) << endl;
    return result;
}


double constant_0(Mffc* testcase, vector<vector<int>>& sim_result_ptr){
    auto sim_result = sim_result_ptr;
    int pisize = testcase->Pi_origin.size();
    int* bound_set = new int[4];
    int* free_set = new int[pisize -4];
    for(int i = 0; i < 4; i++){
        bound_set[i] = i + 1;
    }
    for(int i = 0; i < pisize - 4; i ++){
        free_set[i] = pisize - i;
    }
    double error = 0;
    decomposed_chart* dc = create_decomposed_chart(testcase->ref, bound_set, free_set, pisize);
    Create_Weight(dc, testcase, bound_set, free_set, sim_result);
    for(int i = 0; i < dc->column_size; i ++){
        for(int j = 0; j < dc->row_size; j ++){
            if(dc->decompose_chart[i][j] == 1){
                error += dc->weight[i][j];
            }
        }
    }
    cout<<"For the MFFC whose root is " << Abc_ObjName(testcase->root) << " and ID is "<< testcase->root->Id << ", its error for constant 0 is " << error<<endl;
    return error;
}

inline void swap(int& a, int& b){
    int temp = a;
    a = b;
    b = temp;
}

inline void inst(int *arr, int sz){
    for(int i = 1; i < sz; i++){
        for(int j = 0; j < i; j++){
            if(arr[j] > arr[i])
                swap(arr[j], arr[i]);
        }
    }
}

inline long int factorial(int i){
    long int ans = 1;
    for(int j = 1; j <= i; j++) ans *= j;
    return ans;
}

inline long int NChooseR(int n, int r){
    long int res = factorial(n)/factorial(r)/factorial(n-r);
    return res;

}

inline long int NPR(int n, int r){
    return factorial(n)/factorial(n-r);
}





void print_node(Abc_Ntk_t* pNtk){
    Abc_Obj_t* pObj;
    int i;
    Abc_NtkForEachNode(pNtk, pObj, i){
            cout << "Object with id " << Abc_ObjId(pObj) << " has name " << Abc_ObjName(pObj) << endl;
        }
    Abc_NtkForEachPi(pNtk, pObj, i){
        cout << "Pi with id " << Abc_ObjId(pObj) << " has name " << Abc_ObjName(pObj) << endl;
    }
    Abc_NtkForEachPo(pNtk, pObj, i){
        cout << "Po with id " << Abc_ObjId(pObj) << " has name " << Abc_ObjName(pObj) << endl;
    }

}


void pattern::print_info(bool detail){
    for(int i = 0; i < length; i++) cout << val[i] << " ";
    cout << endl;
    cout << "Error for this pattern is " << error << endl;
    if(detail) {

        for (int i = 0; i < row_size; i++) {
            cout << "Line " << i << "'s pattern type is "
                 << pattern_name[row[i]] << endl;
        }
    }
}
/**Function*************************************************************
  Synopsis    []
  Description [Find the optimized pattern]
***********************************************************************/
//the entry in ith row and jth column is stored in dc[j][i]



inline pattern_type flip(pattern_type a){
    if(a == PATTERN) return COMPLTD_PATTERN;
    else if(a == COMPLTD_PATTERN) return PATTERN;

    assert(0);
}



int* choose_four(const int* arr, int size, int index){
    int* ans = new int[size];
    for(int i = 0; i < size; i++){
        ans[i] = arr[i];
    }
    sort(ans, ans+size);
    int offset = factorial(size)/NChooseR(size, 4);

    for(int i = 0; i < offset*index; i++){
        next_permutation(ans, ans+size);
    }
    reverse(ans, ans+size);
    return ans;
}

string dec2bin(unsigned n, int num_of_digit)
{
    string res;

    for(int i = 0; i != num_of_digit; i++){
        res.push_back((n & 1) + '0');
        n >>= 1;
    }

    if (res.empty())
        res = "0";
    else {
        //std::reverse(res.begin(), res.end());
    }

    return res;
}



void print_array(int* arr, int sz){
    for(int i = 0; i < sz; i++){
        cout << arr[i] << " ";
    }
    cout << endl;
}

void print_array(double* arr, int sz){
    for(int i = 0; i < sz; i++){
        cout << arr << " ";
    }
    cout << endl;
}

void update_id(candidate* cd, int* PiSet, int PiSize, int oldId, int newId, int bSize = 4){
    assert(cd->PiSize == PiSize);
    for(int i = 0; i < bSize; i++) if(cd->bSet[i] == oldId) cd->bSet[i] = newId;
    for(int i = 0; i < PiSize-bSize; i++) if(cd->fSet[i] == oldId) cd->fSet[i] = newId;
    for(int i = 0; i < PiSize; i++) if(PiSet[i] == oldId) PiSet[i] = newId;
}

void consolidate(Abc_Ntk_t* pOld, Abc_Ntk_t* pNew, candidate* cd, int* PiSet, int PiSize, int bSize=4){
    Abc_Obj_t* pObj, *qObj;
    int oldId, newId, i, j;


    int* oldIds = new int[Abc_NtkNodeNum(pOld)];
    int* newIds = new int[Abc_NtkNodeNum(pNew)];
    assert(Abc_NtkNodeNum(pOld) == Abc_NtkNodeNum(pNew));
    int sz = Abc_NtkNodeNum(pOld);


    j = 0;
    Abc_NtkForEachNode(pOld, pObj, i){
            //cout << Abc_ObjId(pObj) << " ";
            oldIds[j] = Abc_ObjId(pObj);
            j++;
        }
    //cout << endl;

    j = 0;
    Abc_NtkForEachNode(pNew, pObj, i){
            //cout << Abc_ObjId(pObj) << " ";
            newIds[j] = Abc_ObjId(pObj);
            j++;
        }
    //cout << endl;

    //print_array(PiSet, PiSize);
    for(i = 0; i < sz; i++){
        if(oldIds[i] != newIds[i]){
            update_id(cd, PiSet, PiSize, oldIds[i], newIds[i], bSize);
        }
    }
    //print_array(PiSet, PiSize);
    delete[] oldIds;
    delete[] newIds;
}

struct cpr_cd{
    bool operator()(candidate* cd, candidate* cd2) const {
        return cd->accu_error <= cd2->accu_error;
    }
}cpr_cd_t;

struct cpr_pt_err{
    bool operator()(pattern* pt1, pattern* pt2) const {
        return pt1->get_error() <= pt2->get_error();
    }
}cpr_pt_err_t;


double compare_truth_table(int** arr0, int** arr1, double** wt, int r, int c){
    double ans = 0;
    for(int i = 0; i < c; i++){
        for(int j = 0; j < r; j++){
            if(arr0[i][j] != arr1[i][j]) ans += wt[i][j];
        }
    }
    return ans;
}




Mffc* local_approx(Abc_Ntk_t* pNtk, Mffc* mffc, int PiSize, vector<vector<int>>& sim_result_ptr, int beamSize, int min_level){
    Mffc* result;
    auto simu_result = sim_result_ptr;

    mffc->root->Level;

    int* PiSet = new int[PiSize];
    int i, j = 0;
    Abc_Obj_t* pNode;
    Abc_NtkForEachPi(mffc->ref, pNode, i){
        pNode->iTemp = 0;
        PiSet[j] = Abc_ObjId(pNode);
        j++;
    }


    pattern *pt = nullptr;
    int *nC4 = nullptr;
    decomposed_chart *dc = nullptr;


    //beamSize = 10;

    //auto errors = new double[beamSize];
    //for(i = 0; i < beamSize; i++) errors[i] = 1;


    Mffc** mffc_set = new Mffc*[beamSize]{nullptr};
    mffc_set[0] = mffc;
    for(i = 1; i < beamSize; i++) mffc_set[i] = mffc->duplicate();

    int** PiSetSet = new int*[beamSize]{nullptr};
    for(i = 0; i < beamSize; i++){
        PiSetSet[i] = new int[PiSize]{-1};
        for(j = 0; j < PiSize; j++) {
            PiSetSet[i][j] = PiSet[j];
        }
        //make five duplicates for the primarily input set
    }

    delete[] PiSet;
    //Update of beam search

    int PiSize_ = PiSize;
    int *bSet = new int[4]{0};
    int *bSet_ = new int[4]{0};
    int *fSet = new int[PiSize - 4]{0};
    int *fSet_ = new int[PiSize-4]{0};
    int *PiSet_copy = new int[PiSize];
    int r = (int)pow(2, PiSize-4);
    int c = (int)pow(2, 4);
    int **t_table = new int*[c];
    auto **wt_table = new double*[c];
    for(i = 0; i < c; i++){
        t_table[i] = new int[r];
        wt_table[i] = new double[r];
    }
    auto candidates= new candidate*[beamSize+1]{nullptr};

    int bSize_ = 4;
    int rounds = (int)ceil(((double)PiSize-4.0)/3.0);

    int residue = (PiSize + 1) % 3;
    int* residues = new int[beamSize];
    for(int ii = 0; ii < beamSize; ii++) residues[ii] = residue;
    bool* use_residues = new bool[beamSize];
    for(int ii = 0; ii < beamSize; ii++) use_residues[ii] = false;

    int residue_level = -1;
    int cnt = 0;
    //rounds = 1;
    for(int ii = 0; ii < rounds; ii++) {
        //primary loop
        int bSize = 4;
        //if(ii == residue_level) bSize = residue;
        long int loop_size = PiSize < 16 ? NChooseR(PiSize, 4) : 2000;


        for(j = 0; j < beamSize; j++){
            candidates[j] = new candidate;
            candidates[j]->accu_error = 10;
        }

        for(int bi = 0; bi < beamSize; bi++) {

            // beam search loop
            for (j = 0; j < PiSize; j++) PiSet_copy[j] = PiSetSet[bi][j];

            for (i = 0; i < loop_size; i++) {
                //Cn4 loop


                for(j = 0; j < factorial(bSize); j++) {
                    next_permutation(PiSet_copy, PiSet_copy + PiSize);
                }

                nC4 = new int[PiSize];
                for (j = 0; j < PiSize; j++) nC4[j] = PiSet_copy[PiSize - j - 1];

                for (j = 0; j < PiSize; j++) {
                    if (j < bSize) bSet[j] = nC4[j];
                    else fSet[j - bSize] = nC4[j];
                }

                delete[] nC4;


                if (i == 0) {
                    dc = create_decomposed_chart(mffc_set[bi]->ref, bSet, fSet, PiSize, min_level, bSize);
                    Create_Weight(dc, mffc_set[bi], bSet, fSet, simu_result);


                    if(ii == 0 && bi == 0){
                        // "remember" the initial truth table
                        assert(cnt == 0);
                        // Save the initial truthe table and bSet fSet
                        for(int jj = 0; jj < bSize; jj++) bSet_[jj] = bSet[jj];
                        for(int jj = 0; jj < PiSize - bSize; jj++) fSet_[jj] = fSet[jj];

                        for(int ic = 0; ic < c; ic++){
                            for(int ir = 0; ir < r; ir++){
                                t_table[ic][ir] = dc->decompose_chart[ic][ir];
                                wt_table[ic][ir] = dc->weight[ic][ir];
                            }
                        }

                        cnt++;
                    }
                } else {
                    assert(dc);
                    dc->change_set(bSet, fSet);
                    // calculate the truth table according to given bound set and free set
                }
                pt = new pattern((int) pow(2, bSize), (int) pow(2, PiSize - bSize));
                //pt = new pattern(16, (int) pow(2, PiSize - 4));
                pt->pattern_cracker2(dc->decompose_chart, dc->weight);
                // solve for optimal decomposition
                //if (ii == rounds - 1) pt->print_info(true);


                double err = pt->get_error() + mffc_set[bi]->error;


                auto cd = new candidate(pt, fSet, bSet, PiSize, err, bi, bSize);
                candidates[beamSize] = cd;
                std::sort(candidates, candidates+beamSize+1, cpr_cd_t);
                // add into the decision tree (beam search) only keep top beamSize branches
                delete candidates[beamSize];

            }
            //end of Cn4 iteration
            delete dc;

            if(ii == 0) break;
            //For the first round, we do not have to go through the beamSearch iteration process
        }
        //end of beam search iteration



        Mffc** mffc_set_ = new Mffc*[beamSize]{nullptr};
        for(j = 0; j < beamSize; j++){
            int index = candidates[j]->bi_index;
            assert(index != -1);
            if(index != j) {
                mffc_set_[j] = mffc_set[index]->duplicate();
                consolidate(mffc_set[index]->ref, mffc_set_[j]->ref, candidates[j], PiSetSet[j], PiSize, bSize);
                //cout << "Consolidating on the iteration ii = " << ii << " , j = " << j << endl;
            }
            else{
                mffc_set_[j] = mffc_set[index];
            }
        }

        for(j = 0; j < beamSize; j++){
            int index = candidates[j]->bi_index;
            if(index != j) delete mffc_set[j];
        }
        delete[] mffc_set;
        mffc_set = mffc_set_;




        bool use_residue;
        for(j = 0; j < beamSize; j++){
            use_residue = false;

            pattern* pt_op = nullptr;
            candidate* cd_op = nullptr;
            if (candidates[j]->pt->get_error() > 0.01 && residues[j] == 1) {

                auto cdd = candidates[j];
                //int bi_index = cdd->bi_index;
                int bi_index = j;

                int fSize = cdd->PiSize-3;
                //for (int k = 0; k < 4; k++) cout << cdd->bSet[k] << " ";
                //cout << endl;
                //for (int k = 0; k < cdd->PiSize - 4; k++) cout << cdd->fSet[k] << " ";
                //cout << endl;

                /*
                Abc_Obj_t* pnode;
                int ij;
                Abc_NtkForEachObj(mffc_set[bi_index]->ref, pnode, ij) {
                        cout << Abc_ObjId(pnode) << " ";
                }
                cout << endl;
                Io_WriteBlifLogic(mffc_set[bi_index]->ref, (char*)"temp_out", 1);
                //print information for debug
                */

                dc = create_decomposed_chart(mffc_set[bi_index]->ref, cdd->bSet, cdd->fSet, cdd->PiSize, min_level);
                Create_Weight(dc, mffc_set[bi_index], cdd->bSet, cdd->fSet, simu_result);
                //cdd->pt->print_info(false);
                //dc->print_decomposed_chart();
                //dc->print_weight();
                pattern** pts = new pattern*[4];

                int four = 4;
                for(int kkk = 0; kkk < four; kkk++){
                    int** aug_dc;
                    double** aug_wt;
                    amend_freeset(dc->decompose_chart, dc->weight,
                                  cdd->PiSize-4, 4, kkk, aug_dc, aug_wt);


                    pattern* niqiu = new pattern(16, (int) pow(2, fSize));
                    pts[kkk] = niqiu;
                    niqiu->copy_val(cdd->pt);
                    niqiu->decide_rowtype(aug_dc, aug_wt);

                    int c_size = (int)pow(2, 4);
                    int r_size = (int)pow(2, cdd->PiSize-3);

                    /*
                    if(niqiu->get_error() > cdd->pt->get_error()+100 && niqiu->get_row_size() <= 16){
                        cdd->pt->print_info(true);
                        dc->print_decomposed_chart();
                        dc->print_weight();
                        niqiu->print_info(true);
                        for(int iii = 0; iii < r_size; iii++){
                            for(int jj = 0; jj < c_size; jj++){
                                cout << aug_wt[jj][iii] << " ";
                            }
                            cout << endl;
                        }
                        for(int iii = 0; iii < r_size; iii++){
                            for(int jj = 0; jj < c_size; jj++){
                                cout << aug_dc[jj][iii] << " ";
                            }
                            cout << endl;
                        }
                        cout << "******************************************" << endl;
                    }
                    */

                    for(int iii = 0; iii < c_size; iii++){
                        delete[] aug_dc[iii];
                        delete[] aug_wt[iii];
                    }
                    delete[] aug_dc;
                    delete[] aug_wt;
                    //pts[kkk] = niqiu;
                }
                //try to add every bound set element as a free set input
                double min = 1;
                int index = -1;
                for(int kkk = 0; kkk < four; kkk++){
                    if(pts[kkk]->get_error() < min){
                        min = pts[kkk]->get_error();
                        index = kkk;
                    }
                }
                assert(index != -1);



                //cdd->pt->print_info(true);
                delete dc;
                std::sort(pts, pts+4, cpr_pt_err_t);
                pt_op = pts[0];
                if(pt_op->get_error() < cdd->pt->get_error()) use_residue = true;


                //cout << "Original error: " << cdd->pt->get_error() << endl;
                //for(int ijk = 0; ijk < four; ijk++) cout << pts[ijk]->get_error() << " ";
                //cout << endl;
                for(int ijk = 1; ijk < four; ijk++) delete pts[ijk];

                if(!use_residue) delete pts[0];
                else{
                    int* bset = new int[4];
                    int* fset = new int[fSize];
                    for(int ijk = 0; ijk < 4; ijk++) bset[ijk] = cdd->bSet[ijk];
                    for(int ijk = 1; ijk < fSize; ijk++) fset[ijk] = cdd->fSet[ijk-1];
                    fset[0] = cdd->bSet[index];
                    //add one bound set input to free set

                    int indexx = cdd->bi_index;
                    double err = pt_op->get_error() + mffc_set[indexx]->error;
                    int sz = cdd->PiSize+1;
                    cd_op = new candidate(pt_op, fset, bset, sz, err, j);
                    assert(cd_op->pt->get_val(0) == 1 || cd_op->pt->get_val(0) == 0);
                }
                delete[] pts;

            }
            //if the error in one round is too large and we haven't use residue,
            //try to apply non-disjoint decomposiiton

            if(!use_residue){
                delete cd_op;
                cd_op = candidates[j];
            }


            bool initial = (ii == 0);

            int newid = modify_mffc(mffc_set[j]->ref, cd_op->pt, cd_op->bSet, cd_op->fSet, cd_op->PiSize - bSize, initial, bSize);
            // Wirte back local transformation to the original circuit
            Io_WriteBlifLogic(mffc_set[j]->ref, (char*)"temp_out.blif", 1);

            /*
            int cst = 0;
            for(int k = 0; k < candidates[j]->pt->get_length(); k++){
                cst += candidates[j]->pt->get_val(k);
            }
            if(cst == 0) mffc_set[j]->hasConstant = true;
            */

            mffc_set[j]->error += cd_op->pt->get_error();

            PiSetSet[j][0] = newid;
            for(int k = 1; k < cd_op->PiSize - 3; k++) PiSetSet[j][k] = cd_op->fSet[k-1];

            auto t_table_new = create_dc(mffc_set[j]->ref, bSet_, fSet_, PiSize_, 0);
            double err_true = compare_truth_table(t_table, t_table_new, wt_table, r, c);
            mffc_set[j]->error = err_true;
            //cout << "At round " << ii << endl;
            //cout << "Real error vs. calculated error: " << err_true << " "
            //     << mffc_set[j]->error << endl;
            for(i = 0; i < c; i++){
                delete[] t_table_new[i];
            }
            delete[] t_table_new;

            if(use_residue) delete cd_op;


            //if (mffc_set[j]->error > 0.15) {
            //    mffc_set[j]->error = 1;
            //}
            use_residues[j] = use_residue;
        }


        for(int jj = 0; jj < beamSize; jj++){
            if(use_residues[jj]) residues[jj] -= 1;
        }

        assert(residue >= 0);
        //modify the mffc after the beam search

        for(j = 0; j < beamSize; j++){

            delete candidates[j];
        }

        PiSize -= (bSize-1);
        min_level++;
        if (PiSize <= 4) break;


    }
    //end of round iteration

    for(j = 0; j < beamSize; j++){
        auto t_table_new = create_dc(mffc_set[j]->ref, bSet_, fSet_, PiSize_, 0);
        double err_t = compare_truth_table(t_table, t_table_new, wt_table, r, c);
        for(i = 0; i < c; i++){
            delete[] t_table_new[i];
        }
        delete[] t_table_new;
        mffc_set[j]->error = err_t;
    }



    for(i = 0; i < beamSize; i++){
        delete[] PiSetSet[i];
    }
    delete[] PiSetSet;
    //delete[] errors;
    delete[] PiSet_copy;
    delete[] bSet;
    delete[] bSet_;
    delete[] fSet;
    delete[] fSet_;
    delete[] candidates;
    delete[] residues;
    delete[] use_residues;
    for(i = 0; i < c; i++){
        delete[] wt_table[i];
        delete[] t_table[i];
        //delete[] t_table_new[i];
    }
    //delete[] t_table_new;
    delete[] t_table;
    delete[] wt_table;
    //Io_WriteBlifLogic(mffc_set[0]->ref, (char*)"temp_out.blif", 1);

    int index = 0;
    for(i = 0; i < beamSize; i++){
        if(mffc_set[i]->error < mffc_set[index]->error) index = i;
    }

    for(i = 0; i < beamSize; i++){
        if(i == index) continue;
        delete mffc_set[i];
    }

    return mffc_set[index];

}

int modify_mffc(Abc_Ntk_t* mffc, pattern* pt, int* bset, int* fset, int fSetSize, bool initial, int bSize){
    //Modify the mffc using bound set and free set, return the ID of the bound node
    Abc_Obj_t* newpi;
    int l;
    int result;
    Abc_Obj_t *pNode, *qNode, *vicNode;
    int j, max = 0;
    if(initial) {
        Abc_NtkForEachNode(mffc, pNode, j) {

                if (Abc_ObjId(pNode) > max) max = Abc_ObjId(pNode);
            }
        Abc_NtkForEachNode(mffc, pNode, j) {
                if (Abc_ObjId(pNode) == max) continue;
                Abc_NtkDeleteObj(pNode);
            }
    }else{
        vicNode = Abc_NtkObj(mffc, Abc_NtkPo(mffc, 0)->vFanins.pArray[0]);
        Abc_NtkDeleteObj(vicNode);
        //
    }
    Abc_NtkForEachObj(mffc, pNode, j) {
            //cout << Abc_ObjId(pNode) << " ";
        }
    //cout << endl;
    //Delete old node

    vector<int> numerical_pattern = pt->decoder();

    bool isConstant1 = false;
    //bool isConstant0 = false;
    if(numerical_pattern.empty()){
        pNode = Abc_NtkCreateNodeConst0(mffc);
        //isConstant0 = true;
    }
    else if(numerical_pattern.size() == (int)pow(2, bSize)){
        pNode = Abc_NtkCreateNodeConst1(mffc);
        isConstant1 = true;
    }
    else {

        pNode = Abc_NtkCreateNode(mffc);
        for (int i = 0; i < bSize; i++) {
            qNode = Abc_NtkObj(mffc, bset[i]);
            Abc_ObjAddFanin(pNode, qNode);
        }
    }

    //Create bound set node
    string temp = std::to_string(Abc_ObjId(pNode));
    Nm_ManStoreIdName(mffc->pManName, Abc_ObjId(pNode), ABC_OBJ_NODE, (char *) "nn", (char *) temp.c_str());
    result = Abc_ObjId(pNode);
    //Asign name
    //cout << "New node has name " << Abc_ObjName(pNode) << endl;


    /*
    for(int i = 0; i < 4; i++) {
        cout << Abc_ObjName(Abc_NtkObj(mffc, bset[i])) << " " << bset[i] << endl;
    }*/

    vector<int>::size_type ix;
    static char buffer[2000] = {0};
    for (int i = 0; i != 2000; i++) {
        buffer[i] = 0;
    }
    string s;

    if(!isConstant1 ) {
        for (ix = 0; ix != numerical_pattern.size(); ix++) {
            //if (isConstant1) break;
            //s = dec2bin((unsigned)numerical_pattern[ix], 4);
            s = dec2bin((unsigned) numerical_pattern[ix], bSize);
            //cout << s << endl;
            sprintf(buffer, "%s%s 1\n", buffer, s.c_str());
        }
        if(numerical_pattern.empty()){
            sprintf(buffer, "%s%s 0\n", buffer, s.c_str());
        }

        char *buffer_use = new char[2000];//Memory leak here
        for (int i = 0; i < 2000; i++) buffer_use[i] = buffer[i];
        pNode->pData = buffer_use;
    }
    //Assign pData

    auto pattern_type = pt->get_pattern_type();
    int zeroRowCnt = 0;
    for(int i = 0; i != pt->get_row_size(); i++){
        if(pattern_type[i] == ALL_ZERO) zeroRowCnt++;
    }


    Abc_Obj_t* endNode;
    if(zeroRowCnt == pt->get_row_size()){
        endNode = Abc_NtkCreateNodeConst0(mffc);
    }
    else{
        endNode = Abc_NtkCreateNode(mffc);
    }
    //Create auxiliary node(used to connected to Po)
    temp = std::to_string(Abc_ObjId(endNode));
    Nm_ManStoreIdName(mffc->pManName, Abc_ObjId(endNode), ABC_OBJ_NODE, (char *) "xx", (char *) temp.c_str());
    //Create name

    if(zeroRowCnt != pt->get_row_size()){


        Abc_ObjAddFanin(endNode, pNode);

        for(int i = 0; i < fSetSize; i++) {
            Abc_ObjAddFanin(endNode, Abc_NtkObj(mffc, fset[i]));
        }
        //Assign free set as endNodes's fanin


        /*
        for(int i = 0; i < fSetSize; i++){
            cout << Abc_ObjName(Abc_NtkObj(mffc, fset[i]))  << " " << fset[i] << endl;
        }*/

        int sz = 1;
        for(int i = 0; i < fSetSize; i++) sz *= 2;
        sz *= 40;
        char* buffer2 = new char[sz]{0};
        //static char buffer2[200000] = {0};
        //For PiSet larger than 15 we need 200000
        //For PiSet larger than 16 we need 2000000
        for(int jk = 0; jk != sz; jk++) {
            buffer2[jk] = 0;
        }
        //auto pattern_type = pt->get_pattern_type();
        for(int i = 0; i != pt->get_row_size(); i++){
            s = dec2bin((unsigned)i, fSetSize);
            switch(pattern_type[i]){
                case(ALL_ONE):{
                    sprintf(buffer2, "%s-%s 1\n", buffer2, s.c_str());
                    break;
                }
                case(ALL_ZERO):{
                    break;
                }
                case(PATTERN):{
                    //s = dec2bin((unsigned)i, fSetSize);
                    sprintf(buffer2, "%s1%s 1\n", buffer2, s.c_str());
                    break;
                }
                case(COMPLTD_PATTERN):{
                    //s = dec2bin((unsigned)i, fSetSize);
                    sprintf(buffer2, "%s0%s 1\n", buffer2, s.c_str());
                    break;
                }
                default: assert(0);
            }
        }

        endNode->pData = buffer2;
        //Endnode's pData
    }

    qNode = Abc_NtkPo(mffc, 0);
    Abc_ObjRemoveFanins(qNode);
    Abc_ObjAddFanin(qNode, endNode);
    //Connect to Po


    //vicNode = Abc_NtkObj(mffc, Abc_NtkPo(mffc, 0)->vFanins.pArray[0]);
    //cout << "Victim has node name: " << Abc_ObjName(vicNode) << endl;
    return result;
}


void write_back(Mffc* mffc, Abc_Ntk_t* original_ntk){
    Abc_Obj_t* pNode, *pFanin, *qNode;
    Abc_Obj_t* victim;
    int i, j, k;

    /*
    vectr<Abc_Obj_t*>::size_type ix;
    for(ix = 0; ix != mffc->Pi_origin.size(); ix++){
        pNode = mffc->Pi_origin[ix];
        Abc_ObjForEachFanout(pNode, victim, i){
            Abc_NtkDeleteObj(victim);
        }
    }
    Abc_ObjForEachFanin(mffc->root, victim, i){
        Abc_NtkDeleteObj(victim);
    }
    */

    Vec_Ptr_t *vCone, *vSupp;
    vCone = Vec_PtrAlloc(100);
    vSupp = Vec_PtrAlloc(100);
    bool null_or_not = false;
    //NodeMffcConeSupp(mffc->root, vCone, vSupp, false);
    for(i = 0; i < mffc->Pi_origin.size(); i++){
        Abc_Obj_t* temp = Abc_NtkObj(original_ntk, mffc->ID_origin[i]);
        if(temp == NULL){
            null_or_not = true;
            cout<<"The node is NULL"<<endl;
            break;
        }
        Vec_PtrPush(vSupp, temp);
    }
    for(i = mffc->Pi_origin.size()+1; i < mffc->ID_origin.size(); i ++){
        Abc_Obj_t* temp = Abc_NtkObj(original_ntk, mffc->ID_origin[i]);
        if(temp == NULL){
            null_or_not = true;
            cout<<"The node is NULL"<<endl;
            break;
        }
        Vec_PtrPush(vCone, temp);
    }


    if(!null_or_not) {
        Vec_PtrForEachEntry(Abc_Obj_t*, vCone, pNode, i) {
            if (Abc_ObjId(pNode) != Abc_ObjId(mffc->root)) {
                Abc_NtkDeleteObj(pNode);
            }
        }

        Vec_PtrFree(vCone);
        Vec_PtrFree(vSupp);

        Abc_NtkForEachPi(mffc->ref, pNode, i) {
            pNode->pCopy = mffc->Pi_origin[i];
        }

        vector<Abc_Obj_t *> internal_nodes;
        Abc_NtkForEachNode(mffc->ref, pNode, i) {
                Abc_NtkDupObj(original_ntk, pNode, 0);
                Abc_ObjForEachFanin(pNode, pFanin, k) {
                    Abc_ObjAddFanin(pNode->pCopy, pFanin->pCopy);
                }
            }
        assert(pNode);
        qNode = mffc->root;
        Abc_ObjRemoveFanins(qNode);
        static char pData[10] = "1 1\n";
        qNode->pData = pData;
        Abc_ObjAddFanin(qNode, pNode->pCopy);
        Abc_ObjReplace(qNode, pNode->pCopy);
    }
    else{
        Vec_PtrFree(vCone);
        Vec_PtrFree(vSupp);
    }

}


void amend_freeset(int** dec_chart, double** wt, int fSetSize, int bSetSize, int aug_ID, int**& aug_dc, double**& aug_wt){
    int new_free_size = fSetSize+1;
    int column_size = pow(2, bSetSize);
    int row_size = pow(2, new_free_size);
    aug_dc = new int*[column_size];
    for(int i = 0; i < column_size; i++){
        aug_dc[i] = new int[row_size];
    }
    aug_wt = new double*[column_size];
    for(int i = 0; i < column_size; i++){
        aug_wt[i] = new double[row_size];
    }
    int* binary_free = new int[new_free_size];
    int* binary_origin = new int[fSetSize];
    int* binary_bound = new int[bSetSize];
    for(int i = 0; i < column_size; i++){
        for(int j = 0; j < row_size/2; j++){
            aug_dc[i][j] = dec_chart[i][j];
            aug_dc[i][j+row_size/2] = dec_chart[i][j];
            //aug_dc[i][2*j+1] = dec_chart[i][j];
            binary(binary_origin, fSetSize, i);
            if(binary_origin[aug_ID]==0){
                aug_wt[i][j] = wt[i][j];
                aug_wt[i][j+row_size/2] = 0;
            }
            else{
                aug_wt[i][j+row_size/2] = wt[i][j];
                aug_wt[i][j] = 0;
            }
        }
    }
    delete[] binary_bound;
    delete[] binary_free;
    delete[] binary_origin;
}

void write_back_input(Abc_Ntk_t* origin_Ntk, int origin_input, int replace_input,
                      int output_id, vector<int>& array){
    auto output_node = Abc_NtkObj(origin_Ntk, output_id);
    auto deleted_node = Abc_NtkObj(origin_Ntk, origin_input);
    auto replace_node = Abc_NtkObj(origin_Ntk, replace_input);
    Abc_ObjDeleteFanin(output_node, deleted_node);
    Abc_ObjAddFanin(output_node, replace_node);
    int gg = 6*pow(2,output_node->vFanins.nSize);
    char* _pdata = new char[gg];
    int count = 0;
    int one_count = 0;
    for(int i = 0; i < pow(2, output_node->vFanins.nSize); i++){
        int* bin = new int[output_node->vFanins.nSize];
        if(array[i] == 1){
            one_count += 1;
            binary(bin, output_node->vFanins.nSize, i);
            for(int j = 0; j < output_node->vFanins.nSize; j++){
                 if(bin[j] == 1) _pdata[count] = '1';
                 else _pdata[count] = '0';
                 count ++;
            }
            _pdata[count] = ' ';
            _pdata[count + 1] = '1';
            _pdata[count + 2] = '\n';
            count += 3;
        }
    }
    if(one_count == 0){
        Abc_ObjRemoveFanins(output_node);
        _pdata[0] = ' ';
        _pdata[1] = '0';
        _pdata[2] = '\n';
    }
    else if(one_count == pow(2, output_node->vFanins.nSize)){
        Abc_ObjRemoveFanins(output_node);
        _pdata[0] = ' ';
        _pdata[1] = '1';
        _pdata[2] = '\n';
    }
    output_node->pData = _pdata;
}
