//
// Created by niqiu on 19-7-29.
//

#include "approx.h"

struct cd_t{
    int index;
    double error;
};

struct cpr0{
    bool operator()(cd_t cd0, cd_t cd1) const {
        return cd0.error < cd1.error;
    }
}cpr0_t;

bool compare_array(int* arr0, int* arr1, int sz){
    for(int i = 0; i < sz; i++){
        if(arr0[i] != arr1[i]) return false;
    }
    return true;
}

inline int inverse_nq(int num){
    assert(num == 0 || num == 1);
    if(num) return 0;
    else return 1;
}

double pattern::error_for_candidate(int **dc, double **wt, int *pt) {
    double error = 0;
    double err1=0, err2=0;
    for(int i = 0; i < row_size; i++){
        if(row[i] == ALL_ONE || row[i] == ALL_ZERO) continue;

        for(int j = 0; j < length; j++){
            err1 += abs(pt[j] - dc[j][i]) * wt[j][i];
            err2 += abs(pt[j] - inverse_nq(dc[j][i])) * wt[j][i];
        }
    }
    assert(err1 <= 2 && err2 <= 2);
    error = err1 < err2? err1: err2;
    return error;
}


double pattern::rowtype_entry_error(int position, int **dc, double **wt) {
    double error = 0;
    int pos = position;
    for(int j = 0; j < length; j++){
        if(row[pos] == PATTERN){
            error += abs(dc[j][pos] - val[j]) * wt[j][pos];
        }else if(row[pos] == COMPLTD_PATTERN){
            int v = inverse_nq(val[j]);
            error += abs(dc[j][pos] - v) * wt[j][pos];
        }else if(row[pos] == ALL_ONE){
            error += abs(dc[j][pos] - 1) * wt[j][pos];
        }else if(row[pos] == ALL_ZERO){
            error += abs(dc[j][pos] - 0) * wt[j][pos];
        }else{
            assert(0);
        }
    }

    return error;
}

double pattern::pt_entry_error(int position, int** dc, double** wt){
    double error = 0;
    int pos = position;
    for(int i = 0; i < row_size; i++){

        if(row[i] == PATTERN) {
            error += abs(dc[pos][i] - val[pos]) * wt[pos][i];
        }else if(row[i] == COMPLTD_PATTERN){
            int v = inverse_nq(val[pos]);
            error += abs(dc[pos][i] - v) * wt[pos][i];
        }else if(row[i] == ALL_ONE){
            error += abs(dc[pos][i] - 1) * wt[pos][i];
        }else if(row[i] == ALL_ZERO){
            error += abs(dc[pos][i] - 0) * wt[pos][i];
        }else{
            assert(0);
        }
    }
    return error;
}

int min_four(double a1_, double a2, double a3, double a4){
    if(a1_ <= a2 && a1_ <= a3 && a1_ <= a4) return 1;
    if(a2 <= a1_ && a2 <= a3 && a2 <= a4) return 2;
    if(a3 <= a1_ && a3 <= a2 && a3 <= a4) return 3;
    if(a4 <= a1_ && a4 <= a2 && a4 <= a3) return 4;
    assert(0);
}

bool pattern::percolate_up(int** dc, double ** wt){
    bool changed = false;
    for(int i = 0; i < row_size; i++){
        if(row[i] == ALL_ONE || row[i] == ALL_ZERO) continue;
        pattern_type row_old = row[i];
        /*
        row[i] = ALL_ZERO;
        double error_0 = rowtype_entry_error(i, dc, wt);
        row[i] = ALL_ONE;
         */
        double error_1 = rowtype_entry_error(i, dc, wt);
        row[i] = PATTERN;
        double error_pt = rowtype_entry_error(i, dc, wt);
        row[i] = COMPLTD_PATTERN;
        double error_cpt = rowtype_entry_error(i, dc, wt);

        int min = min_four(999, 999, error_pt, error_cpt);
        switch(min){
            case 1:{
                row[i] = ALL_ZERO;
                break;
            }
            case 2:{
                row[i] = ALL_ONE;
                break;
            }
            case 3:{
                row[i] = PATTERN;
                break;
            }
            case 4:{
                row[i] = COMPLTD_PATTERN;
                break;
            }
            default: assert(0);
        }

        if(row[i] != row_old) changed = true;
    }
    return changed;
}

bool pattern::percolate_right(int** dc, double** wt){
    bool changed = false;
    for(int j = 0; j < length; j++){
        double error_old = pt_entry_error(j, dc, wt);
        val[j] = inverse_nq(val[j]);
        double error_new = pt_entry_error(j, dc, wt);

        if(error_old <= error_new) val[j] = inverse_nq(val[j]);
        else changed = true;
    }
    return changed;
}

pattern::pattern(int l, int r_size) {
    this->length = l;
    this->row_size = r_size;
    this->val = new int[length];
    this->row = new pattern_type[row_size];
    this->error = 0;
}

pattern::~pattern() {
    delete[] val;
    delete[] row;
}

int** pattern::optimal_apparant_pt(int** dc, double** wt){
    int sz = 5;
    int** ans = new int*[sz];
    for(int i = 0; i < sz; i++) ans[i] = new int[length];
    double error = 1;
    int index = -1;
    int* temp = new int[length];
    cd_t* cd_arr = new cd_t[row_size];
    for(int i = 0; i < row_size; i++){
        for(int j = 0; j < length; j++) temp[j] = dc[j][i];
        cd_arr[i].index = i;

        cd_arr[i].error = error_for_candidate(dc, wt, temp);

    }
    std::sort(cd_arr, cd_arr+row_size, cpr0_t);

    //for(int i = 0; i < row_size; i++) cout << i << " " << cd_arr[i].error << endl;

    index = cd_arr[0].index;
    //cout << "Index for first guy: " << index << endl;
    for(int j = 0; j < length; j++) ans[0][j] = dc[j][index];
    int cnt = 1;
    int i = 1;
    while(cnt < 5){
        if(i == row_size - 1){
            for(int k = cnt; k < 5; k++){
                for(int j = 0; j < length; j++) ans[k][j] = dc[j][index];
            }
            break;
        }


        index = cd_arr[i].index;
        for(int j = 0; j < length; j++) ans[cnt][j] = dc[j][index];
        if(compare_array(ans[cnt-1], ans[cnt], 16)) i++;
        else cnt++;

    }
    //for(i = 0; i < sz; i++) print_array(ans[i], 16);


    delete[] temp;
    delete[] cd_arr;
    return ans;
}


void pattern::pattern_cracker2(int** dc, double** wt){
    int isOne;
    int isZero;

    for(int i = 0; i < row_size; i++){
        isOne = 1;
        isZero = 0;
        row[i] = UNDECIDED;
        for(int j = 0; j < length; j++){
            isOne *= dc[j][i];
            isZero += dc[j][i];
        }
        if(isOne) row[i] = ALL_ONE;
        if(!isZero) row[i] = ALL_ZERO;
    }
    //Pre-compute to determine the case of ALL_ONEs and ALL_ZEROs

    int** pt = optimal_apparant_pt(dc, wt);
    int sz = 5;
    double err = 1;
    int index = -1;
    //cout << "Pattern_cracker2" << endl;
    for(int k = 0; k < sz; k++) {
        //print_array(pt[k], 16);
        this->error = local_error(dc, wt, 0, row_size - 1, pt[k]);
        //Calculate error
        for (int i = 0; i < length; i++) val[i] = pt[k][i];
        //Pass the result to the class attribute

        double error_new;
        for (int i = 0; i < length; i++) {
            val[i] = inverse(val[i]);
            error_new = local_error(dc, wt, 0, row_size - 1, val);
            if (error_new >= this->error) {
                val[i] = inverse(val[i]);
            } else {
                this->error = error_new;
            }
        }
        double err1, err2;
        for (int i = 0; i < row_size; i++) {
            if (row[i] == ALL_ONE || row[i] == ALL_ZERO) continue;
            err1 = err2 = 0;
            for (int j = 0; j < length; j++) {
                err1 += abs(val[j] - dc[j][i]) * wt[j][i];
                err2 += abs(val[j] - inverse_nq(dc[j][i])) * wt[j][i];
            }
            row[i] = err1 < err2 ? PATTERN : COMPLTD_PATTERN;
        }

        assert((this->error - this->calculate_error(dc, wt)) < 0.01);

        for(int z = 0; z < 2; z++) {
            percolate_right(dc, wt);
            percolate_up(dc, wt);
        }
        this->error = this->calculate_error(dc, wt);
        if(this->error < err){
            index = k;
            err = this->error;
        }
        //this->print_info(false);
    }
    //cout << "End of Pattern_cracker2" << endl;
    //cout << "**********************************************" << endl;


    for (int i = 0; i < length; i++) val[i] = pt[index][i];

    for(int i = 0; i < row_size; i++){
        isOne = 1;
        isZero = 0;
        row[i] = UNDECIDED;
        for(int j = 0; j < length; j++){
            isOne *= dc[j][i];
            isZero += dc[j][i];
        }
        if(isOne) row[i] = ALL_ONE;
        if(!isZero) row[i] = ALL_ZERO;
    }

    this->error = local_error(dc, wt, 0, row_size - 1, val);
    double error_new;
    for (int i = 0; i < length; i++) {
        val[i] = inverse(val[i]);
        error_new = local_error(dc, wt, 0, row_size - 1, val);
        if (error_new >= this->error) {
            val[i] = inverse(val[i]);
        } else {
            this->error = error_new;
        }
    }

    //cout << "error is: " << this->error << endl;

    double err1, err2;
    for (int i = 0; i < row_size; i++) {
        if (row[i] == ALL_ONE || row[i] == ALL_ZERO) continue;
        err1 = err2 = 0;
        for (int j = 0; j < length; j++) {
            err1 += abs(val[j] - dc[j][i]) * wt[j][i];
            err2 += abs(val[j] - inverse_nq(dc[j][i])) * wt[j][i];
        }
        row[i] = err1 < err2 ? PATTERN : COMPLTD_PATTERN;
    }
    //Decide the row_pattern_type

    for(int z = 0; z < 2; z++) {
        percolate_right(dc, wt);
        percolate_up(dc, wt);
    }
    this->error = this->calculate_error(dc, wt);

    int zero = 0;
    for(int i = 0; i < 16; i++){
        zero += val[i];
    }
    if(zero == 0){
        for(int i = 0; i < row_size; i++){
            if(row[i] == ALL_ZERO) row[i] = PATTERN;
        }
    }

    for(int i = 0; i < sz; i++) delete[] pt[i];
    delete[] pt;
}

vector<int> pattern::decoder(){
    vector<int> ans;
    for(int i = 0; i < this->length; i++){
        if(this->val[i] == 1) ans.push_back(i);
    }
    return ans;
}

void pattern::decide_rowtype(int** dc, double** wt){
    int isOne;
    int isZero;
    for(int i = 0; i < row_size; i++){
        isOne = 1;
        isZero = 0;
        row[i] = UNDECIDED;
        for(int j = 0; j < length; j++){
            isOne *= dc[j][i];
            isZero += dc[j][i];
        }
        if(isOne) row[i] = ALL_ONE;
        if(!isZero) row[i] = ALL_ZERO;
    }

    double err1, err2, err3, err4;
    for (int i = 0; i < row_size; i++) {
        if (row[i] == ALL_ONE || row[i] == ALL_ZERO) continue;
        err1 = 0;
        err2 = 0;
        err3 = 0;
        err4 = 0;
        for (int j = 0; j < length; j++) {
            err1 += abs(val[j] - dc[j][i]) * wt[j][i];
            err2 += abs(val[j] - inverse_nq(dc[j][i])) * wt[j][i];
            err3 += abs(0 - dc[j][i]) * wt[j][i];
            err4 += abs(1 - inverse_nq(dc[j][i])) * wt[j][i];
        }
        auto candidate1 = err1 < err2 ? PATTERN : COMPLTD_PATTERN;
        double error1 = err1 < err2 ? err1: err2;

        auto candidate2 = err1 < err2 ? ALL_ZERO : ALL_ONE;
        double error2 = err3 < err4 ? err3: err4;

        row[i] = error1 < error2? candidate1: candidate2;


    }


    for(int i = 0; i < row_size; i++){
        bool boom = (row[i] == UNDECIDED);
        assert(!boom);
    }

    this->error = this->calculate_error(dc, wt);
}


double pattern::local_error(int ** dc, double** wt, int row_head, int row_tail, int* candidate){
    double ans = 0;
    for(int i = row_head; i <= row_tail; i++){
        if(row[i] == ALL_ONE || row[i] == ALL_ZERO) continue;

        double error1 = 0;
        double error2 = 0;
        for(int j = 0; j < length; j++){
            error1 += abs(candidate[j] - dc[j][i]) * wt[j][i];
            error2 += abs(candidate[j] - inverse_nq(dc[j][i])) * wt[j][i];
        }
        ans += error1 < error2? error1: error2;
    }

    return ans;
}


int* pattern::merge(int ** dc, double** wt, int row_head, int row_tail, int* candidates[], int candidate_size){
    auto error = new double[candidate_size];
    for(int ctr = 0; ctr < candidate_size; ctr++){
        error[ctr] = local_error(dc, wt, row_head, row_tail, candidates[ctr]);
    }

    int minAt = 0;
    for(int i = 0; i < candidate_size; i++){
        if(error[i] < error[minAt]) minAt = i;
    }

    int* ans = new int[length];
    for(int i = 0; i < length; i++) ans[i] = candidates[minAt][i];
    delete[] error;
    return ans;
}

int* pattern::p_cracker_helper(int ** dc, double** wt, int row_head, int row_tail){
    //assert(row_tail >= row_head);
    if(row_tail <= row_head){
        int* ans = new int[length];
        for(int i = 0; i < length; i++) ans[i] = dc[i][row_head];
        return ans;
    }
    else{

        int mid1 = row_head + (row_tail - row_head)/3;
        int mid2 = row_tail - (row_tail - row_head)/3;
        //int mid2 = mid1 +
        int* candidate1 = pattern::p_cracker_helper(dc, wt, row_head, mid1);
        int* candidate2 = pattern::p_cracker_helper(dc, wt, mid1+1, mid2);
        int* candidate3 = pattern::p_cracker_helper(dc, wt, mid2+1, row_tail);

        int* candidates[3] = {candidate1, candidate2, candidate3};
        const int candidate_size = 3;

        int* ans = merge(dc, wt, row_head, row_tail, candidates, candidate_size);
        for(int i = 0; i != candidate_size; i++) delete[] candidates[i];
        return ans;
    }
}

void pattern::copy_val(pattern* pt){
    for(int i = 0; i < this->length; i++){
        this->set_val(i, pt->get_val(i));
    }
}


double pattern::calculate_error(int **dc, double **wt) {
    double ans = 0;
    //this->print_info(false);
    double check = 0;
    for(int i = 0; i < row_size; i++){
        for(int j = 0; j < length; j++) {
            if (row[i] == ALL_ONE) {
                ans += abs(1 - dc[j][i]) * wt[j][i];
            }
            else if(row[i] == ALL_ZERO) {
                ans += abs(0 - dc[j][i]) * wt[j][i];
            }
            else if(row[i] == PATTERN) {
                ans += abs(val[j] - dc[j][i]) * wt[j][i];
            }
            else{
                ans += abs(val[j] - inverse_nq(dc[j][i])) * wt[j][i];
            }
            check += wt[j][i];
        }
    }
    assert(check < 1.5);
    return ans;
}

void candidate::consolidate() {
    std::sort(this->bSet, this->bSet+4);
    std::sort(this->fSet, this->fSet+PiSize-4);
}


void candidate::print_info(bool isDetail) {
    pt->print_info(isDetail);
    print_array(bSet, 4);
    print_array(fSet, PiSize-4);
    cout << "bi_index: " << bi_index << endl;
    cout << "Accumulate error: " << accu_error << endl;
}

candidate::candidate(pattern *pt_, int *fSet_, int *bSet_, int PiSize_, double accu_error_, int bi_index_, int bSize) {
    pt = pt_;
    PiSize = PiSize_;
    bSet = new int[bSize];
    fSet = new int[PiSize - bSize];
    for (int i = 0; i < bSize; i++) bSet[i] = bSet_[i];
    for (int i = 0; i < PiSize - bSize; i++) fSet[i] = fSet_[i];

    accu_error = accu_error_;
    bi_index = bi_index_;
}

candidate::candidate() {
    pt = nullptr;
    fSet = bSet = nullptr;
    PiSize = 0;
    accu_error= 0;
    bi_index = -1;
}

candidate::~candidate() {

    delete[] bSet;
    delete[] fSet;

    delete pt;

}
