//
// Created by frank on 18-10-17.
//
#include "simulation.h"
#include "approx.h"
#include "error_estimation.h"
#include <map>

extern Pmatrix PM;

extern vector<pair<int, int>> replaced_pair;

extern vector<pair<int, double>> efficience;

set<int> removed;

int replace_count = 0;

int reduced_area = 0;

struct cpr_rl{
    bool operator()(Mffc* pMffc, Mffc* qMffc) const {
        return pMffc->root->Level < qMffc->root->Level;
    }
}cpr_rl_t;

void decomposed_chart::change_set1(int* bound_set, int* free_set) {
    clock_t gg = clock();
    int map[column+row];
    int biggest = 0;
    for(int i = 0; i < row; i ++){
        map[i] = array_row[i];
        if(array_row[i] > biggest) biggest = array_row[i];
    }
    for(int i = 0; i < column; i++){
        map[i+row] = array_column[i];
        if(array_column[i] > biggest) biggest = array_column[i];
    }
    int new_dc[column_size][row_size];
    double new_weight[column_size][row_size];
    for(int n=0; n<column_size; n++) {
        for (int m = 0; m < row_size; m++) {
            new_dc[n][m] = decompose_chart[n][m];
            new_weight[n][m] = weight[n][m];
        }
    }
    int bina[row+column];
    for(int i = 0; i < row+column; i++){
        bina[i] = pow(2,i);
    }
    int array[biggest];
    int bin[biggest];
    /*for(int i = 0; i < column; i ++){
        //if(bound_set[i] > biggest) biggest = bound_set[i];
        bin[bound_set[i]] = pow(2, bound_set[i]);
    }
    for(int i = 0; i < row; i ++){
        //if(free_set[i] > biggest) biggest = free_set[i];
        bin[free_set[i]] = pow(2, free_set[i]);
    }*/
    //int array[column+row];
    int input_column[column];
    int input_row[row];
    int new_column, new_row;
    for(int n=0; n<column_size; n++){
        for(int m=0; m<row_size; m++){
            int sum = row_size * n + m;
            for(int i = 0; i < row+column; i ++){
                if((sum & bina[i]) == 0) array[map[i]] = 0;
                else array[map[i]] = 1;
                //cout<<array[map[i]]<<endl;
            }
            //cout<<"-------------------------"<<endl;
            /*binary(input_column, column, n);
            binary(input_row, row, m);
            for(int i = 0; i < row; i ++){
                cout<< input_row[i] <<endl;
            }
            for(int i = 0; i < column; i ++){
                cout<< input_column[i] <<endl;
            }*/
            //for(int j=0; j<column; j++){
                //array[array_column[j]] = input_column[j];
            //}
            //for(int j=0; j<row; j++){
                //array[array_row[j]] = input_row[j];
            //}
            for(int j=0; j<column; j++){
                input_column[j]=array[bound_set[j]];
            }
            for(int j=0; j<row; j++){
                input_row[j]=array[free_set[j]];
            }
            new_column=binary_inverse(input_column, column);
            new_row=binary_inverse(input_row, row);
            weight[new_column][new_row]= new_weight[n][m];
            decompose_chart[new_column][new_row]=new_dc[n][m];
            //array.clear();
        }
    }
    for(int i = 0; i < row; i ++){
        array_row[i] = free_set[i];
    }
    for(int i = 0; i < column; i ++){
        array_column[i] = bound_set[i];
    }
    gg = clock() - gg;
    //cout << "dec2bi completes in " << (double)gg/CLOCKS_PER_SEC << "s" << endl;

}


void decomposed_chart::change_set(int* bound_set, int* free_set, int new_bound_size, int new_free_size) {
    if(new_bound_size == 0 && new_free_size == 0){
        new_bound_size = column_size;
        new_free_size = row_size;
    }
    clock_t tt = clock();
    clock_t gg = clock();
    int new_dc[column_size][row_size];
    double new_weight[column_size][row_size];
    for(int n=0; n<column_size; n++){
        for(int m=0; m<row_size;m++){
            new_dc[n][m]=decompose_chart[n][m];
            new_weight[n][m]=weight[n][m];
        }
    }
    int biggest = 0;
    for(int i = 0; i < column; i ++){
        if(bound_set[i] > biggest) biggest = bound_set[i];
    }
    for(int i = 0; i < row; i ++){
        if(free_set[i] > biggest) biggest = free_set[i];
    }
    int array[biggest];
    gg = clock() - gg;
    //cout << "dec2bi completes in " << (double)gg/CLOCKS_PER_SEC << "s" << endl;
    //int array[column+row];
    int input_column[column];
    int input_row[row];
    int new_column, new_row;
    for(int n=0; n<column_size; n++){
        for(int m=0; m<row_size; m++){
            gg = clock();
            binary(input_column, column, n);
            binary(input_row, row, m);
            gg = clock() - gg;
            //cout << "dec2bi completes in " << (double)gg/CLOCKS_PER_SEC << "s" << endl;
            gg = clock();
            for(int j=0; j<column; j++){
                array[array_column[j]] = input_column[j];
            }
            for(int j=0; j<row; j++){
                array[array_row[j]] = input_row[j];
            }
            for(int j=0; j<column; j++){
                input_column[j]=array[bound_set[j]];
            }
            for(int j=0; j<row; j++){
                input_row[j]=array[free_set[j]];
            }
            gg = clock() - gg;
            //cout << "dec2bi completes in " << (double)gg/CLOCKS_PER_SEC << "s" << endl;
            gg = clock();
            new_column=binary_inverse(input_column, column);
            new_row=binary_inverse(input_row, row);
            gg = clock() - gg;
            //cout << "dec2bi completes in " << (double)gg/CLOCKS_PER_SEC << "s" << endl;
            weight[new_column][new_row]= new_weight[n][m];
            decompose_chart[new_column][new_row]=new_dc[n][m];
            //array.clear();
        }
    }
    for(int i = 0; i < row; i ++){
        array_row[i] = free_set[i];
    }
    for(int i = 0; i < column; i ++){
        array_column[i] = bound_set[i];
    }
    tt = clock() - tt;
    //cout << "dec2bi completes in " << (double)gg/CLOCKS_PER_SEC << "s" << endl;;
    //cout << "dec2bi completes in " << (double)tt/CLOCKS_PER_SEC << "s" << endl;

}

int inverse(int num){
    if(num==1) return 0;
    else if(num==0) return 1;
    else {
        std::cout<<"the input can only be 0 or 1!"<<std::endl;
        std::cout << num << endl;
    }

    return 0;
}

int binary_inverse(int* array, int size){
    int result=0;
    for(int i=0; i<size; i++){
        result+=array[i]*pow(2,i);
    }
    return result;
}

//base on the size of the bound set or free set, it can change an integer into an array,
// which is the binary form of the integer with this function we can convert
// an integer into the real input case of the bound set or free set. for example,
// if the bound set has 4 elements, they are i1, i2, i3, i4, and now I give a value 7,
// which means i1=i2=i3=1, i4=0 then the parameter size should be 4,
// and num should be 7, and the elements stored in array
// will become array[0]=array[1]=array[2]=1, array[3]=0
void binary(int* array, int size, int num){
    for(int i=0;i<size;i++){
        array[i]=num%2;
        num=num/2;
    }
}

void binary_new(int* array, int size, int value){
    for(int i = size - 1; i >= 0; i --){
        if(value & (1 << i)) array[i] = 1;
        else array[i] = 0;
    }
}

//simulate one specific node, base on the total network, the value of
// the Fanin of the node, and the pData, I can get the value of the output of the node,
// and store it in the dTemp
void simulate_node(Abc_Ntk_t* LUT_network, Abc_Obj_t* pnode){
    char* read=(char*)pnode->pData;
    int inverseornot=0;
    if(strcmp(read," 1\n") == 0){
        pnode->dTemp = 1;
        return;
    }
    else if(strcmp(read," 0\n") == 0){
        pnode->dTemp = 0;
        return;
    }
    int output=0;
    //the value of the output
    int check=0;
    //used to check whether the pData reach its end or not
    while(1){
        int save=1;
        //use to temporarily store the value of one line,
        // cause we know sometimes pData may have two or more lines
        //cout<<"the pnode is "<<pnode->Id<<endl;
        for(int i=0;i<pnode->vFanins.nSize;i++){
            if(check==1){
                break;
            }
            auto faninnode = Abc_NtkObj(LUT_network, pnode->vFanins.pArray[i]);
            //cout<<"ID is "<<faninnode->Id<<" dTemp is "<<faninnode->dTemp<<endl;
            switch((int)*(read+i)){
                case 49: {
                    save *= faninnode->dTemp;
                    break;
                }
                case 48: {
                    save *= inverse(faninnode->dTemp);
                    break;
                }
                case 45:
                    break;
                case 0: {                 //when it's 0, it means the pData is finished
                    check = 1;
                    break;
                }
                default: break;
            }
        }


        if(check==1){
            if(output>=1){
                pnode->dTemp = 1;
            }
            else pnode->dTemp = 0;
            if(inverseornot){
                pnode->dTemp = inverse(pnode->dTemp);
            }
            //pnode->iTemp+=pnode->dTemp;
            break;
        }
        if((int)*(read+pnode->vFanins.nSize+1)==48){
            //if the number represent the output is 0,
            // that means you should inverse the result
            output+=save;
            inverseornot=1;
        }
        else if((int)*(read+pnode->vFanins.nSize+1)==49){
            output+=save;
            inverseornot=0;
        }
        read+=pnode->vFanins.nSize+3;
    }
}

//simulate all the network, input_row and input_column
// is the integer that represent the input case of free set and bound set
int simulation(Abc_Ntk_t* LUT_network, int* bound_set, \
        int* free_set, int bound_size, int free_size, int input_row, int input_column, int min_level){
    int level_max=Abc_NtkLevel(LUT_network);
    //record the whole level of the network
    int output;
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    //base on the level of the network, I use a 2d Vector to store
    // all the nodes with different level
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int j;
    Abc_NtkForEachNode(LUT_network, pNode, j){
        int level=pNode->Level-1;
        if(level>=0) {
            Vec_PtrPush(Nodes[level], pNode);
        }
    }
    int binary_row[free_size];
    //change the integer into the real input case
    binary(binary_row,free_size, input_row);
    int binary_column[bound_size];
    binary(binary_column,bound_size, input_column);
    for(int n=0; n<bound_size; n++){
        //set the dTemp of the PI to be same as the input case
        auto pnode=Abc_NtkObj(LUT_network, bound_set[n]);
        pnode->dTemp=binary_column[n];
        //cout<<"Test: for"<<pnode->Id<<" is "<<binary_column[n]<<endl;
    }
    //cout<<"END"<<endl;
    for(int n=0; n<free_size; n++){
        auto pnode=Abc_NtkObj(LUT_network, free_set[n]);
        pnode->dTemp=binary_row[n];
    }
    //cout<<"The min_level is "<<min_level<<endl;
    for(int n=min_level; n<level_max; n++){
        //simulate all the nodes base on their levels, we should simulate the low-level nodes first, because the high-level nodes may take the low-level nodes as the Fanins
        int size=Nodes[n]->nSize;
        for(int i=0; i<size; i++){
            Abc_Obj_t* output_node=(Abc_Obj_t*)Nodes[n]->pArray[i];
            //cout<<"ID is "<<output_node->vFanins.pArray[0]<<" and dtemp is "<<output_node->dTemp<<endl;
            simulate_node(LUT_network, output_node);
        }
    }
    Abc_Obj_t* tempnode = (Abc_Obj_t*)(LUT_network->vPos->pArray[0]);
    Abc_Obj_t* outnode = Abc_NtkObj(LUT_network, tempnode->vFanins.pArray[0]);
    output = outnode->dTemp;
    //output=((Abc_Obj_t*)Nodes[level_max-1]->pArray[0])->dTemp;
    //since we only concern about he MFFC, so for the MFFC circuit,
    // there will be only one output, that is the node with the highest level,
    // and the dTemp of the node is the output of the whole circuit
    for(int k=0;k<level_max;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
    return output;
}



decomposed_chart* create_decomposed_chart(Abc_Ntk_t* LUT_network, int* bound_set,\
        int* free_set, int num, int min_level, int bound_size){
    int row=pow(2,num-bound_size);
    int column=pow(2,bound_size);
    int** bin_arr=new int*[column];
    for(int n=0; n<column; n++){
        bin_arr[n]=new int[row];
    }
    for(int i=0; i<column; i++){
        for(int j=0; j<row; j++){
            bin_arr[i][j]=simulation(LUT_network, bound_set, free_set, bound_size, num-bound_size, j, i, min_level);
        }
    }
    decomposed_chart* decomposedChart =
             new decomposed_chart(bin_arr, num-bound_size, bound_size, free_set, bound_set);
    return decomposedChart;
}

int** create_dc(Abc_Ntk_t* LUT_network, int* bound_set, int* free_set, int num, int min_level, int bound_size){
    int row=pow(2,num-bound_size);
    int column=pow(2,bound_size);
    int** bin_arr=new int*[column];
    for(int n=0; n<column; n++){
        bin_arr[n]=new int[row];
    }
    for(int i=0; i<column; i++){
        for(int j=0; j<row; j++){
            bin_arr[i][j]=simulation(LUT_network, bound_set, free_set, bound_size, num-bound_size, j, i, min_level);
        }
    }
    return bin_arr;
}

void simulate_weight(decomposed_chart* decomposedChart, Abc_Ntk_t* MFFC_network, vector<Abc_Obj_t*> input){
    int row[decomposedChart->row];
    int column[decomposedChart->column];
    int x, y;
    Abc_Obj_t* PNode;
    for(int m=0; m<decomposedChart->row; m++){
        row[m]=input[decomposedChart->array_row[m]-1]->dTemp;
    }

    for(int m=0; m<decomposedChart->column; m++){
        column[m]=input[decomposedChart->array_column[m]-1]->dTemp;
    }

    x=binary_inverse(column, decomposedChart->column);
    y=binary_inverse(row, decomposedChart->row);

    decomposedChart->weight[x][y]++;

}
//create a weight for the given MFFC circuit and the old LUT circuit,
// remember the weight is stored in the class decomposed_chart,
// the weight is a 2d double array, just like the decomposed chart
//and you can use print_weight to print it
void create_weight(decomposed_chart* decomposedChart, Abc_Ntk_t* MFFC_network, Abc_Ntk_t* LUT_network, vector<Abc_Obj_t*> input){
    srand(time(NULL));
    decomposedChart->clear_weight();
    int inputsize=LUT_network->vPis->nSize;
    int max_level=0;
    int count=100000;
    Abc_Obj_t* pnode;
    Abc_Obj_t* pin;
    for(int i=0; i<input.size(); i++){
        pnode=Abc_NtkObj(LUT_network,input[i]->Id);
        if(pnode->Level>max_level){
            max_level=pnode->Level;
        }
    }
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[max_level];        //base on the level of the network, I use a 2d Vector to store all the nodes with different level
    for(int k=0;k<max_level;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int j;
    Abc_NtkForEachNode(LUT_network, pNode, j){
            int level=pNode->Level-1;
            if(level<max_level && level>=0)
                Vec_PtrPush(Nodes[level], pNode);
        }
    for(int ic=0; ic<count; ic++){
        //if(n%1000 == 0) cout << "Iteration: " << n << endl;
        int *input_condition=new int[inputsize];
        for(int u=0; u<inputsize; u++){
            input_condition[u]=rand()%2;
            //cout<<input_condition[u];
        }
        //cout<<endl;
        for(j=0; j<LUT_network->vPis->nSize; j++){
            pin=(Abc_Obj_t*)(LUT_network->vPis->pArray[j]);
            pin->dTemp=input_condition[j];
            //pin->iTemp+=pin->dTemp;
        }
        for(int n=0; n<max_level; n++){                      //simulate all the nodes base on their levels, we should simulate the low-level nodes first, because the high-level nodes may take the low-level nodes as the Fanins
            int size=Nodes[n]->nSize;
            for(int i=0; i<size; i++){
                Abc_Obj_t* output_node=(Abc_Obj_t*)Nodes[n]->pArray[i];
                simulate_node(LUT_network, output_node);
            }
        }
        simulate_weight(decomposedChart,MFFC_network,input);
        delete[] input_condition;
    }
    for(int o=0; o<decomposedChart->column_size; o++){
        for(int p=0; p<decomposedChart->row_size; p++){
            decomposedChart->weight[o][p]=(double)((decomposedChart->weight[o][p])/100000);
        }
    }
    for(int k=0;k<max_level;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
}

void decomposed_chart::print_weight() {
    for(int i=0; i<column_size; i++){
        std::cout<<" |"<<i;
    }
    std::cout<<" |"<<std::endl;
    for(int i=0; i<row_size; i++){
        for(int j=0; j<=column_size; j++){
            std::cout<<"---";
        }
        std::cout<<std::endl;
        std::cout<<i<<"| ";
        for(int j=0; j<column_size; j++){
            std::cout<<setiosflags(ios::fixed)<<setprecision(4)<<weight[j][i]<<" |";
        }
        std::cout<<std::endl;
    }
}

void decomposed_chart::print_decomposed_chart() {
    for(int i=0; i<column_size; i++){
        std::cout<<" |"<<i;
    }
    std::cout<<" |"<<std::endl;
    for(int i=0; i<row_size; i++){
        for(int j=0; j<=column_size; j++){
            std::cout<<"---";
        }
        std::cout<<std::endl;
        std::cout<<i<<"| ";
        for(int j=0; j<column_size; j++){
            std::cout<<decompose_chart[j][i]<<" |";
        }
        std::cout<<std::endl;
    }
}

void decomposed_chart::clear_weight() {
    for(int i=0; i<this->column; i++){
        for(int j=0; j<this->row; j++){
            this->weight[i][j]=0;
        }
    }
}


/*vector<vector<int>> simulate_whole(Abc_Ntk_t* LUT_network){
    srand(time(NULL));
    vector<vector<int>> array(1);
    cout<<Abc_NtkNodeNum(LUT_network)<<endl;
    cout<<Abc_NtkPiNum(LUT_network)<<endl;
    cout<<Abc_NtkPoNum(LUT_network)<<endl;
    for(int i=0; i<(Abc_NtkNodeNum(LUT_network)+Abc_NtkPiNum(LUT_network)+Abc_NtkPoNum(LUT_network)+1); i++){
        vector<int> temp;
        temp.resize(100000);
        array.push_back(temp);
    }
    //cout<<"!!!"<<array.size()<<endl;
    int level_max=Abc_NtkLevel(LUT_network);
    //record the whole level of the network
    int count=100000;
    Abc_Obj_t* pin;
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    //base on the level of the network, I use a 2d Vector to store all the nodes with different level
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    //..........................
    Abc_Obj_t* pObj;
    int j;
    int min = 99999;
    Abc_NtkForEachNode(LUT_network, pObj, j){
            int id = Abc_ObjId(pObj);
            //cout << id << " " << Abc_ObjName(pObj) << endl;
            if(id < min) min = id;
        }
    PM.offset = min;
    //..................................
    Abc_Obj_t* pNode;
    Abc_NtkForEachNode(LUT_network, pNode, j) {
            int level = pNode->Level - 1;
            if (level >= 0) {
                Vec_PtrPush(Nodes[level], pNode);
            }
        }
    for(int i=0; i<count; i++) {
        for (int k = 0; k < LUT_network->vPis->nSize; k++) {
            pin = (Abc_Obj_t *) (LUT_network->vPis->pArray[k]);
            pin->dTemp = rand() % 2;
            //*((int*)(pin->sample)+i) = pin->dTemp;
            array[pin->Id][i]=pin->dTemp;
        }
        for (int n = 0; n < level_max; n++) {
            //simulate all the nodes base on their levels, we should simulate the low-level nodes first,
            // because the high-level nodes may take the low-level nodes as the Fanins
            int size = Nodes[n]->nSize;
            for (int l = 0; l < size; l++) {
                pin = (Abc_Obj_t *) Nodes[n]->pArray[l];
                simulate_node(LUT_network, pin);
                //*((int*)(pin->sample)+i) = pin->dTemp;
                array[pin->Id][i]=pin->dTemp;
            }
        }

        Abc_NtkForEachNode(LUT_network, pObj, j) {

                int numFanout = Abc_ObjFanoutNum(pObj);
                int nodeid = Abc_ObjId(pObj);
                string name_in = Abc_ObjName(pObj);
                PM.BD[i][nodeid - PM.offset] = new int[numFanout];
                PM.Index[nodeid - PM.offset] = numFanout;
                cout << nodeid << endl;
            }

        for (int n = 0; n < level_max; n++) {
            int size = Nodes[n]->nSize;
            for (int l = 0; l < size; l++) {
                pin = (Abc_Obj_t *) Nodes[n]->pArray[l];
                int k;
                Abc_Obj_t* fanout;
                string name_in = Abc_ObjName(pin);
                Abc_ObjForEachFanout(pin, fanout, k) {
                    string name_out = Abc_ObjName(fanout);
                    if (name_in == name_out) {
                        continue;
                    }
                    int true_val = array[pin->Id][i];

                    pin->dTemp = 0;
                    simulate_node(LUT_network, fanout);
                    int result0 = (int)fanout->dTemp;

                    pin->dTemp = 1;
                    simulate_node(LUT_network, fanout);
                    int result1 = (int)fanout->dTemp;

                    int result = (result0 != result1);
                    PM.BD[i][pin->Id - PM.offset][k] = result;

                    pin->dTemp = true_val;
                }
            }
        }
        Abc_Obj_t* po;
        int k;

        Abc_NtkForEachPo(LUT_network, po, k){
            Abc_Obj_t* node = Abc_ObjFanin0(po);
            int id = Abc_ObjId(node);
            cout << i << " " << id - PM.getoffset() << endl;
            cout << PM.N << endl;
            PM.CPM[i][id-PM.offset][0] = 1;
        }
        cout << k << endl;
        cout << "--------------"<<endl;

        for (int n = level_max - 1; n >= 0; n--){
            int result = 0;
            int size = Nodes[n]->nSize;
            for (int l = 0; l < size; l++) {
                pin = (Abc_Obj_t *) Nodes[n]->pArray[l];
                string name_node = Abc_ObjName(pin);
                int pid = pin->Id;
                if (PM.CPM[i][pid - PM.offset][0] != -1)
                    continue;
                int k;
                Abc_Obj_t* fanout;
                Abc_ObjForEachFanout(pin, fanout, k) {
                    string name_out = Abc_ObjName(fanout);
                    int fanoutid = Abc_ObjId(fanout);
                    if (name_node == name_out) continue;

                    int bd = PM.BD[i][pid-PM.offset][k];
                    int cpm = PM.CPM[i][fanoutid-PM.offset][0];
                    assert(cpm != -1);
                    assert(bd != -1);
                    if(bd == 1 && cpm == 1){
                        result = 1;
                        break;
                    }
                }
                PM.CPM[i][pid-PM.offset][0] = result;
            }
        }
    }
    for(int i = 0; i < PM.M; i++){
        for(int j = 0; j < PM.N; j++){
            if(PM.Index[j] == 0){
                continue;
            }
            else delete[] PM.BD[i][j];
        }
        delete[] PM.BD[i];
    }
    delete[] PM.BD;
    for(int k=0;k<level_max;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
    return array;
}*/
vector<vector<int>> simulate_whole(Abc_Ntk_t* LUT_network){
    srand(time(NULL));
    vector<vector<int>> array(1);
    cout<<Abc_NtkNodeNum(LUT_network)<<endl;
    cout<<Abc_NtkPiNum(LUT_network)<<endl;
    cout<<Abc_NtkPoNum(LUT_network)<<endl;
    for(int i=0; i<(Abc_NtkNodeNum(LUT_network)+Abc_NtkPiNum(LUT_network)+Abc_NtkPoNum(LUT_network)+1); i++){
        vector<int> temp;
        temp.resize(100000);
        array.push_back(temp);
    }
    //cout<<"!!!"<<array.size()<<endl;
    int level_max=Abc_NtkLevel(LUT_network);
    //record the whole level of the network
    int count=100000;
    Abc_Obj_t* pin;
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    //base on the level of the network, I use a 2d Vector to store all the nodes with different level
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int j;
    Abc_NtkForEachNode(LUT_network, pNode, j) {
            int level = pNode->Level - 1;
            if (level >= 0) {
                Vec_PtrPush(Nodes[level], pNode);
            }
        }
    for(int i=0; i<count; i++) {
        for (int k = 0; k < LUT_network->vPis->nSize; k++) {
            pin = (Abc_Obj_t *) (LUT_network->vPis->pArray[k]);
            pin->dTemp = rand() % 2;
            //*((int*)(pin->sample)+i) = pin->dTemp;
            array[pin->Id][i]=pin->dTemp;
        }
        for (int n = 0; n < level_max; n++) {
            //simulate all the nodes base on their levels, we should simulate the low-level nodes first,
            // because the high-level nodes may take the low-level nodes as the Fanins
            int size = Nodes[n]->nSize;
            for (int l = 0; l < size; l++) {
                pin = (Abc_Obj_t *) Nodes[n]->pArray[l];
                simulate_node(LUT_network, pin);
                //*((int*)(pin->sample)+i) = pin->dTemp;
                array[pin->Id][i]=pin->dTemp;
            }
        }
    }
    for(int k=0;k<level_max;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
    return array;
}

/*void Create_Weight(decomposed_chart* decomposedChart, Mffc* MFFC_network, int* bound_set1, int* free_set1, vector<vector<int>>& array){
    int count=100000;
    int x,y;
    int bound_input[decomposedChart->column];
    int free_input[decomposedChart->row];
    Abc_Obj_t* pnode;
    Abc_Ntk_t* MFFC=MFFC_network->ref;
    int level_max=Abc_NtkLevel(MFFC);
    //for(int i=0; i<input.size(); i++){
        //pnode=((Abc_Obj_t*)MFFC->vPis->pArray[i]);
        //pnode->sample=input[i]->sample;
    //}
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    //base on the level of the network, I use a 2d Vector to store all the nodes with different level
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int record_weight = 0;
    int j;
    Abc_NtkForEachNode(MFFC, pNode, j) {
            int level = pNode->Level - 1;
            if (level >= 0) {
                Vec_PtrPush(Nodes[level], pNode);
            }
        }
    for(int i=0; i<count; i++){
        for(int n=0; n<MFFC_network->Pi_origin.size(); n++){
            pnode=((Abc_Obj_t*)MFFC->vPis->pArray[n]);
            //pnode->dTemp=*((int*)MFFC_network.Pi_origin[n]->sample+i);
            pnode->dTemp=array[MFFC_network->Pi_origin[n]->Id][i];
        }
        int check = 0;
        for(int ii = 0; ii < 1; ii++){
            int kkk = MFFC_network->root->Id - PM.getoffset();
            cout << kkk << endl;
            cout << i << endl;
            cout << PM.getN() << endl;
            check = PM.CPM[i][MFFC_network->root->Id - PM.getoffset()][ii];
        }
        for (int n = 0; n < level_max; n++) {
            //simulate all the nodes base on their levels, we should simulate the low-level nodes first,
            // because the high-level nodes may take the low-level nodes as the Fanins
            int size = Nodes[n]->nSize;
            for (int l = 0; l < size; l++) {
                pnode = (Abc_Obj_t *) Nodes[n]->pArray[l];
                simulate_node(MFFC, pnode);
            }
        }
        for(int j=0; j<decomposedChart->column; j++){
            pnode=Abc_NtkObj(MFFC, decomposedChart->array_column[j]);
            bound_input[j]=pnode->dTemp;
        }
        for(int j=0; j<decomposedChart->row; j++){
            pnode=Abc_NtkObj(MFFC, decomposedChart->array_row[j]);
            free_input[j]=pnode->dTemp;
        }
        x=binary_inverse(bound_input, decomposedChart->column);
        y=binary_inverse(free_input, decomposedChart->row);
        //delete[] bound_input;
        //delete[] free_input;
        if(check == 1)
            decomposedChart->weight[x][y]++;
        else record_weight ++;
    }
    cout << "how many weight is invalid" << record_weight << endl;
    if(record_weight > 0){
        cout << "bingo!" << endl;
    }
    for(int i=0; i<decomposedChart->column_size; i++){
        for(int j=0; j<decomposedChart->row_size; j++){
            decomposedChart->weight[i][j]=decomposedChart->weight[i][j]/count;
        }
    }
    for(int k=0;k<level_max;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
}*/

void Create_Weight(decomposed_chart* decomposedChart, Mffc* MFFC_network, int* bound_set1, int* free_set1, vector<vector<int>>& array){
    int count=100000;
    int x,y;
    int bound_input[decomposedChart->column];
    int free_input[decomposedChart->row];
    Abc_Obj_t* pnode;
    Abc_Ntk_t* MFFC=MFFC_network->ref;
    int level_max=Abc_NtkLevel(MFFC);
    /*for(int i=0; i<input.size(); i++){
        pnode=((Abc_Obj_t*)MFFC->vPis->pArray[i]);
        pnode->sample=input[i]->sample;
    }*/
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    //base on the level of the network, I use a 2d Vector to store all the nodes with different level
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int j;
    Abc_NtkForEachNode(MFFC, pNode, j) {
            int level = pNode->Level - 1;
            if (level >= 0) {
                Vec_PtrPush(Nodes[level], pNode);
            }
        }
    for(int i=0; i<count; i++){
        for(int n=0; n<MFFC_network->Pi_origin.size(); n++){
            pnode=((Abc_Obj_t*)MFFC->vPis->pArray[n]);
            //pnode->dTemp=*((int*)MFFC_network.Pi_origin[n]->sample+i);
            pnode->dTemp=array[MFFC_network->Pi_origin[n]->Id][i];
        }
        for (int n = 0; n < level_max; n++) {
            //simulate all the nodes base on their levels, we should simulate the low-level nodes first,
            // because the high-level nodes may take the low-level nodes as the Fanins
            int size = Nodes[n]->nSize;
            for (int l = 0; l < size; l++) {
                pnode = (Abc_Obj_t *) Nodes[n]->pArray[l];
                simulate_node(MFFC, pnode);
            }
        }
        for(int j=0; j<decomposedChart->column; j++){
            pnode=Abc_NtkObj(MFFC, decomposedChart->array_column[j]);
            bound_input[j]=pnode->dTemp;
        }
        for(int j=0; j<decomposedChart->row; j++){
            pnode=Abc_NtkObj(MFFC, decomposedChart->array_row[j]);
            free_input[j]=pnode->dTemp;
        }
        x=binary_inverse(bound_input, decomposedChart->column);
        y=binary_inverse(free_input, decomposedChart->row);
        //delete[] bound_input;
        //delete[] free_input;
        decomposedChart->weight[x][y]++;
    }
    for(int i=0; i<decomposedChart->column_size; i++){
        for(int j=0; j<decomposedChart->row_size; j++){
            decomposedChart->weight[i][j]=decomposedChart->weight[i][j]/count;
        }
    }
    for(int k=0;k<level_max;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
}

void find_new_FFC(Abc_Ntk_t* MFFC, set<int, firstcome>& root_array, set<int, firstcome>& ID_of_FFC,
                 set<int, firstcome> old_root_array, set<int, firstcome>& possible_FFC, set<int, firstcome>& visited,set<int, firstcome>& real_root, int upperbound){
    if(root_array.size() == 0)
        return;
    set<int, firstcome> candidate;
    possible_FFC.clear();
    set<int, firstcome>::iterator it;
    it = root_array.begin();
    int times = random()%(root_array.size());
    for(int i = 0; i < times; i++){
        it++;
    }
    int root_ID = *it;
    root_array.erase(it);
    old_root_array.erase(root_ID);
    //Abc_Obj_t* root = Abc_NtkObj(MFFC, *it);
    Abc_Obj_t* root = Abc_NtkObj(MFFC, root_ID);
    if(Abc_ObjIsPi(root)){
        real_root.insert(root->Id);
        return find_new_FFC(MFFC, root_array, ID_of_FFC, old_root_array, possible_FFC, visited, real_root, upperbound);
    }
    else{
        possible_FFC.insert(root->Id);
        for(int i = 0; i < root->vFanins.nSize; i ++){
            if(real_root.find(root->vFanins.pArray[i]) == real_root.end()) {
                root_array.insert(root->vFanins.pArray[i]);
            }
        }
        for(int i = 0; i < root->vFanouts.nSize; i ++){
            int temp = root->vFanouts.pArray[i];
            if(ID_of_FFC.find(temp) == ID_of_FFC.end() && possible_FFC.find(temp) == possible_FFC.end()){
                candidate.insert(temp);
                root_array.erase(temp);
            }
        }
    }
    while(!candidate.empty()){
        int num = *candidate.begin();
        candidate.erase(candidate.begin());
        possible_FFC.insert(num);
        Abc_Obj_t* pnode = Abc_NtkObj(MFFC, num);
        for(int i = 0; i < pnode->vFanins.nSize; i ++){
            if(possible_FFC.find(pnode->vFanins.pArray[i]) == possible_FFC.end() && candidate.find(pnode->vFanins.pArray[i]) == candidate.end()){
                root_array.insert(pnode->vFanins.pArray[i]);
            }
        }
        for(int i = 0; i < pnode->vFanouts.nSize; i ++){
            int temp1 = pnode->vFanouts.pArray[i];
            if(ID_of_FFC.find(temp1) == ID_of_FFC.end() && possible_FFC.find(temp1) == possible_FFC.end()){
                candidate.insert(temp1);
                root_array.erase(temp1);
            }
        }
    }
    set<int, firstcome> temp_root;
    temp_root.insert(root_array.begin(), root_array.end());
    temp_root.insert(real_root.begin(), real_root.end());
    if(temp_root.size() == upperbound) {
        ID_of_FFC.insert(possible_FFC.begin(), possible_FFC.end());
        real_root.insert(root_array.begin(), root_array.end());
        return;
    }
    else if(temp_root.size() > upperbound) {
        real_root.insert(root_ID);
        root_array.clear();
        root_array.insert(old_root_array.begin(), old_root_array.end());
        return find_new_FFC(MFFC, root_array, ID_of_FFC, old_root_array, possible_FFC, visited, real_root, upperbound);
    }
    else{
        ID_of_FFC.insert(possible_FFC.begin(), possible_FFC.end());
        old_root_array.clear();
        old_root_array.insert(root_array.begin(), root_array.end());
        return find_new_FFC(MFFC, root_array, ID_of_FFC, old_root_array, possible_FFC, visited, real_root, upperbound);
    }
}

Mffc* create_ffc_new(Mffc* origin_MFFC, Abc_Ntk_t* origin_ntk, int upperbound){
    set<int, firstcome> root_array;
    set<int, firstcome> ID_of_FFC;
    set<int, firstcome> possible;
    set<int, firstcome> old_root_array;
    set<int, firstcome> real;
    set<int, firstcome> visited;
    set<int, firstcome>::iterator it;
    int begin = ((Abc_Obj_t*)origin_MFFC->ref->vPos->pArray[0])->vFanins.pArray[0];
    Abc_Obj_t* node = Abc_NtkObj(origin_MFFC->ref, begin);
    for(int i = 0; i < node->vFanins.nSize; i++) {
        old_root_array.insert(node->vFanins.pArray[i]);
        root_array.insert(node->vFanins.pArray[i]);
    }
    ID_of_FFC.insert(begin);
    find_new_FFC(origin_MFFC->ref, root_array, ID_of_FFC, old_root_array, possible, visited, real, upperbound);
    Abc_Ntk_t* new_ffc = Abc_NtkCreateffc(origin_MFFC->ref, node, Abc_ObjName(node), ID_of_FFC, real);
    Mffc* new_FFC = new Mffc;
    Abc_Obj_t* temp;
    new_FFC->ref = new_ffc;
    new_FFC->root = origin_MFFC->root;
    for(it = real.begin(); it != real.end(); it++){
        int ID = origin_MFFC->ID_origin[*it-1];
        temp = Abc_NtkObj(origin_ntk, ID);
        new_FFC->Pi_origin.push_back(temp);
        new_FFC->ID_origin.push_back(ID);
    }
    new_FFC->ID_origin.push_back(0);
    for(it = ID_of_FFC.begin(); it != ID_of_FFC.end(); it++){
        int id = origin_MFFC->ID_origin[*it-1];
        new_FFC->ID_origin.push_back(id);
    }
    return new_FFC;
}

Abc_Ntk_t * Abc_NtkCreateffc( Abc_Ntk_t * pNtk, Abc_Obj_t * pNode, char * pNodeName, set<int, firstcome>& ID_of_FFC,  set<int, firstcome>& root_array )
{
    Abc_Ntk_t * pNtkNew;
    Abc_Obj_t * pObj, * pFanin, * pNodeCoNew;
    Vec_Ptr_t * vCone, * vSupp;
    set<int, weight_compare>::iterator it;
    char Buffer[1000];
    int i, k;

    assert( Abc_NtkIsLogic(pNtk) || Abc_NtkIsStrash(pNtk) );
    assert( Abc_ObjIsNode(pNode) );

    // start the network
    pNtkNew = Abc_NtkAlloc( pNtk->ntkType, pNtk->ntkFunc, 1 );
    // set the name
    sprintf( Buffer, "%s_%s", pNtk->pName, pNodeName );
    pNtkNew->pName = Extra_UtilStrsav(Buffer);

    // establish connection between the constant nodes
    if ( Abc_NtkIsStrash(pNtk) )
        Abc_AigConst1(pNtk)->pCopy = Abc_AigConst1(pNtkNew);

    // collect the nodes in MFFC
    vCone = Vec_PtrAlloc( 100 );
    vSupp = Vec_PtrAlloc( 100 );

    for(it = ID_of_FFC.begin(); it != ID_of_FFC.end(); it++){
        Abc_Obj_t* pnode = Abc_NtkObj(pNtk, *it);
        Vec_PtrPush(vCone, pnode);
        //cout<<"it"<<*it<<endl;
    }

    for(it = root_array.begin(); it != root_array.end(); it++){
        Abc_Obj_t* pnode = Abc_NtkObj(pNtk, *it);
        Vec_PtrPush(vSupp, pnode);
    }
    /*Abc_NodeDeref_rec( pNode );
    Abc_NodeMffcConeSupp( pNode, vCone, vSupp );
    Abc_NodeRef_rec( pNode );*/
    // create the PO
    pNodeCoNew = Abc_NtkCreatePo( pNtkNew );
    Abc_ObjAssignName( pNodeCoNew, pNodeName, NULL );
    // create the PIs
    Vec_PtrForEachEntry( Abc_Obj_t *, vSupp, pObj, i )
    {
        pObj->pCopy = Abc_NtkCreatePi(pNtkNew);
        Abc_ObjAssignName( pObj->pCopy, Abc_ObjName(pObj), NULL );
    }
    // copy the nodes
    Vec_PtrForEachEntry( Abc_Obj_t *, vCone, pObj, i )
    {
        // if it is an AIG, add to the hash table
        if ( Abc_NtkIsStrash(pNtk) )
        {
            pObj->pCopy = Abc_AigAnd( (Abc_Aig_t *)pNtkNew->pManFunc, Abc_ObjChild0Copy(pObj), Abc_ObjChild1Copy(pObj) );
        }
        else
        {
            Abc_NtkDupObj( pNtkNew, pObj, 0 );
            Abc_ObjForEachFanin( pObj, pFanin, k )
                Abc_ObjAddFanin( pObj->pCopy, pFanin->pCopy );
        }
    }
    // connect the topmost node
    Abc_ObjAddFanin( pNodeCoNew, pNode->pCopy );
    Vec_PtrFree( vCone );
    Vec_PtrFree( vSupp );

    if ( !Abc_NtkCheck( pNtkNew ) )
        fprintf( stdout, "Abc_NtkCreateMffc(): Network check has failed.\n" );
    return pNtkNew;
}

bool compare_PI_new(Mffc* one, Mffc* two){
    bool answer = false;
    set<int, firstcome> one_PI;
    set<int, firstcome> two_PI;
    for(int i = 0; i < one->Pi_origin.size(); i ++){
        one_PI.insert(one->Pi_origin[i]->Id);
    }
    for(int i = 0; i < two->Pi_origin.size(); i ++){
        if(one_PI.find(two->Pi_origin[i]->Id) != one_PI.end()){
            answer = true;
        }
        else{
            answer = false;
            break;
        }
    }
    return answer;
}

void create_ffc_set_new(Mffc* origin_FFC, Abc_Ntk_t* origin_ntk, int upperbound, int lowerbound, set<int, firstcome>& changed, vector<vector<Mffc*>>& ffc_set){
    srand(time(NULL));
    bool repeat_or_not = false;
    for(int i = lowerbound; i <= upperbound; i++){
        int round = pow(2, i);
        if(round >= 250) round = 250;
        for(int j = 0; j < round; j++){
            Mffc* temp = create_ffc_new(origin_FFC, origin_ntk, i);
            if(temp->Pi_origin.size() < lowerbound) {
                delete temp;
                continue;
            }
            else{
                Abc_Obj_t* tempnode;
                bool null_or_not = false;
                /*for(int l = 0; l < temp->Pi_origin.size(); l ++){
                    tempnode = temp->Pi_origin[l];
                    if(tempnode == NULL) null_or_not = true;
                }
                for(int l = temp->Pi_origin.size()+1; l < temp->ID_origin.size(); l++){
                    tempnode = Abc_NtkObj(origin_ntk, l);
                    if(tempnode == NULL) null_or_not = true;
                }*/
                for(int l = 0; l < temp->ID_origin.size(); l++){
                    if(changed.find(temp->ID_origin[l]) != changed.end()){
                        if(temp->ID_origin[l] != 0) null_or_not = true;
                    }
                    else{
                        if(temp->ID_origin[l] != 0){
                            Abc_Obj_t* checknode = Abc_NtkObj(origin_ntk, temp->ID_origin[l]);
                            if(checknode == NULL) null_or_not = true;
                        }
                    }
                }
                if(null_or_not == true) {
                    delete temp;
                    continue;
                }
            }
            int ID = temp->Pi_origin.size();
            if(ffc_set[ID-lowerbound].size() == 0){
                ffc_set[ID-lowerbound].push_back(temp);
            }
            else{
                for(int k = 0; k < ffc_set[ID-lowerbound].size(); k ++){
                    bool contemp = compare_PI_new(ffc_set[ID-lowerbound][k],temp);
                    if(contemp){
                        repeat_or_not = true;
                        break;
                    }
                }
                if(repeat_or_not == false){
                    ffc_set[ID-lowerbound].push_back(temp);
                }
                else{
                    delete temp;
                }
                repeat_or_not = false;
            }
        }
    }
    for(int i = 0; i < ffc_set.size(); i++){
        for(int j = 0; j < ffc_set[i].size(); j ++){
            Abc_Ntk_t* temp = ffc_set[i][j]->ref;
            ffc_set[i][j]->ref = Abc_NtkDup(ffc_set[i][j]->ref);
            Abc_NtkDelete(temp);
        }
    }
    return;
}

double constant_0(Mffc* testcase, vector<vector<int>>& sim_result_ptr, int bound_size){
    auto sim_result = sim_result_ptr;
    int pisize = testcase->Pi_origin.size();
    int* bound_set = new int[bound_size];
    int* free_set = new int[pisize - bound_size];
    for(int i = 0; i < bound_size; i++){
        bound_set[i] = i + 1;
    }
    for(int i = 0; i < pisize - bound_size; i ++){
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


double constant_zero(Abc_Ntk_t* origin_Ntk, int upperbound, int lowerbound, vector<vector<int>>& sim_result, set<int, firstcome>& changed, set<int>& critical_path, double& upperbound_error){
    vector<Mffc*> mffc_set = getMffcNtk(origin_Ntk, lowerbound, upperbound);
    std::sort (mffc_set.begin(), mffc_set.end(), cpr_rl_t);
    Mffc* mffc = nullptr;
    Mffc* temp;
    double error = 0;
    double final_error = 0;
    vector<Mffc*>::size_type ix = 0;
    //int ix = mffc_set.size() - 1;
    vector<int> record;
    for(;ix != mffc_set.size(); ix++){
        mffc = mffc_set[ix];
        if(critical_path.find(mffc->root->Id) != critical_path.end()){
            cout<<"the root of MFFC is on the critical path, be careful!" << endl;
        }
        bool null_or_not = false;
        for(int l = 0; l < mffc->ID_origin.size(); l++){
            if(changed.find(mffc->ID_origin[l]) != changed.end()){
                if(mffc->ID_origin[l] != 0) null_or_not = true;
            }
            else{
                if(mffc->ID_origin[l] != 0){
                    Abc_Obj_t* checknode = Abc_NtkObj(origin_Ntk, mffc->ID_origin[l]);
                    if(checknode == NULL) null_or_not = true;
                }
            }
        }
        if(null_or_not) continue;
        vector<vector<Mffc*>> ffc_set;
        for(int m = 0; m <= (12-4); m++){
            vector<Mffc*> candidate;
            ffc_set.push_back(candidate);
        }
        if(mffc->Pi_origin.size() < 12) {
            error = constant_0(mffc, sim_result);
            if(error > upperbound_error) continue;
            cout << "the improvement is " << Abc_NtkNodeNum(mffc->ref) << endl;
            int level_origin = Abc_NtkLevel(mffc->ref);
            double ratio = (double) error/Abc_NtkNodeNum(mffc->ref);
            if (error < 0.05) {
                cout << "do you want to change it? if so, enter 1" << endl;
                int saveornot = 0;
                if(ratio < 0.005) saveornot = 1;
                if (saveornot == 1) {
                    modify_constant0(mffc->ref);
                    for(int i = 0; i < mffc->ID_origin.size(); i ++){
                        changed.insert(mffc->ID_origin[i]);
                    }
                    Io_WriteBlifLogic(mffc->ref, (char*)"new_mffc.blif", 1);
                    write_back(mffc, origin_Ntk);
                    Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
                    upperbound_error -= error;
                    final_error += error;
                }
            }
        }
        else{
            create_ffc_set_new(mffc, origin_Ntk, 12, 4, changed, ffc_set);
            for(int i = 0; i < ffc_set.size(); i ++){
                if (!ffc_set[i].empty()) temp = ffc_set[i][0];
            }
            error = constant_0(temp, sim_result);
            if(error > upperbound_error) continue;
            cout << "the improvement is " << Abc_NtkNodeNum(mffc->ref) << endl;
            int level_origin = Abc_NtkLevel(mffc->ref);
            double ratio = (double) error/Abc_NtkNodeNum(mffc->ref);
            if (error < 0.05) {
                cout << "do you want to change it? if so, enter 1" << endl;
                int saveornot = 0;
                if(ratio < 0.005) saveornot = 1;
                //cin >> saveornot;
                if (saveornot == 1) {
                    modify_constant0(mffc->ref);
                    for(int i = 0; i < mffc->ID_origin.size(); i ++){
                        changed.insert(mffc->ID_origin[i]);
                    }
                    Io_WriteBlifLogic(mffc->ref, (char*)"new_mffc.blif", 1);
                    write_back(mffc, origin_Ntk);
                    Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
                    upperbound_error -= error;
                    final_error += error;
                }
            }
        }
    }
    Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    for(int jj = 0; jj < mffc_set.size(); jj ++){
        delete mffc_set[jj];
    }
    return final_error;
}


int check_level(Mffc* new_mffc){
    Abc_Obj_t* pNode;
    Abc_Obj_t* temp;
    int i, child, level_now;
    level_now = 0;
    map<int, int> level_record;
    map<int, int>::iterator iter;
    for(int j = 0; j < new_mffc->Pi_origin.size(); j ++){
        temp = new_mffc->Pi_origin[j];
        int level_temp = temp->Level;
        level_record.insert(map<int, int>::value_type (j+1, level_temp));
        //cout<<level_temp<<endl;
    }
    queue<int> waiting_list;
    for(int j = 0; j < new_mffc->Pi_origin.size(); j ++){
        waiting_list.push(j + 1);
    }
    while(!waiting_list.empty()){
        int check_point = waiting_list.front();
        waiting_list.pop();
        pNode = Abc_NtkObj(new_mffc->ref, check_point);
        for(int j = 0; j < pNode->vFanouts.nSize; j ++){
            child = pNode->vFanouts.pArray[j];
            iter = level_record.find(child);
            if(iter == level_record.end()) level_record.insert(map<int, int>::value_type (child, level_record[check_point]+1));
            else if((*iter).second < level_record[check_point] + 1) level_record[child] = level_record[check_point] + 1;
            //if(level_record[child-1] < level_record[check_point-1] + 1) level_record[child-1] = level_record[check_point-1] + 1;
            waiting_list.push(child);
            if(level_record[check_point] + 1 > level_now) level_now = level_record[check_point] + 1;
        }
    }
    level_record.clear();
    return level_now - 1;
}

double create_map(Abc_Ntk_t* origin_Ntk, vector<vector<int>>& sim_result, map<int, Mffc*>& best_choise, vector<vector<Mffc *>>& ffc_set, set<int>& critical_path
                ,set<int, firstcome>& changed, double upperbound_error){
    char* old = (char*) "origin_mffc.blif";
    char* newfilename = (char*) "new_mffc.blif";
    Mffc* new_mffc;
    //test code
    Mffc* temp_mffc;
    Mffc* origin_one;
    bool reduce_or_not = false;
    //test code
    for (int h = 0; h < ffc_set.size(); h++) {
        //test code
        if(reduce_or_not) break;
        //test code
        for (int o = 0; o < ffc_set[h].size(); o++) {
            //test code
            if(reduce_or_not) break;
            //test code
            int num = Abc_NtkPiNum(ffc_set[h][o]->ref);
            int optimal = (num - 1) % 3 == 0 ? (num - 1) / 3 : (num - 1) / 3 + 1;
            //reduce_input(ffc_set[h][o], sim_result);
            int* boundset = new int[ffc_set[h][o]->Pi_origin.size()-1];
            int size = pow(2, ffc_set[h][o]->Pi_origin.size()-1);
            int* pattern = new int[size];
            if (Abc_NtkNodeNum(ffc_set[h][o]->ref) <= optimal) {
                cout << "this ffc with input size " << ffc_set[h][o]->Pi_origin.size()<<" is perfect, we cannot change it" << endl;
                //origin_one = ffc_set[h][o]->duplicate();
                temp_mffc = ffc_set[h][o]->duplicate();
                double temp_error, final_error = 0;
                if(reduce_input(temp_mffc, sim_result, pattern, boundset, temp_error) != -1) {
                    int count = 0;
                    while(reduce_input(temp_mffc, sim_result, pattern, boundset, temp_error) != -1) {
                        count ++;
                        int reduce_ID = reduce_input(temp_mffc, sim_result, pattern, boundset, temp_error);
                        final_error += temp_error;
                        temp_mffc->RI_write_back(reduce_ID, pattern);
                        cout<<"now the count equals to "<< count <<endl;
                    }
                    int num_now = temp_mffc->Pi_origin.size();
                    int new_optimal = (num_now - 1) % 3 == 0 ? (num_now - 1) / 3 : (num_now - 1) / 3 + 1;
                    bool better = false;
                    if(new_optimal < optimal) better = true;
                    if(Abc_NtkPiNum(temp_mffc->ref) > 4) {
                        new_mffc = local_approx(origin_Ntk, temp_mffc, Abc_NtkPiNum(temp_mffc->ref), sim_result);
                        final_error += new_mffc->error;
                    }
                    else new_mffc = temp_mffc;
                    cout << new_mffc->error << endl;
                    if (new_mffc->error <= 0.005 && final_error <= 0.01 && final_error < upperbound_error + 0.002) {
                        new_mffc->level_reduce = check_level(ffc_set[h][o]) - check_level(new_mffc);
                        cout << "the level reduce is " << new_mffc->level_reduce << endl;
                        if (new_mffc->level_reduce >= 0 && better) {
                            //test code
                            reduce_or_not = true;
                            for (int l = 0; l < new_mffc->ID_origin.size(); l++) {
                                changed.insert(new_mffc->ID_origin[l]);
                            }
                            //test code
                            Io_WriteBlifLogic(new_mffc->ref, old, 1);
                            write_back(new_mffc, origin_Ntk);
                            Io_WriteBlifLogic(origin_Ntk, newfilename, 1);
                            best_choise.clear();
                            //delete new_mffc;
                            //delete origin_one;
                            delete[] boundset;
                            delete[] pattern;
                            return final_error;
                        }
                    }
                    /*else{
                        Io_WriteBlifLogic(origin_one->ref, old, 1);
                        write_back(origin_one, origin_Ntk);
                        Io_WriteBlifLogic(origin_Ntk, newfilename, 1);
                    }*/
                    //Io_WriteBlifLogic(ffc_set[h][o]->ref, old, 1);
                    //write_back(ffc_set[h][o], origin_Ntk);
                    //Io_WriteBlifLogic(origin_Ntk, newfilename, 1);
                }
                //delete temp_mffc;
                //delete origin_one;
                delete[] boundset;
                delete[] pattern;
                //continue;
            }
            else if(true){
                //origin_one = ffc_set[h][o]->duplicate();
                temp_mffc = ffc_set[h][o]->duplicate();
                double temp_error, final_error = 0;
                //test code
                if(reduce_input(temp_mffc, sim_result, pattern, boundset, temp_error) != -1) {
                    int count = 0;
                    while (reduce_input(temp_mffc, sim_result, pattern, boundset, temp_error) != -1) {
                        count++;
                        int reduce_ID = reduce_input(temp_mffc, sim_result, pattern, boundset, temp_error);
                        final_error += temp_error;
                        temp_mffc->RI_write_back(reduce_ID, pattern);
                        cout << "now the count equals to " << count << endl;
                    }
                    if (Abc_NtkPiNum(temp_mffc->ref) > 4) {
                        new_mffc = local_approx(origin_Ntk, temp_mffc, Abc_NtkPiNum(temp_mffc->ref), sim_result);
                        final_error += new_mffc->error;
                    } else new_mffc = temp_mffc;
                    cout << new_mffc->error << endl;
                    if (new_mffc->error <= 0.005 && final_error <= 0.01 && final_error < upperbound_error + 0.002) {
                        new_mffc->level_reduce = check_level(ffc_set[h][o]) - check_level(new_mffc);
                        cout << "the level reduce is " << new_mffc->level_reduce << endl;
                        if (new_mffc->level_reduce >= 0) {
                            //test code
                            reduce_or_not = true;
                            for (int l = 0; l < new_mffc->ID_origin.size(); l++) {
                                changed.insert(new_mffc->ID_origin[l]);
                            }
                            //test code
                            Io_WriteBlifLogic(new_mffc->ref, old, 1);
                            write_back(new_mffc, origin_Ntk);
                            Io_WriteBlifLogic(origin_Ntk, newfilename, 1);
                            best_choise.clear();
                            //delete new_mffc;
                            //delete origin_one;
                            delete[] boundset;
                            delete[] pattern;
                            return final_error;
                        }
                    }
                    //continue;
                }
                //delete origin_one;
                //delete temp_mffc;
                delete[] boundset;
                delete[] pattern;
                //continue;
            }
            int improve = Abc_NtkNodeNum(ffc_set[h][o]->ref) - optimal;
            Io_WriteBlifLogic(ffc_set[h][o]->ref, old, 1);
            int level_origin = check_level(ffc_set[h][o]);
            int pis = ffc_set[h][o]->Pi_origin.size();
            clock_t time = clock();
            new_mffc = local_approx(origin_Ntk, ffc_set[h][o], num, sim_result);
            time = clock() - time;
            cout << "the PI size is " << pis << " the local approxiamte time is "
                 << (double) time / CLOCKS_PER_SEC << "s" << endl;
            //int level_new = Abc_NtkLevel(new_mffc->ref);
            int level_new = check_level(new_mffc);
            cout<<"the original level is "<< level_origin << " and now the level is " << level_new<< endl;
            new_mffc->level_reduce = level_origin - level_new;
            //Io_WriteBlifLogic(new_mffc->ref, newfilename, 1);
            //int jjj = 0;
            //Io_WriteBlifLogic(ffc_set[h][o]->ref, newfilename, 1);
            map<int, Mffc *>::iterator it;
            //Io_WriteBlifLogic(ffc_set[0][15]->ref, newfilename, 1);
            if (best_choise.find(improve) != best_choise.end()) {
                double error_test = (*best_choise.find(improve)).second->error;
                int level_test = (*best_choise.find(improve)).second->level_reduce;
                if (new_mffc->error < error_test) {
                    it = best_choise.find(improve);
                    delete (*it).second;
                    best_choise.erase(it);
                    Mffc *duplicate = new_mffc->duplicate();
                    best_choise.insert(pair<int, Mffc *>(improve, duplicate));
                    Io_WriteBlifLogic(best_choise[improve]->ref, newfilename, 1);
                } else if (new_mffc->error == error_test && new_mffc->level_reduce > level_test) {
                    it = best_choise.find(improve);
                    delete (*it).second;
                    best_choise.erase(it);
                    Mffc *duplicate = new_mffc->duplicate();
                    best_choise.insert(pair<int, Mffc *>(improve, duplicate));
                    Io_WriteBlifLogic(best_choise[improve]->ref, newfilename, 1);
                } else {
                    delete new_mffc;
                }
            } else {
                Mffc *duplicate = new_mffc->duplicate();
                best_choise.insert(pair<int, Mffc *>(improve, duplicate));
            }
        }
    }
    return 0;
}

double not_zero(Abc_Ntk_t* origin_Ntk, int upperbound, int lowerbound, double upperbound_error, vector<vector<int>>& sim_result, set<int, firstcome>& changed, set<int>& critical_path){
    vector<Mffc*> mffc_set = getMffcNtk(origin_Ntk, lowerbound, upperbound);
    cout<<mffc_set.size()<<endl;
    std::sort (mffc_set.begin(), mffc_set.end(), cpr_rl_t);
    Mffc* mffc = nullptr;
    double error = 0;
    int wb = 0;
    bool critical = false;
    int non_critical = 0;
    //Mffc* new_mffc;
    vector<Mffc*>::size_type ix = 0;
    char* old = (char*) "origin_mffc.blif";
    char* newfilename = (char*) "new_mffc.blif";
    for(;ix != mffc_set.size(); ix++) {
        //test 10/15
        //MFFC_count++;
        //test 10/15
        mffc = mffc_set[ix];
        int rootID = mffc->root->Id;
        if(critical_path.find(mffc->root->Id) != critical_path.end()){
            critical = true;
        }
        vector<vector<Mffc*>> ffc_set;
        if (mffc->Pi_origin.size() >= 12) {
            for(int m = 0; m <= (12-6); m++){
                vector<Mffc*> candidate;
                ffc_set.push_back(candidate);
            }
            create_ffc_set_new(mffc, origin_Ntk, 12, 6, changed, ffc_set);
        } else if (mffc->Pi_origin.size() < 12 && mffc->Pi_origin.size() >= 6) {
            for(int m = 0; m <= (mffc->Pi_origin.size()-6); m++){
                vector<Mffc*> candidate;
                ffc_set.push_back(candidate);
            }
            create_ffc_set_new(mffc, origin_Ntk, mffc->Pi_origin.size(), 6, changed, ffc_set);
            bool null_or_not = false;
            for(int l = 0; l < mffc->ID_origin.size(); l++){
                if(changed.find(mffc->ID_origin[l]) != changed.end()){
                    if(mffc->ID_origin[l] != 0) null_or_not = true;
                }
                /*else{
                    if(mffc->ID_origin[l] != 0){
                        Abc_Obj_t* checknode = Abc_NtkObj(origin_Ntk, mffc->ID_origin[l]);
                        if(checknode == NULL) null_or_not = true;
                    }
                }*/
            }
            if(null_or_not == true) {
                //delete temp;
                continue;
            }
            if(!null_or_not) ffc_set[mffc->Pi_origin.size() - 6].push_back(mffc);
        }
        //testcode 10/15
        /*bool have_tested = 0;
        for(int i = 0; i < ffc_set.size(); i++) {
            FFC_count += ffc_set[i].size();
            if(ffc_set[i].size() != 0 && have_tested == 0){
                for(int j = 0; j < ffc_set[i][0]->ID_origin.size(); j ++){
                    area.insert(ffc_set[i][0]->ID_origin[j]);
                }
            }
            cout<<ffc_set[i].size()<<endl;
        }
        cout<<"one MFFC ends "<<MFFC_count<<" "<<FFC_count<<" "<<area.size()<<endl;
        continue;*/
        //testcode 10/15
        map<int, Mffc *> best_choise;
        error += create_map(origin_Ntk, sim_result, best_choise, ffc_set, critical_path, changed, upperbound_error);
        map<int, Mffc *>::iterator it;
        bool skip = false;
        if(critical) cout<<"the root of MFFC is on the critical path, be careful! root ID is "<< rootID << endl;
        int choisen = 0;
        for (it = best_choise.begin(); it != best_choise.end(); it++) {
            cout << "The improve is " << (*it).first << " and the error is " << (*it).second->error
                 << " and the level reduction is " << (*it).second->level_reduce << endl;
            if ((*it).second->error >= 0.1 || error + (*it).second->error > upperbound + 0.004)
                continue;
                //skip = true;
            else if((*it).second->level_reduce >= 0){
                double ratio = (*it).second->error/(*it).first;
                cout<<"the ratio is "<< ratio << endl;
                if(ratio <= 0.01 && (*it).first > choisen) choisen = (*it).first;
                cout<<"now the choisen is " << choisen << endl;
            }
            else if(!critical && (*it).second->level_reduce >= -2){
                double ratio = (*it).second->error/(*it).first;
                cout<<"the ratio is "<< ratio << endl;
                if(ratio <= 0.01 && (*it).first > choisen) choisen = (*it).first;
                cout<<"now the choisen is " << choisen << endl;
            }
        }
        if (!best_choise.empty())
            cout << "the name of root is " << Abc_ObjName((*best_choise.begin()).second->root) << endl;
        //if (skip) continue;
        if (best_choise.empty()) continue;
        //if(Abc_ObjIsPo((*best_choise.begin()).second->root)) cout<<"the root is PO, be careful"<<endl;
        // int check = 0;
        // while (check != 0) {
        //     cout
        //             << "Please choose the ffc you want to have a look, enter the improve of the ffc, if none or you want to continue, enter 0:"
        //             << endl;
        //     cin >> check;
        //     if (check == 0) break;
        //     Io_WriteBlifLogic(best_choise[check]->ref, newfilename, 1);
        // }
        cout
                << "now please choose the ffc you want to change, please enter the improve of the ffc you want, if none, please enter 0:"
                << endl;
        int choise = choisen;
        cout<<"choisen equal to "<<choisen<<endl;
        if (choise == 0) continue;
        if(error + best_choise[choise]->error > upperbound_error + 0.004) continue;
        if(best_choise[choise]->level_reduce < 0){
            if(non_critical - best_choise[choise]->level_reduce > 2) continue;
            non_critical -= best_choise[choise]->level_reduce;
            cout<<"the root is not on critical path, so we write it back, and now the level increase" << non_critical <<endl;
        }
        for (int l = 0; l < best_choise[choise]->ID_origin.size(); l++) {
            changed.insert(best_choise[choise]->ID_origin[l]);
        }
        Io_WriteBlifLogic(best_choise[choise]->ref, newfilename, 1);
        cout<<check_level(best_choise[choise])<<endl;
        //Io_WriteBlifLogic(fuckyou, old, 1);
        write_back(best_choise[choise], origin_Ntk);
        wb++;
        error += best_choise[choise]->error;
        for (it = best_choise.begin(); it != best_choise.end(); it++) {
            delete (*it).second;
        }
        Io_WriteBlifLogic(origin_Ntk, (char *) "out_whole.blif", 1);
        cout << "do you still want to continue or not? Enter 'N' to stop:" << endl;
        char answer = 'G';
        //cin >> answer;
        //if (answer == 'N') break;
        best_choise.clear();
        if(error >= upperbound_error) break;
    }
    cout << "Write back number: " << wb << endl;
    cout << "Accumulative error: " << error << endl;
    return error;
}

void find_critical_path(Abc_Ntk_t* origin_Ntk, set<int>& critical, int level_max){
    if(level_max == 0) level_max = Abc_NtkLevel(origin_Ntk);
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    Abc_Obj_t* pnode;
    Abc_Obj_t* parent;
    queue<int> candidate;
    int ID_temp;
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int j;
    Abc_NtkForEachNode(origin_Ntk, pNode, j) {
            int level = pNode->Level - 1;
            if (level >= 0) {
                Vec_PtrPush(Nodes[level], pNode);
            }
    }
    int count = level_max;
    for(int i = 0; i < Nodes[count-1]->nSize; i ++){
        pnode = (Abc_Obj_t *) Nodes[count-1]->pArray[i];
        candidate.push(pnode->Id);
    }
    while(!candidate.empty()){
        ID_temp = candidate.front();
        pnode = Abc_NtkObj(origin_Ntk, ID_temp);
        if(pnode->Level == count - 1) count -= 1;
        for(int h = 0; h < pnode->vFanins.nSize; h ++){
            if(critical.find(pnode->vFanins.pArray[h]) != critical.end()) continue;
            parent = Abc_NtkObj(origin_Ntk, pnode->vFanins.pArray[h]);
            if(parent->Level == count - 1) {
                candidate.push(pnode->vFanins.pArray[h]);
                critical.insert(pnode->vFanins.pArray[h]);
            }
        }
        candidate.pop();
    }
}

bool check_path(set<int>& critical_path, int* array, int size){
    for(int i = 0; i < size; i ++){
        if(critical_path.find(array[i]) != critical_path.end()){
            return false;
        }
    }
    return true;

}


/*void reduce_input(Mffc* mffc, decomposed_chart* dc, vector<vector<int>>& sim_result_ptr){
    int input_size = mffc->Pi_origin.size();
    int rows = pow(2, input_size - 1);
    int bSet[1] = {1};
    int fSet[input_size-1];
    for(int i = 0; i < input_size-1; i ++){
        fSet[i] = 2+i;
    }
    double error = 0;
    double min_error = 1;
    int reduced_ID = 0;
    auto simu_result = sim_result_ptr;
    double weight1 = 0;
    //dc->print_weight();
    for(int g = 0; g < dc->column_size; g++){
        for(int h = 0; h < dc->row_size; h++){
            weight1 += dc->weight[g][h];
        }
    }
    cout<<weight1<<endl;
    for(int j = 0; j < input_size; j ++){
        if(bSet[0] != j+1){
            int temp = bSet[0];
            bSet[0] = j+1;
            fSet[j-1] = temp;
            dc->change_set(bSet, fSet);
        }
        double weight = 0;
        for(int g = 0; g < rows; g++){
            weight += dc->weight[0][g];
            weight += dc->weight[1][g];
        }
        if(weight < 0.99){
            cout<<weight<<endl;
        }
        for(int r = 0; r < rows; r ++){
            if(dc->decompose_chart[0][r] != dc->decompose_chart[1][r]){
                if(dc->weight[0][r] <= dc->weight[1][r]) {
                    //dc->decompose_chart[1][r] = dc->decompose_chart[0][r];
                    error += dc->weight[0][r];
                }
                else {
                    //dc->decompose_chart[0][r] = dc->decompose_chart[1][r];
                    error += dc->weight[1][r];
                }
            }
        }

        if(error < 0.0005 && error < min_error){
            min_error = error;
            reduced_ID = j;
        }
        error = 0;
    }
}*/
int reduce_input(Mffc* mffc, vector<vector<int>>& sim_result_ptr, int* pattern, int* boundset, double& origin_error){
    int input_size = mffc->Pi_origin.size();
    if(input_size == 0 || mffc->hasConstant){
        return -1;
    }
    int rows = pow(2, input_size - 1);
    int bSet[1] = {1};
    int fSet[input_size-1];
    for(int i = 0; i < input_size-1; i ++){
        fSet[i] = 2+i;
    }
    double error = 0;
    double min_error = 1;
    int reduced_ID = -1;
    auto sim_result = sim_result_ptr;
    double weight1 = 0;
    decomposed_chart* dc = create_decomposed_chart(mffc->ref, bSet, fSet, mffc->Pi_origin.size(), 0, 1);
    Create_Weight(dc, mffc, bSet, fSet, sim_result);
    //dc->print_weight();
    for(int g = 0; g < dc->column_size; g++){
        for(int h = 0; h < dc->row_size; h++){
            weight1 += dc->weight[g][h];
        }
    }
    cout<<weight1<<endl;
    for(int j = 0; j < input_size; j ++){
        if(bSet[0] != j+1){
            int temp = bSet[0];
            bSet[0] = j+1;
            fSet[j-1] = temp;
            dc->change_set(bSet, fSet);
        }
        double weight = 0;
        for(int g = 0; g < rows; g++){
            weight += dc->weight[0][g];
            weight += dc->weight[1][g];
        }
        if(weight < 0.99||weight > 1.01){
            cout<<weight<<endl;
        }
        for(int r = 0; r < rows; r ++){
            if(dc->decompose_chart[0][r] != dc->decompose_chart[1][r]){
                if(dc->weight[0][r] <= dc->weight[1][r]) {
                    //dc->decompose_chart[1][r] = dc->decompose_chart[0][r];
                    error += dc->weight[0][r];
                }
                else {
                    //dc->decompose_chart[0][r] = dc->decompose_chart[1][r];
                    error += dc->weight[1][r];
                }
            }
        }

        if(error < 0.005 && error < min_error){
            min_error = error;
            reduced_ID = j+1;
        }
        error = 0;
    }
    if(reduced_ID != -1){
        cout<<"The min error is "<<min_error<<endl;
        bSet[0] = reduced_ID;
        int count = 0;
        for(int i = 0; i < input_size; i++){
            if(i+1 == reduced_ID) continue;
            fSet[count] = i+1;
            count++;
        }
        cout<<"The boundset is "<<bSet[0]<<endl;
        cout<<"The freeset is ";
        for(int i = 0; i < input_size-1; i++){
            cout<<fSet[i];
        }
        cout<<endl;
        dc->change_set(bSet, fSet);
        double weight = 0;
        for(int g = 0; g < rows; g++){
            weight += dc->weight[0][g];
            weight += dc->weight[1][g];
        }
        if(weight < 0.99||weight > 1.01){
            cout<<weight<<endl;
        }
        cout<<"rows is "<<rows<<endl;
        for(int r = 0; r < rows; r ++) {
            if (dc->decompose_chart[0][r] != dc->decompose_chart[1][r]) {
                if (dc->weight[0][r] <= dc->weight[1][r]){
                    pattern[r] = dc->decompose_chart[1][r];
                    //dc->decompose_chart[1][r] = dc->decompose_chart[0][r];
                    error += dc->weight[0][r];
                }
                else{
                    pattern[r] = dc->decompose_chart[0][r];
                    error += dc->weight[1][r];
                }
            }
            else {
                pattern[r] = dc->decompose_chart[1][r];
                //dc->decompose_chart[0][r] = dc->decompose_chart[1][r];
                //error += dc->weight[1][r];
            }
            //pattern[r] = dc->decompose_chart[1][r];
        }
        cout<<"Now the error is "<< error << endl;
        origin_error = error;
        for(int i = 0; i < input_size-1; i++){
            boundset[i] = fSet[i];
        }
    }
    cout<<"reduced_ID is "<< reduced_ID<<endl;
    return reduced_ID;
}

void amend_freeset(int** dec_chart, double** wt, int* fSet, int fSetSize, int* bSet, int bSetSize, int aug_ID, int**& aug_dc, double**& aug_wt){
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
            aug_dc[i][2*j] = dec_chart[i][j];
            aug_dc[i][2*j+1] = dec_chart[i][j];
            binary(binary_origin, fSetSize, j);
            if(binary_origin[aug_ID]==0){
                aug_wt[i][2*j] = wt[i][j];
                aug_wt[i][2*j+1] = 0;
            }
            else{
                aug_wt[i][2*j+1] = wt[i][j];
                aug_wt[i][2*j] = 0;
            }
        }
    }
    delete[] binary_bound;
    delete[] binary_free;
    delete[] binary_origin;
}

int find_level(Mffc* mffc, int id, map<int, int>& level_record){
    auto iter = level_record.find(id);
    int max_level = -1;
    if(iter == level_record.end()){
        auto pnode = Abc_NtkObj(mffc->ref, id);
        for(int i = 0; i < pnode->vFanouts.nSize; i++){
            if(max_level == -1) max_level = find_level(mffc, pnode->vFanouts.pArray[i], level_record)+1;
            else{
                if(max_level < find_level(mffc, pnode->vFanouts.pArray[i], level_record)+1){
                    max_level = find_level(mffc, pnode->vFanouts.pArray[i], level_record)+1;
                }
            }
        }
        level_record.insert(map<int, int>::value_type (id, max_level));
        return max_level;
    }
    else{
        return (*iter).second;
    }
}

int max_level(Mffc* new_mffc, int* id_set, int size){
    Abc_Obj_t* pNode;
    Abc_Obj_t* temp;
    int i, child, level_now;
    level_now = 0;
    map<int, int> level_record;
    map<int, int>::iterator iter;
    for(int j = 0; j < new_mffc->Pi_origin.size(); j ++){
        temp = new_mffc->Pi_origin[j];
        int level_temp = temp->Level;
        level_record.insert(map<int, int>::value_type (j+1, level_temp));
        //cout<<level_temp<<endl;
    }
    int max = -1;
    for(int i = 0; i < size; i++){
        if(find_level(new_mffc, id_set[i], level_record) > max){
            max = find_level(new_mffc, id_set[i], level_record);
        }
    }
    return max;
}

double check_error(vector<vector<int>>& array, Abc_Ntk_t* origin_Ntk, vector<int> origin_input_set, int new_input_id, int output_id){
    double round_num = 1000000;
    double error_count = 0;
    int size = origin_input_set.size() + 1;
    int* binary = new int[size];
    vector<vector<int>> decomposed_chart;
    vector<int> pattern;
    for(int i = 0; i < pow(2, size); i++){
        vector<int> temp;
        temp.push_back(0);
        temp.push_back(0);
        decomposed_chart.push_back(temp);
        pattern.push_back(0);
    }
    for(int i = 0; i < round_num; i++){
        for(int j = 0; j < size - 1; j++){
            binary[j] = array[origin_input_set[j]][i];
        }
        binary[size - 1] = array[new_input_id][i];
        int result = binary_inverse(binary, size);
        if(array[output_id][i] == 0){
            decomposed_chart[result][0] += 1;
        }
        else{
            decomposed_chart[result][1] += 1;
        }
    }
    for(int i = 0; i < pow(2, size); i++){
        if(decomposed_chart[i][0] < decomposed_chart[i][1]){
            error_count += decomposed_chart[i][0];
            pattern[i] = 1;
        }
        else{
            error_count += decomposed_chart[i][1];
            pattern[i] = 0;
        }
        if(error_count > 1000){
            return 1;
        }
    }
    return (double)(error_count/round_num);
}

double compare_signal(vector<vector<int>>& array, int first_node, int second_node){
    double round_num = 100000;
    double error_count = 0;
    for(int i = 0; i < round_num; i++){
        if(array[first_node][i] != array[second_node][i]){
            error_count ++;
        }
        if(error_count > 1000){
            return 1;
        }
    }
    return (double) error_count/round_num;
}

void not_zero_sasimi(Abc_Ntk_t* origin_Ntk, vector<vector<int>>& sim_result, double& upperbound_error){
    int level_max=Abc_NtkLevel(origin_Ntk);
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    Abc_Obj_t* pNode;
    int j;
    Abc_NtkForEachNode(origin_Ntk, pNode, j) {
        int level = pNode->Level - 1;
        if (level >= 0) {
            Vec_PtrPush(Nodes[level], pNode);
        }
    }
    for(int i = level_max - 1; i >= 0; i--){
        for(int j = 0; j < Nodes[i]->nSize; j++){
            pNode = (Abc_Obj_t *) Nodes[i]->pArray[j];
            double best_error = 1;
            int best_choice_id = -1;
            for(int k = 0; k <= i; k++){
                for(int l = 0; l < Nodes[k]->nSize; l++){
                    int id = ((Abc_Obj_t*) Nodes[k]->pArray[l])->Id;
                    double error = compare_signal(sim_result, pNode->Id, id);
                    if(error < best_error && id != pNode->Id){
                        best_error = error;
                        best_choice_id = id;
                    }
                }
            }
            if(best_error < 0.05){
                Abc_Ntk_t *pNtkNew = Abc_NtkCreateMffc(origin_Ntk, pNode, Abc_ObjName(pNode));
                replaced_pair.push_back(pair<int, int>(pNode->Id, best_choice_id));
                efficience.push_back(pair<int, double>(Abc_NtkNodeNum(pNtkNew), best_error));
                cout << "the original signal is " << pNode->Id << ", and now the replaced signal is " << best_choice_id << endl;
                cout << "the improvemnet is " << Abc_NtkNodeNum(pNtkNew) << ", and the error is " << best_error << endl;
            }
            else continue;
        }
    }
    for(int i = 0; i < replaced_pair.size(); i++){
        if(removed.find(replaced_pair[i].first) == removed.end() && removed.find(replaced_pair[i].second) == removed.end()){
            if(upperbound_error - efficience[i].second < 0) return;
            cout << replace_count << endl;
            replace_two_signal(origin_Ntk,replaced_pair[i].first, replaced_pair[i].second);
            replace_count += 1;
            reduced_area += efficience[i].first;
            upperbound_error -= efficience[i].second;
            removed.insert(replaced_pair[i].first);
        }
    }
    cout << "here is replaced count: " << replace_count << endl;
    cout << "and here is reduced area: " << reduced_area << endl;
}

void replace_two_signal(Abc_Ntk_t* origin_Ntk, int signal_id_1, int signal_id_2){
    Abc_Obj_t* first_node = Abc_NtkObj(origin_Ntk, signal_id_1);
    Abc_Obj_t* second_node = Abc_NtkObj(origin_Ntk, signal_id_2);
    if(first_node != NULL && second_node != NULL) {
        Abc_ObjReplace(first_node, second_node);
    }
}

pair<double, vector<int>> compare_input(Abc_Ntk_t* origin_Ntk, vector<vector<int>>& array, int origin_id, int replace_id, int output_id){
    Abc_Obj_t* output_node = Abc_NtkObj(origin_Ntk, output_id);
    int input_size = output_node->vFanins.nSize;
    vector<int> pattern;
    vector<vector<int>> record;
    vector<int> input_record;
    int new_input_size = input_size;
    for(int i = 0; i < input_size; i++){
        auto N = output_node->vFanins.pArray[i];
        if(N == replace_id){
            new_input_size -= 1;
            cout << new_input_size << endl;
            vector<int> temp;
            return pair<double, vector<int>>(1, temp);
        }
        else if(N != origin_id) input_record.push_back(N);
    }
    if(new_input_size == input_size) input_record.push_back(replace_id);
    for(int i = 0; i < input_size; i++){
        //cout << input_record[i] << endl;
    }
    //cout << "-----------------------" <<endl;
    //cout << replace_id << endl;
    //cout << origin_id << endl;
    //cout << Abc_ObjName((Abc_Obj_t*)Abc_NtkObj(origin_Ntk, replace_id)) << endl;
    //cout << Abc_ObjName((Abc_Obj_t*)Abc_NtkObj(origin_Ntk, origin_id)) << endl;
    //cout << Abc_ObjName(output_node) <<endl;
    int* bin = new int[new_input_size];
    for(int i = 0; i < pow(2, new_input_size); i++){
        pattern.push_back(0);
        vector<int> temp;
        for(int j = 0; j < 2; j++){
            temp.push_back(0);
        }
        record.push_back(temp);
    }
    for(int i = 0; i < 100000; i++){
        //cout << i << "------------------" << endl;
        for(int j = 0; j < new_input_size; j++){
            bin[j] = array[input_record[j]][i];
            //cout<<input_record[j]<<" "<<Abc_ObjName(Abc_NtkObj(origin_Ntk, input_record[j]))<<" "<<array[input_record[j]][i] << endl;
        }
        int binary_num = binary_inverse(bin, new_input_size);
        //cout << binary_num << endl;
        record[binary_num][array[output_id][i]] += 1;
        //cout << array[output_id][i] << endl;
    }
    int error_count = 0;
    for(int i = 0; i < pow(2, new_input_size); i++){
        if(record[i][0] < record[i][1]){
            pattern[i] = 1;
            error_count += record[i][0];
        }
        else{
            pattern[i] = 0;
            error_count += record[i][1];
        }
        //cout << record[i][1] << " " << record[i][0] << endl;
        if(record[i][1] != 0){
            //cout <<" hahahah "<< endl;
        }
    }
    delete[] bin;
    double error = (double) error_count / (double) 100000;
    cout << "error count is "<< error_count << endl;
    cout << "error is " << error << endl;
    return(pair<double, vector<int>>(error, pattern));
}

double not_zero_input_replace(Abc_Ntk_t* origin_Ntk, vector<vector<int>>& sim_result, double& upper_bound){
    int level_max=Abc_NtkLevel(origin_Ntk);
    Vec_Ptr_t** Nodes=new Vec_Ptr_t*[level_max];
    for(int k=0;k<level_max;k++){
        Nodes[k]=Vec_PtrAlloc(16);
    }
    int write_back_count = 0;
    Abc_Obj_t* pNode;
    double accumulate_error = 0;
    int j;
    Abc_NtkForEachNode(origin_Ntk, pNode, j) {
        int level = pNode->Level - 1;
        if (level >= 0) {
            Vec_PtrPush(Nodes[level], pNode);
        }
    }
    map<int,int> replaced;
    int reduce_LUT = 0;
    for(int i = level_max - 1; i >= 0; i--){
        for(j = 0; j < Nodes[i]->nSize; j++) {
            pNode = (Abc_Obj_t *) Nodes[i]->pArray[j];
            if(pNode->vFanouts.nSize == 0) continue;
            int fanout_replace_count = 0;
            int Fanout_size = pNode->vFanouts.nSize;
            vector<int> Fanout_set;
            for(int k = 0; k < Fanout_size; k++){
                auto fanout = pNode->vFanouts.pArray[k];
                Fanout_set.push_back(fanout);
            }
            Abc_Ntk_t *pNtkNew = Abc_NtkCreateMffc(origin_Ntk, pNode, Abc_ObjName(pNode));
            int possible_reduce = Abc_NtkNodeNum(pNtkNew);
            double possible_error = 0;
            vector<pair<int, vector<int>>> Pattern_set;
            vector<pair<int, int>> Pair_set;
            for(int k = 0; k < Fanout_size; k++){
                auto fanout = Fanout_set[k];
                auto fanout_node = Abc_NtkObj(origin_Ntk, fanout);
                if(fanout_node->vFanouts.nSize == 0) continue;
                double min_error = 1;
                pair<int, int> min_pair;
                vector<int> min_pattern;
                int min_output;
                for(int l = 0; l <= i; l++){
                    for(int ll = 0; ll < Nodes[l]->nSize; ll++) {
                        auto replace_node = (Abc_Obj_t *) Nodes[l]->pArray[ll];
                        if(replaced.find(replace_node->Id) != replaced.end() || replace_node->Id == 0) continue;
                        if(replace_node->Id != pNode->Id){
                            pair<double, vector<int>> P = compare_input(origin_Ntk, sim_result, pNode->Id, replace_node->Id, fanout);
                            if(P.first < min_error && P.first <= 0.02 && P.first < upper_bound){
                                min_error = P.first;
                                min_pattern = P.second;
                                min_pair.first= pNode->Id;
                                min_pair.second = replace_node->Id;
                                min_output = fanout;
                            }
                            //cout << "origin id is " << pNode->Id <<", and replace id is " << replace_node->Id <<", and fanout id is " << fanout << ", and the error is " << P.first << endl;
                        }
                    }
                }
                if(min_error <= 0.02){
                    possible_error += min_error;
                    //upper_bound -= min_error;
                    //write_back_count += 1;
                    //write_back_input(origin_Ntk, min_pair.first, min_pair.second, min_output, min_pattern);
                    Pattern_set.push_back(pair<int, vector<int>>(min_output, min_pattern));
                    Pair_set.push_back(pair<int,int>(min_pair.first, min_pair.second));
                    fanout_replace_count += 1;
                    cout << "the name is "<< Abc_ObjName(pNode) << endl;
                }
                cout << min_pair.first << " " << min_pair.second << " " << min_output << " " << min_error<< endl;
                //Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
                //replaced.insert(pair<int, int>(min_pair.first, 1));
            }
            if(fanout_replace_count == Fanout_size) {
                if(upper_bound - possible_error > 0) {
                    upper_bound -= possible_error;
                    for (int kkk = 0; kkk < Pattern_set.size(); kkk++) {
                        write_back_count += 1;
                        write_back_input(origin_Ntk, Pair_set[kkk].first, Pair_set[kkk].second, Pattern_set[kkk].first, Pattern_set[kkk].second);
                        replaced.insert(pair<int, int>((Pair_set[kkk].first), 1));
                    }
                    Abc_Obj_t *ppnode;
                    int lll;
                    /*Abc_NtkForEachNode(pNtkNew, ppnode, lll) {
                            replaced.insert(pair<int, int>(ppnode->Id, 1));
                        }*/
                    reduce_LUT += possible_reduce;
                    cout << possible_reduce << endl;
                    Pair_set.clear();
                    Pattern_set.clear();
                    Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
                }
            }
        }
    }
    cout << write_back_count << endl;
    cout << reduce_LUT << endl;
    Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    for(int k=0;k<level_max;k++){
        Vec_PtrFree(Nodes[k]);
    }
    delete[] Nodes;
    return 0;
}

