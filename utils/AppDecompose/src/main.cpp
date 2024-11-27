#include <iostream>
#include <cassert>
#include <abc_api.h>
#include <util.h>
#include <simulation.h>
#include <approx.h>
#include <vector>
#include <map>
#include <set>
#include <queue>
#include <algorithm>
#include "cmdline.h"

using namespace std;

//#define DEBUG
void nodetest();
void simu_test();
void print_node_info(Abc_Ntk_t* blif_ntk);
void transform(char* filename);

int FFC_count = 0;
int MFFC_count = 0;
set<int> area;

// Pmatrix PM(1,1,1);

vector<pair<int, int>> replaced_pair;

vector<pair<int, double>> efficience;

void printMffcInfo(Abc_Ntk_t* pNtk){
    cout << "Level = " << Abc_NtkLevel(pNtk) << endl;
    cout << "Number of internal nodes = " << Abc_NtkNodeNum(pNtk) << endl;
    cout << "Number of Pi = " << Abc_NtkPiNum(pNtk) << endl;
    assert(Abc_NtkPoNum(pNtk) == 1);
}

struct cpr_rl1{
    bool operator()(Mffc* pMffc, Mffc* qMffc) const {
        return pMffc->root->Level < qMffc->root->Level;
    }
}cpr_rl_t1;

Abc_Ntk_t* abc_map(Abc_Ntk_t* const in){
    char* fn = (char*) "temp_abc_map.blif";
    Io_WriteBlifLogic(in, fn, 1);
    Abc_Start();
    char Command[1000];
    Abc_Frame_t * pAbc = Abc_FrameGetGlobalFrame();

    sprintf( Command, "read_blif ./%s", fn);
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "strash") ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "rewrite") ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "balance") ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "rewrite") ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "balance") ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "rewrite") ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "if -K 4") ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    sprintf( Command, "write_blif ./%s", fn) ;
    if ( Cmd_CommandExecute( pAbc, Command ) )
    {
        fprintf( stdout, "Cannot execute command \"%s\".\n", Command );
        assert( 0 );
    }
    Abc_Stop();

    return ReadBlif(fn);
}

using namespace cmdline;
parser Cmdline_Parser(int argc, char * argv[])
{
  parser option;
  option.add <string> ("input", 'i', "Input file", false, "circuit/i2c_exact.blif");
  option.add <string> ("error", 'e', "Error Threshold", false, "0.01");
  option.parse_check(argc, argv);
  return option;
}

int main(int argc, char* argv[]) {
    clock_t t = clock();
    // command line parser
    parser option = Cmdline_Parser(argc, argv);
    string filename_str = option.get <string> ("input");
    char* filename = new char[filename_str.length() + 1];
    std::strcpy(filename, filename_str.c_str());

    char* temp_name = (char*) "out_whole.blif";

    string upper_bound_error_str = option.get <string> ("error");
    double upper_bound_error = stod(upper_bound_error_str);

    Abc_Obj_t* pNode;
    int i, j;
    double error = 0;
    set<int, firstcome> changed;
    vector<vector<int>> sim_result;
    Abc_Ntk_t* origin_Ntk = ReadBlif(filename);
    //Io_WriteBlifLogic(origin_Ntk, (char *) "out_whole.blif", 1);
    //test code
    cout << "Total number of nodde: " << Abc_NtkNodeNum(origin_Ntk) << endl;
    Abc_Ntk_t *temp_Ntk;
    /*temp_Ntk = abc_map(origin_Ntk);
    Io_WriteBlifLogic(temp_Ntk, (char *) "out_whole.blif", 1);
    cout << upper_bound_error << endl;
    Abc_NtkDelete(origin_Ntk);
    origin_Ntk = ReadBlif(temp_name);*/
    Abc_Obj_t* tempnode = Abc_NtkObj(origin_Ntk, 1);
    cout<<Abc_ObjName(tempnode)<<endl;
    cout<<tempnode->Id<<endl;
    set<int> critical_path;
    find_critical_path(origin_Ntk, critical_path);
    cout<<"-----------------------"<<endl;
    Abc_Obj_t* tempa;
    int u;
    int o;
    cout << "Total number of nodde: " << Abc_NtkNodeNum(origin_Ntk) << endl;
    cout << "Level is " << Abc_NtkLevel(origin_Ntk) <<endl;
    clock_t t_sim = clock();
    /*Abc_Obj_t* pnode = Abc_NtkObj(origin_Ntk, 320);
    cout << Abc_ObjName(pnode) << endl;
    pnode = Abc_NtkObj(origin_Ntk, 321);
    cout << Abc_ObjName(pnode) << endl;
    pnode = Abc_NtkObj(origin_Ntk, 578);
    cout << Abc_ObjName(pnode) << endl;
    pnode = Abc_NtkObj(origin_Ntk, 405);
    cout << Abc_ObjName(pnode) << endl;*/
    //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
    sim_result = simulate_whole(origin_Ntk);
    t_sim = clock() - t_sim;
    cout << "Simulation completes in " << (double)t_sim/CLOCKS_PER_SEC << "s" << endl;
    // int input_rep;
    // cout <<" do you want to apply input replacement? " <<endl;
    // cin >> input_rep;
    not_zero_input_replace(origin_Ntk, sim_result, upper_bound_error);
    if (upper_bound_error > 0) Io_WriteBlifLogic(origin_Ntk, (char *) "out_whole.blif", 1);
    cout << upper_bound_error << endl;
    Abc_NtkDelete(origin_Ntk);
    origin_Ntk = ReadBlif(temp_name);
    sim_result = simulate_whole(origin_Ntk);
    vector<Mffc*>::size_type ix = 0;
    //PM.test();
    //PM.init(&sim_result, origin_Ntk);
    //PM.calculateCPM();
#ifndef DEBUG
    vector<Mffc*> mffc_set = getMffcNtk(origin_Ntk, 8, 8);
    std::sort (mffc_set.begin(), mffc_set.end(), cpr_rl_t1);
    Mffc* mffc = nullptr;
    int wb = 0;
    Mffc* new_mffc;
    // int choice = 1;
    // cout<<"If you want to test constant 0 case, please input 0:"<<endl;
    // cin >> choice;
    // Mffc* temp;
    // if(choice == 0){
    //     constant_zero(origin_Ntk, 20, 4, sim_result, changed, critical_path,upper_bound_error);
    //     return 0;
    // }
    clock_t t_ = clock();
    double origin_error = upper_bound_error;
    constant_zero(origin_Ntk, 100, 4, sim_result, changed, critical_path, upper_bound_error);
    //if(origin_error != upper_bound_error) {
        /*Abc_Ntk_t *temp_Ntk;
        temp_Ntk = abc_map(origin_Ntk);
        Io_WriteBlifLogic(temp_Ntk, (char *) "out_whole.blif", 1);
        cout << upper_bound_error << endl;
        Abc_NtkDelete(origin_Ntk);
        origin_Ntk = ReadBlif(temp_name);*/
        //sim_result = simulate_whole(origin_Ntk);
    //}
    if(upper_bound_error < 0)
        return 0;
    //upper_bound_error -= not_zero(origin_Ntk, 12, 6, upper_bound_error, sim_result, changed, critical_path);
    //testcode 10/15
    //upper_bound_error -= not_zero(origin_Ntk, 1000, 6, upper_bound_error, sim_result, changed, critical_path);
    //testcode 10/15
    cout<<FFC_count<<" end"<<endl;
    cout<<area.size()<<endl;
    cout<<origin_Ntk->vPis->nSize<<endl;
    t_ = clock() - t_;
    cout << "not zero completes in " << (double)t_/CLOCKS_PER_SEC << "s" << endl;
    if (upper_bound_error > 0) Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    //if(error < upper_bound_error) {
    if(upper_bound_error > 0){
        origin_Ntk = ReadBlif(temp_name);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        sim_result = simulate_whole(origin_Ntk);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        //PM.init(&sim_result, origin_Ntk);
        //PM.calculateCPM();
        t_ = clock();
        upper_bound_error -= not_zero(origin_Ntk, 20, 12, upper_bound_error, sim_result, changed, critical_path);
        t_ = clock() - t_;
        cout << "not zero completes in " << (double) t_ / CLOCKS_PER_SEC << "s" << endl;
        if (upper_bound_error > 0) Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    }
    else cout<<"error reach the upperbound"<<endl;
    //if(error < upper_bound_error) {
    cout << upper_bound_error << endl;
    cout << "Total number of nodde: " << Abc_NtkNodeNum(origin_Ntk) << endl;
    if(upper_bound_error > 0){
        origin_Ntk = ReadBlif(temp_name);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        sim_result = simulate_whole(origin_Ntk);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        //PM.init(&sim_result, origin_Ntk);
        //PM.calculateCPM();
        t_ = clock();
        upper_bound_error -= not_zero(origin_Ntk, 30, 21, upper_bound_error, sim_result, changed, critical_path);
        t_ = clock() - t_;
        cout << "not zero completes in " << (double) t_ / CLOCKS_PER_SEC << "s" << endl;
        if (upper_bound_error > 0) Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    }
    else cout<<"error reach the upperbound"<<endl;
    cout << "Total number of nodde: " << Abc_NtkNodeNum(origin_Ntk) << endl;
    //if(error < upper_bound_error) {
    if(upper_bound_error > 0){
        origin_Ntk = ReadBlif(temp_name);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        sim_result = simulate_whole(origin_Ntk);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        //PM.init(&sim_result, origin_Ntk);
        //PM.calculateCPM();
        t_ = clock();
        upper_bound_error -= not_zero(origin_Ntk, 40, 31, upper_bound_error, sim_result, changed, critical_path);
        t_ = clock() - t_;
        cout << "not zero completes in " << (double) t_ / CLOCKS_PER_SEC << "s" << endl;
        if (upper_bound_error > 0) Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    }
    else cout<<"error reach the upperbound"<<endl;
    cout << "Total number of nodde: " << Abc_NtkNodeNum(origin_Ntk) << endl;
    //if(error < upper_bound_error) {
    if(upper_bound_error > 0){
        origin_Ntk = ReadBlif(temp_name);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        sim_result = simulate_whole(origin_Ntk);
        //PM.resize(100000, Abc_NtkNodeNum(origin_Ntk), Abc_NtkPoNum(origin_Ntk));
        //PM.init(&sim_result, origin_Ntk);
        //PM.calculateCPM();
        t_ = clock();
        upper_bound_error -= not_zero(origin_Ntk, 100, 41, upper_bound_error, sim_result, changed, critical_path);
        t_ = clock() - t_;
        cout << "not zero completes in " << (double) t_ / CLOCKS_PER_SEC << "s" << endl;
        if (upper_bound_error > 0) Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    }
    else cout<<"error reach the upperbound"<<endl;

#endif

    if (upper_bound_error > 0) Io_WriteBlifLogic(origin_Ntk, (char*)"out_whole.blif", 1);
    //assert(mffc);
    Abc_NtkDelete(origin_Ntk);
    delete mffc;
    for(ix = 0; ix != mffc_set.size(); ix++){
        if(mffc_set[ix] == mffc) continue;
        delete mffc_set[ix];
    }

    t = clock() - t;
    cout << "Time taken: " << (double)t/CLOCKS_PER_SEC << endl;
    cout << upper_bound_error << endl;
    delete[] filename;
    return 0;
}
