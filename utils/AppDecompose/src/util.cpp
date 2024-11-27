//
// Created by niqiu on 18-9-24.
//

#include "util.h"
/**Function*************************************************************
  Synopsis    []
  Description [Print information from the vector containing all the Mffcs of pNtk]
***********************************************************************/
void printMffcSetInfo(vector<Mffc*> pNtkSet, bool isNode){
    map<int, int> pNtkSetInfo = isNode?
            getMffcSetNodeInfo(pNtkSet): getMffcSetPiInfo(pNtkSet);
    auto jx = pNtkSetInfo.begin();
    for(jx = jx; jx != pNtkSetInfo.end(); jx++ ){
        if(isNode) cout << "The cone size ";
        else cout << "The Pi size ";
        cout << jx->first << " appears "\
             << jx->second << " Times" << endl;
    }
}
/**Function*************************************************************
  Synopsis    []
  Description [Obtain information from the vector containing all the Mffcs of pNtk]
***********************************************************************/

map<int ,int> getMffcSetNodeInfo(vector<Mffc*> pNtkSet){
    vector<Abc_Ntk_t*>::size_type ix;
    map<int, int> pNtkSetInfo;
    int ConeSize;
    for(ix = 0; ix != pNtkSet.size(); ix++){
        ConeSize = Abc_NtkNodeNum(pNtkSet[ix]->ref);
        pNtkSetInfo[ConeSize]++;
    }
    return pNtkSetInfo;
}
/**Function*************************************************************
  Synopsis    []
  Description [Obtain information from the vector containing all the Mffcs of pNtk]
***********************************************************************/

map<int ,int> getMffcSetPiInfo(vector<Mffc*> pNtkSet){
    vector<Abc_Ntk_t*>::size_type ix;
    map<int, int> pMffcSetInfo;
    int PiSize;
    for(ix = 0; ix != pNtkSet.size(); ix++){
        PiSize = Abc_NtkPiNum(pNtkSet[ix]->ref);
        pMffcSetInfo[PiSize]++;
    }
    return pMffcSetInfo;
}


/**Function*************************************************************
  Synopsis    []
  Description []
***********************************************************************/
void NodeMffcConeSupp(Abc_Obj_t* pNode, Vec_Ptr_t* vCone, Vec_Ptr_t* vSupp, bool print){
    Abc_NodeDeref_rec(pNode);
    Abc_NodeMffcConeSupp(pNode, vCone, vSupp);
    Abc_NodeRef_rec(pNode);

    if(print) {
        printf("Printing the MFFC for node %s\n", Abc_ObjName(pNode));
        Abc_Obj_t *pObj;
        int i;
        printf("Node = %6s : Supp = %3d  Cone = %3d  (",
               Abc_ObjName(pNode), Vec_PtrSize(vSupp), Vec_PtrSize(vCone));
        Vec_PtrForEachEntry(Abc_Obj_t *, vCone, pObj, i)printf(" %s", Abc_ObjName(pObj));
        printf(" )\n");
    }
}
/**Function*************************************************************
  Synopsis    []
  Description []
***********************************************************************/
Abc_Ntk_t* ReadBlif(char* filename){
    Abc_Ntk_t* blif_ntk = Io_ReadBlif(filename,1);
    blif_ntk = Abc_NtkToLogic(blif_ntk);
    return blif_ntk;
}
/**Function*************************************************************
  Synopsis    []
  Description [Retrieve all the Mffcs in pNtk]
***********************************************************************/
vector<Mffc*> getMffcNtk(Abc_Ntk_t* pNtk, int minPiSize, int maxPiSize){
    Vec_Ptr_t *vCone, *vSupp;
    vector<Mffc*> pNtkSet;

    int i = 0;
    Abc_Obj_t* pNode;
    Abc_NtkForEachNode(pNtk, pNode, i){
        auto pMffc = new Mffc;
        pMffc->root = pNode;
        //cout << "GetMffcNtk fron node " << Abc_ObjName(pNode) << endl;
        vCone = Vec_PtrAlloc(100);
        vSupp = Vec_PtrAlloc(100);
        NodeMffcConeSupp(pNode, vCone, vSupp, false);
        Abc_Obj_t* qNode;
        int j;
        //cout << "Mffc size is: " << Abc_NodeMffcSize(pNode) << endl;

        if(vSupp->nSize >= minPiSize && vSupp->nSize <= maxPiSize) {
            auto temp = (Abc_Obj_t *) vCone->pArray[vCone->nSize - 1];
            Abc_Ntk_t *pNtkNew = Abc_NtkCreateMffc(pNtk, temp, Abc_ObjName(pNode));
            //std::cout << Vec_PtrSize(vCone) << " " << Abc_NtkNodeNum(pNtkNew) << endl;
            pMffc->ref = pNtkNew;
            Vec_PtrForEachEntry(Abc_Obj_t*, vSupp, qNode, j){
                pMffc->Pi_origin.push_back(Abc_NtkObj(pNtk,qNode->Id));
            }
            for(int i = 0; i < vSupp->nSize; i++){
                pMffc->ID_origin.push_back(((Abc_Obj_t*)vSupp->pArray[i])->Id);
            }
            pMffc->ID_origin.push_back(0);
            for(int i = 0; i < vCone->nSize; i++){
                pMffc->ID_origin.push_back(((Abc_Obj_t*)vCone->pArray[i])->Id);
            }
            pNtkSet.push_back(pMffc);

        }
        else{
            delete pMffc;
        }

        Vec_PtrFree(vCone);
        Vec_PtrFree(vSupp);
    }
    return pNtkSet;
}



/*void map_LUT(char* old_filename, char* new_filename){
    Abc_Frame_t * pAbc;
    Abc_Start();
    pAbc = Abc_FrameGetGlobalFrame();
    char command[1000];
    sprintf(command,"read_blif %s",old_filename);
    Cmd_CommandExecute(pAbc, command);
    char* command2="if -K 4";
    Cmd_CommandExecute(pAbc, command2);
    sprintf(command,"write_blif %s",new_filename);
    Cmd_CommandExecute(pAbc, command);
    Abc_Stop();
}*/
Mffc::~Mffc(){
    /*if(this->ref != NULL && this->ref->vObjs->nSize > 0 && this->ref->vObjs->nSize < 1000) {
        Abc_NtkDelete(this->ref);
    }*/
}

Mffc* Mffc::duplicate() {
    Mffc* ans = new Mffc;
    ans->error = this->error;
    ans->Pi_origin = this->Pi_origin;
    //unnecessary replicate of memory
    ans->level_reduce = this->level_reduce;
    ans->reduced_num = this->reduced_num;
    ans->root = this->root;
    Abc_Ntk_t* pNtkNew = Abc_NtkDup(this->ref);
    ans->ref = pNtkNew;
    ans->hasConstant = hasConstant;
    for(int i = 0; i < this->ID_origin.size(); i ++){
        ans->ID_origin.push_back(this->ID_origin[i]);
    }
    return ans;
}

void Mffc::RI_write_back(int reduce_ID, int *pattern) {
    vector<Abc_Obj_t*> temp;
    vector<int> temp_int;
    Abc_Obj_t* temp_Obj;
    int temp_ID;
    //vector<Abc_Obj_t*>::iterator it;
    //temp.pop_back();
    //temp_int.pop_back();
    int size = Pi_origin.size();
    for(int i = 0; i < size; i++){
        temp_Obj = *(Pi_origin.end()-1);
        if(size - i != reduce_ID) temp.push_back(temp_Obj);
        Pi_origin.pop_back();
    }
    size = temp.size();
    for(int i = 0; i < size; i++){
        temp_Obj = *(temp.end()-1);
        Pi_origin.push_back(temp_Obj);
        temp.pop_back();
    }
    size = ID_origin.size();
    for(int i = 0; i < size; i++){
        temp_ID = *(ID_origin.end()-1);
        if(size - i != reduce_ID) temp_int.push_back(temp_ID);
        ID_origin.pop_back();
    }
    size = temp_int.size();
    for(int i = 0; i < size; i++){
        temp_ID = *(temp_int.end()-1);
        ID_origin.push_back(temp_ID);
        temp_int.pop_back();
    }
    Abc_Obj_t* pNode, *qNode;
    Abc_Ntk_t* ntk = this->ref;
    int j;
    int max = 0;
    Abc_NtkForEachNode(ntk, pNode, j) {
            if (Abc_ObjId(pNode) > max) max = Abc_ObjId(pNode);
        }
    assert(max != 0);
    Abc_NtkForEachNode(ntk, pNode, j) {
            if (Abc_ObjId(pNode) == max) continue;
            Abc_NtkDeleteObj(pNode);
        }
    //delete every node


    Abc_NtkForEachPi(ntk, pNode, j){
        if(Abc_ObjId(pNode) == reduce_ID) Abc_NtkDeleteObj(pNode);
    }
    //delete input


    int cnt = 0;
    int sz = (int) pow(2, this->Pi_origin.size());
    if(sz == 1){
        cout<<"!!!!"<<endl;
        cout<<"Now the pattern is "<< pattern[0] << endl;
    }
    //assume we have 11 inputs
    for(int i = 0; i < sz; i++){
        if(pattern[i] == 1) cnt++;
    }

    Abc_Obj_t* endNode, *oNode;

    if(cnt == 0) {
        endNode = Abc_NtkCreateNodeConst0(ntk);
        this->hasConstant = true;
    }
    else endNode = Abc_NtkCreateNode(ntk);

    string name = std::to_string(Abc_ObjId(endNode));
    Nm_ManStoreIdName(ntk->pManName, Abc_ObjId(endNode),
                      ABC_OBJ_NODE, (char *) "xx", (char *) name.c_str());

    if(cnt != 0) {
        Abc_NtkForEachPi(ntk, qNode, j) {
            Abc_ObjAddFanin(endNode, qNode);
        }
    }
    //Add all Pis to endnode


    if(cnt != 0) {
        char *buffer = new char[20 * cnt]{0};
        for (int i = 0; i < sz; i++) {
            if (pattern[i] == 1) {
                string s = dec2bin((unsigned) i, this->Pi_origin.size());
                sprintf(buffer, "%s%s 1\n", buffer, s.c_str());
            }
        }
        endNode->pData = buffer;
        //assigning pData (truth table relationship)
    }



    oNode = Abc_NtkPo(ntk, 0);
    Abc_ObjRemoveFanins(oNode);
    Abc_ObjAddFanin(oNode, endNode);
    //add to the primary output (strange rule from abc)

    Abc_Ntk_t* pNtkNew = Abc_NtkDup(ntk);
    Abc_NtkDelete(ntk);
    this->ref = pNtkNew;

    Abc_NtkForEachObj(this->ref, pNode, j){
            cout << Abc_ObjId(pNode) << " ";
        }
    cout << endl;
}

/*Abc_Ntk_t* abc_map(Abc_Ntk_t* const in){
    char* fn = (char*) "temp_abc_map.blif";
    Io_WriteBlifLogic(in, fn, 1);
    //Abc_Start();
    char Command[1000];
    Abc_Frame_t * pAbc = Abc_FrameGetGlobalFrame();

    sprintf( Command, "read_blif ./%s", fn);
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
    //Abc_Stop();

    return ReadBlif(fn);
}*/
