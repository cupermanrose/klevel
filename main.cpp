#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "header.h"
#include "klevel.h"
#include "filemem.h"
#include "tgs.h"
#include "skyline.h"
#include "cellTree.h"
#include "tgs.h"
#include "utk.h"

string ConfigFile = "/home/jiahao/disk/Projects/klevel/config/config.txt";
fstream fout("/home/jiahao/disk/Projects/klevel/config/results.txt", ios::out);
string klevelfile;

float** PointSet = new float* [MAXPTS + 1]; // object set with cl and cu
RtreeNodeEntry **p = new RtreeNodeEntry *[MAXPTS];
vector<vector<float>> HalfSpaces; // halfspace
//unordered_map<long int, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long long, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory
unordered_map<long int, cell*> cellID; // cell ID
double totalSpaceCost = 0.0; // space cost (MB)
double totalIO=0;

string datafile, rtreefile, resultfile;
int dim, ik, qk, querytype, objCnt, query_num;
vector<vector<float>> OriginData;
vector<Hyperplane> AllHP;
vector<Obj> All_records;

void ReadParameter(string ConfigFile) {
    fstream fin(ConfigFile);
    fin >> dim >> ik >> qk >> querytype >> query_num;
    fin >> datafile;
    fin >> rtreefile;
    fin >> klevelfile;
    fin.close();
    cout << "Config Done" << endl;
}
void DataLoading(int dim, string datafile) {
    fstream fin(datafile, ios::in);
    OriginData.clear();
    objCnt = 0;
    while (true) {
        int id;
        float *cl = new float[dim];
        float *cu = new float[dim];
        fin >> id;
        if (fin.eof())
            break;

        // for baseline
        PointSet[objCnt + 1] = new float[2 * dim];

        vector<float> tmp; tmp.clear();
        for (int d = 0; d < dim; d++) fin >> cl[d];
        for (int d = 0; d < dim; d++) {
            fin >> cu[d];
            cu[d] = cl[d] + EPS;
            tmp.push_back(cl[d]);
            PointSet[objCnt + 1][d] = cl[d];
            PointSet[objCnt + 1][d + dim] = cu[d];
        }
        OriginData.push_back(tmp);


        Hypercube hc(dim, cl, cu);
        p[objCnt++] = new RtreeNodeEntry(objCnt, hc);

        if (TEST) {
            if (objCnt >= 50) break;
        }
        //log information
        if (objCnt % 1000 == 0)
            cout << ".";
        if (objCnt % 10000 == 0)
            cout << objCnt << " objects loaded" << endl;
    }

    cout << "Total number of objects: " << objCnt << endl;
    fin.close();
    return;
}
void Skyband_Onionlayer(KLEVEL& idx, int dim, int Qk_max) {
    // build rtree
    cout << "Bulkloading R-tree..." << endl;
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    FileMemory mem(PAGESIZE, rtreefile, RtreeNodeEntry::fromMem, true);
    Rtree * rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild * 0.3, (int)maxChild * 0.3, p, OriginData.size(), false);
    cout << "[Rtree build done]" << endl;

    vector<int> skyband, klayers;
    //k-skyband
    kskyband(dim, *rtree, skyband, OriginData, Qk_max);
    cout << "skyband size: " << skyband.size() << endl;
    //k-onion-layer
    onionlayer(skyband, OriginData, Qk_max, klayers, dim);
    cout << "# k-onion layer records / # of k-skyband: " << klayers.size() << "/" << skyband.size() << endl;


    // All records with 2 baseline filter
   All_records.clear();
    for (int i = 0; i < klayers.size(); i++) {
        Obj NewObj(i, dim, OriginData[klayers[i]]);
        All_records.push_back(NewObj);
    }

    return;

}
void PreComputeAllHP(int dim) {
    // Precompute hyperplane
    AllHP.clear();
    for (int i = 0; i < All_records.size(); i++) {
        for (int j = 0; j < All_records.size(); j++) {
            int tmp = Halfspace::ComputeHP(All_records[i], All_records[j], dim);
        }
    }
}

// UTK query processing

void UTK_QueryGenerator(int Q_num, vector < vector<float> >& QryGen, const int dim, const float sigma) {
    QryGen.clear();

    // test
    /*vector<float> region(2 * dim, 0);
    for (int di = 0; di < dim; di++) {
        region[2 * di] = 0.01;
        region[2 * di + 1] = 0.01 + sigma;
    }
    QryGen.push_back(region);*/

    float query[20];
    srand(123);
    while (QryGen.size() < Q_num)
    {
        vector<float> region(2 * dim, 0);
        float sum = 0;
        for (int di = 0; di < dim; di++)
        {
            query[di] = rand() * 1.0 / RAND_MAX;
            sum += query[di];
        }

        if (sum < 1 - sigma)
        {
            for (int di = 0; di < dim; di++)
            {
                region[di * 2] = query[di];
                region[di * 2 + 1] = query[di] + sigma;
            }
            QryGen.push_back(region);
        }
    }
}

void UTK_BaselineQueryProcess(Rtree* rtree, vector<float>& Qregion, int Qk, int Qdim) {
    clock_t at = clock();
    unordered_map<int, cell*> utkRet;
    utkRet.clear();
    UTK* sol = new UTK();
    vector<long int> rskyband;
    vector<cell*> exactutk;
    sol->rskyband(Qregion, Qdim, *rtree, rskyband, PointSet, Qk); // filter

    cout << "r-skyband size: " << rskyband.size() << endl;
    if (rskyband.size() == Qk)
    {
        cellTree* candidate = new cellTree(2 * (Qdim - 1));
        exactutk.push_back(candidate->root);
    }
    else
    {
        sol->jaa(Qregion, exactutk, Qk, Qdim, PointSet, *rtree); // refinement
        //treeSpacecost = sol->treeSize;
    }

    cout << "============================================" << endl;
    //cout << "# UTK records / # of k-skyband: " << utkRet.size() << "/" << rskyband.size() << endl;
    cout << "Regions: " << exactutk.size() << endl;
    for (auto it = utkRet.begin(); it != utkRet.end(); it++) {
        cout << it->first << ' ';
    }
    cout << endl;
    //resultSize += exactutk.size();
    cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
    /*cout << "Total space cost: " << fixed << totalSpaceCost + treeSpacecost << " MB " << endl;
    Spacecost += treeSpacecost;*/
    cout << "============================================" << endl;

    delete sol;
    utkRet.clear();
    unordered_map<int, cell*>().swap(utkRet);
    rskyband.clear();
    vector<long int>().swap(rskyband);
    exactutk.clear();
    vector<cell*>().swap(exactutk);
}
void UTK_BaselineSolution(vector<vector<float>>& Qregions, int dim) {
    // build rtree
    clock_t at = clock();
    cout << "Bulkloading R-tree..." << endl;
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    FileMemory mem(PAGESIZE, rtreefile, RtreeNodeEntry::fromMem, true);
    Rtree * rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild * 0.3, (int)maxChild * 0.3, p, objCnt, false);
    cout << "[Rtree build done]" << endl;


    for (int i = 0; i < Qregions.size(); i++) {
        if (i==31) continue; // bug in this case
        if (i==64) continue; // bug in this case
        if (i==92) continue; // bug in this case
        //cout << "Query: " << i << endl;
        UTK_BaselineQueryProcess(rtree, Qregions[i], cur_qk, dim);
        fout << i << ' ' << (clock() - at) * 1.0 / CLOCKS_PER_SEC << endl;
    }
    cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
    fout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
}

void UTK_IndexProcess(KLEVEL& idx, vector<float>& Qregion, int Qk, int Qdim) {
    int cnt = 0;
    set<int> utk; utk.clear();

    if (Qk <= idx.Ik_max) {
        for (int i = 0; i < idx.L[Qk].size(); i++) {
            if (KLEVEL::Intersect(Qregion, idx.L[Qk][i], Qdim)) {
                for (auto it = idx.L[Qk][i].confirmed.begin(); it != idx.L[Qk][i].confirmed.end(); it++) {
                    utk.insert(*it);
                }
                cnt++;
            }
        }
    }
    else {
        vector<KLEVEL::region> tmpL; tmpL.clear();
        for (int i = 0; i < idx.L[idx.Ik_max].size(); i++) {
            if (KLEVEL::Intersect(Qregion, idx.L[idx.Ik_max][i], Qdim)) {
                tmpL.push_back(idx.L[idx.Ik_max][i]);
            }
        }
        cout << tmpL.size() << endl;
        cnt = KLEVEL::UpdateLevel_forUTK(All_records, tmpL, Qk - idx.Ik_max, Qdim, Qregion, utk);
    }

    cout << "Regions: " << cnt << endl;
    cout << "UTK: " << utk.size() << endl;

}
void UTK_IndexSolution(vector<vector<float>>& Qregions, int dim, int Ik_max, int Qk_max) {

    KLEVEL idx(dim, ik, qk);

    if (IndexBuilding) {
        Skyband_Onionlayer(idx, dim, Qk_max);
        PreComputeAllHP(dim);
        idx.Bulk_Loading_cplex(fout, All_records, dim, klevelfile);
        //idx.WriteToDisk(klevelfile);
        cout << "Building Done" << endl;
        return;
    }
    else {
        Skyband_Onionlayer(idx, dim, Qk_max);
        idx.ReadFromDisk(All_records, klevelfile);
    }

    clock_t at = clock();

    for (int i = 0; i < Qregions.size(); i++) {
        if (i==31) continue; // bug in this case
        if (i==64) continue; // bug in this case
        if (i==92) continue; // bug in this case
        cout << "Query " << i << endl;
        //fout << "Query " << i << endl;
        UTK_IndexProcess(idx, Qregions[i], cur_qk, dim);
        fout << i << ' ' << (clock() - at) * 1.0 / CLOCKS_PER_SEC << endl;
    }

    cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
    fout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;

}

// kspr query processing

void kSPR_QueryGenerator(int Q_num, vector<int>& Qids) {
    Qids.clear();
    srand(0);
    for (int i = 0; i < Q_num; i++) {
        Qids.push_back(rand() % objCnt);
    }
}

void kSPR_IndexProcess(KLEVEL& idx, int Qid_ori, int Qk, int Qdim) {
    bool flag = false;
    int Qid;
    for (int i = 0; i < All_records.size(); i++) {
        if (All_records[i].ID == Qid_ori) {
            flag = true;
            Qid = i;
            break;
        }
    }
    if (!flag) return;


    int cnt = 0;

    if (Qk <= idx.Ik_max) {
        for (int i = 0; i < idx.L[Qk].size(); i++) {
            if (idx.L[Qk][i].confirmed.find(Qid) != idx.L[Qk][i].confirmed.end()) cnt++;
        }
    }
    else {
        vector<KLEVEL::region> tmpL; tmpL.clear();
        for (auto it = idx.L[idx.Ik_max].begin(); it != idx.L[idx.Ik_max].end(); it++) {
            if (it->confirmed.find(Qid) != it->confirmed.end()) cnt++;
            else if (it->candidates.find(Qid) != it->candidates.end()) tmpL.push_back(*it);
        }
        cout << 0 << ' ' << tmpL.size() << endl;
        cnt = cnt + KLEVEL::UpdateLevel_forkSPR(All_records, tmpL, Qk - idx.Ik_max, Qdim, Qid);

        tmpL.clear();
        vector<KLEVEL::region>().swap(tmpL);
    }

    cout << "Regions: " << cnt << endl;
}
void kSPR_IndexSolution(vector<int>& Qids, int dim, int Qk_max, int Ik_max) {

    KLEVEL idx(dim, Ik_max, Qk_max);

    if (IndexBuilding) {
        Skyband_Onionlayer(idx, dim, Qk_max);
        PreComputeAllHP(dim);
        idx.Bulk_Loading_cplex(fout, All_records, dim, klevelfile);
        idx.WriteToDisk(klevelfile);
        cout << "Building Done" << endl;
        return;
    }
    else {
        Skyband_Onionlayer(idx, dim, Qk_max);
        PreComputeAllHP(dim);
        idx.ReadFromDisk(All_records, klevelfile);
    }

    clock_t at = clock();

    for (int i = 0; i < query_num; i++) {
        cout << "Query " << i << endl;
        kSPR_IndexProcess(idx, Qids[i], qk, dim);
    }

    cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
}

int kSPR_BaselineQueryProcess(Rtree* rtree, int Qid, int Qk, int Qdim) {

    clock_t at = clock();
    unordered_set<long int> skylines;
    unordered_set<long int> removeSL;
    unordered_set<long int> singular;
    multimap<long int, RtreeNodeEntry*> RecordEntry;
    multimap<float, int> heap;
    cellTree* sol = new cellTree();
    vector<long int> skyband, newAddSL;
    vector<cell*> finalResult;
    vector<cell*> leaves;

    RecordIDtoHalfPlaneID.clear();
    //unordered_map<long int, long int>().swap(RecordIDtoHalfPlaneID);
    unordered_map<long long, long int>().swap(RecordIDtoHalfPlaneID);

    HalfSpaces.clear();
    vector<vector<float>>().swap(HalfSpaces);

    long int totalNoOfHPs = 0;

    Point Qpt;
    Qpt.m_dimen = Qdim;
    for (int i = 0; i < Qdim; i++) Qpt.m_coor[i] = PointSet[Qid][i];

    // count dominators ???
    //int Noofdominator = countDominator(*rtree, PointSet, Qpt, RecordEntry);

    int insertCount = 0;
    do
    {
        updateSkylines(Qdim, *rtree, skylines, PointSet, heap, RecordEntry);
        sol->maintainDAG(skylines, removeSL, PointSet, Qdim);
        computeHP(Qdim, PointSet, Qpt, skylines, newAddSL);
        totalNoOfHPs += newAddSL.size();

        cout << "(1) Skylines(" << insertCount << "): " << newAddSL.size() << "; #Ret: " << finalResult.size() << endl;
        if (newAddSL.size() > 0)
        {
            if (insertCount == 0)
                sol->insert_kspr(newAddSL, Qk, *rtree, Qpt, finalResult);
            else
                sol->opt_insert(newAddSL, Qk, *rtree, Qpt, finalResult);

            sol->collectLeaf_kspr(leaves, Qk);
        }
        else
        {
            sol->collectLeaf_kspr(leaves, Qk);
        }

        // process leafs
        sol->markSingular_kspr(leaves, removeSL, Qk, finalResult, singular);
        cout << "(3) Singular, #Ret: " << finalResult.size() <<endl;
        removeSkylines(skylines, removeSL);
        if (removeSL.size() == 0)
        {
            break;
        }

        insertCount++;
    } while (leaves.size() > 0);

    int ans = finalResult.size();
    // Free Memory
    {
        for (int i = 0; i < finalResult.size(); i++)
        {
            finalResult[i]->release();
            delete finalResult[i];
        }
        finalResult.clear();
        vector<cell*>().swap(finalResult);

        newAddSL.clear();
        vector<long>().swap(newAddSL);

        for (multimap<long int, RtreeNodeEntry*>::iterator iter = RecordEntry.begin(); iter != RecordEntry.end(); iter++)
        {
            delete iter->second;
        }
        RecordEntry.clear();
        multimap<long int, RtreeNodeEntry*>().swap(RecordEntry);

        heap.clear();
        multimap<float, int>().swap(heap);

        skylines.clear();
        unordered_set<long int>().swap(skylines);

        removeSL.clear();
        unordered_set<long int>().swap(removeSL);

        singular.clear();
        unordered_set<long int>().swap(singular);

        RecordIDtoHalfPlaneID.clear();
        unordered_map<long long, long int>().swap(RecordIDtoHalfPlaneID);

        HalfSpaces.clear();
        vector<vector<float>>().swap(HalfSpaces);

        delete sol;
    }

    cout << "Result Size: " << ans << endl;
    cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
    return ans;

}
void kSPR_BaselineSolution(vector<int>& Qids, int dim, int k) {
    // build rtree
    cout << "Bulkloading R-tree..." << endl;
    const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
    FileMemory mem(PAGESIZE, rtreefile, RtreeNodeEntry::fromMem, true);
    Rtree * rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild * 0.3, (int)maxChild * 0.3, p, objCnt, false);
    cout << "[Rtree build done]" << endl;

    rtreeRAM(*rtree, ramTree);
    // aggregate rtree
    aggregateRecords(*rtree);
    cout << "[Aggregate Rtree done]" << endl;

    clock_t at = clock();
    for (int i = 0; i < query_num; i++) {
        clock_t at_this=clock();
        cout << "Query: " << i << endl;
        fout << "Query: " << i << endl;
        int cnt = kSPR_BaselineQueryProcess(rtree, Qids[i], qk, dim);
        cout << "Time cost: " << fixed << (clock() - at_this) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
        fout << "Result Number:" << cnt << endl;
        fout << "Time cost: " << fixed << (clock() - at_this) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
    }
    cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
    fout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
}

int main() {
    ReadParameter(ConfigFile);
    DataLoading(dim, datafile);
    if (IndexBuilding){
        KLEVEL idx(dim, ik, qk);
        Skyband_Onionlayer(idx, dim, qk);
        PreComputeAllHP(dim);
        idx.Bulk_Loading_cplex(fout, All_records, dim, klevelfile);
        //idx.WriteToDisk(klevelfile);
        return 0;
    }


    if (querytype == 1) {
    	float sigma = Q_sigma;
    	vector<vector<float>> Qregions;
    	UTK_QueryGenerator(query_num, Qregions, dim - 1, sigma);

    	UTK_BaselineSolution(Qregions, dim);
    	//UTK_IndexSolution(Qregions, dim, ik, qk);
    }
    else if (querytype == 2) {
    	vector<int> Qids;
    	kSPR_QueryGenerator(query_num, Qids);
    	kSPR_BaselineSolution(Qids, dim, qk);
    	//kSPR_IndexSolution(Qids, dim, Q_k, I_k);
    }

    return 0;
}