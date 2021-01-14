#include "header.h"
#include "object.h"
#include "klevel.h"
#include "hypercube.h"
#include "rentry.h"
#include "rnode.h"
#include "rtree.h"
#include "tgs.h"
#include "filemem.h"
#include "skyline.h"
#include "qhull_klevel.h"
#include "cellTree.h"
#include "utk.h"
#include "global.h"
#include "ilcplex/ilocplex.h"

// global variables
int objCnt = 0; // # of data objects
double totalIO = 0;
//double totalSpaceCost = 0.0; // space cost (MB)
//unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory
RtreeNodeEntry** p = new RtreeNodeEntry * [MAXPTS];
RtreeNodeEntry** Rp = new RtreeNodeEntry * [MAXREG];
float** PointSet = new float* [MAXPTS + 1]; // object set with cl and cu
float** TotalObj = new float* [MAXPTS + 1]; // (cl+cu)/2.0
vector<Hyperplane> AllHP;
vector<Obj> All_records, tmp_records;

vector<vector<float>> HalfSpaces; // halfspace 
//unordered_map<long int, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long long, long int> RecordIDtoHalfPlaneID;  //  record ID to halfplane ID
unordered_map<long int, RtreeNode*> ramTree; // load Rtree to main-memory
unordered_map<long int, cell*> cellID; // cell ID
double totalSpaceCost = 0.0; // space cost (MB)

//int kspr_Qid[10] = { 377378,91783,1351,72472,145257,287445,44387,9610,69497,103675 }; // 20
int kspr_Qid[10] = { 103675,297607,15825,258564,201384,302803,283662,264839,189965,166527 };

//fstream ResultOut("E:\\klevel\\Project\\klevel_index\\results.txt", ios::out);

string ReadParameter(const int a_argc, const char** a_argv, const char* a_param, const char* a_def)
{
	for (int i = 0; i < a_argc; i++)
	{
		if (strcmp(a_argv[i], a_param) == 0)
		{
			if (i + 1 == a_argc)
				return "";
			else
				return a_argv[i + 1];
		}
	}
}

void DataLoading(int dim, string datafile) {
	// data loading
	cout << "Load data points from file" << endl;
	fstream fpdata;
	fpdata.open(datafile, ios::in);
	while (true)
	{
		int id;
		float* cl = new float[dim];
		float* cu = new float[dim];
		fpdata >> id;
		if (fpdata.eof())
			break;

		PointSet[objCnt + 1] = new float[2 * dim];

		for (int d = 0; d < dim; d++)
		{
			fpdata >> cl[d];
			PointSet[objCnt + 1][d] = cl[d];
		}

		for (int d = 0; d < dim; d++)
		{
			fpdata >> cu[d];
			cu[d] = cl[d] + eps_klevel;
			PointSet[objCnt + 1][d + dim] = cu[d];
		}

		Hypercube hc(dim, cl, cu);
		p[objCnt++] = new RtreeNodeEntry(id, hc);

		if (TEST) {
			if (objCnt >= D_size) break;
		}
		//log information
		if (objCnt % 1000 == 0)
			cout << ".";
		if (objCnt % 10000 == 0)
			cout << objCnt << " objects loaded" << endl;
	}

	cout << "Total number of objects: " << objCnt << endl;
	//double rawSize = (objCnt * sizeof(float) * 2 * dim) / MB;
	//totalSpaceCost += rawSize;
	fpdata.close();
	return;
}
void UTK_QueryGenerator(vector < vector<float> >& QryGen, const int dim, const float sigma)
{
	QryGen.clear();

	// test
	/*vector<float> region(2 * dim, 0);
	for (int di = 0; di < dim; di++) {
		region[2 * di] = 0.01;
		region[2 * di + 1] = 0.01 + sigma;
	}
	QryGen.push_back(region);*/

	float query[20];
	srand(0);
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
void kSPR_QueryGenerator(vector<int>& Qids) {
	Qids.clear();
	srand(0);
	for (int i = 0; i < Q_num; i++) {
		Qids.push_back(rand() % objCnt);
	}
}

void Skyband_Onionlayer(KLEVEL& idx, int dim, int Qk_max) {
	// build rtree
	cout << "Bulkloading R-tree..." << endl;
	const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
	FileMemory mem(PAGESIZE, rtreefile, RtreeNodeEntry::fromMem, true);
	Rtree * rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild * 0.3, (int)maxChild * 0.3, p, objCnt, false);
	cout << "[Rtree build done]" << endl;

	////rtreeRAM(*rtree, ramTree);
	//totalSpaceCost += ramTree.size() * 4096.00 / MB;
	vector<long int> skyband, klayers, min_level;
	kskyband(dim, *rtree, skyband, PointSet, Qk_max);
	cout << "skyband size: " << skyband.size() << endl;

	onionlayer(skyband, PointSet, Qk_max, klayers, min_level, dim);
	cout << "# k-onion layer records / # of k-skyband: " << klayers.size() << "/" << skyband.size() << endl;


	for (int i = 0; i < klayers.size(); i++) {
		TotalObj[i] = new float[dim];
		for (int j = 0; j < dim; j++) {
			TotalObj[i][j] = (PointSet[klayers[i]][j] + PointSet[klayers[i]][j + dim]) / 2.0;
		}
	}

	All_records.clear();
	for (int i = 0; i < klayers.size(); i++) {
		//if (i >= 138) break;
		//Obj NewObj(klayers[i], dim, TotalObj[i]);
		Obj NewObj(i, dim, TotalObj[i]);
		All_records.push_back(NewObj);
	}

	fstream fout(recordfile, ios::out);
	for (int i = 0; i < All_records.size(); i++) {
		//if (i >= 138) break;
		fout << All_records[i].ID << ' ' << min_level[i] << endl;
	}
	fout.close();

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

void UTK_BaselineQueryProcess(Rtree* rtree, vector<float>& Qregion, int Qk, int Qdim) {
	clock_t at = clock();
	unordered_map<int, cell*> utkRet;
	utkRet.clear();
	UTK* sol = new UTK();
	vector<long int> rskyband;
	vector<cell*> exactutk;
	sol->rskyband(Qregion, Qdim, *rtree, rskyband, PointSet, Qk); // filter
	/*for (int i = 0; i < rskyband.size(); i++) {
		for (int d = 0; d < Qdim; d++) {
			cout << PointSet[rskyband[i]][d] << ' ';
		}
		cout << endl;
	}*/

	//cout << "k vs. r-skyband size: " << skyband.size() << "; " << rskyband.size() << endl;
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

}
void UTK_BaselineSolution(vector<vector<float>>& Qregions, int dim) {
	// build rtree
	cout << "Bulkloading R-tree..." << endl;
	const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
	FileMemory mem(PAGESIZE, rtreefile, RtreeNodeEntry::fromMem, true);
	Rtree * rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild * 0.3, (int)maxChild * 0.3, p, objCnt, false);
	cout << "[Rtree build done]" << endl;

	clock_t at = clock();
	for (int i = 0; i < Qregions.size(); i++) {
		cout << "Query: " << i << endl;
		//UTK_BaselineQueryProcess(rtree, Qregions[i], Q_k, dim);
		UTK_BaselineQueryProcess(rtree, Qregions[i], 30, dim);
	}
	cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
}


//void UTK_GIndexProcess(KLEVEL& idx, vector<float>& Qregion, int Qk, int Qdim) {
//	
//	int cnt = 0;
//	set<int> utk; utk.clear();
//
//	for (int i = 0; i < idx.L[Qk].size(); i++) {
//		if (KLEVEL::Intersect(Qregion, idx.L[Qk][i], Qdim)) {
//
//			for (auto it = idx.L[Qk][i].confirmed.begin(); it != idx.L[Qk][i].confirmed.end(); it++) {
//				utk.insert(*it);
//			}
//
//			cnt++;
//		}
//	}
//
//	for (int l = Qk + 1; l <= idx.Ik_max; l++) {
//		for (int i = 0; i < idx.L[l].size(); i++) {
//			if (idx.L[l][i].objID != -1) continue;
//			if (idx.L[l][i].confirmed.size() >= Qk) continue;
//			if (KLEVEL::Intersect(Qregion, idx.L[l][i], Qdim)) {
//
//				cnt = cnt + KLEVEL::GroupingProcess(All_records, Qk, Qdim, Qregion, idx.L[l][i], utk);
//			}
//		}
//	}
//	cout << "Regions: " << cnt << endl;
//	cout << "UTK: " << utk.size() << endl;
//
//}
void UTK_IndexSolution(vector<vector<float>>& Qregions, int dim, int Qk_max, int Ik_max) {
	fstream fout("E:\\klevel\\Project\\klevel_index\\info.txt", ios::out);
	KLEVEL idx(dim, Ik_max, Qk_max);

	if (IndexBuilding) {
		Skyband_Onionlayer(idx, dim, Qk_max);
		PreComputeAllHP(dim);
		idx.Bulk_Loading_singleObj(fout, All_records, dim, klevelfile);
		//idx.Bulk_Loading_cplex(fout, All_records, dim, klevelfile);
		//idx.WriteToDisk(klevelfile);
		cout << "Building Done" << endl;
		return;
	}
	else {
		Skyband_Onionlayer(idx, dim, Qk_max);
		//PreComputeAllHP(dim);
		idx.ReadFromDisk(All_records, klevelfile);
	}

	clock_t at = clock();

	for (int i = 0; i < Qregions.size(); i++) {
		cout << "Query " << i << endl;
		//UTK_IndexProcess(idx, Qregions[i], Q_k, dim);
		UTK_IndexProcess(idx, Qregions[i], 30, dim);
	}

	cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;

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
	fstream fout("E:\\klevel\\Project\\klevel_index\\info.txt", ios::out);
	KLEVEL idx(dim, Ik_max, Qk_max);

	if (IndexBuilding) {
		Skyband_Onionlayer(idx, dim, Qk_max);
		PreComputeAllHP(dim);
		idx.Bulk_Loading(fout, All_records, dim, klevelfile);
		//idx.WriteToDisk(klevelfile);
		cout << "Building Done" << endl;
		return;
	}
	else {
		Skyband_Onionlayer(idx, dim, Qk_max);
		//PreComputeAllHP(dim);
		idx.ReadFromDisk(All_records, klevelfile);
	}

	clock_t at = clock();

	for (int i = 0; i < Q_num; i++) {
		cout << "Query " << i << endl;
		kSPR_IndexProcess(idx, Qids[i], 50, dim);
	}

	cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
}

int kSPR_BaselineQueryProcess(Rtree* rtree, int Qid, int Qk, int Qdim) {

	//clock_t at = clock();
	//vector<long int> skyband, newAddSL;
	//vector<cell*> finalResult;

	//Point Qpt;
	//Qpt.m_dimen = Qdim;
	//for (int i = 0; i < Qdim; i++) Qpt.m_coor[i] = PointSet[Qid][i];

	//int insertCount = 0;
	//Getkskyband(Qdim, *rtree, skyband, Qpt, PointSet, Qk);
	////totalNoOfHPs += kskyband.size();
	//computeHPforkskyband(Qdim, PointSet, Qpt, skyband, newAddSL);
	//cout << "k-skyband: " << skyband.size() << "; precentage: " << 1.0 * skyband.size() / objCnt << endl;
	//cellTree* sol = new cellTree();
	//sol->insert_kspr(newAddSL, Qk, *rtree, Qpt, finalResult);
	//sol->collectLeaf_kspr(finalResult, Qk);
	//printf("Result Size: %d\n", finalResult.size());
	//cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;

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

	// count dominators
	int Noofdominator = countDominator(*rtree, PointSet, Qpt, RecordEntry);

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
	fstream fout("E:\\klevel\\Project\\klevel_index\\info.txt", ios::out);
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
	for (int i = 0; i < Q_num; i++) {
		clock_t at_this=clock();
		cout << "Query: " << i << endl;
		fout << "Query: " << i << endl;
		//int cnt = kSPR_BaselineQueryProcess(rtree, Qids[i], Q_k, dim);
		int cnt = kSPR_BaselineQueryProcess(rtree, Qids[i], 50, dim);
		cout << "Time cost: " << fixed << (clock() - at_this) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
		fout << "Result Number:" << cnt << endl;
		fout << "Time cost: " << fixed << (clock() - at_this) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
	}
	cout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
	fout << "Total time cost: " << fixed << (clock() - at) * 1.0 / CLOCKS_PER_SEC << " SEC " << endl;
}

int main(const int argc, const char** argv) {

	cout << "Parse Parameters" << endl;

	int dim = stoi(ReadParameter(argc, argv, "-d", ""));
	string datafile = ReadParameter(argc, argv, "-f", "");
	
	DataLoading(dim, datafile);

	fstream fout("E:\\klevel\\Project\\klevel_index\\info.txt", ios::out);
	KLEVEL idx(dim, 50, 50);
	Skyband_Onionlayer(idx, dim, 50);
	PreComputeAllHP(dim);
	idx.Bulk_Loading_cplex(fout, All_records, dim, klevelfile);

	//if (Q_TYPE == "UTK") {
	//	float sigma = Q_sigma;
	//	vector<vector<float>> Qregions;
	//	UTK_QueryGenerator(Qregions, dim - 1, sigma);

	//	//UTK_BaselineSolution(Qregions, dim);
	//	UTK_IndexSolution(Qregions, dim, Q_k, I_k);
	//}
	//else if (Q_TYPE == "kSPR") {
	//	vector<int> Qids;
	//	kSPR_QueryGenerator(Qids);
	//	kSPR_BaselineSolution(Qids, dim, Q_k);
	//	//kSPR_IndexSolution(Qids, dim, Q_k, I_k);
	//}
	
	system("pause");
}

//void qhulltest() {
//	int dim = DIM;             /* dimension of points */
//	int numpoints;            /* number of points */
//	coordT points[(DIM + 1) * TOTpoints]; /* array of coordinates for each point */
//	coordT* rows[TOTpoints];
//	boolT ismalloc = False;    /* True if qhull should free points in qh_freeqhull() or reallocation */
//	char flags[250];          /* option flags for qhull, see qh-quick.htm */
//	FILE* outfile = stdout;    /* output from qh_produce_output()
//								 use NULL to skip qh_produce_output() */
//	FILE* errfile = stderr;    /* error messages from qhull code */
//	int exitcode;             /* 0 if no error from qhull */
//	facetT* facet;            /* set by FORALLfacets */
//	int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
//	int i;
//
//	qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
//	qhT* qh = &qh_qh;
//
//	QHULL_LIB_CHECK
//
//		qh_zero(qh, errfile);
//
//
//
//
//
//	/*
//	  Run 3: halfspace intersection about the origin
//	*/
//	printf("\n========\ncompute halfspace intersection about the origin for a diamond\n");
//	//sprintf(flags, "qhull H0 s Tcv %s", argc >= 4 ? argv[3] : "Fp");
//	numpoints = SIZEcube;
//	makehalf(points, numpoints, dim);
//	for (i = numpoints; i--; )
//		rows[i] = points + (dim + 1) * i;
//	qh_printmatrix(qh, outfile, "input as halfspace coefficients + offsets", rows, numpoints, dim + 1);
//	fflush(NULL);
//	/* use qh_sethalfspace_all to transform the halfspaces yourself.
//	   If so, set 'qh->feasible_point and do not use option 'Hn,...' [it would retransform the halfspaces]
//	*/
//	sprintf(flags, "qhull H0 s Tcv %s", "Fp");
//	coordT * feasible_point = new coordT[dim - 1];
//	for (int i = 0; i < dim - 1; i++) feasible_point[i] = 0.001;
//	exitcode = qh_new_qhull_jiahao(qh, dim + 1, numpoints, points, ismalloc,
//		flags, feasible_point, outfile, errfile);
//
//	fflush(NULL);
//	if (!exitcode)
//		print_summary(qh);
//#ifdef qh_NOmem
//	qh_freeqhull(qh, qh_ALL);
//#else
//	qh_freeqhull(qh, !qh_ALL);
//	qh_memfreeshort(qh, &curlong, &totlong);
//	if (curlong || totlong)  /* could also check previous runs */
//		fprintf(stderr, "qhull internal warning (user_eg, #3): did not free %d bytes of long memory (%d pieces)\n",
//			totlong, curlong);
//#endif
//}
//void IndexBuilding_BB(int dim, int k, string indexfile) {
//	vector<long int> klayers, min_level;
//	kskylinekonionlayer_Filter(klayers, min_level, dim, k, indexfile);
//
//	KLEVEL klevel_index(dim, k);
//
//	clock_t at = clock();
//	for (int i = 0; i < All_records.size(); i++) {
//		cout << "Insert: " << i << endl;
//		klevel_index.Insert(All_records[i], min_level[i]);
//	}
//	cout << "BB Insert Time Cost: " << (clock() - at) / (float)CLOCKS_PER_SEC << endl;
//	system("pause");
//
//	// 2D Example
//	/*for (int i = 0; i < 5; i++)skybandObj[i] = new float[2];
//	skybandObj[0][0] = 0.4; skybandObj[0][1] = 0.8;
//	skybandObj[1][0] = 0.6; skybandObj[1][1] = 0.3;
//	skybandObj[2][0] = 0.2; skybandObj[2][1] = 0.9;
//	skybandObj[3][0] = 1.0; skybandObj[3][1] = 0.2;
//	skybandObj[4][0] = 0.6; skybandObj[4][1] = 0.2;
//	dim = 2; k = 4;*/
//
//	// 3D Example
//	/*for (int i = 0; i < 3; i++)skybandObj[i] = new float[3];
//	skybandObj[0][0] = 0.8; skybandObj[0][1] = 0.3; skybandObj[0][2] = 0.3;
//	skybandObj[1][0] = 0.4; skybandObj[1][1] = 0.8; skybandObj[1][2] = 0.7;
//	skybandObj[2][0] = 0.1; skybandObj[2][1] = 0.9; skybandObj[2][2] = 0.1;
//	dim = 3; k = 3;*/
//
//}