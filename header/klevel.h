#ifndef _KLEVEL_H
#define _KLEVEL_H

#include <vector>
#include <set>


#define max_object 10000;

//#include "object.h"
#include "hyperplane.h"
//#include "hypercube.h"
//#include "skyline.h"
//#include "rentry.h"
//#include "rnode.h"
//#include "rtree.h"
//#include "tgs.h"
//#include "filemem.h"

using namespace std;

class KLEVEL {

public:

	struct region {
		int level, objID;
		Halfspace r;
		set<int> confirmed; // the top-kth in this region
		set<int> candidates; // candidates.first -> objID, candidates.second -> min level in this region
		set<int> grouping; // grouping options
		~region(){
		    r.vertices.clear();
            vector<vector<float>>().swap(r.vertices);
            confirmed.clear();
            candidates.clear();
            grouping.clear();
		}
	};

public:

	int dim, Ik_max, ave_size, Qk_max;
	vector<vector<region>> L;


public:
	KLEVEL();
	KLEVEL(int dimension, int Ikmax, int Qkmax);
	~KLEVEL();

    void WriteToDisk(string outfile);
    void WriteToDisk(string outfile, int i);
    void ReadFromDisk(vector<Obj>& All_records, string infile);

    // utk query
    static bool Intersect(vector<float>& Qregion, region& cur, int dim);
    static bool isIn(vector<float>& v, vector<float>& Qregion, int dim);
    static bool isIn(vector<float>& v, unordered_map<int, bool>& HPs, int dim);
    static bool CheckRepeatTopk(vector<region>& L, region& pre_r, int kth_obj, vector<int>& k_rskyband);
    static int UpdateLevel_forUTK(vector<Obj>& All_records, vector<region>& L, int k, int dim, vector<float>& Qregion, set<int>& utk);
    static void AddHP(vector<Obj>& All_records, region& cur, vector<pair<int, bool>>& QregionHP);

    //kspr query
    static int UpdateLevel_forkSPR(vector<Obj>& All_records, vector<region>& L, int k, int dim, int Qid);


    bool XdominateY(Obj& o, Obj& q, int dim);
    int ComputeMinLevel(int q, vector<Obj>& All_records);
    void GetMinLevel(vector<Obj>& All_records, vector<int>& min_level);

    void Bulk_Loading_Init(vector<Obj>& All_records);
    static void krskyband(Halfspace& r, int a_k, vector<Obj>& All_records, set<int>& candidates, vector<vector<int>>& skyband_list, vector<vector<int>>& dominate_list, int dim);
    static int isRdominated(vector<vector<float>>& vertices, float focal[], float entry[], int dim); //X is dominated by Y in R;

    // cplex

    void Bulk_Loading_cplex(fstream& fout, vector<Obj>& All_records, int dim, string indexfile);
    static void CreateNewRegion(vector<Obj>& All_records, region& tmp, int level, int cur, region& pre_r, vector<int>& one_rskyband, vector<int>& k_rskyband);
    static void AddHP(vector<Obj>& All_records, region& cur);
    bool CheckRepeatTopk(int level, region& pre_r, int kth_obj, vector<int>& k_rskyband);

    void AddHP_grouping(vector<Obj>& All_records, region& cur);
    //static void AddHP_cplex(vector<Obj>& All_records, region& cur);
    //bool CheckRepeatTopk_cplex(int level, region& pre_r, int kth_obj, vector<int>& k_rskyband);
    //static void CreateNewRegion_cplex(vector<Obj>& All_records, region& tmp, int level, int cur, region& pre_r, vector<int>& one_rskyband, vector<int>& k_rskyband);


    //bool CheckRepeatTopk(int level, region& cur);


	//void Bulk_Loading(fstream& fout, vector<Obj>& All_records, int dim, string indexfile);

	//int LevelGrouping_step(int k, vector<int>& candidates, vector<int>& total_candidates, region& this_r, vector<Obj>& All_records, vector<int>& dominators);
	//int ComputeMinLevel(int q, vector<Obj>& All_records);


	//static int GroupingProcess(vector<Obj>& All_records, int k, int dim, vector<float>& Qregion, region& r_st, set<int>& results);
	
	//static void rskyband_singlek(Halfspace& r, int a_k, vector<Obj>& All_records, set<int>& candidates, vector<int>& skyband, int dim);



	//bool countRegionDominator(int k, float pt[], vector<int>& rskyband_set, float* PG[], float r[], unordered_set<int>& dominators);
	//bool countRegionDominator(float r[], vector<Obj>& All_records, Obj& cur_obj, vector<int>& rskyband_set, int k, unordered_set<int>& dominators);
	//void rskyband(float** PG, float r[], vector<Obj>& All_records, set<int>& candidates, int a_k, vector<int>& rskyband_set);
	//void rskyband(Rtree* tmp_rtree, float** PG, float r[], vector<Obj>& All_records, int a_k, vector<int>& rskyband_set);
	//void rskyband(float r[], vector<Obj>& All_records, int a_k, vector<int>& rskyband_set);
	//void rskyband(float r[], vector<Obj>& All_records, set<int>& candidates, int a_k, vector<int>& rskyband_set);
	
	
	//void rskyband(Halfspace& r, vector<Obj>& All_records, set<int>& candidates, int a_k, vector<int>& rskyband_set);
	//void rskyband(Halfspace& r, vector<Obj>& All_records, set<int>& candidates, int a_k, vector<int>& rskyband_set, vector<int>& dominators);
	
	
	//void One_rskyband_Filter(Halfspace& r, vector<Obj>& All_records, vector<int>& one_rskyband);
	
	//bool MergeLevel(region_rd& r, vector<Obj>& All_records, set<int>& utk_results);
	//int LevelGrouping(int level, vector<region_rd>& this_level, region_rd& r, vector<Obj>& All_records, set<int>& utk_results, vector<int>& k_rskyband);
	//int LevelGrouping(int k, vector<int>& candidates, vector<int>& dominators, region_rd& this_r, vector<Obj>& All_records);
	
	// Full-level
	/*int CreateNode_full(int ID, Halfspace& r);
	void Insert_full(Obj& q);
	void Update_SplitNode_full(int cur, Halfspace& r_LE, Halfspace& r_GE);
	void Update_NewNode_full(int cur, vector<int>& aboveNodes, vector<int>& belowNodes);
	void FindKlevel_full(vector<int>& Klevel, int k);
	int FindFirstNode_in_Klevel_full(int k);
	void FindVoronoiWithRank_full(Obj& q, vector<Halfspace>& Rq, vector<int>& rank_q, vector<int>& HP_q);*/

	// Max_k-level
	//bool CheckRepeatTopk(region& r1, region& r2, int kth_obj);
	//void Bulk_Loading(vector<Obj>& All_records, vector<long int>& min_level);
	//
	

	//float orderScore(vector<float>& pivot, float entry[]);
	//float orderScore(float pivot[], float entry[]);
	//void pivotRegion(float r[], vector<float>& pivot);
	//bool isRdominated(float r[], float focal[], vector<float>& entry, bool& fDe);
	////bool isRdominated(vector<vector<float>>& vertices, float focal[], float entry[]); //X is dominated by Y in R;
	
	//void SpaceSplit(vector<Obj>& All_records, int split_size, int st[], float offset, int dim, int mk);
	//void AllSubSpaceBuilding(vector<Obj>& All_records, int split_size, int dim, int mk);
	//void SingleSubSpaceBuilding(vector<Obj>& All_records, SubSpace& space, int dim, float offset, int st[], int mk);
	//void testing(vector<Obj>& All_records, vector<long int>& min_level, int dim);
	
	
	//int CreateNode(int ID, Halfspace& r);
	//void Insert(Obj& q, int min_level);
	////void Insert(Obj& q);
	//int Update_SplitNode(int cur, int HP);
	//void Update_SplitNode_lastLevel(int cur, int HP);
	//void Update_NewNode(int cur, vector<int>& aboveNodes, vector<int>& belowNodes);
	//void FindKlevel(int level_limit); // return false if inserted less than Max_k
	//void FindVoronoiWithRank(Obj& q, vector<Halfspace>& Rq, vector<int>& rank_q, vector<int>& HP_q);
	//void DeleteNode(int cur);
	//void CheckLinkERROR();
	//float CompRadius(Obj& o);
};

#endif // !_KLEVEL_H

