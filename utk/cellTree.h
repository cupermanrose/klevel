/* ----------------------------------------------------------------------------
This header file includes class Advanced RegionTree declaration.
It is data structure to maintain the cells in query plane with hyperplanes
---------------------------------------------------------------------------- */

#ifndef _CELL_TREE_H_
#define _CELL_TREE_H_

#include "header.h"
#include "lp_lib.h"
#include "rtree.h"
#include "rentry.h"
#include "rnode.h"
#include "filemem.h"
#include "global.h"

extern int objCnt;

typedef struct rBound
{
	int minimum;
	int maximum;

	rBound()
	{
		minimum = 0;
		maximum = objCnt - 1;
	}
};

typedef struct cell
{
	unordered_map<long int, bool> IntersectHP;
	unordered_set<long int> AboveHP;
	unordered_set<long int> BelowHP;

	rBound lu;
	int rank;
	bool isPruned;

	unordered_set<long int> focaled;
	unordered_set<long int> ignoreset;

	cell* left;
	cell* right;

	cell()
	{
		left = NULL;
		right = NULL;
		isPruned = false;
		rank = 0;
	}

	cell(cell* copy)
	{
		if (copy->isPruned == false)
		{
			left = NULL;
			right = NULL;
			rank = copy->rank;
			isPruned = copy->isPruned;
		}
	}

	void appendto(cell* node)
	{
		for (auto iter = IntersectHP.begin(); iter != IntersectHP.end(); iter++)
		{
			node->IntersectHP[iter->first] = iter->second;
		}

		for (auto iter = AboveHP.begin(); iter != AboveHP.end(); iter++)
		{
			node->AboveHP.insert(*iter);
		}

		for (auto iter = BelowHP.begin(); iter != BelowHP.end(); iter++)
		{
			node->BelowHP.insert(*iter);
		}
		node->rank = rank;
		node->isPruned = isPruned;
	}

	void copyleaf(cell* node)
	{
		for (auto iter = IntersectHP.begin(); iter != IntersectHP.end(); iter++)
		{
			node->IntersectHP[iter->first] = iter->second;
		}
		for (auto iter = BelowHP.begin(); iter != BelowHP.end(); iter++)
		{
			node->BelowHP.insert(*iter);
		}
		for (auto iter = focaled.begin(); iter != focaled.end(); iter++)
		{
			node->focaled.insert(*iter);
		}
		for (auto iter = ignoreset.begin(); iter != ignoreset.end(); iter++)
		{
			node->ignoreset.insert(*iter);
		}
		node->rank = rank;
		node->isPruned = isPruned;
	}

	void release()
	{
		IntersectHP.clear();
		unordered_map<long int, bool>().swap(IntersectHP);
		AboveHP.clear();
		unordered_set<long int>().swap(AboveHP);
		BelowHP.clear();
		unordered_set<long int>().swap(BelowHP);
	}

	double cellSize()
	{
		double ret = 0;
		ret += IntersectHP.size() * 5;
		ret += AboveHP.size() * 4;
		ret += BelowHP.size() * 4;
		if (left != NULL)
			ret += 4;
		else if (right != NULL)
			ret + 4;
		return ret*1.0 / MB;
	}
};

extern unordered_map<long int, cell*> cellID;
extern unordered_map<long int, RtreeNode*> ramTree;

typedef struct gNode
{
	long int rID;
	unordered_set<long int> rDominator;
	gNode(long int& recordID) :rID(recordID){};
};


struct cellCompare
{
	bool operator()(const cell* a, const cell* b) const
	{
		return a->rank < b->rank;
	}
};

struct cellpairCompare
{
	bool operator()(const pair<cell*, long int> &a, const pair<cell*, long int> &b) const
	{
		return a.first->rank < b.first->rank;
	}
};

typedef struct Score
{
	float min;
	float max;
	Score() 
	{
		min = 0;
		max = 1;
	}
};

class cellTree
{
public:
	cellTree();
	cellTree(cell* node);
	cellTree(int size);
	~cellTree();

	void releaseCell(cell* node);

	void insert_kspr(vector<long int>& hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult);
	void insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult); // for utk
	bool isFeasible(unordered_map<long int, bool>& touchhs, long int hpid, bool sideindicator);
	void lpModel(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs);
	void addHP(lprec* model, long int hpid, bool sideindicator);

	void updateRank(cell* node, const int mink);
	void inserthp(long int & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs);

	void collectLeaf(vector<cell*>& leaves, const int& mink);
	void dfsTraversal(cell* node, cell* leaf, vector<cell*>& leaves);


	void opt_insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult);
	void opt_inserthp(long int & hpid, const int mink, cell* node, cell* all);
	int isDominatorInserted(unordered_set<long int>& rdominators, cell* all);

	void maintainDAG(unordered_set<long int>& skylines, unordered_set<long int>& removeSL, float* PG[], const int dimen);
	void markSingular_kspr(vector<cell*>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular);
	void markSingular(vector<cell*>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular);
	void markSingular(vector<pair<cell*, long int>>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular);

	void updateCellTree(cell* node);

	int Lemma2(vector<cell*>& cells);

	void lpModelwithSIDELEN(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs);
	void addHPwithSIDELEN(lprec* model, long int hpid, bool sideindicator);
	void halfspace2polytope(cell* ret, string outfile, string hsfile);

	// look-ahead techniques

	void collectLeaf(vector<pair<cell*, long int>>& leaves, const int& mink);
	void collectLeaf_kspr(vector<cell*>& leaves, const int& mink);

	void dfsTraversal(cell* node, cell* leaf, vector<pair<cell*, long int>>& leaves);
	void updateSkyline(unordered_set<long int>& skylines, vector<pair<cell*, long int>>& leaves, unordered_set<long int>& singular);

	void rankBound(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult);
	void cellBound(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink);
	void recordScore(lprec* lp, RtreeNodeEntry* e, Score& rScore, bool& isLeafNode);
	void focalScore(lprec* lp, Point& pt, Score& pScore);

	void findCellMBR(lprec* lp, vector<float>& cl, vector<float>& cu);
	void approxiRecordScore(RtreeNodeEntry* e, Score& entryR, vector<float>& cl, vector<float>& cu, bool& isLeafNode);
	void rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);


	// optimized look-ahead techniques
	void rankBoundOpt(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult);
	void cellBoundOpt(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink);

	double lpOptObj(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax);
	void fastOptObj(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore);

	void rankBoundOptFULL(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult);
	void cellBoundOptFULL(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink);
	double lpOptObjFULL(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax);
	void fastOptObjFULL(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore);

	cell* root;
	double treeSize;

private:
	unordered_map<long int, gNode*> dagNode;
};
#endif