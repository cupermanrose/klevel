#ifndef _HYPERPLANE_H
#define _HYPERPLANE_H

#include "object.h"
#include "lp_lib.h"
#include "qhull_ra.h"
#include <fstream>
#include <vector>
#include <unordered_map>
#include "ilcplex/ilocplex.h"


using namespace std;

struct Hyperplane {
	vector<float> w;
	int o1, o2;
};

extern vector<Hyperplane> AllHP;

class Halfspace {

public:
	unordered_map<int, bool> HPs; // bool = false, o1<=o2; bool = true, o1>=o2
	vector<vector<float>> vertices;
	float innerPoint[Max_Dimension];

public:

	Halfspace();
	Halfspace(vector<float>& Qregion, int dim);
	~Halfspace();

	void lpModel(lprec* lp, int dim); // the lp model of this region
	void addHP(lprec* lp, int dim, int HP, bool sideindicator); // add a hyperplane into the lp model;
	//static bool is_Feasible(Halfspace& r, int dim, float range[]); // store the dimension range
	static bool is_Feasible(Halfspace& r, int dim); // r is non-empty
	//static bool is_Feasible(Halfspace& r, int dim, int HP, bool sideindicator); // r is non-empty after add HP
	//static bool is_Feasible(Halfspace& r_a, Halfspace& r_b, int dim); // the intersection of r_a and r_b is non-empty
	//static bool is_intersected(Halfspace& r, int dim, int HP); // HP and r intersect
	static int ComputeHP(Obj& o1, Obj& o2, int dim);
	//static void Split(Halfspace& r, int HP, Halfspace& r_LE, Halfspace& r_GE); // o1 <= o2 in r_LE; o1 >= o2 in r_GE
	//static void Intersect(Halfspace& r_s, Halfspace& r_a, Halfspace& r_b);
	//static bool Subset(Halfspace& r_a, Halfspace& r_b); // ture: r_b is a subset of r_a;

	static void ComputeVertexByQhull(int dim, Halfspace& r);

	//static bool is_Feasible_cplex(Halfspace& r, int dim);
};

#endif // !_HYPERPLANE_H