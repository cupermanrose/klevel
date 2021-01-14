#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "time.h"
#include <unordered_map>
#include <random>
#include <queue>
#include <unordered_set>
#include <map>
#include <set>
#include "string.h"
#include <sstream>
#include <iomanip>
#include <bitset>

#include "collection.h"

// Constant
#ifdef DEBUG
#define PAGESIZE 144 // for debugging
#else
#define PAGESIZE 4096   // macro define the pagesize (KB)
#endif
#define MAXPTS 10000001  // macro define the maximum number of points
#define MB 1048576      // macro define the size of MB

#define EPS 1e-6
#define TEST 0
#define cur_qk 15
//#define IndexBuilding TRUE
#define IndexBuilding FALSE

#define MAXREG 10000001  // macro define the maximum number of regions
#define eps_klevel 1e-6

// UTK parameters
#define LEAFNODE -1
#define NONLEAFNODE 1
#define ABOVE	1
#define BELOW	2
#define ON 3
#define HP_OVERLAPPED	4
#define SIDELEN 0.0001
#define CELLING 1.000

// utk query
#define Q_sigma 0.1

//#define TEST TRUE



//#define rtreefile "E:\\klevel\\Project\\klevel_index\\record.idx"
#define recordfile "E:\\klevel\\Project\\klevel_index\\All_records.txt"


// Prune Rules

// Check whether query record is in <=k-level from k-level to 1-level; 
#define Pruning1 TRUE
// r_q and G[cur].r do not have intersection if they have opposite halfspaces
// r_q and G[cur].r have intersection if they have two same halfspaces (when all halfspaces are minimum borders)
#define Pruning2 TRUE

#define QHULL TRUE
//#define GROUPING FALSE
#define Delta 2


using namespace std;
