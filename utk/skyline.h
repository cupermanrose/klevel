#ifndef _SKYLINE_H_
#define _SKYLINE_H_

#include "header.h"
#include "rtree.h"
#include "rnode.h"
#include "rentry.h"
#include "virtualRNode.h"
#include "filemem.h"

#define MAXPAGEID 49999999

//define pageid and virtual rtree node pair
typedef multimap<long int, VirtualRNode*>::value_type PintVRN;

//define float and int vaule pair
typedef multimap<float, long int>::value_type PfltINT;

typedef map<int, int>::value_type PintInt;

//define cell range for two dim case
typedef multimap<float, float>::value_type PfltFLT;


float minDist(float p1[], float p2[], int dimen);


bool IsDominatedBy(const int dimen, const float pt[], vector<long> a_skylines, float* PG[]);


void updateSkylines(const int dimen, Rtree& a_rtree, unordered_set<long int>& a_skylines, float* PG[], multimap<float, int>& heap, multimap<long int, RtreeNodeEntry*>& RecordEntry);

int countkDominator(const int dimen, const float pt[], vector<long> kskyband, float* PG[]);

void computeHPforkskyband(const int dimen, float* PG[], Point& a_pt, vector<long int>& kskyband, vector<long int>& newAddSL);

void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k);

bool isDominateByFocal(const int dimen, const float pt[], Point& focal);

void Getkskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, Point& a_pt, float* PG[], const int k);

void onionlayer(vector<long int>& skyband, float* PG[], const int k, vector<long int >& klayers, vector<long int >& min_level, int& dim);

void computeHP(const int dimen, float* PG[], Point& a_pt, unordered_set<long int>& initSkyline, vector<long int>& newAddSL);

void initHS(const int dimen, vector<float>& regions);

void rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree);

int countRecords(Rtree& a_rtree, int pageID);

void aggregateRecords(Rtree& a_rtree);

int countDominator(Rtree& a_rtree, float* PG[], Point& a_pt, multimap<long int, RtreeNodeEntry*>& RecordEntry);

// for klevel_index
void kskyband(const int dimen, Rtree& a_rtree, vector<int>& kskyband, vector<vector<float>>& OriginData, const int k);

int countkDominator(const int dimen, const float pt[], vector<int> kskyband, vector<vector<float>>& OriginData);

void onionlayer(vector<int>& skyband, vector<vector<float>>& OriginData, const int k, vector<int>& klayers, int& dim);

#endif
