#include "skyline.h"

extern vector<vector<float>> HalfSpaces;
//extern unordered_map<long int, long int> RecordIDtoHalfPlaneID;
extern unordered_map<long long, long int> RecordIDtoHalfPlaneID;
extern double totalIO;
extern unordered_map<long int, RtreeNode*> ramTree;

float minDist(float p1[], float p2[], int dimen)
{
	float mindist = 0;

	for (int i = 0; i < dimen; i++)
	{
		float dist = p1[i] - p2[i];
		mindist += (dist * dist);
	}
	return (float)sqrt(mindist);
}

bool IsDominatedBy(const int dimen, const float pt[], vector<long> a_skylines, float* PG[])
{
	vector<long>::iterator iter;
	if (a_skylines.size() == 0)
		return false;

	for (iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		long pid = *iter;
		float s[MAXDIMEN];
		bool dominated = true;
		for (int i = 0; i < dimen; i++)
		{
			if (PG[pid][i] + eps_klevel < pt[i])
			{
				dominated = false;
				break;
			}
		}
		if (dominated)
			return dominated;
	}
	return false;
}



int countkDominator(const int dimen, const float pt[], vector<long> kskyband, float* PG[])
{
	vector<long>::iterator iter;
	if (kskyband.size() == 0)
		return false;

	int count = 0;
	for (iter = kskyband.begin(); iter != kskyband.end(); iter++)
	{
		long pid = *iter;
		float s[MAXDIMEN];
		bool dominated = true;
		for (int i = 0; i < dimen; i++)
		{
			if (PG[pid][i] + eps_klevel < pt[i])
			{
				dominated = false;
				break;
			}
		}
		if (dominated)
			count++;
	}
	return count;
}
int countDominator(Rtree& a_rtree, float* PG[], Point& a_pt, multimap<long int, RtreeNodeEntry*>& RecordEntry)
{
	queue<long int> H;
	float pt[MAXDIMEN];
	float rd[MAXDIMEN];
	int dimen = a_rtree.m_dimen;
	int NoOfDominators = 0;
	int NoofDominees = 0;
	RtreeNode* node;
	RtreeNodeEntry* e0;
	long int NegPageid, pageid;

	float cl[MAXDIMEN], cu[MAXDIMEN];
	for (int d = 0; d < dimen; d++)
	{
		cl[d] = 0;
		cu[d] = a_pt[d];
	}
	Hypercube Dominee_hc(dimen, cl, cu);
	for (int d = 0; d < dimen; d++)
	{
		cl[d] = a_pt[d];
		cu[d] = 1;
	}
	Hypercube Dominator_hc(dimen, cl, cu);

	H.push(a_rtree.m_memory.m_rootPageID);
	while (!H.empty())
	{
		pageid = H.front();
		H.pop();

		node = a_rtree.m_memory.loadPage(pageid);
		totalIO++;
		if (node->isLeaf())
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				for (int d = 0; d < dimen; d++)
					rd[d] = (PG[node->m_entry[i]->m_id][d] + PG[node->m_entry[i]->m_id][dimen + d]) / 2.0;
				Point tmpPt(dimen, rd);
				if (tmpPt == a_pt)
					continue;
				else if (Dominator_hc.enclose(tmpPt) == true)   //current point lies in the dominator window
					NoOfDominators++;
				else if (Dominee_hc.enclose(tmpPt) == false)
				{
					NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
					RtreeNodeEntry* Nentry = node->m_entry[i]->clone();
					RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
				}
				else if (Dominee_hc.enclose(tmpPt))
				{
					NoofDominees++;
				}
			}
		}
		else
		{
			e0 = node->genNodeEntry();
			if (Dominator_hc.enclose(e0->m_hc) == true)
			{
				for (int i = 0; i < node->m_usedspace; i++)
					NoOfDominators += node->m_entry[i]->num_records;
			}
			else if (Dominee_hc.enclose(e0->m_hc) == false)
			{
				for (int i = 0; i < node->m_usedspace; i++)
					H.push(node->m_entry[i]->m_id);
			}
			else if (Dominee_hc.enclose(e0->m_hc) == true)
			{
				for (int i = 0; i < node->m_usedspace; i++)
					NoofDominees += node->m_entry[i]->num_records;
			}
			delete e0;
		}
		delete node;
	}
	return NoOfDominators;
}

void updateSkylines(const int dimen, Rtree& a_rtree, unordered_set<long int>& a_skylines, float* PG[], multimap<float, int>& heap, multimap<long int, RtreeNodeEntry*>& RecordEntry)
{
	vector<long int> newskyline;
	RtreeNode* node;
	//RtreeNodeEntry* e0;
	multimap<long int, RtreeNodeEntry*>::iterator recordIter;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;
	multimap<float, int>::iterator heapIter;

	for (auto iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		newskyline.push_back(*iter);
	}

	node = a_rtree.m_memory.loadPage(a_rtree.m_memory.m_rootPageID);
	totalIO++;
	if (node->isLeaf())
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
				//pt[j] = (node->m_entry[i]->m_hc.getLower()[j] + node->m_entry[i]->m_hc.getUpper()[j]) / 2;
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
		}
	}
	else
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
		}
	}
	delete node;

	while (heap.size() != 0)
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			recordIter = RecordEntry.find(pageID);
			if (recordIter != RecordEntry.end())
			{
				for (int d = 0; d < dimen; d++)
				{
					pt[d] = recordIter->second->m_hc.getLower()[d] + SIDELEN;
					//pt[d] = (recordIter->second->m_hc.getLower()[d] + recordIter->second->m_hc.getUpper()[d]) / 2;
				}
				if (IsDominatedBy(dimen, pt, newskyline, PG) == false)
				{
					newskyline.push_back(pageID - MAXPAGEID);
				}
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			totalIO++;
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
					{
						pt[d] = node->m_entry[i]->m_hc.getLower()[d] + SIDELEN;
						//pt[d] = (node->m_entry[i]->m_hc.getLower()[d] + node->m_entry[i]->m_hc.getUpper()[d]) / 2;
					}
					if (IsDominatedBy(dimen, pt, newskyline, PG) == false)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
					if (IsDominatedBy(dimen, pt, newskyline, PG) == false)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
					}
				}
			}
			delete node;
		}
	}


	a_skylines.clear();
	for (int i = 0; i < newskyline.size(); i++)
	{
		recordIter = RecordEntry.find(newskyline[i] + MAXPAGEID);
		if (recordIter != RecordEntry.end())
		{
			RecordEntry.erase(recordIter);
		}
		a_skylines.insert(newskyline[i]);
	}
}
void computeHP(const int dimen, float* PG[], Point& a_pt, unordered_set<long int>& initSkyline, vector<long int>& newAddSL)
{
	newAddSL.clear();
	for (auto sIter = initSkyline.begin(); sIter != initSkyline.end(); sIter++)
	{
		long int id = (*sIter);
		if (RecordIDtoHalfPlaneID.find(id) == RecordIDtoHalfPlaneID.end()) // if this record has hyperplane in halfspaces, then we skip it.
		{
			float rd[MAXDIMEN];
			int tmpDim = 0;
			for (int d = 0; d < dimen; d++)
			{
				rd[d] = (PG[id][d] + PG[id][dimen + d]) / 2.0;
				if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
			}
			if (tmpDim == dimen)
			{
				continue;
			}

			// Reduce d dimension record to d-1 halfspace;
			vector<float> tmpHS;
			float rd_d = rd[dimen - 1];
			float p_d = a_pt[dimen - 1];
			for (int d = 0; d < dimen - 1; d++)
			{
				tmpHS.push_back((rd[d] - rd_d) - (a_pt[d] - p_d));
			}
			tmpHS.push_back(p_d - rd_d);
			tmpHS.push_back(id);            //store the ID of incomparable record in HalfSpace
			HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p	
			RecordIDtoHalfPlaneID.insert(PintInt(id, HalfSpaces.size()));
			newAddSL.push_back(id);
		}
	}
}

void initHS(const int dimen, vector<float>& regions)
{
	for (int i = 0; i < regions.size(); i++)
	{
		vector<float> tmpHS;
		int loc = i / 2;
		int flag = i % 2;

		for (int di = 0; di < dimen - 1; di++)
		{
			if (di == loc)
			{
				if (flag == 1)
				{
					tmpHS.push_back(1.0);
				}
				else
				{
					tmpHS.push_back(-1.0);
				}
				
			}
			else
			{
				tmpHS.push_back(0.0);
			}
		}
		if (flag == 1)
		{
			tmpHS.push_back(regions[i]);
		}
		else
		{
			tmpHS.push_back(-regions[i]);
		}
		tmpHS.push_back(-(1 + i)); // init R extent as -1
		HalfSpaces.push_back(tmpHS); // store to Halfspaces
		RecordIDtoHalfPlaneID.insert(PintInt(-(1 + i), HalfSpaces.size()));
	}
}

void computeHPforkskyband(const int dimen, float* PG[], Point& a_pt, vector<long int>& kskyband, vector<long int>& newAddSL)
{
	newAddSL.clear();
	int docount;
	int decount;

	for (int id = 0; id < kskyband.size(); id++)
	{
		docount = 0;
		decount = 0;

		float rd[MAXDIMEN];
		int tmpDim = 0;
		for (int d = 0; d < dimen; d++)
		{
			rd[d] = (PG[kskyband[id]][d] + PG[kskyband[id]][dimen + d]) / 2.0;
			if (fabs(rd[d] - a_pt[d]) <= 1e-6) tmpDim++;
		}
		if (tmpDim == dimen)
		{
			continue;
		}


		for (int di = 0; di < dimen; di++)
		{
			if (rd[di] < a_pt.m_coor[di])
			{
				docount++;
			}
			else if (rd[di] > a_pt.m_coor[di])
			{
				decount++;
			}
		}

		if (docount != dimen && decount != dimen)
		{
			// Reduce d dimension record to d-1 halfspace;
			vector<float> tmpHS;
			float rd_d = rd[dimen - 1];
			float p_d = a_pt[dimen - 1];
			for (int d = 0; d < dimen - 1; d++)
			{
				tmpHS.push_back((rd[d] - rd_d) - (a_pt[d] - p_d));
			}
			tmpHS.push_back(p_d - rd_d);
			tmpHS.push_back(kskyband[id]);            //store the ID of incomparable record in HalfSpace
			HalfSpaces.push_back(tmpHS);          //form the half-space defined by the incomparable record and p	
			RecordIDtoHalfPlaneID.insert(PintInt(kskyband[id], HalfSpaces.size()));
			newAddSL.push_back(kskyband[id]);
		}
	}
}

void kskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, float* PG[], const int k)
{
	RtreeNode* node;
	multimap<float, int> heap;
	multimap<float, int>::iterator heapIter;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;

	heap.insert(PfltINT(INFINITY, a_rtree.m_memory.m_rootPageID));

	while (!heap.empty())
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			for (int d = 0; d < dimen; d++)
				pt[d] = (PG[pageID - MAXPAGEID][d] + PG[pageID - MAXPAGEID][d + dimen]) / 2;
			if (countkDominator(dimen, pt, kskyband, PG) <= k)
			{
				kskyband.push_back(pageID - MAXPAGEID);
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
					{
						pt[d] = node->m_entry[i]->m_hc.getLower()[d] + eps_klevel;
					}
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
					}
				}
			}
		}
	}
}

void rtreeRAM(Rtree & rt, unordered_map<long int, RtreeNode*> & ramTree)
{
	ramTree.clear();
	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = rt.m_memory.loadPage(pageID);
		ramTree[pageID] = node;
		if (node->isLeaf() == false)
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				H.push(node->m_entry[i]->m_id);
			}
		}
	}
}

int countRecords(Rtree& a_rtree, int pageID)
{
	int sumRecords = 0;
	RtreeNode* node = a_rtree.m_memory.loadPage(pageID);
	ramTree[pageID] = node;
	if (node->isLeaf())
	{
		sumRecords = node->m_usedspace;
	}
	else
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			sumRecords += countRecords(a_rtree, node->m_entry[i]->m_id);
		}
	}
	delete node;
	return sumRecords;
}

void aggregateRecords(Rtree& a_rtree)
{
	queue<long int> H;
	RtreeNode* node;
	H.push(a_rtree.m_memory.m_rootPageID);
	long int pageID;

	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = a_rtree.m_memory.loadPage(pageID);

		if (node->isLeaf() == false)
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				node->m_entry[i]->num_records = countRecords(a_rtree, node->m_entry[i]->m_id);
				H.push(node->m_entry[i]->m_id);
			}
		}
		a_rtree.m_memory.writePage(pageID, node);
		ramTree[pageID] = node;
		//delete node;
	}
}

void onionlayer(vector<long int> & skyband, float* PG[], const int k, vector<long int > & klayers, vector<long int> & min_level, int& dim) // Jiahao
{
	const char* tmpfile = "E:\\klevel\\Project\\klevel_index\\qhull\\tmp_qhull.txt";

	vector<long int> records;
	vector<long int> dataset = skyband;

	for (int ki = 0; ki < k; ki++) {
		if (dataset.size() == 0)
		{
			break;
		}
		if (dataset.size() == 1)
		{
			klayers.push_back(dataset[0]);
			min_level.push_back(ki + 1);
			break;
		}
		FILE* fp;
		if ((fp = fopen(tmpfile, "w")) != NULL)
		{
			fprintf(fp, "%d\n", dim);
			fprintf(fp, "%d\n", dataset.size() + dim + 1);

			vector<int> boundary_records; boundary_records.clear();
			for (int i = 0; i < dim; i++) boundary_records.push_back(dataset[0]);

			for (int i = 0; i < dataset.size(); i++) {
				for (int di = 0; di < dim; di++)
				{
					float value = (PG[dataset[i]][di] + PG[dataset[i]][di + dim]) / 2;
					fprintf(fp, "%.4f ", value);
					float d_max_value = (PG[boundary_records[di]][di] + PG[boundary_records[di]][di + dim]) / 2;
					if (value > d_max_value) boundary_records[di] = dataset[i];
				}
				fprintf(fp, "\n");
			}

			for (int i = 0; i < dim; i++) {
				for (int di = 0; di < dim; di++) {
					if (i == di) fprintf(fp, "%.4f ", (PG[boundary_records[i]][di] + PG[boundary_records[i]][di + dim]) / 2);
					else fprintf(fp, "%.4f ", 0.0);
				}
				fprintf(fp, "\n");
			}

			for (int i = 0; i < dim; i++) {
				fprintf(fp, "%.4f ", 0.0);
			}
			fprintf(fp, "\n");

			fclose(fp);
		}

		char command[2048] = "E:\\klevel\\Project\\qhull\\Release\\qconvex.exe";
		strcat(command, " < ");
		strcat(command, tmpfile);
		strcat(command, " Fx > E:\\klevel\\Project\\klevel_index\\qhull\\tmp.ret");
		system(command);

		vector<long int> layerrecords;
		fstream fpdata;
		fpdata.open("E:\\klevel\\Project\\klevel_index\\qhull\\tmp.ret", ios::in);

		int count;
		fpdata >> count;
		int id;
		for (int i = 0; i < count; i++)
		{
			fpdata >> id;
			if (id >= dataset.size()) break;
			layerrecords.push_back(dataset[id]);
		}
		fpdata.close();

		vector<bool> check; check.clear();

		for (int i = 0; i < layerrecords.size(); i++) {
			check.push_back(true);
			bool isDominated = false;
			for (int j = 0; j < layerrecords.size(); j++) {
				if (i != j) {
					int dimcount = 0;
					for (int di = 0; di < dim; di++) {
						//cout << layerrecords[i] << ' ' << layerrecords[j] << endl;
						if (PG[layerrecords[i]][di] < PG[layerrecords[j]][di]) {
							dimcount++;
						}
					}
					if (dimcount == dim) {
						isDominated = true;
						break;
					}
				}
			}
			if (isDominated == false) {
				klayers.push_back(layerrecords[i]);
				min_level.push_back(ki + 1);
			}
			else check[i] = false;
		}

		vector<long int> temp; temp.clear();
		for (int i = 0; i < dataset.size(); i++) {
			bool flag = true;
			for (int j = 0; j < layerrecords.size(); j++) {
				//if (!check[j]) continue;
				if (dataset[i] == layerrecords[j]) {
					flag = false;
					break;
				}
			}
			if (flag) temp.push_back(dataset[i]);
		}
		dataset.clear();
		for (int i = 0; i < temp.size(); i++) dataset.push_back(temp[i]);
	}
}

bool isDominateByFocal(const int dimen, const float pt[], Point& focal)
{
	bool dominated = true;
	for (int i = 0; i < dimen; i++)
	{
		if (focal.m_coor[i] + SIDELEN < pt[i])
		{
			dominated = false;
			break;
		}
	}
	return dominated;
}

void Getkskyband(const int dimen, Rtree& a_rtree, vector<long int>& kskyband, Point& a_pt, float* PG[], const int k)
{
	RtreeNode* node;
	multimap<long int, RtreeNodeEntry*> RecordEntry;
	multimap<long int, RtreeNodeEntry*>::iterator recordIter;
	multimap<float, int> heap;
	int NegPageid;

	float pt[MAXDIMEN];
	float ORIGNIN[MAXDIMEN];
	float mindist;
	for (int i = 0; i < dimen; i++)
		ORIGNIN[i] = 1;

	int pageID;
	float dist_tmp;
	multimap<float, int>::iterator heapIter;

	node = a_rtree.m_memory.loadPage(a_rtree.m_memory.m_rootPageID);
	totalIO++;
	if (node->isLeaf())
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));

			NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
			RtreeNodeEntry * Nentry = node->m_entry[i]->clone();
			RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
		}
	}
	else
	{
		for (int i = 0; i < node->m_usedspace; i++)
		{
			for (int j = 0; j < dimen; j++)
			{
				pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
			}
			mindist = minDist(pt, ORIGNIN, dimen);
			heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
		}
	}
	delete node;

	while (heap.size() != 0)
	{
		heapIter = heap.begin();
		dist_tmp = heapIter->first;
		pageID = heapIter->second;
		heap.erase(heapIter);

		if (pageID > MAXPAGEID)
		{
			recordIter = RecordEntry.find(pageID);
			if (recordIter != RecordEntry.end())
			{
				for (int d = 0; d < dimen; d++)
				{
					pt[d] = recordIter->second->m_hc.getLower()[d] + SIDELEN;
				}
				if (countkDominator(dimen, pt, kskyband, PG) <= k && isDominateByFocal(dimen, pt, a_pt) == false)
				{
					kskyband.push_back(pageID - MAXPAGEID);
				}
			}
			else
			{
				cout << "an error incured in Getkskyband" << endl;
				exit(0);
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
					{
						pt[d] = node->m_entry[i]->m_hc.getLower()[d] + SIDELEN;
					}
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));

						NegPageid = node->m_entry[i]->m_id + MAXPAGEID;
						RtreeNodeEntry * Nentry = node->m_entry[i]->clone();
						RecordEntry.insert(pair<long int, RtreeNodeEntry* >(NegPageid, Nentry));
					}
				}
			}
			else
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int d = 0; d < dimen; d++)
						pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
					if (countkDominator(dimen, pt, kskyband, PG) <= k)
					{
						mindist = minDist(pt, ORIGNIN, dimen);
						heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
					}
				}
			}
			delete node;
		}
	}
}

// for klevel_index

void kskyband(const int dimen, Rtree& a_rtree, vector<int>& kskyband, vector<vector<float>>& OriginData, const int k)
{
    RtreeNode* node;
    multimap<float, int> heap;
    multimap<float, int>::iterator heapIter;

    float pt[MAXDIMEN];
    float ORIGNIN[MAXDIMEN];
    float mindist;
    for (int i = 0; i < dimen; i++)
        ORIGNIN[i] = 1;

    int pageID;
    float dist_tmp;

    heap.insert(PfltINT(INFINITY, a_rtree.m_memory.m_rootPageID));

    while (!heap.empty())
    {
        heapIter = heap.begin();
        dist_tmp = heapIter->first;
        pageID = heapIter->second;
        heap.erase(heapIter);

        if (pageID > MAXPAGEID)
        {
            for (int d = 0; d < dimen; d++)
                pt[d] = OriginData[pageID - MAXPAGEID][d];
            if (countkDominator(dimen, pt, kskyband, OriginData) <= k)
            {
                kskyband.push_back(pageID - MAXPAGEID);
            }
        }
        else
        {
            node = a_rtree.m_memory.loadPage(pageID);
            if (node->isLeaf())
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    for (int d = 0; d < dimen; d++)
                    {
                        pt[d] = node->m_entry[i]->m_hc.getLower()[d] + eps_klevel;
                    }
                    if (countkDominator(dimen, pt, kskyband, OriginData) <= k)
                    {
                        mindist = minDist(pt, ORIGNIN, dimen);
                        heap.insert(PfltINT(mindist, node->m_entry[i]->m_id + MAXPAGEID));
                    }
                }
            }
            else
            {
                for (int i = 0; i < node->m_usedspace; i++)
                {
                    for (int d = 0; d < dimen; d++)
                        pt[d] = node->m_entry[i]->m_hc.getUpper()[d];
                    if (countkDominator(dimen, pt, kskyband, OriginData) <= k)
                    {
                        mindist = minDist(pt, ORIGNIN, dimen);
                        heap.insert(PfltINT(mindist, node->m_entry[i]->m_id));
                    }
                }
            }
        }
    }
}

int countkDominator(const int dimen, const float pt[], vector<int> kskyband, vector<vector<float>>& OriginData)
{
    vector<int>::iterator iter;
    if (kskyband.size() == 0)
        return false;

    int count = 0;
    for (iter = kskyband.begin(); iter != kskyband.end(); iter++)
    {
        long pid = *iter;
        float s[MAXDIMEN];
        bool dominated = true;
        for (int i = 0; i < dimen; i++)
        {
            if (OriginData[pid][i] + eps_klevel < pt[i])
            {
                dominated = false;
                break;
            }
        }
        if (dominated)
            count++;
    }
    return count;
}

void onionlayer(vector<int>& skyband, vector<vector<float>>& OriginData, const int k, vector<int> & klayers, int& dim) {// Jiahao {
    const char* tmpfile = "/home/jiahao/disk/Projects/klevel/config/tmp_qhull.txt";

    vector<int> records;
    vector<int> dataset = skyband;

    for (int ki = 0; ki < k; ki++) {
        if (dataset.size() == 0)
        {
            break;
        }
        if (dataset.size() == 1)
        {
            klayers.push_back(dataset[0]);
            //min_level.push_back(ki + 1);
            break;
        }
        FILE* fp;
        if ((fp = fopen(tmpfile, "w")) != NULL)
        {
            fprintf(fp, "%d\n", dim);
            fprintf(fp, "%d\n", dataset.size() + dim + 1);

            vector<int> boundary_records; boundary_records.clear();
            for (int i = 0; i < dim; i++) boundary_records.push_back(dataset[0]);

            for (int i = 0; i < dataset.size(); i++) {
                for (int di = 0; di < dim; di++)
                {
                    float value = OriginData[dataset[i]][di];
                    fprintf(fp, "%.4f ", value);
                    float d_max_value = OriginData[boundary_records[di]][di];
                    if (value > d_max_value) boundary_records[di] = dataset[i];
                }
                fprintf(fp, "\n");
            }

            for (int i = 0; i < dim; i++) {
                for (int di = 0; di < dim; di++) {
                    if (i == di) fprintf(fp, "%.4f ", OriginData[boundary_records[i]][di]);
                    else fprintf(fp, "%.4f ", 0.0);
                }
                fprintf(fp, "\n");
            }

            for (int i = 0; i < dim; i++) {
                fprintf(fp, "%.4f ", 0.0);
            }
            fprintf(fp, "\n");

            fclose(fp);
        }

        char command[2048] = "/home/jiahao/library/qhull-2020.2/bin/qconvex";
        strcat(command, " < ");
        strcat(command, tmpfile);
        strcat(command, " Fx > /home/jiahao/disk/Projects/klevel/config//tmp.ret");
        system(command);

        vector<long int> layerrecords;
        fstream fpdata;
        fpdata.open("/home/jiahao/disk/Projects/klevel/config//tmp.ret", ios::in);

        int count;
        fpdata >> count;
        int id;
        for (int i = 0; i < count; i++)
        {
            fpdata >> id;
            if (id >= dataset.size()) break;
            layerrecords.push_back(dataset[id]);
        }
        fpdata.close();

        vector<bool> check; check.clear();

        for (int i = 0; i < layerrecords.size(); i++) {
            check.push_back(true);
            bool isDominated = false;
            for (int j = 0; j < layerrecords.size(); j++) {
                if (i != j) {
                    int dimcount = 0;
                    for (int di = 0; di < dim; di++) {
                        //cout << layerrecords[i] << ' ' << layerrecords[j] << endl;
                        if (OriginData[layerrecords[i]][di] < OriginData[layerrecords[j]][di]) {
                            dimcount++;
                        }
                    }
                    if (dimcount == dim) {
                        isDominated = true;
                        break;
                    }
                }
            }
            if (isDominated == false) {
                klayers.push_back(layerrecords[i]);
            }
            else check[i] = false;
        }

        vector<long int> temp; temp.clear();
        for (int i = 0; i < dataset.size(); i++) {
            bool flag = true;
            for (int j = 0; j < layerrecords.size(); j++) {
                //if (!check[j]) continue;
                if (dataset[i] == layerrecords[j]) {
                    flag = false;
                    break;
                }
            }
            if (flag) temp.push_back(dataset[i]);
        }
        dataset.clear();
        for (int i = 0; i < temp.size(); i++) dataset.push_back(temp[i]);
    }
}
