#include "utk.h"

UTK::UTK()
{
	treeSize = 0.0;
}

UTK::~UTK()
{
	layers.clear();
    unordered_map<long int, vector<long int>>().swap(layers);
    daGraph.clear();
    unordered_map<long int, unordered_set<long int>>().swap(daGraph);
	dominationship.clear();
    unordered_map<long int, unordered_set<long int>>().swap(dominationship);
}

void UTK::pivotRegion(vector<float>& R, vector<float>& pivot)
{
	int dim = R.size() / 2;
	pivot.clear(); 
	for (int i = 0; i < dim; i++)
		pivot.push_back((R[2 * i] + R[2 * i + 1]) / 2);
}

bool UTK::isRdominated(const int  dim, vector<float>& R, float focal[], vector<float>& entry, bool& fDe)
{
	// generate HP;
	vector<float> tmpHS;
	float entry_d = entry[dim - 1];
	float focal_d = focal[dim - 1];
	for (int d = 0; d < dim - 1; d++)
	{
		tmpHS.push_back((focal[d] - focal_d) - (entry[d] - entry_d));
	}
	tmpHS.push_back(entry_d - focal_d);


	// verify each vertic
	int posCount = 0;
	int negCount = 0;
	int totVertics = pow(2, dim - 1);
	for (int i = 0; i < totVertics; i++)
	{
		stringstream ss;
		string tmpS = bitset<MAXDIMEN>(i).to_string();
		tmpS = tmpS.substr(tmpS.size() - (dim - 1), dim - 1);

		float sum = 0;
		for (int si = 0; si < tmpS.size(); si++)
		{
			if (tmpS[si] == '0')
				sum += R[si * 2] * tmpHS[si];
			else if (tmpS[si] == '1')
				sum += R[si * 2 + 1] * tmpHS[si];
			else
				cout << "bug here!!!" << endl;
		}

		if (sum > tmpHS[dim - 1])
			negCount++;
		else
			posCount++;
	}

	// dominationship
	if (negCount == totVertics)
	{
		fDe = true;
		return true;
	}
	else if (posCount == totVertics)
	{
		fDe = false;
		return true;
	}
	else
		return false;
}

bool UTK::countRegionDominator(int dimen, float pt[], vector<long int>& rskyband, float* PG[], vector<float>& R, const int k, unordered_set<long int>& dominators)
{
	vector<float> record(dimen,0);
	dominators.clear();
	bool fDe;
	int count = 0;
	/* // should comment this for large k
	if (rskyband.size() < k)
	{
		return true;
	}
	*/
	for (int i = 0; i < rskyband.size(); i++)
	{
		for (int di = 0; di < dimen; di++)
			record[di] = (PG[rskyband[i]][di] + PG[rskyband[i]][di + dimen]) / 2;
		fDe = false;
		if (isRdominated(dimen, R, pt, record, fDe))
		{
			if (fDe == false)
			{
				count++;
				dominators.insert(rskyband[i]);
			}
		}
	}
	if (count < k)
	{
		return true;
	}
	else
	{
		return false;
	}
}

float UTK::orderScore(vector<float>& pivot, float entry[], int dim)
{
	float ret=0;
	float weight = 0;
	for (int i = 0; i < dim - 1; i++)
	{
		weight += pivot[i];
		ret += pivot[i] * entry[i];
	}
	ret += (1 - weight)*entry[dim-1];
	return ret;
}

// okay, dominance graph could be maintain when count dominators
// I did not delete it in case I need it in futher.
void UTK::maintainGraph(const int recordID, float* PG[], vector<float>& region, int dim)
{
	long int tmpID = recordID;
	vector<float> record(dim, 0);
	float* focal = new float[dim];
	bool fDe = true;

	for (int i = 0; i < dim; i++)
		focal[i] = (PG[recordID][i] + PG[recordID][i + dim]) / 2;

	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		for (int i = 0; i < dim; i++)
			record[i] = (PG[iter->first][i] + PG[iter->first][i + dim]) / 2;
		if (isRdominated(dim, region, focal, record, fDe) && fDe==false)
		{
			daGraph[recordID].insert(iter->first);
		}
	}
	unordered_set<long int> newSet;
	daGraph[tmpID] = newSet;
}

void UTK::rskyband(vector<float>& region, const int dimen, Rtree& a_rtree, vector<long int>& rskyband, float* PG[], const int k)
{
	vector<float> pivot;
	pivotRegion(region, pivot);

	RtreeNode* node;
	priority_queue<pair<float, int>> heap;
	int NegPageid;

	float pt[MAXDIMEN];
	float maxscore;
	int pageID;
	float tmpScore;
	unordered_set<long int> dominators;

	heap.push(make_pair(INFINITY, a_rtree.m_memory.m_rootPageID));

	while (!heap.empty())
	{
		tmpScore = heap.top().first;
		pageID = heap.top().second;
		heap.pop();

		if (pageID > MAXPAGEID)
		{
			for (int j = 0; j < dimen; j++)
				pt[j] = (PG[pageID - MAXPAGEID][j] + PG[pageID - MAXPAGEID][j + dimen]) / 2;
			if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
			{
				rskyband.push_back(pageID - MAXPAGEID);
				daGraph[pageID - MAXPAGEID] = dominators;
			}
		}
		else
		{
			node = a_rtree.m_memory.loadPage(pageID);
			//node = ramTree[pageID];
			if (node->isLeaf())
			{
				for (int i = 0; i < node->m_usedspace; i++)
				{
					for (int j = 0; j < dimen; j++)
						pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
					
					if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
					{
						maxscore = orderScore(pivot, pt, dimen);
						heap.push(make_pair(maxscore, node->m_entry[i]->m_id + MAXPAGEID));
					}
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
					if (countRegionDominator(dimen, pt, rskyband, PG, region, k, dominators))
					{
						maxscore = orderScore(pivot, pt, dimen);
						heap.push(make_pair(maxscore, node->m_entry[i]->m_id));
					}
				}
			}
		}
	}
}

void UTK::rskybandLayer()
{
	unordered_set<long int> initSet;

	for (auto iter = layers.begin(); iter != layers.end(); iter++)
	{
		iter->second.clear();
	}

	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		layers[iter->second.size()].push_back(iter->first);
		if (dominationship.find(iter->first) == dominationship.end())
		{
			dominationship[iter->first] = initSet;
		}
		for (auto second = iter->second.begin(); second != iter->second.end(); second++)
		{
			dominationship[*second].insert(iter->first);
		}
	}
}

bool UTK::drill(Point& focal, vector<cell*>& leaves, vector<cell*>& finalResult, const int k, float* PG[])
{
	vector<float> pivot;
	int dim = focal.m_dimen;
	multimap<float, long int> anchors;
	float ret = 0;
	float weight = 0;

	for (int i = 0; i < leaves.size(); i++)
	{
		findPivot(leaves[i], pivot);
		
		for (int i = 0; i < layers.size(); i++)
		{
			for (auto ai = layers[i].begin(); ai != layers[i].end(); ai++)
			{
				ret = 0;
				weight = 0;
				for (int di = 0; di < dim - 1; di++)
				{
					weight += pivot[di];
					ret += pivot[di] * (PG[*ai][di] + PG[*ai][di + dim]) / 2;
				}
				ret += (1 - weight)*(PG[*ai][dim - 1] + PG[*ai][dim - 1 + dim]) / 2;
				
				if (anchors.size() < k)
				{
					anchors.insert(make_pair(ret, *ai));
				}
				else
				{
					if (ret>anchors.begin()->first)
					{
						anchors.erase(anchors.begin());
						anchors.insert(make_pair(ret, *ai));
					}
				}
			}
		}
		
		ret = 0;
		weight = 0;
		for (int di = 0; di < dim - 1; di++)
		{
			weight += pivot[di];
			ret += pivot[di] * focal.m_coor[di];
		}
		ret += (1 - weight)*focal.m_coor[dim-1];

		if (ret>anchors.begin()->first)
		{
			finalResult.push_back(leaves[i]);
			return true;
		}
	}
	return false;
}

//void UTK::computeHS(Point& pt, long int b, float* PG[], int dim)
void UTK::computeHS(Point& pt, long long b, float* PG[], int dim)
{
	if (RecordIDtoHalfPlaneID.find(b) == RecordIDtoHalfPlaneID.end())
	{
		vector<float> entry(dim, 0);
		for (int i = 0; i < dim; i++)
		{
			entry[i] = (PG[b%objCnt][i] + PG[b%objCnt][i + dim]) / 2;
		}

		vector<float> tmpHS;
		float entry_d = entry[dim - 1];
		float focal_d = pt.m_coor[dim - 1];
		for (int d = 0; d < dim - 1; d++)
		{
			tmpHS.push_back((entry[d] - entry_d) - (pt.m_coor[d] - focal_d));
		}
		tmpHS.push_back(focal_d - entry_d);
		tmpHS.push_back(b);            
		HalfSpaces.push_back(tmpHS);   
		RecordIDtoHalfPlaneID.insert(PintInt(b, HalfSpaces.size()));
	}
}

void UTK::rsa(vector<float>& region, unordered_map<int, cell*>& utkRet, const int k, const int dim, float* PG[], Rtree& rtree)
{
	vector<pair<long int, unordered_set<long int>>> orderedgNode;
	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		orderedgNode.push_back(make_pair(iter->first, iter->second));
	}
	sort(orderedgNode.begin(), orderedgNode.end(), daGraphCompare());
	rskybandLayer();

	for (int nodei = 0; nodei < orderedgNode.size() && !daGraph.empty(); nodei++)
	{

		cout << "============================================" << endl;
		cout << "No. " << nodei + 1 << ", Record ID: " << orderedgNode[nodei].first << endl;

		if (daGraph.find(orderedgNode[nodei].first)==daGraph.end())
		{
			continue;
		}

		vector<cell*> finalResult;
		vector<cell*> leaves;
		unordered_set<long int> removeSL;
		unordered_set<long int> singular;
		vector<long int> processRecords;

		HalfSpaces.clear();
		RecordIDtoHalfPlaneID.clear();

		// init focal record
		pair<long int, unordered_set<long int>> iter = orderedgNode[nodei];
		Point pt;
		pt.m_dimen = dim;
		for (int j = 0; j < dim; j++)
		{
			pt.m_coor[j] = (PG[iter.first][j] + PG[iter.first][j + dim]) / 2.0;
		}
		int updatek = k - iter.second.size() - 1;

		initHS(dim, region);
		cellTree* sol = new cellTree(2 * (dim - 1));

		for (int ki = 0; ki < k; ki++)
		{	
			if (ki == 0)
			{
				for (int i = 0; i < layers[ki].size(); i++)
				{
					if (iter.first != layers[ki][i] && iter.second.find(layers[ki][i]) == iter.second.end())
					{
						processRecords.push_back(layers[0][i]);
						computeHS(pt, layers[0][i], PG, dim);
					}
				}
				cout << "(1) Skylines(" << k << "): " << processRecords.size() << endl;
				sol->insert(processRecords, updatek, rtree, pt, finalResult);
			}
			else
			{
				sol->opt_insert(processRecords, updatek, rtree, pt, finalResult);
			}
			sol->collectLeaf(leaves, updatek);
			if (treeSize < sol->treeSize)
			{
				treeSize = sol->treeSize;
			}

			if (leaves.empty())
			{
				cout << "(3) Singular, #Ret: " << finalResult.size() << endl;
				break;
			}

			sol->markSingular(leaves, removeSL, updatek, finalResult, singular);
			if (finalResult.size() != 0)
			{
				cout << "(3) Singular, #Ret: " << finalResult.size() << endl;
				long int utkID = iter.first;
				utkRet[utkID] = finalResult[0]; 
				for (auto idIter = daGraph[utkID].begin(); idIter != daGraph[utkID].end(); idIter++)
				{
					if (daGraph.find(*idIter) != daGraph.end())
					{
						cout << "remove: " << *idIter << endl;
						utkRet[*idIter] = finalResult[0]; // need update
						daGraph.erase(*idIter);
					}
				}
				daGraph.erase(utkID);
				break;
			}
			else if (ki == k)
			{
				cout << "please check why this happens" << endl;
				break;
			}
			else // obtain next batch
			{
				//if (drill(pt, leaves, finalResult, k, PG))
				//{
				//	cout << "(3) Singular, #Ret: " << finalResult.size() << endl;
				//	long int utkID = iter.first;
				//	utkRet[utkID] = finalResult[0];
				//	for (auto idIter = daGraph[utkID].begin(); idIter != daGraph[utkID].end(); idIter++)
				//	{
				//		if (daGraph.find(*idIter) != daGraph.end())
				//		{
				//			cout << "remove: " << *idIter << endl;
				//			utkRet[*idIter] = finalResult[0]; // need update
				//			daGraph.erase(*idIter);
				//		}
				//	}
				//	daGraph.erase(utkID);
				//	break;
				//}
				//else
				//{
					processRecords.clear();
					for (int i = 0; i < layers[ki + 1].size(); i++)
					{
						bool isNext = false;
						for (auto newIter = removeSL.begin(); newIter != removeSL.end(); newIter++)
						{
							if (daGraph[layers[ki + 1][i]].find(*newIter) != daGraph[layers[ki + 1][i]].end())
							{
								isNext = true;
								break;
							}
						}
						if (isNext)
						{
							processRecords.push_back(layers[ki + 1][i]);
							computeHS(pt, layers[ki + 1][i], PG, dim);
						}
					}
				//}
			}
		}
		sol->~cellTree();
	}
}

void UTK::findAnchor(vector<float>& pivot, unordered_set<long int>& ignoreset, long int& focal, unordered_set<long int>& focaled, const int dim, float* PG[], const int k)
{
	multimap<float, long int> anchors;
	for (int i = 0; i < layers.size(); i++)
	{
		for (auto ai = layers[i].begin(); ai != layers[i].end(); ai++)
		{
			if (ignoreset.find(*ai) == ignoreset.end() && focaled.find(*ai) == focaled.end()) // ignore these dominators
			{
				float ret = 0;
				float weight = 0;

				for (int di = 0; di < dim - 1; di++)
				{
					weight += pivot[di];
					ret += pivot[di] * (PG[*ai][di] + PG[*ai][di + dim]) / 2;
				}
				ret += (1 - weight)*(PG[*ai][dim - 1] + PG[*ai][dim - 1 + dim]) / 2;

				if (anchors.size() < k)
				{
					anchors.insert(make_pair(ret, *ai));
				}
				else
				{
					if (ret>anchors.begin()->first)
					{
						anchors.erase(anchors.begin());
						anchors.insert(make_pair(ret, *ai));
					}
				}
			}
		}
	}
	focal = anchors.begin()->second;
	
	//float anchorScore = 0.0;
	//for (int i = 0; i < layers.size(); i++)
	//{
	//	for (auto ai = layers[i].begin(); ai != layers[i].end(); ai++)
	//	{
	//		if (ignoreset.find(*ai) == ignoreset.end() && focaled.find(*ai) == focaled.end()) // ignore these dominators
	//		{
	//			float ret = 0;
	//			float weight = 0;

	//			for (int di = 0; di < dim - 1; di++)
	//			{
	//				weight += pivot[di];
	//				ret += pivot[di] * (PG[*ai][di] + PG[*ai][di + dim]) / 2;
	//			}
	//			ret += (1 - weight)*(PG[*ai][dim - 1] + PG[*ai][dim - 1 + dim]) / 2;

	//			if (ret > anchorScore)
	//			{
	//				anchorScore = ret;
	//				focal = *ai;
	//			}
	//		}
	//	}
	//}
}

void UTK::findPivot(cell* subreg, vector<float>& pivot)
{
	if (subreg->isPruned == false) // update ignoreset iff subreg is less-than region
	{
		// identify ignore set
		for (auto iter = subreg->IntersectHP.begin(); iter != subreg->IntersectHP.end(); iter++)
		{
			if (iter->second == true)
			{
				if (subreg->ignoreset.find(iter->first%objCnt) == subreg->ignoreset.end())
				{
					subreg->ignoreset.insert(iter->first%objCnt);
				}
			}
		}
		for (auto iter = subreg->BelowHP.begin(); iter != subreg->BelowHP.end(); iter++)
		{
			if (subreg->ignoreset.find(*iter%objCnt) == subreg->ignoreset.end())
			{
				subreg->ignoreset.insert(*iter%objCnt);
			}
		}

		//// this is unnecessory for whole program, i keep it as a bug checker!
		//if (subreg->ignoreset.size() != subreg->rank)
		//{
		//	cout << "cell error!, please check" << endl;
		//}
	}

	// identify bound of subregion
	int dimen = HalfSpaces[0].size() - 2;
	double row[MAXDIMEN];
	// lp model;
	cellTree* canTree = new cellTree(subreg);
	lprec *lp = make_lp(0, dimen);
	canTree->lpModel(lp, dimen, subreg->IntersectHP);
	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = -1;
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int ret = solve(lp);
	get_variables(lp, row);

	for (int dii = 0; dii < dimen; dii++)
	{
		pivot.push_back(row[dii]);
	}

	canTree->~cellTree();
	delete_lp(lp);
}

void UTK::insertingRecords(cellTree* subreg, vector<long int>& processRecords, const long int focal, float* PG[], Point& pt, long int& epoch, const int dim)
{
	long int recordID;
	processRecords.clear();
	unordered_set<long int> invalidRecords;
	invalidRecords.insert(focal);
	
	for (auto iter = subreg->root->focaled.begin(); iter != subreg->root->focaled.end(); iter++) // ignore focaled records
	{
		invalidRecords.insert(*iter);
	}
	for (auto iiter = dominationship[focal].begin(); iiter != dominationship[focal].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
	{
		invalidRecords.insert(*iiter);
	}
	// There is another opimization, but I need think how to express it correctly.
	for (auto iter = subreg->root->ignoreset.begin(); iter != subreg->root->ignoreset.end(); iter++)
	{
		invalidRecords.insert(*iter);
		for (auto iiter = daGraph[*iter].begin(); iiter != daGraph[*iter].end(); iiter++)
		{
			invalidRecords.insert(*iiter);
		}
	}

	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		if (invalidRecords.find(iter->first) == invalidRecords.end())
		{
			recordID = epoch*objCnt + iter->first;
			processRecords.push_back(recordID);
			computeHS(pt, recordID, PG, dim);
		}
	}
}

void UTK::InsertRecords(cellTree* subreg, vector<long int>& processRecords, const long int focal, float* PG[], Point& pt, long int& epoch, const int dim)
{
	//long int recordID;
	long long recordID;
	processRecords.clear();
	unordered_set<long int> invalidRecords;
	invalidRecords.insert(focal);
	if (subreg->root->isPruned == true)
	{
		subreg->root->rank = subreg->root->ignoreset.size();
		subreg->root->isPruned = false; 
		for (auto iter = subreg->root->focaled.begin(); iter != subreg->root->focaled.end(); iter++) // ignore focaled records
		{
			invalidRecords.insert(*iter);
			for (auto iiter = dominationship[*iter].begin(); iiter != dominationship[*iter].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
			{
				invalidRecords.insert(*iiter);
			}
		}
		for (auto iter = subreg->root->ignoreset.begin(); iter != subreg->root->ignoreset.end(); iter++) // skip ignored records
		{
			invalidRecords.insert(*iter);
		}
		for (auto iiter = dominationship[focal].begin(); iiter != dominationship[focal].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
		{
			invalidRecords.insert(*iiter);
		}
	}
	else
	{
		for (auto iter = subreg->root->focaled.begin(); iter != subreg->root->focaled.end(); iter++) // ignore focaled records
		{
			invalidRecords.insert(*iter);
		}
		for (auto iiter = dominationship[focal].begin(); iiter != dominationship[focal].end(); iiter++) // Lemma 1: ignore records which are r-dominated by focal record
		{
			invalidRecords.insert(*iiter);
		}
		// There is another opimization, but I need think how to express it correctly.
		for (auto iter = subreg->root->ignoreset.begin(); iter != subreg->root->ignoreset.end(); iter++)
		{
			invalidRecords.insert(*iter);
		}
	}
	for (auto iter = daGraph.begin(); iter != daGraph.end(); iter++)
	{
		if (invalidRecords.find(iter->first) == invalidRecords.end())
		{
			recordID =(long long) epoch*objCnt + iter->first;
			processRecords.push_back(recordID);
			if ((recordID < 0) || (iter->first < 0) || (iter->first >= objCnt)) {
				cout << recordID << ' ' << epoch << ' ' << iter->first << endl;
				system("pause");
			}
			computeHS(pt, recordID, PG, dim);
		}
	}
}

void UTK::jointArrangement(long int& epoch, cellTree* sol, const long int focal, vector<cell*>& exactutk, unordered_set<long int>& focaled,
	const int k, const int dim, unordered_set<long int>& ignoreset, float* PG[], Rtree& rtree, queue<cell*>& subregions)
{
	vector<cell*> finalResult;
	vector<cell*> leaves;
	unordered_set<long int> removeSL;
	unordered_set<long int> singular;
	vector<long int> processRecords;

	// init focal record
	Point pt;
	pt.m_dimen = dim;
	for (int j = 0; j < dim; j++)
	{
		pt.m_coor[j] = (PG[focal][j] + PG[focal][j + dim]) / 2.0;
	}
	InsertRecords(sol, processRecords, focal, PG, pt, epoch, dim);

	if (processRecords.size() == 0)
	{
		exactutk.push_back(sol->root);
	}
	else
	{
		if (processRecords.size() > 0)
		{
			//cout << "(1) Insert Records: " << processRecords.size() << endl;
			sol->insert(processRecords, k, rtree, pt, finalResult);
		}
		sol->collectLeaf(leaves, INFINITY);
		if (treeSize < sol->treeSize)
		{
			treeSize = sol->treeSize;
		}
		for (int ci = 0; ci < leaves.size(); ci++)
		{
			if (leaves[ci]->rank == k - 1)
			{
				exactutk.push_back(leaves[ci]);
			}
			else if (leaves[ci]->rank < k - 1)
			{
				leaves[ci]->ignoreset.insert(focal);
				leaves[ci]->rank++;
				subregions.push(leaves[ci]);
			}
			else
			{
				leaves[ci]->focaled.insert(focal);
				leaves[ci]->BelowHP.clear();
				leaves[ci]->AboveHP.clear();
				subregions.push(leaves[ci]);
			}
		}
	}
	sol->~cellTree();
}

// unfortunately, this function does not work well for JAA
void UTK::initHalfspace(cell* subreg)
{
	vector<vector<float>> newHalfSpaces; 
	//unordered_map<long int, long int> newRecordIDtoHalfPlaneID; 
	unordered_map<long long, long int> newRecordIDtoHalfPlaneID;
	pair<long int, bool> tmpPair;

	unordered_map<long int, bool> tmpIntersectedHP = subreg->IntersectHP;

	int dim = HalfSpaces[0].size();
	int i = 0;
	for (auto iter = tmpIntersectedHP.begin(); iter != tmpIntersectedHP.end(); iter++)
	{
		if (iter->first < 0)
			newHalfSpaces.push_back(HalfSpaces[RecordIDtoHalfPlaneID[iter->first] - 1]);
		else
		{
			HalfSpaces[RecordIDtoHalfPlaneID[iter->first] - 1][dim - 1] = -(i + 1);
			newHalfSpaces.push_back(HalfSpaces[RecordIDtoHalfPlaneID[iter->first] - 1]);
			subreg->IntersectHP.erase(iter->first);
			subreg->IntersectHP[-(i + 1)] = iter->second;
		}
		newRecordIDtoHalfPlaneID.insert(make_pair(-(i + 1), newHalfSpaces.size()));

		i++;
	}

	HalfSpaces.clear();
	RecordIDtoHalfPlaneID.clear();
	HalfSpaces = newHalfSpaces;
	RecordIDtoHalfPlaneID = newRecordIDtoHalfPlaneID;
}

void UTK::jaa(vector<float>& region, vector<cell*>& exactutk, int& k, const int dim, float* PG[], Rtree& a_rtree)
{
	rskybandLayer(); // update layer

	// obtain pivot
	vector<float> pivot;
	pivotRegion(region, pivot);
	// obtain anchor
	unordered_set<long int> ignoreset;
	unordered_set<long int> focaled;
	long int focal = 0;
	findAnchor(pivot, ignoreset, focal, focaled, dim, PG, k);
	// process focal
	queue<cell*> subregions;

	HalfSpaces.clear();
	RecordIDtoHalfPlaneID.clear();
	initHS(dim, region);
	cellTree* candidate = new cellTree(2 * (dim - 1));
	long int epoch = 0;

	jointArrangement(epoch, candidate, focal, exactutk, focaled, k, dim, ignoreset, PG, a_rtree, subregions);
	cell* subreg = NULL;

	
	while (!subregions.empty())
	{
		//cout << epoch << ", result size: " << exactutk.size() << "; remain subregions: " << subregions.size() << endl;
		epoch++;
		subreg = subregions.front();
		subregions.pop();

		pivot.clear(); 
		findPivot(subreg, pivot);
		if (subreg->ignoreset.size() >=k )
		{
			exactutk.push_back(subreg);
		}
		else
		{
			findAnchor(pivot, subreg->ignoreset, focal, subreg->focaled, dim, PG, k - subreg->ignoreset.size());
			//initHalfspace(subreg);
			candidate = new cellTree(subreg);
			jointArrangement(epoch, candidate, focal, exactutk, subreg->focaled, k, dim, subreg->ignoreset, PG, a_rtree, subregions);
		}
	}
}