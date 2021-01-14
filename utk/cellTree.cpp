#include "cellTree.h"

extern vector<vector<float>> HalfSpaces;
//extern unordered_map<long int, long int> RecordIDtoHalfPlaneID;
extern unordered_map<long long, long int> RecordIDtoHalfPlaneID;
extern double totalSpaceCost;
extern unordered_map<long int, RtreeNode*> ramTree;

cellTree::cellTree()
{
	root = new cell();
	treeSize = 0;
}

cellTree::cellTree(cell* node)
{
	root = new cell();
	node->copyleaf(root);
	treeSize = 0;
}

cellTree::cellTree(int size)
{
	root = new cell();
	for (int i = 0; i < size; i++)
	{
		root->IntersectHP[-(i+1)] = false;
	}
	treeSize = 0;
}

cellTree::~cellTree()
{
	for (auto iter = dagNode.begin(); iter != dagNode.end(); iter++)
	{
		iter->second->rDominator.clear();
		unordered_set<long int>().swap(iter->second->rDominator);
		delete iter->second;
	}
	dagNode.clear();
	unordered_map<long int, gNode*>().swap(dagNode);
	releaseCell(root);
}

void cellTree::releaseCell(cell* node)
{
	cell* tmp;
	for (tmp = node; tmp != nullptr && tmp->left != nullptr; tmp = tmp->left);

	while (node != nullptr)
	{
		for (tmp = node; tmp != nullptr && tmp->left != nullptr; tmp = tmp->left);
		cell* old = node;
		node = node->left;
		old->IntersectHP.clear();
		unordered_map<long int, bool>().swap(old->IntersectHP);
		old->BelowHP.clear();
		unordered_set<long int>().swap(old->BelowHP);
		old->AboveHP.clear();
		unordered_set<long int>().swap(old->AboveHP);
		delete old;
	}
}

void cellTree::lpModel(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs)
{
	double row[MAXDIMEN];
	if (lp == NULL)
	{
		fprintf(stderr, "Unable to create new LP model\n");
		exit(0);
	}

	// constraints in each dimension
	for (int di = 0; di < dimen; di++)
	{
		row[0] = 0;
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			if (dii - 1 == di)
			{
				row[dii] = 1;
			}
			else
			{
				row[dii] = 0;
			}
		}

		add_constraint(lp, row, GE, 0);
		add_constraint(lp, row, LE, 1);
	}

	// in reduced space, sum_{q_i} should less than 1
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = 1;
	}
	add_constraint(lp, row, GE, 0);
	add_constraint(lp, row, LE, 1);

	// constraints in intersected hyperplanes
	for (unordered_map<long int, bool>::iterator iter = touchhs.begin(); iter != touchhs.end(); iter++)
	{
		addHP(lp, iter->first, iter->second);
	}

	// set scale 
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

void cellTree::lpModelwithSIDELEN(lprec* lp, int& dimen, unordered_map<long int, bool>& touchhs)
{
	double row[MAXDIMEN];
	if (lp == NULL)
	{
		fprintf(stderr, "Unable to create new LP model\n");
		exit(0);
	}

	// constraints in each dimension
	for (int di = 0; di < dimen; di++)
	{
		row[0] = 0;
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			if (dii - 1 == di)
			{
				row[dii] = 1;
			}
			else
			{
				row[dii] = 0;
			}
		}
		add_constraint(lp, row, GE, SIDELEN);
		add_constraint(lp, row, LE, 1-SIDELEN);
	}

	// in reduced space, sum_{q_i} should less than 1
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = 1;
	}
	add_constraint(lp, row, GE, SIDELEN);
	add_constraint(lp, row, LE, 1-SIDELEN);

	// constraints in intersected hyperplanes
	for (unordered_map<long int, bool>::iterator iter = touchhs.begin(); iter != touchhs.end(); iter++)
	{
		addHPwithSIDELEN(lp, iter->first, iter->second);
	}

	// set scale 
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

void cellTree::addHP(lprec* model, long int hpid, bool sideindicator)
{
	int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	int dimen = HalfSpaces[0].size() - 2;
	double row[MAXDIMEN];

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = HalfSpaces[hsID][dii - 1];
	}
	if (sideindicator == false)
	{
		add_constraint(model, row, LE, HalfSpaces[hsID][dimen]);
	}
	else if (sideindicator == true)
	{
		add_constraint(model, row, GE, HalfSpaces[hsID][dimen]);
	}
	else
	{
		std::cout << "Unable to detect half plane direction!!!" << endl;
	}
}

void cellTree::addHPwithSIDELEN(lprec* model, long int hpid, bool sideindicator)
{
	int hsID = RecordIDtoHalfPlaneID[hpid] - 1;
	int dimen = HalfSpaces[0].size() - 2;
	double row[MAXDIMEN];

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = HalfSpaces[hsID][dii - 1];
	}
	if (sideindicator == false)
	{
		add_constraint(model, row, LE, HalfSpaces[hsID][dimen]-SIDELEN);
	}
	else if (sideindicator == true)
	{
		add_constraint(model, row, GE, HalfSpaces[hsID][dimen]+SIDELEN);
	}
	else
	{
		std::cout << "Unable to detect half plane direction!!!" << endl;
	}
}

bool cellTree::isFeasible(unordered_map<long int, bool>& touchhs, long int hpid, bool sideindicator)
{
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;
	lprec *lp = make_lp(0, dimen);
	
	lpModel(lp, dimen, touchhs);
	addHP(lp, hpid, sideindicator);

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

	// for reduced space
	delete_lp(lp);
	if (ret == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

void cellTree::updateRank(cell* node, const int mink)
{
	node->rank++;
	if (node->rank >= mink)
	{
		node->isPruned = true;
		if (node->left != NULL)
		{
			releaseCell(node->left);
			node->left = NULL;
		}
		if (node->right != NULL)
		{
			releaseCell(node->right);
			node->right = NULL;
		}

	}
	else
	{
		if (node->left != NULL)
			updateRank(node->left, mink);
		if (node->right != NULL)
			updateRank(node->right, mink);
	}
}

void cellTree::inserthp(long int & hpid, const int mink, cell* node, unordered_map<long int, bool>& touchhs)
{
	if (node->isPruned == false)
	{
		for (auto iter = node->IntersectHP.begin(); iter != node->IntersectHP.end(); iter++)
		{
			touchhs[iter->first] = iter->second;
		}
		if (isFeasible(touchhs, hpid, false) == false)
		{
			node->BelowHP.insert(hpid);
			updateRank(node, mink);
		}
		else if (isFeasible(touchhs, hpid, true) == false)
		{
			node->AboveHP.insert(hpid);
		}
		else
		{
			if (node->left == NULL&&node->right == NULL)
			{
				if (node->rank > mink)
				{
					node->isPruned = true;
				}
				else
				{
					node->left = new cell(node);
					node->left->IntersectHP[hpid] = false;

					node->right = new cell(node);
					node->right->IntersectHP[hpid] = true;
					node->right->rank++;

					if (node->right->rank >= mink)
					{
						node->right->isPruned = true;
					}
				}
			}
			else if (node->left != NULL&&node->right != NULL)
			{
				if (node->left->isPruned == false)
				{
					unordered_map<long int, bool> leftths = touchhs;
					inserthp(hpid, mink, node->left, leftths);
				}
				if (node->right->isPruned == false)
				{
					unordered_map<long int, bool> rightths = touchhs;
					inserthp(hpid, mink, node->right, rightths);
				}
				else if (node->left->isPruned == true && node->right->isPruned == true)
				{
					node->isPruned = true;
				}
			}
		}
	}
}

void cellTree::insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult)
{
	for (int pos = 0; pos < hps.size(); pos++)
	{
		if (root->isPruned == false)
		{
			unordered_map<long int, bool> hs;
			inserthp(hps[pos], mink, root, hs);
		}
	}
	updateCellTree(root);
}

void cellTree::insert_kspr(vector<long int>& hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult)
{
	vector<cell*> leaves;
	if (root->left == NULL && root->right == NULL)
	{
		root->left = new cell(root);
		root->left->IntersectHP[hps[0]] = false;
		totalSpaceCost += 38;

		root->right = new cell(root);
		root->right->IntersectHP[hps[0]] = true;
		root->right->rank++;
		totalSpaceCost += 38;
	}
	for (int pos = 1; pos < hps.size(); pos++)
	{
		if (pos % 20 == 0)
		{
			cout << pos << endl;
		}
		if (root->isPruned == false)
		{
			unordered_map<long int, bool> leftths;
			inserthp(hps[pos], mink, root->left, leftths);
			unordered_map<long int, bool> rightths;
			inserthp(hps[pos], mink, root->right, rightths);
		}
	}
	updateCellTree(root);
}

void cellTree::collectLeaf_kspr(vector<cell*>& leaves, const int& mink)
{
	leaves.clear();
	if (root->isPruned == false)
	{
		if (root->left != NULL && root->left->isPruned == false)
		{
			cell* leaf = new cell();
			dfsTraversal(root->left, leaf, leaves);
		}
		if (root->right != NULL && root->right->isPruned == false)
		{
			cell* leaf = new cell();
			dfsTraversal(root->right, leaf, leaves);
		}
	}
	sort(leaves.begin(), leaves.end(), cellCompare());
}


void cellTree::collectLeaf(vector<cell*>& leaves, const int& mink)
{
	leaves.clear();
	if (root->isPruned == false)
	{
		treeSize += root->cellSize();
		//if (root->left != NULL && root->left->isPruned == false)
		if (root->left != NULL)
		{
			cell* leaf = new cell();
			root->copyleaf(leaf);
			dfsTraversal(root->left, leaf, leaves);
		}
		//if (root->right != NULL && root->right->isPruned == false)
		if (root->right != NULL)
		{
			cell* leaf = new cell();
			root->copyleaf(leaf);
			dfsTraversal(root->right, leaf, leaves);
		}
		if (root->left == NULL&& root->right == NULL)
		{
			cell* leaf = new cell();
			root->copyleaf(leaf);
			leaves.push_back(leaf);
		}
	}
	sort(leaves.begin(), leaves.end(), cellCompare());
}

void cellTree::collectLeaf(vector<pair<cell*, long int>>& leaves, const int& mink)
{
	leaves.clear();
	if (root->isPruned == false)
	{
		if (root->left != NULL && root->left->isPruned == false)
		{
			cell* leaf = new cell();
			dfsTraversal(root->left, leaf, leaves);
		}
		if (root->right != NULL && root->right->isPruned == false)
		{
			cell* leaf = new cell();
			dfsTraversal(root->right, leaf, leaves);
		}
	}
	sort(leaves.begin(), leaves.end(), cellpairCompare());
}

void cellTree::dfsTraversal(cell* node, cell* leaf, vector<cell*>& leaves)
{
	node->copyleaf(leaf);
	treeSize += node->cellSize();
	if (node->left == NULL && node->right == NULL)
	{
		leaves.push_back(leaf);
	}
	else
	{
//		if (node->left->isPruned == false)
		if (node->left != NULL)
		{
			cell* left = new cell();
			leaf->copyleaf(left);
			dfsTraversal(node->left, left, leaves);
		}
		//if (node->right->isPruned == false)
		if (node->right != NULL)
		{
			dfsTraversal(node->right, leaf, leaves);
		}
		else
		{
			leaf->release();
			delete leaf;
		}
	}
}

void cellTree::dfsTraversal(cell* node, cell* leaf, vector<pair<cell*, long int>>& leaves)
{
	node->copyleaf(leaf);

	if (node->left == NULL && node->right == NULL)
	{
		leaves.push_back(make_pair(leaf,leaves.size()));
		cellID[leaves.size() - 1] = node;
	}
	else
	{
		if (node->left->isPruned == false)
		{
			cell* left = new cell();
			leaf->copyleaf(left);
			dfsTraversal(node->left, left, leaves);
		}
		if (node->right->isPruned == false)
		{
			dfsTraversal(node->right, leaf, leaves);
		}
		else
		{
			leaf->release();
			delete leaf;
		}
	}
}

void cellTree::opt_insert(vector<long int> &hps, int mink, Rtree& rt, Point& a_pt, vector<cell*>& finalResult)
{
	for (int pos = 0; pos < hps.size(); pos++)
	{
		if (root->isPruned == false)
		{
			unordered_map<long int, bool> leftths;
			inserthp(hps[pos], mink, root->left, leftths);
			unordered_map<long int, bool> rightths;
			inserthp(hps[pos], mink, root->right, rightths);
		}
	}
	updateCellTree(root);
}

void cellTree::opt_inserthp(long int & hpid, const int mink, cell* node, cell* all)
{
	if (node->isPruned == false)
	{
		node->appendto(all);
		int indicator = isDominatorInserted(dagNode[hpid]->rDominator, all);
		if (indicator == -1) // found negative halfspace from hpid's dominators
		{
			node->AboveHP.insert(hpid);
			totalSpaceCost += 4;
			releaseCell(all);
			delete all;
		}
		else if (indicator == 1) // found positive halfspace from hpid's dominators
		{
			inserthp(hpid, mink, node, all->IntersectHP);
		}
		else if (indicator == 0) // go to its children
		{
			if (node->left->isPruned == false)
			{
				cell* left = new cell();
				all->appendto(left);
				opt_inserthp(hpid, mink, node->left, left);
			}
			if (node->right->isPruned == false)
			{
				opt_inserthp(hpid, mink, node->left, all);
			}
			else 
			{
				node->isPruned = true;
				all->release();
				delete all;
			}
		}
		else
		{
			cout << "there is a bug in opt_inserthp" << endl;
		}
	}
}

void cellTree::maintainDAG(unordered_set<long int>& skylines, unordered_set<long int>& removeSL, float* PG[], const int dimen)
{
	if (removeSL.size() == 0)
	{
		for (auto iter = skylines.begin(); iter != skylines.end(); iter++)
		{
			gNode* tmp = new gNode((long int&)*iter);
			dagNode[*iter] = tmp;
		}
	}
	else
	{
		float oldR[MAXDIMEN], newR[MAXDIMEN];
		for (auto iter = skylines.begin(); iter != skylines.end(); iter++)
		{
			if (dagNode.find(*iter) == dagNode.end())
			{
				gNode* tmp = new gNode((long int&)*iter);
				dagNode[*iter] = tmp;
				totalSpaceCost += 5;
				for (auto riter = removeSL.begin(); riter != removeSL.end(); riter++)
				{
					for (int i = 0; i < dimen; i++)
					{
						newR[i] = (PG[*iter][i] + PG[*iter][i + dimen]) / 2;
						oldR[i] = (PG[*riter][i] + PG[*riter][i + dimen]) / 2;
						if (dominateRecords(oldR, newR, dimen))
						{
							tmp->rDominator.insert(*riter);
						}
					}
				}
			}
			else
			{
				// this node has been inserted in dag
			}
		}
	}
}

int cellTree::isDominatorInserted(unordered_set<long int>& rdominators, cell* all)
{
	for (auto iter = rdominators.begin(); iter != rdominators.end(); iter++)
	{
		if (all->AboveHP.find(*iter) != all->AboveHP.end())
		{
			return -1;
		}
		if (all->BelowHP.find(*iter) != all->BelowHP.end())
		{
			return 1;
		}
		if (all->IntersectHP.find(*iter) != all->IntersectHP.end())
		{
			return all->IntersectHP[*iter] == true ? 1 : -1;
		}
	}
	return 0;
}

void cellTree::markSingular_kspr(vector<cell*>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular)
{
	long int recordID;
	a_skylines.clear();
	bool flag = false;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i]->isPruned == false)
		{
			if (leaves[i]->rank < mink)
			{
				for (auto iter = leaves[i]->BelowHP.begin(); iter != leaves[i]->BelowHP.end(); iter++)
				{
					recordID = *iter;
					if (singular.find(recordID) == singular.end())
					{
						flag = true;
						if (a_skylines.find(recordID) == a_skylines.end())
						{
							a_skylines.insert(recordID);
						}
					}
				}

				for (auto iter = leaves[i]->IntersectHP.begin(); iter != leaves[i]->IntersectHP.end(); iter++)
				{
					if (iter->second == true)
					{
						recordID = iter->first;
						if (singular.find(recordID) == singular.end())
						{
							flag = true;
							if (a_skylines.find(recordID) == a_skylines.end())
							{
								a_skylines.insert(recordID);
							}
						}
					}
				}

				if (flag == false)
				{
					leaves[i]->isPruned = true;
					cell* tmp = new cell(leaves[i]);
					//leaves[i]->appendto(tmp);
					finalResult.push_back(tmp);
				}
				leaves[i]->release();
				delete leaves[i];
			}
			else
			{
				leaves[i]->release();
				delete leaves[i];
			}
		}
	}
	for (auto iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		singular.insert(*iter);
	}
}

void cellTree::markSingular(vector<cell*>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular)
{
	long int recordID;
	a_skylines.clear();
	//finalResult.clear();
	bool flag = false;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i]->isPruned == false)
		{
			if (leaves[i]->rank <= mink)
			{
				for (auto iter = leaves[i]->BelowHP.begin(); iter != leaves[i]->BelowHP.end(); iter++)
				{
					recordID = *iter;
					if (singular.find(recordID) == singular.end())
					{
						flag = true;
						if (a_skylines.find(recordID) == a_skylines.end())
						{
							a_skylines.insert(recordID);
						}
					}
				}

				for (auto iter = leaves[i]->IntersectHP.begin(); iter != leaves[i]->IntersectHP.end(); iter++)
				{
					if (iter->second == true)
					{
						recordID = iter->first;
						if (singular.find(recordID) == singular.end())
						{
							flag = true;
							if (a_skylines.find(recordID) == a_skylines.end())
							{
								a_skylines.insert(recordID);
							}
						}
					}
				}

				if (flag == false)
				{
					leaves[i]->isPruned = true;
					//cell* tmp = new cell(leaves[i]);
					cell* tmp = new cell();
					leaves[i]->copyleaf(tmp);
					//leaves[i]->appendto(tmp);
					finalResult.push_back(tmp);
					return; // return as we are sure that record r will be utk;
				}
				leaves[i]->release();
				delete leaves[i];
			}
			else
			{
				leaves[i]->release();
				delete leaves[i];
			}
		}
	}
	for (auto iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		singular.insert(*iter);
	}
}

void cellTree::markSingular(vector<pair<cell*, long int>>& leaves, unordered_set<long int>& a_skylines, const int& mink, vector<cell*>& finalResult, unordered_set<long int>& singular)
{
	long int recordID;
	//finalResult.clear();
	bool flag = false;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (leaves[i].first->rank < mink)
			{
				for (auto iter = leaves[i].first->BelowHP.begin(); iter != leaves[i].first->BelowHP.end(); iter++)
				{
					recordID = *iter;
					if (singular.find(recordID) == singular.end())
					{
						flag = true;
						if (a_skylines.find(recordID) == a_skylines.end())
						{
							a_skylines.insert(recordID);
						}
						break;
					}
				}

				for (auto iter = leaves[i].first->IntersectHP.begin(); iter != leaves[i].first->IntersectHP.end(); iter++)
				{
					if (iter->second == true)
					{
						recordID = iter->first;
						if (singular.find(recordID) == singular.end())
						{
							flag = true;
							if (a_skylines.find(recordID) == a_skylines.end())
							{
								a_skylines.insert(recordID);
							}
							break;
						}
					}
				}

				if (flag == false)
				{
					leaves[i].first->isPruned = true;
					cellID[leaves[i].second]->isPruned = true;
					cell* tmp = new cell();
					leaves[i].first->copyleaf(tmp);
					//leaves[i]->appendto(tmp);
					finalResult.push_back(tmp);
				}
				//leaves[i].first->release();
				//delete leaves[i];
			}
			else
			{
				leaves[i].first->release();
				//delete leaves[i];
			}
		}
	}
	/*for (auto iter = a_skylines.begin(); iter != a_skylines.end(); iter++)
	{
		singular.insert(*iter);
	}*/
}

void cellTree::updateSkyline(unordered_set<long int>& skylines, vector<pair<cell*, long int>>& leaves, unordered_set<long int>& singular)
{
	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			for (auto iter = leaves[i].first->BelowHP.begin(); iter != leaves[i].first->BelowHP.end(); iter++)
			{
				if (singular.find(*iter) == singular.end())
				{
					singular.insert(*iter);
					if (skylines.find(*iter) != skylines.end())
					{
						skylines.erase(*iter);
					}
					else
					{
						cout << "please check what happened!" << endl;
					}
				}
			}

			for (auto iter = leaves[i].first->IntersectHP.begin(); iter != leaves[i].first->IntersectHP.end(); iter++)
			{
				if (iter->second == true && singular.find(iter->first) == singular.end())
				{
					singular.insert(iter->first);
					if (skylines.find(iter->first) != skylines.end())
					{
						skylines.erase(iter->first);
					}
					else
					{
						cout << "please check what happened!" << endl;
					}
				}
			}
		}
	}
}

void cellTree::updateCellTree(cell* node)
{
	if (node->isPruned == true || node == NULL)
		return ;
	else if (node->left == NULL&&node->right == NULL)
		return;
	else if (node->left->isPruned == true && node->right->isPruned == true)
	{
		node->isPruned = true;
		releaseCell(node->left);
		releaseCell(node->right);
		node->left = NULL;
		node->right = NULL;
		return;
	}
	else
	{
		updateCellTree(node->left);
		updateCellTree(node->right);
		if (node->left->isPruned == true && node->right->isPruned == true)
		{
			node->isPruned = true;
			releaseCell(node->left);
			releaseCell(node->right);
			node->left = NULL;
			node->right = NULL;
			return;
		}
	}
}

int cellTree::Lemma2(vector<cell*>& cells)
{
	int touching = 0;
	int count = 0;
	for (int i = 0; i < cells.size() && count<=100; i++)
	{
		double row[MAXDIMEN];
		int dimen = HalfSpaces[0].size() - 2;
		lprec *lp = make_lp(0, dimen);
		lpModel(lp, dimen, cells[i]->IntersectHP);
		touching += cells[i]->IntersectHP.size();
		//addHP(lp, hpid, sideindicator);

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
		delete_lp(lp);
		count++;
	}
	return touching;
}

void cellTree::halfspace2polytope(cell* ret, string outfile, string hsfile)
{
	double row[MAXDIMEN];

	int dimen = HalfSpaces[0].size() - 2;
	lprec *lp = make_lp(0, dimen);

	lpModelwithSIDELEN(lp, dimen, ret->IntersectHP);

	row[0] = 0;
	for (int dii = 1; dii < dimen + 1; dii++)
	{
		row[dii] = 1;
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int tmp = solve(lp);
	get_variables(lp, row);

	if (tmp == 0)
	{
		// qhull
		FILE *fout1 = fopen(hsfile.c_str(), "w+");
		if (fout1 == NULL)
		{
			cout << "Error in opening file innerpoint.txt" << endl;
			getchar();
			exit(0);
		}

		fprintf(fout1, "%d 1\n", dimen);   //the dimensionality
		for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i]); //the feasible point found by lp above   *****************
		fprintf(fout1, "\n");
		fprintf(fout1, "%d\n", dimen + 1);   //dimensionality + 1
		fprintf(fout1, "%ld\n", ret->IntersectHP.size() + 2 * (dimen + 1));  //the total number of hyperplanes for intersection

		// constraints in each dimension
		for (int di = 0; di < dimen; di++)
		{
			row[0] = 0;
			for (int dii = 1; dii < dimen + 1; dii++)
			{
				if (dii - 1 == di)
				{
					row[dii] = 1;
				}
				else
				{
					row[dii] = 0;
				}
			}

			//add_constraint(lp, row, GE, 0);
			for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -row[i + 1]); // dimension constraint  *****************
			fprintf(fout1, "0 \n ");
			for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i + 1]); // dimension constraint *****************
			fprintf(fout1, "-1 \n ");
		}

		// in reduced space, sum_{q_i} should less than 1
		for (int dii = 1; dii < dimen + 1; dii++)
		{
			row[dii] = 1;
		}
		for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -1 * row[i + 1]); // dimension constraint  *****************
		fprintf(fout1, "0 \n ");
		for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i + 1]); // dimension constraint *****************
		fprintf(fout1, "-1 \n ");

		// constraints in intersected hyperplanes
		int hpid, hsID;
		bool sideindicator;
		for (unordered_map<long int, bool>::iterator iter = ret->IntersectHP.begin(); iter != ret->IntersectHP.end(); iter++)
		{
			hpid = iter->first;
			sideindicator = iter->second;
			hsID = RecordIDtoHalfPlaneID[hpid] - 1;

			row[0] = 0;
			for (int dii = 1; dii < dimen + 1; dii++)
			{
				row[dii] = HalfSpaces[hsID][dii - 1];
			}
			if (sideindicator == false)
			{
				for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", row[i + 1]); // dimension constraint *****************
				fprintf(fout1, "%f \n ", -1 * HalfSpaces[hsID][dimen]);
			}
			else if (sideindicator == true)
			{
				for (int i = 0; i < dimen; i++) fprintf(fout1, "%f ", -row[i + 1]); // dimension constraint  *****************
				fprintf(fout1, "%f \n ", HalfSpaces[hsID][dimen]);
			}
			else
			{
				std::cout << "Unable to detect half plane direction!!!" << endl;
			}
		}
		fclose(fout1);

		char sys_string[4096] = "D:\\MonoTopK\\src\\qhull\\qhalf.exe Fp < ";
		strcat(sys_string, hsfile.c_str());
		strcat(sys_string, " | D:\\MonoTopK\\src\\qhull\\qconvex.exe FS > ");
		strcat(sys_string, outfile.c_str());

		system(sys_string);
	}
	return;
}

void cellTree::rankBound(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult)
{
	rBound itemR;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (i % 100 == 0)
				cout << i << endl;
			cellBound(leaves[i].first, rt, a_pt, itemR, mink);
			if (itemR.minimum > mink)
			{
				leaves[i].first->isPruned = true;
			}
			else if (itemR.maximum <= mink)
			{
				leaves[i].first->isPruned = true;
				cellID[leaves[i].second]->isPruned = true;

				cell* tmp = new cell();
				leaves[i].first->copyleaf(tmp);
				tmp->rank = itemR.maximum;
				finalResult.push_back(tmp);
				
			}
			//releaseCell(leaves[i]);
		}
	}
	//leaves.clear();
	//updateCellTree(root);
}

void cellTree::rtreeRAM(Rtree& rt, unordered_map<long int, RtreeNode*>& ramTree)
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
		totalSpaceCost += 4104;
		if (node->isLeaf() == false)
		{
			for (int i = 0; i < node->m_usedspace; i++)
			{
				H.push(node->m_entry[i]->m_id);
			}
		}
	}
}

void cellTree::cellBound(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink)
{
	Score pScore, eScore;
	itemR.minimum = 0;
	itemR.maximum = objCnt - 1;
	int dimen = HalfSpaces[0].size() - 2;
	lprec* lp = make_lp(0, dimen);
	int isOptimal;
	lpModel(lp, dimen, item->IntersectHP);
	focalScore(lp, a_pt, pScore);

	vector<float> cl, cu;
	float clsum = 0, cusum = 0;
	findCellMBR(lp, cl, cu);
	for (int di = 0; di < dimen; di++)
	{
		clsum += cl[di];
		cusum += cu[di];
	}

	cl.push_back(max(1.0 - cusum, 0.0));
	cu.push_back(min(1.0 - clsum, 1.0));

	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	bool isLeafNode;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = ramTree[pageID];
		if (node->isLeaf() == false)
		{
			isLeafNode = false;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				approxiRecordScore(node->m_entry[i], eScore, cl, cu, isLeafNode);
				if (eScore.min >= pScore.max)
				{
					itemR.minimum += node->m_entry[i]->num_records;
					if (itemR.minimum > mink)
						return;
				}
				else if (eScore.max < pScore.min)
				{
					// the records in this rtree node will affect the rank of given cell "item"
					itemR.maximum -= node->m_entry[i]->num_records;
					if (itemR.maximum < mink)
						return;
				}
				else
				{
					H.push(node->m_entry[i]->m_id);
				}
			}
		}
		else
		{
			isLeafNode = true;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				approxiRecordScore(node->m_entry[i], eScore, cl, cu, isLeafNode);
				if (eScore.min >= pScore.max)
				{
					itemR.minimum += 1;
					if (itemR.minimum > mink)
						return;
				}
				else if (eScore.max < pScore.min)
				{
					itemR.maximum -= 1;
					if (itemR.maximum < mink)
						return;
				}
				else
				{
					/*recordScore(lp, node->m_entry[i], eScore, isLeafNode);
					if (eScore.min >= pScore.max)
					{
						itemR.minimum += 1;
						if (itemR.minimum > mink)
							return;
					}
					else if (eScore.max < pScore.min)
					{
						itemR.maximum -= 1;
						if (itemR.maximum < mink)
							return;
					}*/
				}
			}
		}
	}
	delete_lp(lp);
}

void cellTree::focalScore(lprec* lp, Point& a_pt, Score& pScore)
{
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;
	int isOptimal;

	// obtain min score for focal record
	row[0] = 0;
	for (int dii = 0; dii < dimen; dii++)
	{
		//row[dii + 1] = a_pt[dii];
		row[dii + 1] = a_pt[dii] - a_pt[dimen];
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, row);

	if (isOptimal == 2)
	{
		cout << "focalScore LP_solver: please check what happened!!!" << endl;
	}

	pScore.min = a_pt[dimen];
	for (int dii = 0; dii < dimen; dii++)
	{
		pScore.min += (a_pt[dii] - a_pt[dimen]) * row[dii];
	}

	//double tmp = a_pt[dimen] + get_objective(lp);

	// obtain max score for focal record
	row[0] = 0;
	for (int dii = 0; dii < dimen; dii++)
	{
		//row[dii + 1] = -a_pt[dii];
		row[dii + 1] = -(a_pt[dii] - a_pt[dimen]);
	}
	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, row);

	if (isOptimal == 2)
	{
		cout << "focalScore LP_solver: please check what happened!!!" << endl;
	}

	pScore.max = a_pt[dimen];
	for (int dii = 0; dii < dimen; dii++)
	{
		pScore.max += (a_pt[dii] - a_pt[dimen]) * row[dii];
	}
}

void cellTree::recordScore(lprec* lp, RtreeNodeEntry* e, Score& entryR, bool& isLeafNode)
{
	int dimen = HalfSpaces[0].size() - 2;
	double r_lower[MAXDIMEN];
	double r_upper[MAXDIMEN];
	int isOptimal;

	r_lower[0] = 0.0;
	r_upper[0] = 0.0;

	if (isLeafNode == false)
	{
		for (int di = 0; di < dimen; di++)
		{
			r_upper[di + 1] = -(e->m_hc.getUpper()[di] - e->m_hc.getUpper()[dimen]);
			r_lower[di + 1] = e->m_hc.getLower()[di] - e->m_hc.getLower()[dimen];
		}
	}
	else
	{
		for (int di = 0; di < dimen; di++)
		{
			r_upper[di + 1] = -((e->m_hc.getLower()[di] + e->m_hc.getUpper()[di]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2);
			r_lower[di + 1] = ((e->m_hc.getLower()[di] + e->m_hc.getUpper()[di]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2);
		}
	}

	set_obj_fn(lp, r_upper);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, r_upper);

	if (isOptimal == 2)
	{
		cout << "recordScore LP_solver: please check what happened!!!" << endl;
	}

	if (isLeafNode == false)
	{
		entryR.max = e->m_hc.getUpper()[dimen];
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.max += (e->m_hc.getUpper()[dii] - e->m_hc.getUpper()[dimen]) * r_upper[dii];
		}
	}
	else
	{
		entryR.max = (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2;
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.max += ((e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2) * r_upper[dii];
		}
	}


	set_obj_fn(lp, r_lower);
	set_verbose(lp, IMPORTANT);
	isOptimal = solve(lp);
	get_variables(lp, r_lower);

	if (isOptimal == 2)
	{
		cout << "recordScore LP_solver: please check what happened!!!" << endl;
	}

	if (isLeafNode == false)
	{
		entryR.min = e->m_hc.getLower()[dimen];
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.min += (e->m_hc.getLower()[dii] - e->m_hc.getLower()[dimen]) * r_lower[dii];
		}
	}
	else
	{
		entryR.min = (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2;
		for (int dii = 0; dii < dimen; dii++)
		{
			entryR.min += ((e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 - (e->m_hc.getUpper()[dimen] + e->m_hc.getLower()[dimen]) / 2) * r_lower[dii];
		}
	}
}

void cellTree::findCellMBR(lprec* lp, vector<float>& cl, vector<float>& cu)
{
	double row[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;
	int isOptimal;
	cl.clear();
	cu.clear();

	// obtain min score for focal record
	for (int i = 0; i < dimen; i++)
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (dii == i)
				row[dii + 1] = 1;
			else
				row[dii + 1] = 0;
		}
		set_obj_fn(lp, row);
		set_verbose(lp, IMPORTANT);
		isOptimal = solve(lp);
		get_variables(lp, row);
		if (isOptimal == 2)
		{
			cout << "focalScore LP_solver: please check what happened!!!" << endl;
		}
		else
		{
			cl.push_back(row[i]);
		}


		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (dii == i)
				row[dii + 1] = -1;
			else
				row[dii + 1] = 0;
		}
		set_obj_fn(lp, row);
		set_verbose(lp, IMPORTANT);
		isOptimal = solve(lp);
		get_variables(lp, row);
		if (isOptimal == 2)
		{
			cout << "focalScore LP_solver: please check what happened!!!" << endl;
		}
		else
		{
			cu.push_back(row[i]);
		}
	}
}

void cellTree::approxiRecordScore(RtreeNodeEntry* e, Score& entryR, vector<float>& cl, vector<float>& cu, bool& isLeafNode)
{
	int dimen = HalfSpaces[0].size() - 2;
	entryR.max = 0;
	entryR.min = 0;

	if (isLeafNode == false)
	{
		for (int dii = 0; dii < dimen + 1; dii++)
		{
			entryR.max += e->m_hc.getUpper()[dii] * cu[dii];
			entryR.min += e->m_hc.getLower()[dii] * cl[dii];
		}
	}
	else
	{
		for (int dii = 0; dii < dimen + 1; dii++)
		{
			entryR.max += (e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 * cu[dii];
			entryR.min += (e->m_hc.getUpper()[dii] + e->m_hc.getLower()[dii]) / 2 * cl[dii];
		}
	}
}

void cellTree::rankBoundOpt(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult)
{
	rBound itemR;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (i % 100 == 0)
				cout << i << endl;
			cellBoundOpt(leaves[i].first, rt, a_pt, itemR, mink);
			if (itemR.minimum > mink)
			{
				leaves[i].first->isPruned = true;
			}
			else if (itemR.maximum <= mink)
			{
				leaves[i].first->isPruned = true;
				cellID[leaves[i].second]->isPruned = true;

				cell* tmp = new cell();
				leaves[i].first->copyleaf(tmp);
				tmp->rank = itemR.maximum;
				finalResult.push_back(tmp);
			}
		}
	}
}

void cellTree::cellBoundOpt(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink)
{
	itemR.minimum = 0;
	itemR.maximum = 0;
	int dimen = HalfSpaces[0].size() - 2;
	lprec* lp = make_lp(0, dimen);
	int isOptimal;
	lpModel(lp, dimen, item->IntersectHP);
	double obj = 0;

	vector<float> cl, cu;
	float clsum = 0, cusum = 0;
	findCellMBR(lp, cl, cu);
	
	for (int di = 0; di < dimen; di++)
	{
		clsum += cl[di];
		cusum += cu[di];
	}

	cl.push_back(max(1.0 - cusum, 0.0));
	cu.push_back(min(1.0 - clsum, 1.0));
	Score objScore;

	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	bool isLeafNode;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = ramTree[pageID];
		if (node->isLeaf() == false)
		{
			isLeafNode = false;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObj(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min > 0)
				{
					itemR.minimum += node->m_entry[i]->num_records;
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else
				{
					if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += node->m_entry[i]->num_records;
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else
					{
						H.push(node->m_entry[i]->m_id);
					}
				}
			}
		}
		else
		{
			isLeafNode = true;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObj(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min >= 0)
				{
					itemR.minimum += 1;
					itemR.maximum += 1;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += 1;
				}
				else
				{
					if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += 1;
						itemR.maximum += 1;
					}
					else if ((obj = lpOptObj(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += 1;
					}
				}
			}
		}
	}
	delete_lp(lp);
}

double cellTree::lpOptObj(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax)
{
	double ret;
	double row[MAXDIMEN];
	double optobj[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;

	if (isMax == false)
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = (e->m_hc.getLower()[dii] + -e->m_hc.getLower()[dimen]) - (a_pt[dii] - a_pt[dimen]);
			else
				row[dii + 1] = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2) - (a_pt[dii] - a_pt[dimen]);
		}
		if (isLeafNode == false)
		{
			ret = e->m_hc.getLower()[dimen] - a_pt[dimen];
		}
		else
		{
			ret = (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2 - a_pt[dimen];
		}
		
	}
	else
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = -((e->m_hc.getUpper()[dii] - e->m_hc.getUpper()[dimen]) - (a_pt[dii] - a_pt[dimen]));
			else
				row[dii + 1] = -(((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2) - (a_pt[dii] - a_pt[dimen]));
		}
		
		if (isLeafNode == false)
		{
			ret = (e->m_hc.getUpper()[dimen] - a_pt[dimen]);
		}
		else
		{
			ret = (e->m_hc.getLower()[dimen] + e->m_hc.getUpper()[dimen]) / 2 - a_pt[dimen];
		}
	}

	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int tmp = solve(lp);
	get_variables(lp, optobj);
	if (tmp != 0)
	{
		cout << "optimized objective has errors, please check!" << endl;
	}
	else
	{
		double tmp = ret + get_objective(lp);
		if (isMax)
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += -1 * row[dii + 1] * optobj[dii];
			}
		}
		else
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += row[dii + 1] * optobj[dii];
			}
		}
	}
	return ret;
}

void cellTree::fastOptObj(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore)
{
	double weight = 0.0;
	double maxw = 0.0;
	double minw = 0.0;
	int dimen = HalfSpaces[0].size() - 1;
	
	for (int dii = 0; dii < dimen; dii++)
	{
		if (isLeafNode == false)
		{
			minw = (e->m_hc.getLower()[dii] - a_pt[dii]);
			maxw = (e->m_hc.getUpper()[dii] - a_pt[dii]);

			if (maxw < 0)
			{
				objScore.max += minw*cl[dii];
				objScore.min += maxw*cu[dii];
			}
			else if (minw <= 0)
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cu[dii];
			}
			else
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cl[dii];
			}
		}
		else
		{
			weight = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
			if (weight >= 0)
			{
				objScore.max += weight*cu[dii];
				objScore.min += weight*cl[dii];
			}
			else
			{
				objScore.max += weight*cl[dii];
				objScore.min += weight*cu[dii];
			}
		}
	}
	
}

void cellTree::cellBoundOptFULL(cell* item, Rtree& rt, Point& a_pt, rBound& itemR, const int& mink)
{
	itemR.minimum = 0;
	itemR.maximum = 0;
	int dimen = HalfSpaces[0].size() - 2;
	lprec* lp = make_lp(0, dimen);
	int isOptimal;
	lpModel(lp, dimen, item->IntersectHP);
	double obj = 0;

	vector<float> cl, cu;
	float clsum = 0, cusum = 0;
	findCellMBR(lp, cl, cu);
	Score objScore;

	queue<long int> H;
	RtreeNode* node;
	H.push(rt.m_memory.m_rootPageID);
	long int pageID;
	bool isLeafNode;
	while (!H.empty())
	{
		pageID = H.front();
		H.pop();
		node = ramTree[pageID];
		if (node->isLeaf() == false)
		{
			isLeafNode = false;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObjFULL(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min > 0)
				{
					itemR.minimum += node->m_entry[i]->num_records;
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += node->m_entry[i]->num_records;
				}
				else
				{
					if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += node->m_entry[i]->num_records;
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += node->m_entry[i]->num_records;
					}
					else
					{
						H.push(node->m_entry[i]->m_id);
					}
				}
			}
		}
		else
		{
			isLeafNode = true;
			for (int i = 0; i < node->m_usedspace; i++)
			{
				fastOptObjFULL(cu, cl, node->m_entry[i], a_pt, isLeafNode, objScore);
				if (objScore.min >= 0)
				{
					itemR.minimum += 1;
					itemR.maximum += 1;
				}
				else if (objScore.max > 0)
				{
					itemR.maximum += 1;
				}
				else
				{
					if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, false)) >= 0)
					{
						itemR.minimum += 1;
						itemR.maximum += 1;
					}
					else if ((obj = lpOptObjFULL(lp, node->m_entry[i], a_pt, isLeafNode, true)) > 0)
					{
						itemR.maximum += 1;
					}
				}
			}
		}
	}
	delete_lp(lp);
}

double cellTree::lpOptObjFULL(lprec* lp, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, const bool& isMax)
{
	double ret = 0.0;
	double row[MAXDIMEN];
	double optobj[MAXDIMEN];
	int dimen = HalfSpaces[0].size() - 2;

	if (isMax == false)
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = (e->m_hc.getLower()[dii] - a_pt[dii]);
			else
				row[dii + 1] = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
		}
	}
	else
	{
		row[0] = 0;
		for (int dii = 0; dii < dimen; dii++)
		{
			if (isLeafNode == false)
				row[dii + 1] = -(e->m_hc.getUpper()[dii] - a_pt[dii]);
			else
				row[dii + 1] = -((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
		}

	}

	set_obj_fn(lp, row);
	set_verbose(lp, IMPORTANT);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);

	int tmp = solve(lp);
	get_variables(lp, optobj);
	if (tmp != 0)
	{
		cout << "optimized objective has errors, please check!" << endl;
	}
	else
	{
		double tmp = ret + get_objective(lp);
		if (isMax)
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += -1 * row[dii + 1] * optobj[dii];
			}
		}
		else
		{
			for (int dii = 0; dii < dimen; dii++)
			{
				ret += row[dii + 1] * optobj[dii];
			}
		}
	}
	return ret;
}

void cellTree::rankBoundOptFULL(vector<pair<cell*, long int>>& leaves, Rtree& rt, Point& a_pt, const int& mink, vector<cell*>& finalResult)
{
	rBound itemR;

	for (int i = 0; i < leaves.size(); i++)
	{
		if (leaves[i].first->isPruned == false)
		{
			if (i % 100 == 0)
				cout << i << endl;
			cellBoundOptFULL(leaves[i].first, rt, a_pt, itemR, mink);
			if (itemR.minimum > mink)
			{
				leaves[i].first->isPruned = true;
			}
			else if (itemR.maximum <= mink)
			{
				leaves[i].first->isPruned = true;
				cellID[leaves[i].second]->isPruned = true;

				cell* tmp = new cell();
				leaves[i].first->copyleaf(tmp);
				tmp->rank = itemR.maximum;
				finalResult.push_back(tmp);
			}
		}
	}
}

void cellTree::fastOptObjFULL(vector<float>& cu, vector<float>& cl, RtreeNodeEntry* e, Point& a_pt, bool& isLeafNode, Score& objScore)
{
	double weight = 0.0;
	double maxw = 0.0;
	double minw = 0.0;
	int dimen = HalfSpaces[0].size() - 2;

	for (int dii = 0; dii < dimen; dii++)
	{
		if (isLeafNode == false)
		{
			minw = (e->m_hc.getLower()[dii] - a_pt[dii]);
			maxw = (e->m_hc.getUpper()[dii] - a_pt[dii]);

			if (maxw < 0)
			{
				objScore.max += minw*cl[dii];
				objScore.min += maxw*cu[dii];
			}
			else if (minw <= 0)
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cu[dii];
			}
			else
			{
				objScore.max += maxw*cu[dii];
				objScore.min += minw*cl[dii];
			}
		}
		else
		{
			weight = ((e->m_hc.getLower()[dii] + e->m_hc.getUpper()[dii]) / 2 - a_pt[dii]);
			if (weight >= 0)
			{
				objScore.max += weight*cu[dii];
				objScore.min += weight*cl[dii];
			}
			else
			{
				objScore.max += weight*cl[dii];
				objScore.min += weight*cu[dii];
			}
		}
	}

}