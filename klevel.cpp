#include "klevel.h"
#include "header.h"
#include <algorithm>
#include <sstream>
#include <bitset>

KLEVEL::KLEVEL() {
}

KLEVEL::KLEVEL(int d, int Ikmax, int Qkmax) {
	Ik_max = Ikmax;
	dim = d;
	Qk_max = Qkmax;
}

KLEVEL::~KLEVEL() {
}

void KLEVEL::WriteToDisk(string outfile) {
    ofstream fout(outfile.c_str(), std::ofstream::binary);
    fout.write((char*)& Ik_max, sizeof(int));
    for (int i = 1; i <= Ik_max; i++) {
        int level_size = L[i].size();
        fout.write((char*)& level_size, sizeof(int));
        for (int j = 0; j < level_size; j++) {
            // objID
            fout.write((char*)& L[i][j].objID, sizeof(int));

            // topkset
            int confirmed_size = L[i][j].confirmed.size();
            fout.write((char*)& confirmed_size, sizeof(int));
            for (auto it = L[i][j].confirmed.begin(); it != L[i][j].confirmed.end(); it++) {
                fout.write((char*) & *it, sizeof(int));
            }

            // grouping set;
            int g_size = L[i][j].grouping.size();
            fout.write((char*)& g_size, sizeof(int));
            for (auto it = L[i][j].grouping.begin(); it != L[i][j].grouping.end(); it++) {
                fout.write((char*) & *it, sizeof(int));
            }

            // V
            int v_size = L[i][j].r.vertices.size();
            fout.write((char*)& v_size, sizeof(int));
            for (int t = 0; t < v_size; t++) {
                for (int d = 0; d < dim - 1; d++) {
                    fout.write((char*)& L[i][j].r.vertices[t][d], sizeof(float));
                }
            }

            // H
            int HP_size = L[i][j].r.HPs.size();
            fout.write((char*)& HP_size, sizeof(int));
            for (auto it = L[i][j].r.HPs.begin(); it != L[i][j].r.HPs.end(); it++) {
                fout.write((char*)& (it->first), sizeof(int));
                fout.write((char*)& (it->second), sizeof(bool));
            }

            if (i == Ik_max) {
                // candidates;
                int cand_size = L[i][j].candidates.size();
                fout.write((char*)& cand_size, sizeof(int));
                for (auto it = L[i][j].candidates.begin(); it != L[i][j].candidates.end(); it++) {
                    fout.write((char*) & (*it), sizeof(int));
                }
            }
        }
    }
    fout.close();
}

// need fix
void KLEVEL::WriteToDisk(string outfile, int i) {
    if (i==0){
        ofstream fout(outfile.c_str(), std::ofstream::binary);
        fout.write((char*)& Ik_max, sizeof(int));
        fout.close();
        return;
    }

    ofstream fout(outfile.c_str(), std::ofstream::binary | std::ofstream::app);

    int level_size = L[i].size();
    fout.write((char*)& level_size, sizeof(int));
    for (int j = 0; j < level_size; j++) {
        // objID
        fout.write((char*)& L[i][j].objID, sizeof(int));

        // topkset
        int confirmed_size = L[i][j].confirmed.size();
        fout.write((char*)& confirmed_size, sizeof(int));
        for (auto it = L[i][j].confirmed.begin(); it != L[i][j].confirmed.end(); it++) {
            fout.write((char*) & *it, sizeof(int));
        }

        // grouping set;
        int g_size = L[i][j].grouping.size();
        fout.write((char*)& g_size, sizeof(int));
        for (auto it = L[i][j].grouping.begin(); it != L[i][j].grouping.end(); it++) {
            fout.write((char*) & *it, sizeof(int));
        }

        // V
        int v_size = L[i][j].r.vertices.size();
        fout.write((char*)& v_size, sizeof(int));
        for (int t = 0; t < v_size; t++) {
            for (int d = 0; d < dim - 1; d++) {
                fout.write((char*)& L[i][j].r.vertices[t][d], sizeof(float));
            }
        }

        // H
        int HP_size = L[i][j].r.HPs.size();
        fout.write((char*)& HP_size, sizeof(int));
        for (auto it = L[i][j].r.HPs.begin(); it != L[i][j].r.HPs.end(); it++) {
            fout.write((char*) & (it->first), sizeof(int));
            fout.write((char*) & (it->second), sizeof(bool));
        }

        if (i == Ik_max) {
            // candidates;
            int cand_size = L[i][j].candidates.size();
            fout.write((char*)& cand_size, sizeof(int));
            for (auto it = L[i][j].candidates.begin(); it != L[i][j].candidates.end(); it++) {
                fout.write((char*) & (*it), sizeof(int));
            }
        }
    }

    fout.close();
}

void KLEVEL::ReadFromDisk(vector<Obj>& All_records, string infile) {
    // Precompute hyperplane
    AllHP.clear();
    for (int i = 0; i < All_records.size(); i++) {
        for (int j = 0; j < All_records.size(); j++) {
            int tmp = Halfspace::ComputeHP(All_records[i], All_records[j], dim);
        }
    }

    L.clear();
    ifstream fin(infile.c_str(), std::ifstream::binary);
    int k_max;
    fin.read((char*)& k_max, sizeof(int));
    for (int i = 0; i <= k_max; i++) {
        vector<region> tmp; tmp.clear();
        if (i == 0) {
            L.push_back(tmp);
            continue;
        }
        int level_size;
        fin.read((char*)& level_size, sizeof(int));
        for (int j = 0; j < level_size; j++) {
            region cur;
            // objID
            fin.read((char*)& cur.objID, sizeof(int));

            // topkset
            cur.confirmed.clear();
            int confirmed_size;
            fin.read((char*)& confirmed_size, sizeof(int));
            for (int t = 0; t < confirmed_size; t++) {
                int p;
                fin.read((char*)& p, sizeof(int));
                cur.confirmed.insert(p);
            }

            // grouping set;
            cur.grouping.clear();
            int g_size;
            fin.read((char*)& g_size, sizeof(int));
            for (int t = 0; t < g_size; t++) {
                int p;
                fin.read((char*)& p, sizeof(int));
                cur.grouping.insert(p);
            }

            // V
            cur.r.vertices.clear();
            int v_size;
            fin.read((char*)& v_size, sizeof(int));
            for (int t = 0; t < v_size; t++) {
                vector<float> v; v.clear();
                for (int d = 0; d < dim - 1; d++) {
                    float v_d;
                    fin.read((char*)& v_d, sizeof(float));
                    v.push_back(v_d);
                }
                cur.r.vertices.emplace_back(v);
            }

            // H
            cur.r.HPs.clear();
            int HP_size;
            fin.read((char*)& HP_size, sizeof(int));
            for (int t = 0; t < HP_size; t++) {
                int HP;
                bool side;
                fin.read((char*)& HP, sizeof(int));
                fin.read((char*)& side, sizeof(bool));
                cur.r.HPs.insert({ HP,side });
            }

            if (i == Ik_max) {
                // candidates;
                cur.candidates.clear();
                int cand_size;
                fin.read((char*)& cand_size, sizeof(int));
                for (int t = 0; t < cand_size; t++) {
                    int cand;
                    fin.read((char*)& cand, sizeof(int));
                    cur.candidates.insert(cand);
                }
            }

            tmp.emplace_back(cur);
        }
        L.emplace_back(tmp);
    }
    fin.close();
}

bool KLEVEL::Intersect(vector<float>& Qregion, region& cur, int dim) {
    if (cur.r.vertices.size() == 0) return false;
    for (int i = 0; i < cur.r.vertices.size(); i++) {
        if (isIn(cur.r.vertices[i], Qregion, dim)) return true;
    }
    int bitset = 1 << (dim - 1);
    for (int i = 0; i < bitset; i++) {
        vector<float> Qv; Qv.clear();
        int tmp_i = i;
        for (int d = 0; d < dim - 1; d++) {
            if (tmp_i % 2 == 0) Qv.push_back(Qregion[d * 2]);
            else Qv.push_back(Qregion[d * 2 + 1]);
            tmp_i = tmp_i >> 1;
        }
        if (isIn(Qv, cur.r.HPs, dim)) return true;
    }
    return false;
}
bool KLEVEL::isIn(vector<float>& v, vector<float>& Qregion, int dim) {
    for (int i = 0; i < dim - 1; i++) {
        if ((v[i] + eps_klevel >= Qregion[2 * i]) && (v[i] - eps_klevel <= Qregion[2 * i + 1])) continue;
        return false;
    }
    return true;
}
bool KLEVEL::isIn(vector<float>& v, unordered_map<int, bool>& HPs, int dim) {
    for (auto it = HPs.begin(); it != HPs.end(); it++) {
        float sum = 0.0;
        for (int i = 0; i < dim - 1; i++) {
            sum = sum + v[i] * AllHP[it->first].w[i];
        }
        if (it->second == false) {
            if (sum <= AllHP[it->first].w[dim - 1]) continue;
            return false;
        }
        else {
            if (sum >= AllHP[it->first].w[dim - 1]) continue;
            return false;
        }
    }
    return true;
}
void KLEVEL::AddHP(vector<Obj>& All_records, region& cur, vector<pair<int, bool>>& QregionHP) {
    // Re add for compute the vertices of new region(merged)
    cur.r.HPs.clear();
    for (auto it = cur.confirmed.begin(); it != cur.confirmed.end(); it++) {
        if (*it == cur.objID) continue;
        int HP = cur.objID * All_records.size() + *it;
        cur.r.HPs.insert({ HP,false }); // the halfspace S(cur)<S(*it)
    }

    for (auto it = cur.candidates.begin(); it != cur.candidates.end(); it++) {
        int HP = cur.objID * All_records.size() + *it;
        cur.r.HPs.insert({ HP,true }); // the halfspace S(cur)>S(*it)
    }

    for (auto it = QregionHP.begin(); it != QregionHP.end(); it++) {
        cur.r.HPs.insert({it->first,it->second});
    }
}
bool KLEVEL::CheckRepeatTopk(vector<region>& L, region& pre_r, int kth_obj, vector<int>& k_rskyband) {
    bool flag = false;

    for (int i = 0; i < L.size(); i++) {
        if (L[i].objID != kth_obj) continue;
        if (L[i].confirmed.find(kth_obj) == L[i].confirmed.end()) continue;
        bool found = true;
        for (auto it = pre_r.confirmed.begin(); it != pre_r.confirmed.end(); it++) {
            if (L[i].confirmed.find(*it) == L[i].confirmed.end()) {
                found = false;
                break;
            }
        }
        if (found) {
            flag = true;
            for (auto it = k_rskyband.begin(); it != k_rskyband.end(); it++) {
                if (*it == kth_obj) continue;
                L[i].candidates.insert(*it);
            }
        }
    }

    return flag;
}
int KLEVEL::UpdateLevel_forUTK(vector<Obj>& All_records, vector<region>& L, int k, int dim, vector<float>& Qregion, set<int>& utk) {
    // add boundary of query regions
    vector<pair<int, bool>> QregionHP; QregionHP.clear();
    for (int i = 0; i < dim - 1; i++) {
        Hyperplane tmpHP;
        tmpHP.w.clear(); tmpHP.o1 = -1; tmpHP.o2 = -1;
        for (int j = 0; j < dim - 1; j++) {
            if (i == j) tmpHP.w.push_back(1.0);
            else tmpHP.w.push_back(0.0);
        }

        tmpHP.w.push_back(Qregion[2 * i]);
        AllHP.push_back(tmpHP);
        QregionHP.push_back({ AllHP.size() - 1, true });

        tmpHP.w.pop_back();

        tmpHP.w.push_back(Qregion[2 * i + 1]);
        AllHP.push_back(tmpHP);
        QregionHP.push_back({ AllHP.size() - 1, false });
    }

    // add all related cells into L_ans
    vector<vector<region>> L_ans; L_ans.clear();

    L_ans.push_back(L);
    vector<region> tmp; tmp.clear();
    for (int level = 0; level <= k; level++) {
        L_ans.push_back(tmp);
    }

    for (int level = 1; level <= k; level++) {

        for (int i = 0; i < L_ans[level - 1].size(); i++) {
            AddHP(All_records, L_ans[level - 1][i], QregionHP);
            if (level == 1) {
                if (!Halfspace::is_Feasible(L_ans[level - 1][i].r, dim)) continue;
            }
            Halfspace::ComputeVertexByQhull(dim, L_ans[level - 1][i].r);
            if (L_ans[level - 1][i].r.vertices.size() == 0) continue;

            vector<int> one_rskyband, k_rskyband;
            vector<vector<int>> skyband_list, dominators_list;
            krskyband(L_ans[level - 1][i].r, k - level + 1, All_records, L_ans[level - 1][i].candidates, skyband_list, dominators_list, dim);
            one_rskyband = skyband_list[0];
            k_rskyband = skyband_list[k - level];

            for (int j = 0; j < one_rskyband.size(); j++) {

                int cur = one_rskyband[j];
                if (CheckRepeatTopk(L_ans[level], L_ans[level - 1][i], cur, k_rskyband)) continue;

                region tmp;
                CreateNewRegion(All_records, tmp, level, cur, L_ans[level - 1][i], one_rskyband, k_rskyband);

                if (Halfspace::is_Feasible(tmp.r, dim)) {
                    L_ans[level].push_back(tmp);
                }

            }
        }
    }

    for (int i = 0; i < L_ans[k].size(); i++) {
        AddHP(All_records, L_ans[k][i], QregionHP);
        Halfspace::ComputeVertexByQhull(dim, L_ans[k][i].r);
    }

    utk.clear();
    for (int i = 0; i < L_ans[k].size(); i++) {

        for (auto it = L_ans[k][i].confirmed.begin(); it != L_ans[k][i].confirmed.end();it++) {
            utk.insert(*it);
        }
    }

    for (int i = 0; i < dim - 1; i++) {
        AllHP.pop_back();
        AllHP.pop_back();
    }

    return L_ans[k].size();
}

// non-fixed
int KLEVEL::UpdateLevel_forkSPR(vector<Obj>& All_records, vector<region>& L, int k, int dim, int Qid) {
    /*int cnt = 0;

    vector<vector<region>> L_ans; L_ans.clear();
    L_ans.push_back(L);
    vector<region> tmp; tmp.clear();
    for (int level = 1; level <= k; level++) {
        L_ans.push_back(tmp);
    }

    for (int level = 1; level <= k; level++) {
        for (int i = 0; i < L_ans[level - 1].size(); i++) {

            if (level > 1) {
                AddHP_cplex(All_records, L_ans[level - 1][i]);
                Halfspace::ComputeVertexByQhull(dim, L_ans[level - 1][i].r);
            }

            if (L_ans[level - 1][i].r.vertices.size() == 0) {
                //cout << "wrong" << endl;
                continue;
            }


            AddHP(All_records, L_ans[level - 1][i]);
			if (level == 1) {
				if (!Halfspace::is_Feasible_cplex(L_ans[level - 1][i].r, dim)) continue;
			}
			Halfspace::ComputeVertexByQhull(dim, L_ans[level - 1][i].r);

			if (L_ans[level - 1][i].r.vertices.size() == 0) continue;


			vector<int> one_rskyband, k_rskyband;
			vector<vector<int>> skyband_list, dominators_list;
			krskyband(L_ans[level - 1][i].r, k - level + 1, All_records, L_ans[level - 1][i].candidates, skyband_list, dominators_list, dim);
			one_rskyband = skyband_list[0];
			k_rskyband = skyband_list[k - level];

			for (int j = 0; j < one_rskyband.size(); j++) {

				int cur = one_rskyband[j];
				if (CheckRepeatTopk(L_ans[level], L_ans[level - 1][i], cur, k_rskyband)) continue;

				region tmp;
				CreateNewRegion(All_records, tmp, level, cur, L_ans[level - 1][i], one_rskyband, k_rskyband);

				if (Halfspace::is_Feasible_cplex(tmp.r, dim)) {
					if (tmp.confirmed.find(Qid) != tmp.confirmed.end()) cnt++;
					else if (tmp.candidates.find(Qid)!=tmp.candidates.end()) L_ans[level].push_back(tmp);
				}

			}
		}
		cout << level << ' ' << L_ans[level].size() << endl;
		L_ans[level - 1].clear();
		vector<region>().swap(L_ans[level - 1]);
	}

	return cnt;*/
}


bool KLEVEL::XdominateY(Obj& o, Obj& q, int dim) {
    for (int i = 0; i < dim; i++) {
        if (o.w[i] < q.w[i]) return false;
    }
    return true;
}

int KLEVEL::ComputeMinLevel(int q, vector<Obj> & All_records) {
    vector<float> max_diff; max_diff.clear();
    for (int i = 0; i < All_records.size(); i++) {
        float dif = 0.0;
        for (int d = 0; d < dim; d++) {
            dif = max(dif, All_records[q].w[d] - All_records[i].w[d]);
        }
        max_diff.push_back(dif);
    }

    set<int> candidates; candidates.clear();
    int dc = 0;
    for (int i = 0; i < All_records.size(); i++) {
        if (i == q) continue;
        if (XdominateY(All_records[q], All_records[i], dim)) continue;
        if (XdominateY(All_records[i], All_records[q], dim)) {
            dc++;
            continue;
        }
        candidates.insert(i);
    }



    //int remain = Ik_max - dc;
    int remain = (Qk_max - dc) * 1;

    vector<int> temp; temp.clear();

    while (remain > 0) {
        for (int d = 0; d < dim; d++) {
            int i = -1; float i_value = 0.0;
            for (auto it = candidates.begin(); it != candidates.end(); it++) {
                int cur = *it;
                if (All_records[q].w[d] > All_records[cur].w[d]) continue;
                if (((All_records[cur].w[d] - All_records[q].w[d]) / max_diff[cur]) > i_value) {
                    i = cur;
                    i_value = (All_records[cur].w[d] - All_records[q].w[d]) / max_diff[cur];
                }
            }
            temp.push_back(i); candidates.erase(i);
        }
        Halfspace r;
        for (int i = 0; i < temp.size(); i++) {
            if (temp[i] == -1) continue;
            int HP = q * All_records.size() + temp[i];
            r.HPs.insert({ HP,true });
        }
        if (!Halfspace::is_Feasible(r, dim)) {
            dc++;
            temp.clear();
        }
        if (dc >= Qk_max) break;
        if (candidates.size() == 0) break;
        //if (!Halfspace::is_Feasible(r, dim)) dc++;
        remain--;
    }

    return dc + 1;
}

void KLEVEL::GetMinLevel(vector<Obj>& All_records, vector<int>& min_level) {
    clock_t at = clock();

	for (int i = 0; i < All_records.size(); i++) {
		if (i % 100 == 0) cout << "Minlevel Processing: " << i << endl;
		min_level.push_back(ComputeMinLevel(i, All_records));
	}

	cout << "MinLevel Time Cost: " << (clock() - at) / (float)CLOCKS_PER_SEC << endl;
}

void KLEVEL::AddHP(vector<Obj>& All_records, region& cur) {
    // Re add for compute the vertices of new region(merged)
    cur.r.HPs.clear();
    for (auto it = cur.confirmed.begin(); it != cur.confirmed.end(); it++) {
        if (*it == cur.objID) continue;
        int HP = cur.objID * All_records.size() + *it;
        cur.r.HPs.insert({ HP,false }); // the halfspace S(cur)<S(*it)
    }

    for (auto it = cur.candidates.begin(); it != cur.candidates.end(); it++) {
        int HP = cur.objID * All_records.size() + *it;
        cur.r.HPs.insert({ HP,true }); // the halfspace S(cur)>S(*it)
    }
}

void KLEVEL::AddHP_grouping(vector<Obj>& All_records, region& cur) {
    // Re add for compute the vertices of new region(merged)
    cur.r.HPs.clear();
    for (auto it_g = cur.grouping.begin(); it_g != cur.grouping.end(); it_g++) {
        for (auto it = cur.confirmed.begin(); it != cur.confirmed.end(); it++) {
            int HP = *it_g * All_records.size() + *it;
            cur.r.HPs.insert({ HP,false }); // the halfspace S(cur)<S(*it)
        }

        for (auto it = cur.candidates.begin(); it != cur.candidates.end(); it++) {
            int HP = *it_g * All_records.size() + *it;
            cur.r.HPs.insert({ HP,true }); // the halfspace S(cur)>S(*it)
        }
    }
}


void KLEVEL::Bulk_Loading_Init(vector<Obj>& All_records) {
    vector<int> min_level;
    GetMinLevel(All_records, min_level);

    // Init Level-0
    region r_ori;
    Halfspace hs_ori;
    hs_ori.vertices.clear();

    //Generate vertices for the whole space
    vector<float> origin; origin.clear();
    for (int d = 0; d < dim - 1; d++) origin.push_back(0);
    hs_ori.vertices.push_back(origin);
    for (int d = 0; d < dim - 1; d++) {
        origin[d] = 1.0;
        hs_ori.vertices.push_back(origin);
        origin[d] = 0.0;
    }
    for (int d = 0; d < dim - 1; d++) hs_ori.innerPoint[d] = 1.0 / (float)dim;

    r_ori.level = 0; r_ori.objID = -1; r_ori.r = hs_ori;
    r_ori.confirmed.clear(); r_ori.candidates.clear(); r_ori.grouping.clear();

    for (int i = 0; i < All_records.size(); i++) {
        //if (min_level[i] > Ik_max) continue;
        if (min_level[i] > Qk_max) continue;
        r_ori.candidates.insert(i);
    }

    cout << "MinLevel Size: " << r_ori.candidates.size() << endl;

    vector<region> L_zero; L_zero.clear(); L_zero.push_back(r_ori);
    L.push_back(L_zero);

    for (int level = 1; level <= Ik_max; level++) {
        vector<region> tmp_level; tmp_level.clear();
        L.push_back(tmp_level);
    }

}
void KLEVEL::krskyband(Halfspace& r, int a_k, vector<Obj>& All_records, set<int>& candidates, vector<vector<int>>& skyband_list, vector<vector<int>>& dominate_list, int dim) {
    vector<int> dominate, dominators;
    dominate.clear();
    dominators.clear();
    for (auto it = candidates.begin(); it != candidates.end(); it++) {
        dominate.push_back(0);
        dominators.push_back(0);
    }

    int i = 0;
    for (auto it = candidates.begin(); it != candidates.end(); it++) {
        int j = i + 1;
        for (auto itt = next(it); itt != candidates.end(); itt++) {
            int res = isRdominated(r.vertices, All_records[*it].w, All_records[*itt].w, dim);
            if (res == 0) {
                //it<itt
                dominators[i]++;
                dominate[j]++;
            }
            else if (res == 1) {
                // it>itt
                dominate[i]++;
                dominators[j]++;
            }
            else {
                // do thinig;
            }
            j++;
        }
        i++;
    }

    skyband_list.clear();
    for (int i = 1; i <= a_k; i++) {
        vector<int> i_skyband, i_dominate; i_skyband.clear(); i_dominate.clear();
        int j = 0;
        for (auto it = candidates.begin(); it != candidates.end(); it++) {
            if (dominators[j] < i) {
                i_skyband.push_back(*it);
                i_dominate.push_back(dominate[j]);
            }
            j++;
        }
        skyband_list.emplace_back(i_skyband);
        dominate_list.emplace_back(i_dominate);
    }

}
int KLEVEL::isRdominated(vector<vector<float>>& vertices, float focal[], float entry[], int dim) {
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
    int cnt = 0;
    for (int i = 0; i < vertices.size(); i++)
    {
        float sum = 0;
        for (int d = 0; d < dim - 1; d++)
        {
            sum = sum + vertices[i][d] * tmpHS[d];
        }

        if (sum > tmpHS[dim - 1]) cnt++;
        if ((cnt != 0) && (i + 1 > cnt)) return 2;
    }

    if (cnt == 0) return 0;
    if (cnt == vertices.size()) return 1;
    return 2;
}

void KLEVEL::Bulk_Loading_cplex(fstream& fout, vector<Obj>& All_records, int dim, string indexfile) {

    Bulk_Loading_Init(All_records);

    set<int> utk_results; utk_results.clear();

    clock_t start, now;
    start = now = clock();

    // Compute Klevel from level-1 to level-Max_k
    int total_region = 0;

    for (int level = 1; level <= Ik_max; level++) {

        if (level >1) {
            int cnt = 0;
            for (int i = 0; i < L[level - 1].size(); i++) {
                bool flag = true;
                for (auto it = L[level - 1][i].candidates.begin(); it != L[level - 1][i].candidates.end(); it++) {
                    if (utk_results.find(*it)==utk_results.end()) flag = false;
                }
                if (!flag) cnt++;
            }
            cout << L[level - 1].size() << " " << cnt << endl;
            //system("pause");
        }

        float one_candidate = 0;
        float spliting = 0;
        float ave_candidate = 0;
        float ave_hull_vertex = 0;
        float grouping_node = 0;
        float subgraph_pruning = 0;

        vector<int> one_rskyband, k_rskyband;
        vector<vector<int>> skyband_list, dominators_list;

        for (int i = 0; i < L[level - 1].size(); i++) {

            if (level > 1) {
                AddHP(All_records, L[level - 1][i]);
                Halfspace::ComputeVertexByQhull(dim, L[level - 1][i].r);
            }

            if (L[level - 1][i].r.vertices.size() == 0) {
                //cout << "wrong" << endl;
                continue;
            }
            ave_hull_vertex = ave_hull_vertex + L[level - 1][i].r.vertices.size();

            //vector<int> one_rskyband, k_rskyband;
            //vector<vector<int>> skyband_list, dominators_list;
            krskyband(L[level - 1][i].r, Qk_max - level + 1, All_records, L[level - 1][i].candidates, skyband_list, dominators_list, dim);
            one_rskyband = skyband_list[0];
            k_rskyband = skyband_list[Qk_max - level];

            one_candidate = one_candidate + one_rskyband.size();
            ave_candidate = ave_candidate + k_rskyband.size();

			for (int j = 0; j < one_rskyband.size(); j++) {

                int cur = one_rskyband[j];

                if (CheckRepeatTopk(level, L[level - 1][i], cur, k_rskyband)) continue;
                region tmp;
                CreateNewRegion(All_records, tmp, level, cur, L[level - 1][i], one_rskyband, k_rskyband);

                if (Halfspace::is_Feasible(tmp.r, dim)) {
                //if (Halfspace::is_Feasible_cplex(tmp.r, dim)) {
                    spliting++;
                    utk_results.insert(cur);
                    L[level].push_back(tmp);
                }
                //if (i==27046) cout << "lpsolver done" << endl;
            }
        }

        total_region = total_region + L[level].size();

        cout << "LEVEL: " << level << endl;
        cout << "The region size of LEVEL " << level << ": " << L[level].size() << endl;
        cout << "The utk size of LEVEL " << level << ": " << utk_results.size() << endl;
        cout << "Time Cost of LEVEL " << level << ": " << (clock() - now) / (float)CLOCKS_PER_SEC << endl;
        cout << "Average one_skyband size of LEVEL" << ": " << one_candidate / (float)L[level - 1].size() << endl;
        cout << "Average spliting number of LEVEL" << ": " << spliting / (float)L[level - 1].size() << endl;
        cout << "Average candidates of LEVEL" << ": " << ave_candidate / (float)L[level - 1].size() << endl;
        cout << "Average verterices of region in LEVEL" << ": " << ave_hull_vertex / (float)L[level - 1].size() << endl;
        cout << "subgraph_pruning in LEVEL" << ": " << subgraph_pruning / (float)L[level - 1].size() << endl;


        fout << "LEVEL: " << level << endl;
        fout << "The region size of LEVEL " << level << ": " << L[level].size() << endl;
        fout << "The utk size of LEVEL " << level << ": " << utk_results.size() << endl;
        fout << "Time Cost of LEVEL " << level << ": " << (clock() - now) / (float)CLOCKS_PER_SEC << endl;
        fout << "Average one_skyband size of LEVEL" << ": " << one_candidate / (float)L[level - 1].size() << endl;
        fout << "Average spliting number of LEVEL" << ": " << spliting / (float)L[level - 1].size() << endl;
        fout << "Average candidates of LEVEL" << ": " << ave_candidate / (float)L[level - 1].size() << endl;
        fout << "Average verterices of region in LEVEL" << ": " << ave_hull_vertex / (float)L[level - 1].size() << endl;
        fout << "subgraph_pruning in LEVEL" << ": " << subgraph_pruning / (float)L[level - 1].size() << endl;
        now = clock();
        WriteToDisk(indexfile, level-1);
        L[level-1].clear();
        vector<region>().swap(L[level-1]);
    }

	for (int i = 0; i < L[Ik_max].size(); i++) {
		if (L[Ik_max][i].objID != -1) {
			AddHP(All_records, L[Ik_max][i]);
			Halfspace::ComputeVertexByQhull(dim, L[Ik_max][i].r);
		}
		/*else {
			AddHP_grouping(All_records, L[Ik_max][i]);
			Halfspace::ComputeVertexByQhull(dim, L[Ik_max][i].r);
		}*/
	}

	WriteToDisk(indexfile,Ik_max);
	L[Ik_max].clear();

	cout << "total_region: " << total_region << endl;
	cout << "Total Time Cost: " << (clock() - start) / (float)CLOCKS_PER_SEC << endl;
	fout << "total_region: " << total_region << endl;
	fout << "Total Time Cost: " << (clock() - start) / (float)CLOCKS_PER_SEC << endl;

	AllHP.clear();
}

void KLEVEL::CreateNewRegion(vector<Obj>& All_records, region & tmp, int level, int cur, region& pre_r, vector<int>& one_rskyband, vector<int> & k_rskyband) {
    tmp.level = level; tmp.objID = cur;
    tmp.grouping.clear();
    tmp.confirmed = pre_r.confirmed; tmp.confirmed.insert(cur);
    tmp.candidates.clear();
    for (auto it = k_rskyband.begin(); it != k_rskyband.end(); it++) {
        if (*it == cur) continue;
        tmp.candidates.insert(*it);
    }

    // HPs just for verification
    tmp.r.vertices.clear();
    tmp.r.HPs = pre_r.r.HPs;
    for (auto it = one_rskyband.begin(); it != one_rskyband.end(); it++) {
        if (*it == tmp.objID) continue;
        int HP = tmp.objID * All_records.size() + *it;
        tmp.r.HPs.insert({ HP,true }); // the halfspace S(cur)>S(*it)
    }

}

bool KLEVEL::CheckRepeatTopk(int level, region& pre_r, int kth_obj, vector<int>& k_rskyband) {
    bool flag = false;

    for (int i = 0; i < L[level].size(); i++) {
        if (L[level][i].objID != kth_obj) continue;
        if (L[level][i].confirmed.find(kth_obj) == L[level][i].confirmed.end()) continue;
        bool found = true;
        for (auto it = pre_r.confirmed.begin(); it != pre_r.confirmed.end(); it++) {
            if (L[level][i].confirmed.find(*it) == L[level][i].confirmed.end()) {
                found = false;
                break;
            }
        }
        if (found) {
            flag = true;
            for (auto it = k_rskyband.begin(); it != k_rskyband.end(); it++) {
                if (*it == kth_obj) continue;
                L[level][i].candidates.insert(*it);
            }
        }
    }

    return flag;
}

/*bool KLEVEL::CheckRepeatTopk_cplex(int level, region& pre_r, int kth_obj, vector<int>& k_rskyband) {
	bool flag = false;

	for (int i = 0; i < L[level].size(); i++) {
		if (L[level][i].confirmed.find(kth_obj) == L[level][i].confirmed.end()) continue;
		bool found = true;
		for (auto it = pre_r.confirmed.begin(); it != pre_r.confirmed.end(); it++) {
			if (L[level][i].confirmed.find(*it) == L[level][i].confirmed.end()) {
				found = false;
				break;
			}
		}
		if (found) {
			flag = true;
			for (auto it = k_rskyband.begin(); it != k_rskyband.end(); it++) {
				if (*it == kth_obj) continue;
				L[level][i].candidates.insert(*it);
			}
		}
	}

	return flag;
}*/

/*
void KLEVEL::CreateNewRegion_cplex(vector<Obj>& All_records, region& tmp, int level, int cur, region& pre_r, vector<int>& one_rskyband, vector<int>& k_rskyband) {
	tmp.level = level; tmp.objID = -1;
	tmp.grouping.clear();
	tmp.confirmed = pre_r.confirmed; tmp.confirmed.insert(cur);
	tmp.candidates.clear();
	for (auto it = k_rskyband.begin(); it != k_rskyband.end(); it++) {
		if (*it == cur) continue;
		tmp.candidates.insert(*it);
	}

	// HPs just for verification
	tmp.r.HPs = pre_r.r.HPs;
	for (auto it = one_rskyband.begin(); it != one_rskyband.end(); it++) {
		if (*it == tmp.objID) continue;
		int HP = tmp.objID * All_records.size() + *it;
		tmp.r.HPs.insert({ HP,true }); // the halfspace S(cur)>S(*it)
	}

}
*/

/*












bool KLEVEL::CheckRepeatTopk(int level, region& cur) {
	bool flag = false;
	for (int i = 0; i < L[level].size(); i++) {
		if (L[level][i].objID != -1) continue;
		bool found = true;
		for (auto it = cur.confirmed.begin(); it != cur.confirmed.end(); it++) {
			if (L[level][i].confirmed.find(*it) == L[level][i].confirmed.end()) {
				found = false;
				break;
			}
		}
		for (auto it = cur.grouping.begin(); it != cur.grouping.end(); it++) {
			if (L[level][i].grouping.find(*it) == L[level][i].grouping.end()) {
				found = false;
				break;
			}
		}
		if (found) {
			flag = true;
			for (auto it =cur.candidates.begin(); it != cur.candidates.end(); it++) {
				L[level][i].candidates.insert(*it);
			}
		}
	}
	return flag;
}












*/

/*
int KLEVEL::LevelGrouping_step(int offset, vector<int>& candidates, vector<int>& total_candidates, region& this_r, vector<Obj>& All_records, vector<int>& dominators) {
	region tmp;
	tmp.level = this_r.level + offset; tmp.objID = -1; tmp.confirmed = this_r.confirmed; tmp.candidates.clear();
	for (auto it = total_candidates.begin(); it != total_candidates.end(); it++) tmp.candidates.insert(*it);
	for (auto it = candidates.begin(); it != candidates.end(); it++) tmp.grouping.insert(*it);
	int cnt = 0;
	string bitmask(offset, 1);
	bitmask.resize(candidates.size(), 0);
	do {
		bool flag = true;
		for (int i = 0; i < candidates.size(); i++) { // candidates[i] cannot be last (size-k) records
			if ((bitmask[i] == '0') && (dominators[i] >= candidates.size() - offset)) flag = false;
		}
		if (!flag) continue;

		for (int i = 0; i < candidates.size(); i++) {
			if (bitmask[i] == '0') {
				tmp.grouping.erase(candidates[i]);
			}
		}
		for (int i = 0; i < candidates.size(); i++) {
			if ((bitmask[i] == '1')) tmp.candidates.erase(candidates[i]);
		}

		if (!CheckRepeatTopk(tmp.level, tmp)) {
			Halfspace tmp_r = this_r.r;
			// update tmp_r.r
			for (int i = 0; i < candidates.size(); i++) {
				if (bitmask[i] == '0') {
					for (int j = 0; j < candidates.size(); j++) {
						if (bitmask[i] == '1') {
							int HP = candidates[i] * All_records.size() + candidates[j];
							tmp_r.HPs.insert({ HP,false }); // the halfspace S(i)<S(j)
						}
					}
				}
			}

			if (Halfspace::is_Feasible(tmp_r, dim)) {
				L[this_r.level + offset].push_back(tmp);
				cnt++;
			}
		}

		for (int i = 0; i < candidates.size(); i++) {
			if (bitmask[i] == '0') {
				tmp.grouping.insert(candidates[i]);
			}
		}
		for (int i = 0; i < candidates.size(); i++) {
			if ((bitmask[i] == '1')) tmp.candidates.insert(candidates[i]);
		}

	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
	return cnt;
}
*/

/*
void KLEVEL::Bulk_Loading(fstream& fout, vector<Obj>& All_records, int dim, string indexfile) {
	
	Bulk_Loading_Init(All_records);

	set<int> utk_results; utk_results.clear();

	clock_t start, now;
	start = now = clock();

	// Compute Klevel from level-1 to level-Max_k
	//int cnt = 0;
	int total_region = 0;

	for (int level = 1; level <= Ik_max; level++) {

		float one_candidate = 0;
		float spliting = 0;
		float ave_candidate = 0;
		float ave_hull_vertex = 0;
		float grouping_node = 0;

		for (int i = 0; i < L[level - 1].size(); i++) {

			if (L[level - 1][i].objID != -1) {
				AddHP(All_records, L[level - 1][i]);
				Halfspace::ComputeVertexByQhull(dim, L[level - 1][i].r);
			}
			else if (level > 1) {
				AddHP_grouping(All_records, L[level - 1][i]);
				Halfspace::ComputeVertexByQhull(dim, L[level - 1][i].r);
			}

			if (L[level - 1][i].r.vertices.size() == 0) continue;
			ave_hull_vertex = ave_hull_vertex + L[level - 1][i].r.vertices.size();

			*/
/*vector<int> one_rskyband, k_rskyband;
			vector<vector<int>> skyband_list, dominators_list;
			krskyband(L[level - 1][i].r, Ik_max - level + 1, All_records, L[level - 1][i].candidates, skyband_list, dominators_list, dim);
			one_rskyband = skyband_list[0];
			k_rskyband = skyband_list[Ik_max - level];*//*


			vector<int> one_rskyband, k_rskyband;
			vector<vector<int>> skyband_list, dominators_list;
			krskyband(L[level - 1][i].r, Qk_max - level + 1, All_records, L[level - 1][i].candidates, skyband_list, dominators_list, dim);
			one_rskyband = skyband_list[0];
			k_rskyband = skyband_list[Qk_max - level];

			if (GROUPING) {
				// If grouping, CheckRepeat and add candidates!!!
				bool grouping = false;
				for (int offset = Ik_max - level + 1; offset > 1; offset--) {
					if (skyband_list[offset - 1].size() <= offset + Delta) {
						int node_inc = LevelGrouping_step(offset, skyband_list[offset - 1], k_rskyband, L[level - 1][i], All_records, dominators_list[offset - 1]);
						grouping_node = grouping_node + node_inc;
						one_candidate =one_candidate + one_rskyband.size();
						ave_candidate = ave_candidate + k_rskyband.size();
						grouping = true;
						break;
					}
				}
				if (grouping) continue;
			}

			

			one_candidate = one_candidate + one_rskyband.size();
			ave_candidate = ave_candidate + k_rskyband.size();

			for (int j = 0; j < one_rskyband.size(); j++) {
				int cur = one_rskyband[j];
				if (CheckRepeatTopk(level, L[level - 1][i], cur, k_rskyband)) continue;

				region tmp;
				CreateNewRegion(All_records, tmp, level, cur, L[level - 1][i], one_rskyband, k_rskyband);

				if (Halfspace::is_Feasible(tmp.r, dim)) {
					spliting = spliting + 1;
					utk_results.insert(cur);
					L[level].push_back(tmp);
				}
			}
		}

		total_region = total_region + L[level].size();

		cout << "LEVEL: " << level << endl;
		cout << "The region size of LEVEL " << level << ": " << L[level].size() << endl;
		cout << "The utk size of LEVEL " << level << ": " << utk_results.size() << endl;
		cout << "Time Cost of LEVEL " << level << ": " << (clock() - now) / (float)CLOCKS_PER_SEC << endl;
		cout << "Average one_skyband size of LEVEL" << ": " <<  one_candidate / (float)L[level - 1].size() << endl;
		cout << "Average spliting number of LEVEL" << ": " << spliting / (float)L[level - 1].size() << endl;
		cout << "Average candidates of LEVEL" << ": " << ave_candidate / (float)L[level - 1].size() << endl;
		cout << "Average verterices of region in LEVEL" << ": " << ave_hull_vertex / (float)L[level - 1].size() << endl;
		cout << "grouping nodes in LEVEL" << ": " << grouping_node << endl;


		fout << "LEVEL: " << level << endl;
		fout << "The region size of LEVEL " << level << ": " << L[level].size() << endl;
		fout << "The utk size of LEVEL " << level << ": " << utk_results.size() << endl;
		fout << "Time Cost of LEVEL " << level << ": " << (clock() - now) / (float)CLOCKS_PER_SEC << endl;
		fout << "Average one_skyband size of LEVEL" << ": " << one_candidate / (float)L[level - 1].size() << endl;
		fout << "Average spliting number of LEVEL" << ": " << spliting / (float)L[level - 1].size() << endl;
		fout << "Average candidates of LEVEL" << ": " << ave_candidate / (float)L[level - 1].size() << endl;
		fout << "Average verterices of region in LEVEL" << ": " << ave_hull_vertex / (float)L[level - 1].size() << endl;
		fout << "grouping nodes in LEVEL" << ": " << grouping_node << endl;
		now = clock();

		WriteToDisk(indexfile, level - 1);
		L[level - 1].clear();
	}

	for (int i = 0; i < L[Ik_max].size(); i++) {
		if (L[Ik_max][i].objID != -1) {
			AddHP(All_records, L[Ik_max][i]);
			Halfspace::ComputeVertexByQhull(dim, L[Ik_max][i].r);
		}
		else {
			AddHP_grouping(All_records, L[Ik_max][i]);
			Halfspace::ComputeVertexByQhull(dim, L[Ik_max][i].r);
		}
	}

	WriteToDisk(indexfile, Ik_max);
	L[Ik_max].clear();
	
	cout << "total_region: " << total_region << endl;
	cout << "Total Time Cost: " << (clock() - start) / (float)CLOCKS_PER_SEC << endl;
	fout << "total_region: " << total_region << endl;
	fout << "Total Time Cost: " << (clock() - start) / (float)CLOCKS_PER_SEC << endl;
	AllHP.clear();
}




int KLEVEL::GroupingProcess(vector<Obj>& All_records, int k, int dim, vector<float>& Qregion, region& r_st, set<int>& results) {
	region tmp;
	int k_st = r_st.confirmed.size();
	tmp.level = k_st; tmp.r = r_st.r;
	tmp.candidates = r_st.grouping; tmp.confirmed = r_st.confirmed; tmp.grouping.clear();

	double Qr[Max_Dimension + 1];
	for (int i = 0; i < dim - 1; i++) {
		Hyperplane tmpHP;
		tmpHP.w.clear(); tmpHP.o1 = -1; tmpHP.o2 = -1;
		for (int j = 0; j < dim - 1; j++) {
			if (i == j) tmpHP.w.push_back(1.0);
			else tmpHP.w.push_back(0.0);
		}

		tmpHP.w.push_back(Qregion[2 * i]);
		AllHP.push_back(tmpHP);
		tmp.r.HPs.insert({ AllHP.size() - 1, true });

		tmpHP.w.pop_back();

		tmpHP.w.push_back(Qregion[2 * i + 1]);
		AllHP.push_back(tmpHP);
		tmp.r.HPs.insert({ AllHP.size() - 1, false });
		
	}
	Halfspace::ComputeVertexByQhull(dim, tmp.r);

	vector<int> one_rskyband, k_rskyband;
	vector<vector<int>> skyband_list, dominators_list;

	list<region> q; q.emplace_back(tmp); int last = q.size();
	for (int offset = 0; offset < k - tmp.level; offset++) {
		for (int i = 0; i < last; i++) {
			auto pre = q.begin();

			krskyband(pre->r, k - k_st, All_records, pre->candidates, skyband_list, dominators_list, dim);
			one_rskyband = skyband_list[0];
			k_rskyband = skyband_list[k - k_st];

			for (int j = 0; j < one_rskyband.size(); j++) {
				int cur = one_rskyband[j];
				region new_r;
				CreateNewRegion(All_records, new_r, pre->level+1, cur, *pre, one_rskyband, k_rskyband);

				if (Halfspace::is_Feasible(new_r.r, dim)) {
					Halfspace::ComputeVertexByQhull(dim, new_r.r);
					q.emplace_back(new_r);
				}
			}

			q.pop_front();
		}
		last = q.size();
	}
	for (int i = 0; i < dim - 1; i++) {	
		AllHP.pop_back();
		AllHP.pop_back();
	}
	return q.size();
}
*/


// large K

/*
void KLEVEL::rskyband_singlek(Halfspace& r, int a_k, vector<Obj>& All_records, set<int>& candidates, vector<int>& skyband, int dim) {

	skyband.clear();
	priority_queue<pair<float, int>> heap;
	for (auto it = candidates.begin(); it != candidates.end(); it++) {
		float score = 0.0;
		float weight = 0.0;
		for (int d = 0; d < dim-1; d++) {
			score = score + r.innerPoint[d] * All_records[*it].w[d];
			weight = weight + All_records[*it].w[d];
		}
		score = score + (1 - weight) * All_records[*it].w[dim - 1];
		heap.push(make_pair(score, *it));
	}

	while (!heap.empty()) {
		int id = heap.top().second;
		heap.pop();

		int dominate = 0;
		for (auto it = skyband.begin(); it != skyband.end(); it++) {
			int res = isRdominated(r.vertices, All_records[id].w, All_records[*it].w, dim);
			if (res == 0) dominate++;
			if (dominate >= a_k) break;
		}
		if (dominate < a_k) skyband.push_back(id);
	}

}











void KLEVEL::Bulk_Loading_singleObj(fstream& fout, vector<Obj>& All_records, int dim, string indexfile) {
	vector<int> min_level;
	GetMinLevel(All_records, min_level);

	for (int Oq = 1000; Oq < All_records.size(); Oq++) {
		if (min_level[Oq] > Ik_max) continue;
		clock_t now = clock();

		vector<region> L; L.clear();

		// Init Level-0
		region r_ori;
		Halfspace hs_ori;
		hs_ori.HPs.clear();
		hs_ori.vertices.clear();

		//Generate vertices for the whole space
		vector<float> origin; origin.clear();
		for (int d = 0; d < dim - 1; d++) origin.push_back(0);
		hs_ori.vertices.push_back(origin);
		for (int d = 0; d < dim - 1; d++) {
			origin[d] = 1.0;
			hs_ori.vertices.push_back(origin);
			origin[d] = 0.0;
		}
		for (int d = 0; d < dim - 1; d++) hs_ori.innerPoint[d] = 1.0 / (float)dim;

		r_ori.level = 0; r_ori.objID = -1; r_ori.r = hs_ori;
		r_ori.confirmed.clear(); r_ori.candidates.clear(); r_ori.grouping.clear();

		for (int i = 0; i < All_records.size(); i++) {
			//if (min_level[i] > Ik_max) continue;
			if (min_level[i] > Qk_max) continue;
			r_ori.candidates.insert(i);
		}
		L.push_back(r_ori);

		int rq_size = UpdateLevel_forkSPR(All_records, L, Ik_max, dim, Oq);
		cout << Oq << ": " << rq_size << endl;
		cout << "Time Cost of " << Oq << ": " << (clock() - now) / (float)CLOCKS_PER_SEC << endl;

		fout << Oq << ": " << rq_size << endl;
		fout << "Time Cost of " << Oq << ": " << (clock() - now) / (float)CLOCKS_PER_SEC << endl;
	}
}
*/


//int KLEVEL::CreateNode_full(int ID, Halfspace& r) {
//	Node NewNode;
//	NewNode.pID = ID;
//	NewNode.r = r;
//	NewNode.plusEdge.clear();
//	NewNode.minusEdge.clear();
//	G.push_back(NewNode);
//	if (ID != -1) node_list[ID].push_back(G.size() - 1);
//	return G.size() - 1;
//}
//void KLEVEL::Insert_full(Obj & q) {
//	vector<int> q_list; q_list.clear(); node_list.push_back(q_list);
//	vector<int> HP_q; HP_q.clear();
//
//	for (int i = 0; i < AllObj.size(); i++) {
//		HP_q.push_back(Halfspace::ComputeHP(q, AllObj[i], dim));
//	}
//	// Split every region that intersects with HP(q,i)
//	for (int i = 0; i < AllObj.size(); i++) {
//		int HP = HP_q[i];
//		int node_size = node_list[i].size();
//		for (auto j = 0; j < node_size; j++) {
//			if (Halfspace::is_intersected(G[node_list[i][j]].r, dim, HP)) {
//				Halfspace r_LE, r_GE; // q <= o_i in r_LE; q >= o_i in r_GE;
//				Halfspace::Split(G[node_list[i][j]].r, HP, r_LE, r_GE); // ???
//				Update_SplitNode_full(node_list[i][j], r_LE, r_GE); // ???
//			}
//		}
//	}
//
//	// PreCompute Klevel before insertion
//	vector<Halfspace> Rq;
//	vector<int> rank_q, abovelevel, belowlevel;
//	FindVoronoiWithRank_full(q, Rq, rank_q, HP_q);
//
//	set<int> Ks; Ks.clear(); abovelevel.clear(); belowlevel.clear();
//	for (int i = 0; i < rank_q.size(); i++) {
//		if (Ks.find(rank_q[i] - 1) == Ks.end()) Ks.insert(rank_q[i] - 1);
//		if (Ks.find(rank_q[i]) == Ks.end()) Ks.insert(rank_q[i]);
//		abovelevel.push_back(-1);
//		belowlevel.push_back(-1);
//	}
//
//	Klevel_buffer.clear();
//	level temp;
//	for (auto it = Ks.begin(); it != Ks.end(); it++) {
//		temp.k = *it;
//		FindKlevel_full(temp.Nodes, temp.k);
//		Klevel_buffer.push_back(temp);
//		for (int i = 0; i < rank_q.size(); i++) {
//			if (rank_q[i] - 1 == temp.k) abovelevel[i] = Klevel_buffer.size() - 1; // pre store above level
//			if (rank_q[i] == temp.k) belowlevel[i] = Klevel_buffer.size() - 1; // pre stroe below level
//		}
//	}
//
//	// Insert all regions generated from q
//	vector<int> aboveNodes, belowNodes, klevel;
//
//	for (int i = 0; i < Rq.size(); i++) {
//
//		aboveNodes.clear();
//		klevel = Klevel_buffer[abovelevel[i]].Nodes;
//		for (int j = 0; j < klevel.size(); j++) {
//			if (Halfspace::is_Feasible(G[klevel[j]].r, Rq[i], dim)) aboveNodes.push_back(klevel[j]);
//		}
//
//		belowNodes.clear();
//		klevel = Klevel_buffer[belowlevel[i]].Nodes;
//		for (int j = 0; j < klevel.size(); j++) {
//			if (Halfspace::is_Feasible(G[klevel[j]].r, Rq[i], dim)) belowNodes.push_back(klevel[j]);
//		}
//
//		int NewNode = CreateNode_full(q.ID, Rq[i]);
//		Update_NewNode_full(NewNode, aboveNodes, belowNodes);
//	}
//
//	AllObj.push_back(q);
//}
//void KLEVEL::Update_SplitNode_full(int cur, Halfspace & r_LE, Halfspace & r_GE) {
//	G[cur].r = r_LE;
//	int NewNode = CreateNode_full(G[cur].pID, r_GE);
//
//	for (set<int>::iterator it = G[cur].minusEdge.begin(); it != G[cur].minusEdge.end();) {
//		if (Halfspace::is_Feasible(G[*it].r, G[NewNode].r, dim)) { // connect it and NewNode
//			G[NewNode].minusEdge.insert(*it);
//			G[*it].plusEdge.insert(NewNode);
//		}
//		if (!Halfspace::is_Feasible(G[*it].r, G[cur].r, dim)) { // disconnect it and cur
//			G[*it].plusEdge.erase(G[*it].plusEdge.find(cur));
//			it = G[cur].minusEdge.erase(it);
//		}
//		else it++;
//	}
//
//	for (set<int>::iterator it = G[cur].plusEdge.begin(); it != G[cur].plusEdge.end();) {
//		if (Halfspace::is_Feasible(G[*it].r, G[NewNode].r, dim)) {
//			G[NewNode].plusEdge.insert(*it);
//			G[*it].minusEdge.insert(NewNode);
//		}
//		if (!Halfspace::is_Feasible(G[*it].r, G[cur].r, dim)) {
//			G[*it].minusEdge.erase(G[*it].minusEdge.find(cur));
//			it = G[cur].plusEdge.erase(it);
//		}
//		else it++;
//	}
//}
//void KLEVEL::Update_NewNode_full(int cur, vector<int> & aboveNodes, vector<int> & belowNodes) {
//
//	for (int i = 0; i < aboveNodes.size(); i++) {
//		G[cur].minusEdge.insert(aboveNodes[i]);
//		G[aboveNodes[i]].plusEdge.insert(cur);
//	}
//	for (int i = 0; i < belowNodes.size(); i++) {
//		G[cur].plusEdge.insert(belowNodes[i]);
//		G[belowNodes[i]].minusEdge.insert(cur);
//	}
//
//	//Halfspace r_cur, r_ori;
//	// ??? delete some plusEdges of aboveNodes[i]
//	//for (auto i = aboveNodes.begin(); i != aboveNodes.end(); i++) {
//	//	for (auto j = G[*i].plusEdge.begin(); j != G[*i].plusEdge.end();) {
//	//		Halfspace::Intersect(r_cur, G[*i].r, G[cur].r);
//	//		Halfspace::Intersect(r_ori, G[*i].r, G[*j].r);
//	//		if (Halfspace::Subset(r_cur, r_ori)) {
//	//			G[*j].minusEdge.erase(*i);
//	//			G[*i].plusEdge.erase(*j);
//	//		}
//	//		else j++;
//	//	}
//	//}
//	//// ??? delete some minusEdge of belowNodes[i], is it useless because we already delete plusEdges of aboveNodes[i]
//	//for (auto i = belowNodes.begin(); i != belowNodes.end(); i++) {
//	//	for (auto j = G[*i].minusEdge.begin(); j != G[*i].minusEdge.end();) {
//	//		Halfspace::Intersect(r_cur, G[*i].r, G[cur].r);
//	//		Halfspace::Intersect(r_ori, G[*i].r, G[*j].r);
//	//		if (Halfspace::Subset(r_cur, r_ori)) {
//	//			G[*j].plusEdge.erase(*i);
//	//			G[*i].minusEdge.erase(*j);
//	//		}
//	//		else j++;
//	//	}
//	//}
//
//	// ??? this trick is correct or not?
//
//	for (auto i = aboveNodes.begin(); i != aboveNodes.end(); i++) {
//		for (auto j = belowNodes.begin(); j != belowNodes.end(); j++) {
//			if (G[*i].plusEdge.find(*j) != G[*i].plusEdge.end()) {
//				G[*i].plusEdge.erase(*j);
//				G[*j].minusEdge.erase(*i);
//			}
//		}
//	}
//
//}
//
//void KLEVEL::FindKlevel_full(vector<int> & Klevel, int k) {
//	Klevel.clear();
//	set<int> visit_node; visit_node.clear();
//	int cur = FindFirstNode_in_Klevel_full(k);
//	Klevel.push_back(cur); visit_node.insert(cur);
//
//	int cnt = 0;
//	while (cnt < Klevel.size()) {
//		cur = Klevel[cnt];
//		// minus then plus, connect based on above nodes
//		for (auto above_level = G[cur].minusEdge.begin(); above_level != G[cur].minusEdge.end(); above_level++) {
//			if (visit_node.find(*above_level) != visit_node.end()) continue;
//			visit_node.insert(*above_level);
//			for (auto this_again = G[*above_level].plusEdge.begin(); this_again != G[*above_level].plusEdge.end(); this_again++) {
//				if (visit_node.find(*this_again) != visit_node.end()) continue;
//				visit_node.insert(*this_again);
//				Klevel.push_back(*this_again);
//			}
//		}
//		// plus then minus, connect based on below nodes
//		for (auto below_level = G[cur].plusEdge.begin(); below_level != G[cur].plusEdge.end(); below_level++) {
//			if (visit_node.find(*below_level) != visit_node.end()) continue;
//			visit_node.insert(*below_level);
//			for (auto this_again = G[*below_level].minusEdge.begin(); this_again != G[*below_level].minusEdge.end(); this_again++) {
//				if (visit_node.find(*this_again) != visit_node.end()) continue;
//				visit_node.insert(*this_again);
//				Klevel.push_back(*this_again);
//			}
//		}
//		cnt++;
//	}
//	return;
//}
//
//int KLEVEL::FindFirstNode_in_Klevel_full(int k) {
//	int cur = top;
//	for (int i = 0; i < k; i++) {
//		cur = *G[cur].plusEdge.begin();
//		if (cur == bot) return -1;
//	}
//	return cur;
//}
//
//void KLEVEL::FindVoronoiWithRank_full(Obj & q, vector<Halfspace> & Rq, vector<int> & rank_q, vector<int> & HP_q) {
//	Halfspace origin_r;
//	Rq.clear(); rank_q.clear();
//	Rq.push_back(origin_r); rank_q.push_back(1);
//
//	for (int i = 0; i < AllObj.size(); i++) {
//		int tmpHP = HP_q[i];
//		int Rq_size = Rq.size();
//		for (int j = 0; j < Rq_size; j++) {
//			// rank(o_q) <= rank(o_i) // ???
//			if (Halfspace::is_Feasible(Rq[j], dim, tmpHP, false)) rank_q[j]++;
//			if (Halfspace::is_intersected(Rq[j], dim, tmpHP)) {
//				Halfspace r_LE, r_GE;
//				Halfspace::Split(Rq[j], tmpHP, r_LE, r_GE);
//				Rq[j] = r_LE;
//				Rq.push_back(r_GE); rank_q.push_back(rank_q[j] - 1); // rank(o_q)>=rank(o_i) in r_GE
//			}
//		}
//	}
//
//}
//float KLEVEL::CompRadius(Obj & o) {
//	float sum = 0;
//	for (int i = 0; i < dim; i++) {
//		sum = sum + o.w[i] * o.w[i];
//	}
//	return sqrt(sum);
//}
//void KLEVEL::Insert(Obj& q, int min_level) {
//
//	vector<int> HP_q; HP_q.clear();
//	for (int i = 0; i < AllObj.size(); i++) {
//		HP_q.push_back(Halfspace::ComputeHP(q, AllObj[i], dim));
//	}
//
//	// Splitting Method 1:
//	// Split every region that intersects with HP(q,i) in 1-level to Max_k-level
//	if (Pruning1) {
//		bool flag = false;
//		for (int i = Max_k; i >= 1; i--) {
//			if (min_level > i) continue; // q cannot be in <=i-level
//			bool flag_this_level = false;
//			int level_size = Klevel_buffer[i].Nodes.size();
//			for (int j = 0; j < level_size; j++) {
//				int cur = Klevel_buffer[i].Nodes[j];
//				int HP = HP_q[G[cur].pID];
//				// {HP,true}, q>=o_i
//				// {HP,false}, q<=o_i
//				if (XdominateY(AllObj[G[cur].pID], q, dim)) continue; 
//				if (!Halfspace::is_Feasible(G[cur].r, dim, HP, true)) continue; // q<=o_i in this region
//				if (!Halfspace::is_Feasible(G[cur].r, dim, HP, false)) continue; // HP doesn't intersect with G[cur].r
//				flag_this_level = true;
//				// discuss!!!
//				if (i < Max_k) {
//					int NewNode = Update_SplitNode(cur, HP);
//					Klevel_buffer[i].Nodes.push_back(NewNode);
//				}
//				else Update_SplitNode_lastLevel(cur, HP);
//			}
//			flag = flag || flag_this_level;
//			if ((!flag_this_level) && (AllObj.size() >= Max_k)) break; // q less than all o_i in this_level
//		}
//		if ((!flag) && (AllObj.size() >= Max_k)) {
//			AllObj.push_back(q);
//			return; // q less than all o_i in <=k-level
//		}
//	}
//	else {
//		for (int i = 1; i <= Max_k; i++) {
//			int level_size = Klevel_buffer[i].Nodes.size();
//			for (int j = 0; j < level_size; j++) {
//				int cur = Klevel_buffer[i].Nodes[j];
//				int HP = HP_q[G[cur].pID];
//				if (Halfspace::is_intersected(G[cur].r, dim, HP)) {
//					// discuss!!!
//					if (i < Max_k) {
//						int NewNode = Update_SplitNode(cur, HP);
//						Klevel_buffer[i].Nodes.push_back(NewNode);
//					}
//					else Update_SplitNode_lastLevel(cur, HP);
//
//					//
//					/*int NewNode = Update_SplitNode(cur, HP);
//					Klevel_buffer[i].Nodes.push_back(NewNode);*/
//				}
//			}
//		}
//	}
//	
//
//	// PreCompute Klevel before insertion
//	vector<Halfspace> Rq;
//	vector<int> rank_q;
//	FindVoronoiWithRank(q, Rq, rank_q, HP_q);
//
//	// Insert all regions generated from q
//	vector<int> aboveNodes, belowNodes, klevel;
//
//	for (int i = 0; i < Rq.size(); i++) {
//		if (rank_q[i] > Max_k) continue;
//		if (rank_q[i] > (AllObj.size() + 1)) cout << "rank error" << endl;
//
//		aboveNodes.clear();
//		klevel = Klevel_buffer[rank_q[i] - 1].Nodes;
//		for (int j = 0; j < klevel.size(); j++) {
//			if (Halfspace::is_Feasible(G[klevel[j]].r, Rq[i], dim)) aboveNodes.push_back(klevel[j]);
//		}
//
//		belowNodes.clear();
//		if (AllObj.size() >= rank_q[i]) {
//			// discuss!!!
//			if (rank_q[i] < Max_k) {
//				klevel = Klevel_buffer[rank_q[i]].Nodes;	
//				for (int j = 0; j < klevel.size(); j++) {
//					if (Halfspace::is_Feasible(G[klevel[j]].r, Rq[i], dim)) belowNodes.push_back(klevel[j]);
//				}
//			}
//			else belowNodes.push_back(bot);
//		}
//		else belowNodes.push_back(bot);
//		
//		int NewNode = CreateNode(q.ID, Rq[i]);
//		Update_NewNode(NewNode, aboveNodes, belowNodes);
//	}
//
//	AllObj.push_back(q);
//
//	if (AllObj.size() >= Max_k) FindKlevel(Max_k);
//	else FindKlevel(AllObj.size());
//
//	// delete Max_k+1 level
//	if (AllObj.size() >= Max_k) {
//		G[bot].minusEdge.clear();
//		for (int i = 0; i < Klevel_buffer[Max_k].Nodes.size(); i++) {
//			int cur = Klevel_buffer[Max_k].Nodes[i];
//			G[bot].minusEdge.insert(cur);
//			G[cur].plusEdge.clear();
//			G[cur].plusEdge.insert(bot);
//		}
//	}
//
//	//CheckLinkERROR();
//
//	Klevel_Obj.clear();
//	for (int i = 1; i <= Max_k; i++) {
//		cout << i << "-level size: " << Klevel_buffer[i].Nodes.size() << endl;
//		for (int j = 0; j < Klevel_buffer[i].Nodes.size(); j++) {
//			Klevel_Obj.insert(G[Klevel_buffer[i].Nodes[j]].pID);
//		}
//	}
//	cout << "Objects in <=Klevel: " << Klevel_Obj.size() << endl;
//}
//void KLEVEL::testing(vector<Obj>& All_records, vector<long int>& min_level, int dim) {
//	AllHP.clear();
//	for (int i = 0; i < All_records.size(); i++) {
//		if (min_level[i] > 1) continue;
//		Halfspace r;
//		for (int j = 0; j < All_records.size(); j++) {
//			if (i == j) continue;
//			int HP = Halfspace::ComputeHP(All_records[i], All_records[j], dim);
//			r.HPs.insert({ HP,true });
//		}
//
//		for (int j = 0; j < All_records.size(); j++) {
//			if (i == j) continue;
//			int min_k = 1;
//			for (int t = 0; t < All_records.size(); t++) {
//				if ((t == i) || (t == j)) continue;
//				int HP = Halfspace::ComputeHP(All_records[j], All_records[t], dim);
//				// j is dominated by t in this region
//				if (!Halfspace::is_Feasible(r, dim, HP, true)) min_k++;
//			}
//			if (min_k > min_level[j]) min_level[j] = min_k;
//		}
//	}
//}
//void KLEVEL::SpaceSplit(vector<Obj>& All_records, int split_size, int st[], float offset, int dim, int mk) {
//	SubSpace space;
//	SingleSubSpaceBuilding(All_records, space, dim, offset, st, mk);
//	subspace_num++;
//	ave_size = ave_size + space.candidate.size();
//	//cout << space.candidate.size() << endl;
//	//system("pause");
//	for (auto it = space.candidate.begin(); it != space.candidate.end(); it++) All_candidates.insert(*it);
//	int sum = 0;
//	for (int i = 0; i < dim - 1; i++) sum = sum + st[i];
//	if (sum >= split_size - 1) return;
//	for (int i = 0; i < dim - 1; i++) {
//		if (st[i] == split_size - 1) continue;
//		st[i]++;
//		SpaceSplit(All_records, split_size, st, offset, dim, mk);
//		st[i]--;
//	}
//}
//
//void KLEVEL::AllSubSpaceBuilding(vector<Obj>& All_records, int split_size, int dim, int mk) {
//	float offset = 1.0 / (float)split_size;
//	int st[Max_Dimension];
//	for (int d = 0; d < dim - 1; d++) st[d] = 0;
//	ave_size = 0; subspace_num = 0;
//	SpaceSplit(All_records, split_size, st, offset, dim, mk);
//	cout << All_candidates.size() << endl;
//	cout << subspace_num << endl;
//	cout << (float) ave_size/subspace_num << endl;
//	return;
//}
//
//bool tempfunc(pair<float, int> A, pair<float, int> B) {
//	return (A.first < B.first);
//}
//
//void KLEVEL::SingleSubSpaceBuilding(vector<Obj>& All_records, SubSpace& space, int dim, float offset, int st[], int mk) {
//	space.Max_k = mk;
//	for (int i = 0; i < dim - 1; i++) {
//		space.L[i] = st[i] * offset;
//		space.U[i] = space.L[i] + offset;
//	}
//	space.candidate.clear();
//	vector<float> ls, us; 
//	ls.clear(); us.clear();
//	for (int i = 0; i < All_records.size(); i++) {
//		float vl[Max_Dimension], vu[Max_Dimension];
//		Obj cur = All_records[i];
//		vector<pair<float, int>> temp;
//		float resvl = 1.0, resvu = 1.0;
//		for (int d = 0; d < dim - 1; d++) {
//			vl[d] = space.L[d]; resvl = resvl - space.L[d];
//			vu[d] = space.L[d]; resvu = resvu - space.L[d];
//			temp.push_back({ cur.w[d],d });
//		}
//		temp.push_back({ cur.w[dim - 1],dim - 1 });
//		sort(temp.begin(), temp.end(), tempfunc);
//		// compute vl
//		for (int td = dim-1; td <= 0; td--) {
//			int d = temp[td].second;
//			if (d != dim - 1) {
//				if (resvl > offset) {
//					resvl = resvl - offset;
//					vl[d] = vl[d] + offset;
//				}
//				else {
//					resvl = 0;
//					vl[d] = vl[d] + resvl;
//					break;
//				}
//			}
//			else {
//				vl[d] = resvl;
//				break;
//			}
//		}
//		// compute vu
//		for (int td = 0; td < dim-1; td++) {
//			int d = temp[td].second;
//			if (d != dim - 1) {
//				if (resvu > offset) {
//					resvu = resvu - offset;
//					vu[d] = vu[d] + offset;
//				}
//				else {
//					resvu = 0;
//					vu[d] = vu[d] + resvu;
//					break;
//				}
//			}
//			else {
//				vu[d] = resvu;
//				break;
//			}
//		}
//		
//		ls.push_back(0.0); us.push_back(0.0);
//		for (int d = 0; d < dim; d++) {
//			ls[i] = ls[i] + cur.w[d] * vl[d];
//			us[i] = us[i] + cur.w[d] * vu[d];
//		}
//
//	}
//	sort(ls.begin(), ls.end());
//	for (int i = 0; i < All_records.size(); i++) {
//		int dominate = 0;
//		for (auto j=All_records.size()-1; j >=0; j--) {
//			//cout << ls[j] << ' ' << us[i] << endl;
//			if (ls[j] > us[i]) dominate++;
//			if (dominate >= space.Max_k) break;
//			if (ls[j] < us[i]) break;
//		}
//		//cout << dominate << endl;
//		if (dominate < space.Max_k) space.candidate.insert(i);
//		//cout << space.candidate.size() << endl;
//	}
//}
//
//bool KLEVEL::CheckRepeatTopk(region& r1, region& r2, int kth_obj) {
//	if (r1.objID != kth_obj) return false;
//	if (r1.records.find(kth_obj) == r1.records.end()) return false;
//	for (auto it = r2.records.begin(); it != r2.records.end(); it++) {
//		if (r1.records.find(*it) == r1.records.end()) return false;
//	}
//	return true;
//}
//
//void KLEVEL::Bulk_Loading(vector<Obj>& All_records, vector<long int>& min_level) {
//	// Precompute hypero=plane
//	AllHP.clear();
//	for (int i = 0; i < All_records.size(); i++) {
//		for (int j = 0; j < All_records.size(); j++) {
//			int tmp = Halfspace::ComputeHP(All_records[i], All_records[j], dim);
//		}
//	}
//	
//	// Init Level-0
//	region r_ori; Halfspace hs_ori;
//	r_ori.k = 0; r_ori.objID = -1; r_ori.r = hs_ori; r_ori.records.clear(); r_ori.candidates.clear();
//	for (int i = 0; i < All_records.size(); i++) {
//		if (min_level[i] > Max_k) continue;
//		r_ori.candidates.push_back({ i,min_level[i] });
//	}
//
//	vector<region> L_zero; L_zero.push_back(r_ori);
//	L.push_back(L_zero);
//	
//	set<int> utk_results; utk_results.clear();
//
//	// Compute Klevel from level-1 to level-Max_k
//	for (int i = 1; i <= Max_k; i++) {
//		// Init
//		vector<region> tmp_level; tmp_level.clear();
//		L.push_back(tmp_level);
//
//		for (int j = 0; j < L[i - 1].size(); j++) {
//			//bool flag = true;
//			//for (int q = 0; q < L[i - 1][j].candidates.size(); q++) {
//			//	if (L[i - 1][j].candidates[q].second != 1) flag = false;
//			//}
//			//if (flag) continue;
//			//vector<int> candidates; candidates.clear();
//			//for (int q = 0; q < All_records.size(); q++)
//			//	if ((min_level[q] <= i) && (L[i - 1][j].records.find(q) == L[i - 1][j].records.end())) candidates.push_back(q);
//
//			for (int q = 0; q < L[i - 1][j].candidates.size(); q++) {
//				if (L[i - 1][j].candidates[q].second > 1) continue;
//				int cur = L[i - 1][j].candidates[q].first;
//				bool flag = false;
//				for (int ii = 0; ii < L[i].size(); ii++) {
//					if (CheckRepeatTopk(L[i][ii], L[i - 1][j], cur)) flag = true;
//					if (flag) break;
//				}
//				if (flag) continue;
//
//				Halfspace r = L[i - 1][j].r;
//				for (int t = 0; t < L[i - 1][j].candidates.size(); t++) {
//					if (q == t) continue;
//					int HP = cur * All_records.size() + L[i - 1][j].candidates[t].first;
//					r.HPs.insert({ HP,true });// the halfspace S(q)>S(t)
//				}
//				if (Halfspace::is_Feasible(r, dim)) {
//					if (utk_results.find(cur) == utk_results.end()) utk_results.insert(cur);
//
//					region tmp;
//					tmp.k = i; tmp.objID = cur; 
//					tmp.records = L[i - 1][j].records;
//					
//					tmp.r.HPs.clear();
//					for (auto it = tmp.records.begin(); it != tmp.records.end(); it++) {
//						int HP = cur * All_records.size() + *it;
//						tmp.r.HPs.insert({ HP,false }); // the halfspace S(cur)<S(*it)
//					}
//					for (int t = 0; t < L[i - 1][j].candidates.size(); t++) {
//						if (q == t) continue;
//						int HP = cur * All_records.size() + L[i - 1][j].candidates[t].first;
//						tmp.r.HPs.insert({ HP,true });// the halfspace S(cur)>S(t)
//					}
//					//tmp.r = r; 
//					
//					tmp.records.insert(tmp.objID);
//
//					if (i != Max_k) { // Compute candidates for next level
//						tmp.candidates.clear();
//						for (int q_next = 0; q_next < L[i - 1][j].candidates.size(); q_next++) {
//							if (q == q_next) continue;
//							
//							//int cand = L[i - 1][j].candidates[q_next].first;
//							//tmp.candidates.push_back({ cand,min_level[cand]-i });
//
//							int rdominated = 0;
//							for (int t_next = 0; t_next < L[i - 1][j].candidates.size(); t_next++) {
//								if (q == t_next) continue;
//								if (q_next == t_next) continue;
//								if (rdominated >= Max_k - i) continue;
//								if (L[i - 1][j].candidates[t_next].second > L[i - 1][j].candidates[q_next].second) continue;
//								int HP = L[i - 1][j].candidates[q_next].first * All_records.size() + L[i - 1][j].candidates[t_next].first;
//								if (!Halfspace::is_Feasible(r, dim, HP, true)) { // S(q_next) > S(t_next) is impossible
//									rdominated++;
//								}
//							}
//							if (rdominated < (Max_k - i)) {
//								tmp.candidates.push_back({ L[i - 1][j].candidates[q_next].first, rdominated + 1 });
//							}
//						}
//					}
//
//					L[i].push_back(tmp);
//				}
//			}
//		}
//
//		cout << "The region size of LEVEL " << i << ": " << L[i].size() << endl;
//		cout << "The utk size of LEVEL " << i << ": " << utk_results.size() << endl;
//	}
//	AllHP.clear();
//}
//
//float KLEVEL::orderScore(vector<float>& pivot, float entry[]) {
//	float ret = 0;
//	float weight = 0;
//	for (int i = 0; i < dim - 1; i++)
//	{
//		weight += pivot[i];
//		ret += pivot[i] * entry[i];
//	}
//	ret += (1 - weight) * entry[dim - 1];
//	return ret;
//}
//float KLEVEL::orderScore(float pivot[], float entry[]) {
//	float ret = 0;
//	float weight = 0;
//	for (int i = 0; i < dim - 1; i++)
//	{
//		weight += pivot[i];
//		ret += pivot[i] * entry[i];
//	}
//	ret += (1 - weight) * entry[dim - 1];
//	return ret;
//}
//
//bool KLEVEL::isRdominated(float r[], float focal[], vector<float>& entry, bool& fDe) {
//	// generate HP;
//	vector<float> tmpHS;
//	float entry_d = entry[dim - 1];
//	float focal_d = focal[dim - 1];
//	for (int d = 0; d < dim - 1; d++)
//	{
//		tmpHS.push_back((focal[d] - focal_d) - (entry[d] - entry_d));
//	}
//	tmpHS.push_back(entry_d - focal_d);
//
//
//	// verify each vertic
//	int posCount = 0;
//	int negCount = 0;
//	int totVertics = pow(2, dim - 1);
//	for (int i = 0; i < totVertics; i++)
//	{
//		//stringstream ss;
//		string tmpS = bitset<MAXDIMEN>(i).to_string();
//		tmpS = tmpS.substr(tmpS.size() - (dim - 1), dim - 1);
//
//		float sum = 0;
//		for (int si = 0; si < tmpS.size(); si++)
//		{
//			if (tmpS[si] == '0')
//				sum += r[si * 2] * tmpHS[si];
//			else if (tmpS[si] == '1')
//				sum += r[si * 2 + 1] * tmpHS[si];
//			else
//				cout << "bug here!!!" << endl;
//		}
//
//		if (sum > tmpHS[dim - 1]) return false;
//		/*else
//			posCount++;*/
//	}
//
//	return true;
//
//	// dominationship
//	/*if (negCount == totVertics)
//	{
//		fDe = true;
//		return true;
//	}
//	else if (posCount == totVertics)
//	{
//		fDe = false;
//		return true;
//	}
//	else
//		return false;*/
//}
//
//bool KLEVEL::countRegionDominator(int k, float pt[], vector<int>& rskyband_set, float* PG[], float r[], unordered_set<int>& dominators) {
//	vector<float> record(dim, 0);
//	dominators.clear();
//	bool fDe;
//	int count = 0;
//	/* // should comment this for large k
//	if (rskyband.size() < k)
//	{
//		return true;
//	}
//	*/
//	for (int i = 0; i < rskyband_set.size(); i++)
//	{
//		for (int di = 0; di < dim; di++)
//			record[di] = (PG[rskyband_set[i]][di] + PG[rskyband_set[i]][di + dim]) / 2;
//		fDe = false;
//		if (isRdominated(r, pt, record, fDe))
//		{
//			if (fDe == false)
//			{
//				count++;
//				dominators.insert(rskyband_set[i]);
//			}
//		}
//	}
//	if (count < k)
//	{
//		return true;
//	}
//	else
//	{
//		return false;
//	}
//}
//void KLEVEL::pivotRegion(float r[], vector<float>& pivot) {
//	pivot.clear();
//	for (int d = 0; d < dim; d++)
//		pivot.push_back((r[2 * d] + r[2 * d + 1]) / 2);
//	return;
//}
//void KLEVEL::rskyband(float** PG, float r[], vector<Obj>& All_records, set<int>& candidates, int a_k, vector<int>& rskyband_set) {
//	rskyband_set.clear();
//	if (a_k == 0) {
//		for (auto it = candidates.begin(); it != candidates.end(); it++) {
//			rskyband_set.push_back(*it);
//		}
//		return;
//	}
//	
//	// build rtree
//	RtreeNodeEntry** tmp = new RtreeNodeEntry * [candidates.size()];
//	int cnt = 0;
//	for (auto it = candidates.begin(); it != candidates.end(); it++) {
//		Hypercube hc = MBRs[*it];
//		tmp[cnt] = new RtreeNodeEntry(*it, hc);
//		cnt++;
//	}
//
//	string indexfile = "E:\\klevel\\Project\\klevel_index\\tmp_rtree.idx";
//	const int maxChild = (PAGESIZE - RtreeNode::size()) / RtreeNodeEntry::size(dim);
//	FileMemory mem(PAGESIZE, indexfile, RtreeNodeEntry::fromMem, true);
//	Rtree * tmp_rtree = TGS::bulkload(mem, dim, maxChild, maxChild, (int)maxChild * 0.3, (int)maxChild * 0.3, tmp, cnt, false);
//
//
//	vector<float> pivot;
//	pivotRegion(r, pivot);
//
//	RtreeNode* node;
//	priority_queue<pair<float, int>> heap;
//	int NegPageid;
//
//	float pt[MAXDIMEN];
//	float maxscore;
//	int pageID;
//	float tmpScore;
//	unordered_set<int> dominators;
//
//	heap.push(make_pair(INFINITY, tmp_rtree->m_memory.m_rootPageID));
//
//	while (!heap.empty())
//	{
//		tmpScore = heap.top().first;
//		pageID = heap.top().second;
//		heap.pop();
//
//		if (pageID >= MAXPAGEID)
//		{
//			//cout << pageID - MAXPAGEID << endl;
//			for (int j = 0; j < dim; j++)
//				pt[j] = (PG[pageID - MAXPAGEID][j] + PG[pageID - MAXPAGEID][j + dim]) / 2;
//			if (countRegionDominator(a_k, pt, rskyband_set, PG, r, dominators))
//			{
//				rskyband_set.push_back(pageID - MAXPAGEID);
//				//daGraph[pageID - MAXPAGEID] = dominators;
//			}
//		}
//		else
//		{
//			node = tmp_rtree->m_memory.loadPage(pageID);
//			//node = ramTree[pageID];
//			if (node->isLeaf())
//			{
//				for (int i = 0; i < node->m_usedspace; i++)
//				{
//					for (int j = 0; j < dim; j++)
//						pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
//
//					if (countRegionDominator(a_k, pt, rskyband_set, PG, r, dominators))
//					{
//						maxscore = orderScore(pivot, pt);
//						heap.push(make_pair(maxscore, node->m_entry[i]->m_id + MAXPAGEID));
//					}
//				}
//			}
//			else
//			{
//				for (int i = 0; i < node->m_usedspace; i++)
//				{
//					for (int j = 0; j < dim; j++)
//					{
//						pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
//					}
//					if (countRegionDominator(a_k, pt, rskyband_set, PG, r, dominators))
//					{
//						maxscore = orderScore(pivot, pt);
//						heap.push(make_pair(maxscore, node->m_entry[i]->m_id));
//					}
//				}
//			}
//		}
//	}
//	delete(tmp);
//}
//void KLEVEL::rskyband(Rtree* tmp_rtree, float** PG, float r[], vector<Obj>& All_records, int a_k, vector<int>& rskyband_set) {
//	rskyband_set.clear();
//
//	vector<float> pivot;
//	pivotRegion(r, pivot);
//
//	RtreeNode * node;
//	priority_queue<pair<float, int>> heap;
//	int NegPageid;
//
//	float pt[MAXDIMEN];
//	float maxscore;
//	int pageID;
//	float tmpScore;
//	unordered_set<int> dominators;
//
//	heap.push(make_pair(INFINITY, tmp_rtree->m_memory.m_rootPageID));
//
//	while (!heap.empty())
//	{
//		tmpScore = heap.top().first;
//		pageID = heap.top().second;
//		heap.pop();
//
//		if (pageID >= MAXPAGEID)
//		{
//			//cout << pageID - MAXPAGEID << endl;
//			for (int j = 0; j < dim; j++)
//				pt[j] = (PG[pageID - MAXPAGEID][j] + PG[pageID - MAXPAGEID][j + dim]) / 2;
//			if (countRegionDominator(a_k, pt, rskyband_set, PG, r, dominators))
//			{
//				rskyband_set.push_back(pageID - MAXPAGEID);
//				//daGraph[pageID - MAXPAGEID] = dominators;
//			}
//		}
//		else
//		{
//			node = tmp_rtree->m_memory.loadPage(pageID);
//			//node = ramTree[pageID];
//			if (node->isLeaf())
//			{
//				for (int i = 0; i < node->m_usedspace; i++)
//				{
//					for (int j = 0; j < dim; j++)
//						pt[j] = node->m_entry[i]->m_hc.getLower()[j] + SIDELEN;
//
//					if (countRegionDominator(a_k, pt, rskyband_set, PG, r, dominators))
//					{
//						maxscore = orderScore(pivot, pt);
//						heap.push(make_pair(maxscore, node->m_entry[i]->m_id + MAXPAGEID));
//					}
//				}
//			}
//			else
//			{
//				for (int i = 0; i < node->m_usedspace; i++)
//				{
//					for (int j = 0; j < dim; j++)
//					{
//						pt[j] = node->m_entry[i]->m_hc.getUpper()[j];
//					}
//					if (countRegionDominator(a_k, pt, rskyband_set, PG, r, dominators))
//					{
//						maxscore = orderScore(pivot, pt);
//						heap.push(make_pair(maxscore, node->m_entry[i]->m_id));
//					}
//				}
//			}
//		}
//	}
//}
//
//bool KLEVEL::countRegionDominator(float r[], vector<Obj>& All_records, Obj& cur_obj, vector<int>& rskyband_set, int k, unordered_set<int>& dominators) {
//	float pt[MAXDIMEN];
//	for (int i = 0; i < dim; i++) pt[i] = cur_obj.w[i];
//	
//	vector<float> record(dim, 0);
//	//dominators.clear();
//	bool fDe;
//	int count = 0;
//	for (int i = 0; i < rskyband_set.size(); i++)
//	{
//		for (int di = 0; di < dim; di++)
//			record[di] = All_records[rskyband_set[i]].w[di];
//
//		fDe = false;
//		if (isRdominated(r, pt, record, fDe))
//		{
//			if (fDe == false)
//			{
//				count++;
//				//dominators.insert(rskyband_set[i]);
//			}
//		}
//		if (count >= k) break;
//	}
//	if (count < k)
//	{
//		return true;
//	}
//	else
//	{
//		return false;
//	}
//}
//void KLEVEL::rskyband(float r[], vector<Obj>& All_records, set<int>& candidates, int a_k, vector<int>& rskyband_set) {
//	rskyband_set.clear();
//
//	float maxscore;
//	vector<float> pivot;
//	pivotRegion(r, pivot);
//	
//	priority_queue<pair<float, int>> heap;
//
//	for (auto it = candidates.begin(); it != candidates.end(); it++) {
//		maxscore = orderScore(pivot, All_records[*it].w);
//		heap.push(make_pair(maxscore, *it));
//	}
//
//	unordered_set<int> dominators;
//
//	float dr[Max_Dimension * 2];
//	while (!heap.empty()) {
//		int cur = heap.top().second;
//		heap.pop();
//
//		if (countRegionDominator(r, All_records, All_records[cur], rskyband_set, a_k, dominators)) {
//			//if (a_k == 1) {
//			//	// divide multiple region
//			//	bool flag = true;
//			//	int region_divide = pow(2, dim - 1);
//			//	for (int i = 0; i < region_divide; i++)
//			//	{
//			//		//stringstream ss;
//			//		string tmpS = bitset<MAXDIMEN>(i).to_string();
//			//		tmpS = tmpS.substr(tmpS.size() - (dim - 1), dim - 1);
//
//			//		for (int j = 0; j < dim - 1; j++) {
//			//			if (tmpS[j] == '0') {
//			//				dr[2 * j] = r[2 * j];
//			//				dr[2 * j + 1] = (r[2 * j] + r[2 * j + 1]) / 2.0;
//			//			}
//			//			else {
//			//				dr[2 * j] = (r[2 * j] + r[2 * j + 1]) / 2.0;
//			//				dr[2 * j + 1] = r[2 * j + 1];
//			//			}
//			//		}
//
//			//		// cur in rskyband if it is in top-k within any divide space
//			//		if (countRegionDominator(dr, All_records, All_records[cur], rskyband_set, a_k, dominators)) {
//			//			flag = false;
//			//			break;
//			//		}
//			//	}
//			//	if (!flag) rskyband_set.push_back(cur);
//			//}
//			//else 
//			rskyband_set.push_back(cur);
//		}
//	}
//}
//bool KLEVEL::isRdominated(vector<vector<float>>& vertices, float focal[], float entry[]) {
//	// generate HP;
//	vector<float> tmpHS;
//	float entry_d = entry[dim - 1];
//	float focal_d = focal[dim - 1];
//	for (int d = 0; d < dim - 1; d++)
//	{
//		tmpHS.push_back((focal[d] - focal_d) - (entry[d] - entry_d));
//	}
//	tmpHS.push_back(entry_d - focal_d);
//
//
//	// verify each vertic
//	for (int i = 0; i < vertices.size(); i++)
//	{
//		float sum = 0;
//		for (int d = 0; d < dim-1; d++)
//		{
//			sum = sum + vertices[i][d] * tmpHS[d];
//		}
//
//		if (sum > tmpHS[dim - 1]) return false;
//	}
//
//	return true;
//
//}
//
//
//void KLEVEL::rskyband(Halfspace& r, vector<Obj>& All_records, set<int>& candidates, int a_k, vector<int>& rskyband_set, vector<int>& dominators) {
//	rskyband_set.clear();
//
//	priority_queue<pair<float, int>> heap;
//
//	for (auto it = candidates.begin(); it != candidates.end(); it++) {
//		heap.push(make_pair(orderScore(r.innerPoint,All_records[*it].w), *it));
//	}
//
//	dominators.clear();
//	
//	float dr[Max_Dimension * 2];
//	while (!heap.empty()) {
//		int cur = heap.top().second;
//		heap.pop();
//
//		int cnt = 0;
//		int dominate = 0;
//		vector<int> tmp_cnt(rskyband_set.size(), 0);
//		for (int i = 0; i < rskyband_set.size(); i++) {
//			//if (isRdominated(r.vertices, All_records[cur].w, All_records[rskyband_set[i]].w)) cnt++;
//			int res = isRdominated(r.vertices, All_records[cur].w, All_records[rskyband_set[i]].w);
//			if (res == 0) {
//				//cur<rskyband_set[i]
//				cnt++;
//				tmp_cnt[i] = 1;
//			}
//			else if (res==1) {
//				// cur>rskyband_set[i]
//				dominate++;
//			}
//			else {
//				// do thinig;
//			}
//		}
//		if (cnt < a_k) {
//			for (int i = 0; i < rskyband_set.size(); i++) {
//				dominators[i] = dominators[i] + tmp_cnt[i];
//			}
//			rskyband_set.push_back(cur);
//			dominators.push_back(dominate);
//		}
//		/*for (auto it = rskyband_set.begin(); it != rskyband_set.end();it++) {
//			if (isRdominated(r.vertices, All_records[cur].w, All_records[*it].w)) cnt++;
//		}*/
//		//if (cnt<a_k) rskyband_set.push_back(cur);
//	}
//}
//bool sfunc(pair<int, float> x, pair<int, float> y) {
//	return x.second > y.second;
//}
//MBR version
//bool KLEVEL::MergeLevel(region_rd& r, vector<Obj>& All_records, set<int>& utk_results) {
//	vector<vector<pair<int, float>>> rank;
//	rank.clear();
//	vector<unordered_set<int>> topkset; 
//	topkset.clear();
//
//	int vertex = 1 << (dim - 1);
//	for (int i = 0; i < vertex; i++) {
//		float wi[Max_Dimension];
//		wi[dim - 1] = 1.0;
//		int tmp = i;
//		for (int d = 0; d < dim - 1; d++) {
//			if (tmp % 2 == 0) wi[d] = r.boundary[2 * d];
//			else wi[d] = r.boundary[2 * d + 1];
//			tmp = tmp >> 1;
//			wi[dim - 1] = wi[dim - 1] - wi[d];
//		}
//
//		vector<pair<int, float>> scores; scores.clear();
//		for (auto it = r.candidates.begin(); it != r.candidates.end(); it++) {
//			int cur = *it;
//			float tmpS = 0.0;
//			for (int d = 0; d < dim; d++) {
//				tmpS = tmpS + wi[d] * All_records[cur].w[d];
//			}
//			scores.push_back({ cur,tmpS });
//		}
//		sort(scores.begin(), scores.end(), sfunc);
//		rank.push_back(scores);
//
//		unordered_set<int> topk; topk.clear();
//		topkset.push_back(topk);
//	}
//	for (int l = r.level + 1; l <= Max_k; l++) {
//		for (int i = 0; i < vertex; i++) {
//			if (rank[i].size() == 0) cout << 123 << endl;
//			topkset[i].insert(rank[i][l - r.level-1].first);
//		}
//		bool flag = true;
//		for (int i = 1; i < vertex; i++) {
//			if (topkset[i] != topkset[0]) {
//				flag = false;
//				break;
//			}
//		}
//		if ((flag)&&(l-r.level>1)) {
//			region_rd NewRegion;
//			NewRegion = r;
//			if (l - r.level == 1) NewRegion.objID = *topkset[0].begin();
//			else NewRegion.objID = -1;
//			NewRegion.level = l;
//			for (auto it = topkset[0].begin(); it != topkset[0].end(); it++) {
//				NewRegion.candidates.erase(*it);
//				NewRegion.confirmed.insert(*it);
//			}
//			L_rd[l].push_back(NewRegion);
//			for (auto it = topkset[0].begin(); it != topkset[0].end(); it++) {
//				if (utk_results.find(*it) == utk_results.end()) utk_results.insert(*it);
//			}
//			return true;
//		}
//	}
//	return false;
//}
//bool KLEVEL::MergeLevel(region_rd& r, vector<Obj>& All_records, set<int>& utk_results) {
//	vector<vector<pair<int, float>>> rank;
//	rank.clear();
//	vector<unordered_set<int>> topkset;
//	topkset.clear();
//	if (r.r.vertices.size() == 0) return false;
//	for (int i = 0; i < r.r.vertices.size(); i++) {
//		float wi[Max_Dimension];
//		wi[dim - 1] = 1.0;
//		for (int d = 0; d < dim - 1; d++) {
//			wi[d] = r.r.vertices[i][d];
//			wi[dim - 1] = wi[dim - 1] - wi[d];
//		}
//
//		vector<pair<int, float>> scores; scores.clear();
//		for (auto it = r.candidates.begin(); it != r.candidates.end(); it++) {
//			int cur = *it;
//			float tmpS = 0.0;
//			for (int d = 0; d < dim; d++) {
//				tmpS = tmpS + wi[d] * All_records[cur].w[d];
//			}
//			scores.push_back({ cur,tmpS });
//		}
//		sort(scores.begin(), scores.end(), sfunc);
//		rank.push_back(scores);
//
//		unordered_set<int> topk; topk.clear();
//		topkset.push_back(topk);
//	}
//
//	for (int l = r.level + 1; l <= Max_k; l++) {
//		for (int i = 0; i < r.r.vertices.size(); i++) {
//			if (rank[i].size() == 0) cout << 123 << endl;
//			topkset[i].insert(rank[i][l - r.level - 1].first);
//		}
//		bool flag = true;
//		for (int i = 1; i < r.r.vertices.size(); i++) {
//			if (topkset[i] != topkset[0]) {
//				flag = false;
//				break;
//			}
//		}
//		if ((flag) && (l - r.level > 2)) {
//			region_rd NewRegion;
//			NewRegion = r;
//			if (l - r.level == 1) NewRegion.objID = *topkset[0].begin();
//			else NewRegion.objID = -1;
//			NewRegion.level = l;
//			for (auto it = topkset[0].begin(); it != topkset[0].end(); it++) {
//				if (NewRegion.candidates.find(*it) != NewRegion.candidates.end()) NewRegion.candidates.erase(*it);
//				NewRegion.confirmed.insert(*it);
//			}
//			L_rd[l].push_back(NewRegion);
//			for (auto it = topkset[0].begin(); it != topkset[0].end(); it++) {
//				if (utk_results.find(*it) == utk_results.end()) utk_results.insert(*it);
//			}
//			return true;
//		}
//	}
//	return false;
//}
//
//int KLEVEL::LevelGrouping(int k, vector<int>& candidates, vector<int>& dominators, region_rd& this_r, vector<Obj>& All_records) {
//	region_rd tmp;
//	tmp.level = Max_k; tmp.objID = -1; tmp.confirmed = this_r.confirmed; tmp.candidates.clear();
//	for (auto it = candidates.begin(); it != candidates.end(); it++) tmp.confirmed.insert(*it);
//	int cnt = 0;
//	string bitmask(k, 1);
//	bitmask.resize(candidates.size(), 0);
//	do {
//		bool flag = true;
//		for (int i = 0; i < candidates.size(); i++) { // candidates[i] cannot be last (size-k) records
//			if ((bitmask[i] == '0') && (dominators[i] >= candidates.size() - k)) flag = false;
//		}
//		if (!flag) continue;
//
//		for (int i = 0; i < candidates.size(); i++) {
//			if (bitmask[i] == '0') {
//				tmp.confirmed.erase(candidates[i]);
//			}
//		}
//		flag = false;
//		for (int i = 0; i < L_rd[Max_k].size(); i++) {
//			if (CheckRepeatTopk(L_rd[Max_k][i], tmp)) {
//				flag = true;
//				break;
//			}
//		}
//
//		if (!flag) {
//			Halfspace tmp_r = this_r.r;
//			// update tmp_r.r
//			for (int i = 0; i < candidates.size(); i++) {
//				if (bitmask[i] == '0') {
//					for (int j = 0; j < candidates.size(); j++) {
//						if (bitmask[i] == '1') {
//							int HP = candidates[i] * All_records.size() + candidates[j];
//							tmp_r.HPs.insert({ HP,false }); // the halfspace S(i)<S(j)
//						}
//					}
//				}
//			}
//
//			if (Halfspace::is_Feasible(tmp_r, dim)) {
//				L_rd[Max_k].push_back(tmp);
//				cnt++;
//			}
//		}
//
//		for (int i = 0; i < candidates.size(); i++) {
//			if (bitmask[i] == '0') {
//				tmp.confirmed.insert(candidates[i]);
//			}
//		}
//	} while (std::prev_permutation(bitmask.begin(), bitmask.end()));
//	return cnt;
//}
//int KLEVEL::LevelGrouping(int level, vector<region_rd>& this_level, region_rd& r, vector<Obj>& All_records, set<int>& utk_results, vector<int>& k_rskyband) {
//	region_rd tmp;
//	tmp.level = level; tmp.objID = -1;
//	tmp.confirmed = r.confirmed;
//	for (auto it = k_rskyband.begin(); it != k_rskyband.end(); it++) tmp.confirmed.insert(*it);
//	tmp.candidates.clear();
//	
//	if (k_rskyband.size() == Max_k - level + 1) {
//		tmp.r = r.r;
//		this_level.push_back(tmp);
//		return 1;
//	}
//
//	int cnt = 0;
//	for (auto it = k_rskyband.begin(); it != k_rskyband.end(); it++) {
//		tmp.confirmed.erase(*it);
//		bool flag = false;
//		for (int i = 0; i < L_rd[Max_k].size(); i++) {
//			if (CheckRepeatTopk(L_rd[Max_k][i], tmp)) {
//				flag = true;
//				break;
//			}
//		}
//		if (flag) {
//			tmp.confirmed.insert(*it);
//			continue;
//		}
//		Halfspace tmp_r = r.r;
//		// (*it) is the (k+1)^th record
//		for (auto itt = k_rskyband.begin(); itt != k_rskyband.end(); itt++) {
//			if (*it == *itt) continue;
//			int HP = (*it) * All_records.size() + *itt;
//			tmp_r.HPs.insert({ HP,false }); // the halfspace S(it)<S(itt)
//		}
//		if (Halfspace::is_Feasible(tmp_r, dim)) {
//			L_rd[Max_k].push_back(tmp);
//			cnt++;
//		}
//		tmp.confirmed.insert(*it);
//	}
//	return cnt;
//}
//void KLEVEL::Insert(Obj & q, int min_level) {
//	vector<int> HP_q; HP_q.clear();
//	for (int i = 0; i < AllObj.size(); i++) {
//		HP_q.push_back(Halfspace::ComputeHP(q, AllObj[i], dim));
//	}
//
//	int RegionSet_size = RegionSet.size();
//	bool flag = false;
//	for (int i = 0; i < RegionSet_size; i++) {
//		if (RegionSet[i].k > Max_k) continue;
//		if (RegionSet[i].k < min_level) continue;
//		int o_i = RegionSet[i].objID;
//		int HP = HP_q[o_i];
//		// {HP,true},  q>=o_i
//		// {HP,false}, q<=o_i
//		if (XdominateY(AllObj[o_i], q, dim)) continue;
//		if (!Halfspace::is_Feasible(RegionSet[i].r, dim, HP, true)) continue;
//		flag = true;
//		if (!Halfspace::is_Feasible(RegionSet[i].r, dim, HP, false)) {
//			// HP doesn't intersect with G[cur].r
//			RegionSet[i].k++;
//			continue;
//		}
//		else {
//			Halfspace r_LE, r_GE; // q <= o_i in r_LE; q >= o_i in r_GE;
//			Halfspace::Split(RegionSet[i].r, HP, r_LE, r_GE); // {HP,false}->r_LE, {HP,true}->r_GE
//			RegionSet[i].r = r_LE;
//
//			if (RegionSet[i].k < Max_k) {
//				kregion NewRegion;
//				NewRegion.k = RegionSet[i].k + 1; // q>=o_i
//				NewRegion.objID = RegionSet[i].objID;
//				NewRegion.r = r_GE;
//				RegionSet.push_back(NewRegion);
//			}
//		}
//	}
//	if ((!flag) && (AllObj.size() >= Max_k)) {
//		AllObj.push_back(q);
//		return;
//	}
//
//	// PreCompute Klevel before insertion
//	vector<Halfspace> Rq;
//	vector<int> rank_q;
//	FindVoronoiWithRank(q, Rq, rank_q, HP_q);
//
//	// Insert all regions generated from q
//
//	for (int i = 0; i < Rq.size(); i++) {
//		if (rank_q[i] > Max_k) continue;
//		if (rank_q[i] > (AllObj.size() + 1)) cout << "rank error" << endl;
//
//		kregion NewRegion;
//		NewRegion.k = rank_q[i];
//		NewRegion.objID = q.ID;
//		NewRegion.r = Rq[i];
//
//		RegionSet.push_back(NewRegion);
//	}
//
//	AllObj.push_back(q);
//	Klevel_Obj.clear();
//	kth_Obj.clear();
//
//	for (int i = 0; i <= Max_k; i++) Klevel_buffer[i].Nodes.clear();
//
//	for (int i = 0; i < RegionSet.size(); i++) {
//		if (RegionSet[i].k > Max_k) continue;
//		Klevel_buffer[RegionSet[i].k].Nodes.push_back(i);
//	}
//
//	for (int i = 1; i <= Max_k; i++) {
//		cout << i << "-level size: " << Klevel_buffer[i].Nodes.size() << endl;
//		for (int j = 0; j < Klevel_buffer[i].Nodes.size(); j++) {
//			if (i == Max_k) {
//				kth_Obj.insert(RegionSet[Klevel_buffer[i].Nodes[j]].objID);
//			}
//			Klevel_Obj.insert(RegionSet[Klevel_buffer[i].Nodes[j]].objID);
//		}
//	}
//	cout << "Objects in <=Klevel: " << Klevel_Obj.size() << endl;
//}
//
//void KLEVEL::FindVoronoiWithRank(Obj & q, vector<Halfspace> & Rq, vector<int> & rank_q, vector<int> & HP_q) {
//	Halfspace origin_r;
//	Rq.clear(); rank_q.clear();
//	Rq.push_back(origin_r); rank_q.push_back(1);
//
//	//for (int i = 0; i < AllObj.size(); i++) {
//	//	int tmpHP = HP_q[i];
//	for (auto it = Klevel_Obj.begin(); it != Klevel_Obj.end(); it++) {
//		int tmpHP = HP_q[*it];
//		int Rq_size = Rq.size();
//		for (int j = 0; j < Rq_size; j++) {
//			if (rank_q[j] > Max_k) continue;
//			if (rank_q[j] < Max_k) {
//				// rank(o_q) <= rank(o_i) 
//				if (Halfspace::is_Feasible(Rq[j], dim, tmpHP, false)) rank_q[j]++;
//				if (Halfspace::is_intersected(Rq[j], dim, tmpHP)) {
//					Halfspace r_LE, r_GE;
//					Halfspace::Split(Rq[j], tmpHP, r_LE, r_GE);
//					Rq[j] = r_LE;
//					Rq.push_back(r_GE); rank_q.push_back(rank_q[j] - 1); // rank(o_q)>=rank(o_i) in r_GE
//				}
//			}
//			else { // rank_q[j]=Max_k;
//				// rank(o_q) <= rank(o_i) 
//				if (Halfspace::is_Feasible(Rq[j], dim, tmpHP, false)) rank_q[j]++;
//				if (Halfspace::is_intersected(Rq[j], dim, tmpHP)) {
//					Halfspace r_LE, r_GE;
//					Halfspace::Split(Rq[j], tmpHP, r_LE, r_GE);
//					// ranking of o_q in r_LE is larger than Max_k, so ignore r_LE
//					Rq[j] = r_GE; rank_q[j]--; // rank_q[j] still equal to Max_k
//				}
//			}
//		}
//	}
//
//	//ignore rank_q[i]>Max_k 
//
//}
//
//int KLEVEL::Update_SplitNode(int cur, int HP) {
//	Halfspace r_LE, r_GE; // q <= o_i in r_LE; q >= o_i in r_GE;
//	Halfspace::Split(G[cur].r, HP, r_LE, r_GE); // {HP,false}->r_LE, {HP,true}->r_GE
//	G[cur].r = r_LE;
//	int NewNode = CreateNode(G[cur].pID, r_GE);
//
//	for (set<int>::iterator it = G[cur].minusEdge.begin(); it != G[cur].minusEdge.end();) {
//		if (Halfspace::is_Feasible(G[*it].r, dim, HP, true)) { // connect it and NewNode;
//		//if (Halfspace::is_Feasible(G[*it].r, G[NewNode].r, dim)) { discuss!!!
//			G[NewNode].minusEdge.insert(*it);
//			G[*it].plusEdge.insert(NewNode);
//		}
//		if (!Halfspace::is_Feasible(G[*it].r, dim, HP, false)) { // disconnect it and cur
//			auto itt = G[*it].plusEdge.find(cur);
//			if (itt != G[*it].plusEdge.end())G[*it].plusEdge.erase(itt);
//			it = G[cur].minusEdge.erase(it);
//		}
//		else it++;
//	}
//
//	for (set<int>::iterator it = G[cur].plusEdge.begin(); it != G[cur].plusEdge.end();) {
//		if (Halfspace::is_Feasible(G[*it].r, dim, HP, true)) {
//			G[NewNode].plusEdge.insert(*it);
//			G[*it].minusEdge.insert(NewNode);
//		}
//		if (!Halfspace::is_Feasible(G[*it].r, dim, HP, false)) {
//			auto itt = G[*it].minusEdge.find(cur);
//			if (itt != G[*it].minusEdge.end()) G[*it].minusEdge.erase(G[*it].minusEdge.find(cur));
//			it = G[cur].plusEdge.erase(it);
//		}
//		else it++;
//	}
//
//	return NewNode;
//}
//
//void KLEVEL::Update_SplitNode_lastLevel(int cur, int HP) {
//	Halfspace r_LE, r_GE; // q <= o_i in r_LE; q >= o_i in r_GE;
//	Halfspace::Split(G[cur].r, HP, r_LE, r_GE); // {HP,false}->r_LE, {HP,true}->r_GE
//	G[cur].r = r_LE;
//
//	for (set<int>::iterator it = G[cur].minusEdge.begin(); it != G[cur].minusEdge.end();) {
//		if (!Halfspace::is_Feasible(G[*it].r, dim, HP, false)) { // disconnect it and cur
//			auto itt = G[*it].plusEdge.find(cur);
//			if (itt != G[*it].plusEdge.end())G[*it].plusEdge.erase(itt);
//			it = G[cur].minusEdge.erase(it);
//		}
//		else it++;
//	}
//
//	for (set<int>::iterator it = G[cur].plusEdge.begin(); it != G[cur].plusEdge.end();) {
//		if (!Halfspace::is_Feasible(G[*it].r, dim, HP, false)) {
//			auto itt = G[*it].minusEdge.find(cur);
//			if (itt != G[*it].minusEdge.end()) G[*it].minusEdge.erase(G[*it].minusEdge.find(cur));
//			it = G[cur].plusEdge.erase(it);
//		}
//		else it++;
//	}
//
//	return;
//}
//
//void KLEVEL::FindKlevel(int level_limit) { // buffer 1-level to k-level
//	set<int> temp;
//	for (int i = 1; i <= level_limit; i++) { // bfs
//		temp.clear();
//		for (auto it = Klevel_buffer[i - 1].Nodes.begin(); it != Klevel_buffer[i - 1].Nodes.end(); it++) {
//			for (auto itt = G[*it].plusEdge.begin(); itt != G[*it].plusEdge.end(); itt++) {
//				temp.insert(*itt);
//			}
//		}
//
//		Klevel_buffer[i].Nodes.clear();
//		for (auto it = temp.begin(); it != temp.end(); it++) {
//			if (*it != bot) Klevel_buffer[i].Nodes.push_back(*it);
//			//else cout << "Find K level ERROR!!!" << endl; // discuss!!!
//		}
//	}
//	return;
//}
//
//void KLEVEL::Update_NewNode(int cur, vector<int> & aboveNodes, vector<int> & belowNodes) {
//
//	for (int i = 0; i < aboveNodes.size(); i++) {
//		G[cur].minusEdge.insert(aboveNodes[i]);
//		G[aboveNodes[i]].plusEdge.insert(cur);
//	}
//	for (int i = 0; i < belowNodes.size(); i++) {
//		G[cur].plusEdge.insert(belowNodes[i]);
//		G[belowNodes[i]].minusEdge.insert(cur);
//	}
//
//	for (auto i = aboveNodes.begin(); i != aboveNodes.end(); i++) {
//		for (auto j = belowNodes.begin(); j != belowNodes.end(); j++) {
//			if (G[*i].plusEdge.find(*j) != G[*i].plusEdge.end()) {
//				G[*i].plusEdge.erase(*j);
//				G[*j].minusEdge.erase(*i);
//			}
//		}
//	}
//	return;
//}
//
//void KLEVEL::DeleteNode(int cur) {
//	// cur is in Max_k+1 level
//	if (G[bot].minusEdge.find(cur) == G[bot].minusEdge.end()) cout << "Delete Error!!!" << endl;
//	G[bot].minusEdge.erase(cur);
//	for (auto it = G[cur].minusEdge.begin(); it != G[cur].minusEdge.end(); it++) {
//		G[*it].plusEdge.clear();
//		G[*it].plusEdge.insert(bot);
//		G[bot].minusEdge.insert(*it);
//	}
//	return;
//}
//
//int KLEVEL::CreateNode(int ID, Halfspace & r) {
//	Node NewNode;
//	NewNode.pID = ID;
//	NewNode.r = r;
//	NewNode.plusEdge.clear();
//	NewNode.minusEdge.clear();
//	G.push_back(NewNode);
//	return G.size() - 1;
//}
//
//void KLEVEL::CheckLinkERROR() {
//	for (int i = 0; i <= Max_k + 1; i++) {
//		for (int j = 0; j < Klevel_buffer[i].Nodes.size(); j++) {
//			int cur = Klevel_buffer[i].Nodes[j];
//			for (auto it = G[cur].minusEdge.begin(); it != G[cur].minusEdge.end(); it++) {
//				if (G[*it].plusEdge.find(cur) == G[*it].plusEdge.end()) {
//					cout << i << ' ' << cur << ' ' << *it << endl;
//					cout << "plus Link ERROR" << endl;
//					system("pause");
//				}
//			}
//			for (auto it = G[cur].plusEdge.begin(); it != G[cur].plusEdge.end(); it++) {
//				if (G[*it].minusEdge.find(cur) == G[*it].minusEdge.end()) {
//					cout << i << ' ' << cur << ' ' << *it << endl;
//					cout << "minus Link ERROR" << endl;
//					system("pause");
//				}
//			}
//		}
//	}
//	return;
//}
