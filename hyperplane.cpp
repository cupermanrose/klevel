#include "hyperplane.h"

Halfspace::Halfspace() {
	HPs.clear(); // Halfspace is the whole space
}

Halfspace::~Halfspace() {
    HPs.clear();
    unordered_map<int, bool>().swap(HPs);
    vertices.clear();
    vector<vector<float>>().swap(vertices);
}

Halfspace::Halfspace(vector<float>& Qregion, int dim) {
	HPs.clear();
	vertices.clear();
	for (int i = 0; i < (1 << (dim - 1)); i++) {
		vector<float> p; p.clear();
		int tmp = i;
		for (int d = 0; d < dim - 1; d++) {
			if (tmp % 2 == 0) p.push_back(Qregion[2 * d]);
			else p.push_back(Qregion[2 * d + 1]);
			tmp >> 1;
		}
		vertices.emplace_back(p);
	}

	for (int d = 0; d < dim - 1; d++) innerPoint[d] = (Qregion[2 * d] + Qregion[2 * d + 1]) / 2.0;
}

/*
bool Halfspace::is_Feasible(Halfspace& r, int dimm, float range[]) {

    double row[Max_Dimension];

    int dim = dimm - 1;

    lprec* lp = make_lp(0, dim);

    r.lpModel(lp, dim);

    for (unordered_map<int, bool>::iterator it = r.HPs.begin(); it != r.HPs.end(); it++) {
        r.addHP(lp, dim, it->first, it->second);
    }

    //row[0] = 0;
    //for (int d = 1; d < dim + 1; d++)
    //{
    //	row[d] = -1;
    //}
    //set_obj_fn(lp, row); //discuss!!!
    //set_verbose(lp, IMPORTANT);
    set_verbose(lp, SEVERE);
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
    //set_presolve(lp, PRESOLVE_ROWDOMINATE, get_presolveloops(lp));
    int ret = solve(lp);


    if (ret == 0)
    {

        for (int i = 1; i <= dim; i++) {
            double var[Max_Dimension];
            row[0] = 0;
            for (int j = 1; j <= dim; j++) {
                if (i == j) row[j] = 1;
                else row[j] = 0;
            }

            set_minim(lp);
            set_obj_fn(lp, row);
            solve(lp);
            get_variables(lp, var);
            range[(i - 1) * 2] = var[i - 1];

            set_maxim(lp);
            set_obj_fn(lp, row);
            solve(lp);
            get_variables(lp, var);
            range[(i - 1) * 2 + 1] = var[i - 1];
        }
        // for reduced space
        delete_lp(lp);
        return true;
    }
    else
    {
        // for reduced space
        delete_lp(lp);
        return false;
    }
}




bool Halfspace::is_Feasible(Halfspace& r, int dimm, float range[]) {

	double row[Max_Dimension];

	int dim = dimm - 1;

	lprec* lp = make_lp(0, dim);

	r.lpModel(lp, dim);

	for (unordered_map<int, bool>::iterator it = r.HPs.begin(); it != r.HPs.end(); it++) {
		r.addHP(lp, dim, it->first, it->second);
	}

	//row[0] = 0;
	//for (int d = 1; d < dim + 1; d++)
	//{
	//	row[d] = -1;
	//}
	//set_obj_fn(lp, row); //discuss!!!
	//set_verbose(lp, IMPORTANT);
	set_verbose(lp, SEVERE);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
	//set_presolve(lp, PRESOLVE_ROWDOMINATE, get_presolveloops(lp));
	int ret = solve(lp);
	

	if (ret == 0)
	{
		
		for (int i = 1; i <= dim; i++) {
			double var[Max_Dimension];
			row[0] = 0;
			for (int j = 1; j <= dim; j++) {
				if (i == j) row[j] = 1;
				else row[j] = 0;
			}

			set_minim(lp);
			set_obj_fn(lp, row);
			solve(lp);
			get_variables(lp, var);
			range[(i - 1) * 2] = var[i - 1];

			set_maxim(lp);
			set_obj_fn(lp, row);
			solve(lp);
			get_variables(lp, var);
			range[(i - 1) * 2 + 1] = var[i - 1];
		}
		// for reduced space
		delete_lp(lp);
		return true;
	}
	else
	{
		// for reduced space
		delete_lp(lp);
		return false;
	}
}

bool Halfspace::is_Feasible(Halfspace & r, int dimm, int HP, bool sideindicator) {
	// check the same HP and opposite sideindicator
	for (auto it = r.HPs.begin(); it != r.HPs.end(); it++) {
		if ((it->first == HP) && (it->second != sideindicator)) return false;
	}

	double row[Max_Dimension];

	int dim = dimm - 1;

	lprec * lp = make_lp(0, dim);

	r.lpModel(lp, dim);

	for (unordered_map<int, bool>::iterator it = r.HPs.begin(); it != r.HPs.end(); it++) {
		r.addHP(lp, dim, it->first, it->second);
	}

	r.addHP(lp, dim, HP, sideindicator);

	row[0] = 0;
	for (int d = 1; d < dim + 1; d++)
	{
		row[d] = -1;
	}
	set_obj_fn(lp, row); //discuss!!!
	//set_verbose(lp, IMPORTANT);
	set_verbose(lp, SEVERE);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
	//set_presolve(lp, PRESOLVE_ROWS, get_presolveloops(lp));
	//set_presolve(lp, PRESOLVE_ROWS | PRESOLVE_COLS | PRESOLVE_LINDEP, get_presolveloops(lp));

	int ret = solve(lp);
	//get_variables(lp, row);

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

bool Halfspace::is_Feasible(Halfspace & r_a, Halfspace & r_b, int dimm) {
	// check the same HP and opposite sideindicator
	for (auto it_a = r_a.HPs.begin(); it_a != r_a.HPs.end(); it_a++) {
		for (auto it_b = r_b.HPs.begin(); it_b != r_b.HPs.end(); it_b++) {
			if ((it_a->first == it_b->first) && (it_a->second != it_b->second)) return false;
		}
	}

	double row[Max_Dimension];

	int dim = dimm - 1;

	lprec* lp = make_lp(0, dim);
	r_a.lpModel(lp, dim);

	for (unordered_map<int, bool>::iterator it = r_a.HPs.begin(); it != r_a.HPs.end(); it++) {
		r_a.addHP(lp, dim, it->first, it->second);
	}

	for (unordered_map<int, bool>::iterator it = r_b.HPs.begin(); it != r_b.HPs.end(); it++) {
		r_a.addHP(lp, dim, it->first, it->second);
	}

	row[0] = 0;
	for (int d = 1; d < dim + 1; d++)
	{
		row[d] = -1;
	}
	set_obj_fn(lp, row); //discuss!!!
	//set_verbose(lp, IMPORTANT);
	set_verbose(lp, SEVERE);
	set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
	//set_presolve(lp, PRESOLVE_ROWDOMINATE, get_presolveloops(lp));
	int ret = solve(lp);
	//get_variables(lp, row);

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

bool Halfspace::is_intersected(Halfspace & r, int dim, int HP) {
	if (!is_Feasible(r, dim, HP, false)) return false;
	if (!is_Feasible(r, dim, HP, true)) return false;
	return true;
}



void Halfspace::Split(Halfspace & r, int HP, Halfspace & r_LE, Halfspace & r_GE) {
	r_LE = r;
	r_LE.HPs.insert(pair<int, bool>{HP, false});
	r_GE = r;
	r_GE.HPs.insert(pair<int, bool>{HP, true});
}
*/


void Halfspace::addHP(lprec* lp, int dim, int HP, bool sideindicator) {
    // EPS
    float EPS_control = 0.000001; // discuss!!!

    double row[Max_Dimension];
    row[0] = 0;
    for (int d = 1; d <= dim; d++)
    {
        row[d] = AllHP[HP].w[d - 1];
    }
    if (sideindicator == false) // o1 <= o2
    {
        add_constraint(lp, row, LE, AllHP[HP].w[dim] - EPS_control);
    }
    else if (sideindicator == true) /// o1>=o2
    {
        add_constraint(lp, row, GE, AllHP[HP].w[dim] + EPS_control);
    }
    else
    {
        std::cout << "Unable to detect half plane direction!!!" << endl;
    }
}

void Halfspace::lpModel(lprec * lp, int dim) {
    double row[Max_Dimension];

    if (lp == NULL)
    {
        fprintf(stderr, "Unable to create new LP model\n");
        exit(0);
    }

    // constraints in each dimension //???
    for (int d = 1; d <= dim; d++)
    {
        row[0] = 0;
        for (int d_s = 1; d_s <= dim; d_s++)
        {
            if (d_s == d)
            {
                row[d_s] = 1.0;
            }
            else
            {
                row[d_s] = 0.0;
            }
        }

        add_constraint(lp, row, GE, 0.0);
        add_constraint(lp, row, LE, 1.0);
    }

    // in reduced space, sum_{w_i} should less than 1 // add it later discuss!!!
    row[0] = 0.0;
    for (int d = 1; d <= dim; d++)
    {
        row[d] = 1.0;
    }
    add_constraint(lp, row, GE, 0);
    add_constraint(lp, row, LE, 1);

    // constraints in intersected hyperplanes
    for (unordered_map<int, bool>::iterator iter = HPs.begin(); iter != HPs.end(); iter++)
    {
        addHP(lp, dim, iter->first, iter->second);
    }

    // set scale
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
}

bool Halfspace::is_Feasible(Halfspace & r, int dimm) {

    double row[Max_Dimension];

    int dim = dimm - 1;

    lprec* lp = make_lp(0, dim);

    r.lpModel(lp, dim);

    /*for (unordered_map<int, bool>::iterator it = r.HPs.begin(); it != r.HPs.end(); it++) {
        r.addHP(lp, dim, it->first, it->second);
    }*/


    //set_verbose(lp, IMPORTANT);
    set_verbose(lp, SEVERE);
    set_scaling(lp, SCALE_GEOMETRIC + SCALE_EQUILIBRATE + SCALE_INTEGERS);
    //set_presolve(lp, PRESOLVE_ROWDOMINATE, get_presolveloops(lp));

    double var[Max_Dimension], var1[Max_Dimension];
    for (int i = 0; i < dim + 1; i++) row[i] = 0.0;
    row[1] = 1.0;
    set_maxim(lp);
    set_obj_fn(lp, row);
    set_timeout(lp,1);
    int ret = solve(lp);

    get_variables(lp, var);

    // for reduced space

    if (ret == 0)
    {
        if (dim > 1) {
            row[1] = 0.0; row[2] = 1.0;
            set_obj_fn(lp, row);
            ret = solve(lp);
            get_variables(lp, var1);
            for (int i = 0; i < dim; i++) var[i] = (var[i] + var1[i]) / 2.0;
        }
        if (dim > 2) {
            row[2] = 0.0; row[3] = 1.0;
            set_obj_fn(lp, row);
            ret = solve(lp);
            get_variables(lp, var1);
            for (int i = 0; i < dim; i++) var[i] = (var[i] + var1[i]) / 2.0;
        }
        for (int i = 0; i < dim; i++) r.innerPoint[i] = var[i];
        delete_lp(lp);
        return true;
    }
    else
    {
        if (ret==7) cout <<"lpsolver time out" << endl;
        delete_lp(lp);
        return false;
    }
}

int Halfspace::ComputeHP(Obj & o1, Obj & o2, int dim) {
    Hyperplane tmpHP;
    tmpHP.w.clear(); tmpHP.o1 = o1.ID; tmpHP.o2 = o2.ID;
    float o1_d = o1.w[dim - 1];
    float o2_d = o2.w[dim - 1];
    for (int d = 0; d < dim - 1; d++) {
        tmpHP.w.push_back((o1.w[d] - o1_d) - (o2.w[d] - o2_d));
    }
    tmpHP.w.push_back(o2_d - o1_d);
    AllHP.push_back(tmpHP);
    return AllHP.size() - 1;
}

void Halfspace::ComputeVertexByQhull(int dim, Halfspace& r) {
	//int dim = DIM;             /* dimension of points */
	const int numpoints = r.HPs.size() + (dim - 1) * 2 + 1;            /* number of points */
	//coordT* halfspaces = (coordT*)malloc(dim * numpoints * sizeof(coordT)); /* array of coordinates for each point */
	coordT* halfspaces = new coordT[dim*numpoints];

	//const int numpoints = (dim - 1) * 2;            /* number of points */
	//coordT* halfspaces = (coordT*)malloc(dim * numpoints * sizeof(coordT)); /* array of coordinates for each point */

	// Ax+B<=0
	int res = 0;
	for (int d = 0; d < dim; d++) {
		if (d < dim-1) {
			halfspaces[res] = 1.0;
			res++;
		}
		else {
			halfspaces[res] = -1.0;
			res++;
		}
	}
	for (int i = 0; i < dim-1; i++) {
		for (int d = 0; d < dim - 1; d++) {
			if (i == d) halfspaces[res] = -1.0;
			else halfspaces[res] = 0.0;
			res++;
		}
		halfspaces[res] = 0.0; res++;
		for (int d = 0; d < dim - 1; d++) {
			if (i == d) halfspaces[res] = 1.0;
			else halfspaces[res] = 0.0;
			res++;
		}
		halfspaces[res] = -1.0; res++;
	}
	for (auto it = r.HPs.begin(); it != r.HPs.end(); it++) {
		for (int d = 0; d < dim - 1; d++) {
			if (it->second) halfspaces[res] = -AllHP[it->first].w[d];
			else halfspaces[res] = AllHP[it->first].w[d];
			res++;
		}
		if (it->second) halfspaces[res] = AllHP[it->first].w[dim - 1];
		else halfspaces[res] = -AllHP[it->first].w[dim - 1];
		res++;
	}

	//coordT* rows[TOTpoints];
	boolT ismalloc = True;    /* True if qhull should free points in qh_freeqhull() or reallocation */
	char flags[250];          /* option flags for qhull, see qh-quick.htm */
	//FILE* outfile = stdout;    /* output from qh_produce_output()
								// use NULL to skip qh_produce_output() */
	//FILE* errfile = stderr;    /* error messages from qhull code */
	
	int exitcode;             /* 0 if no error from qhull */
	facetT* facet;            /* set by FORALLfacets */
	int curlong, totlong;     /* memory remaining after qh_memfreeshort, used if !qh_NOmem  */
	//int i;

	qhT qh_qh;                /* Qhull's data structure.  First argument of most calls */
	qhT* qh = &qh_qh;

	QHULL_LIB_CHECK
	qh_zero(qh, NULL);
	fflush(NULL);

	/* use qh_sethalfspace_all to transform the halfspaces yourself.
	   If so, set 'qh->feasible_point and do not use option 'Hn,...' [it would retransform the halfspaces]
	*/
	//char qhull_cmd[] = "qhull H0 s Tcv Fp";	
	char qhull_cmd[] = "qhull H0 Fp";

	coordT* feasible_point = new coordT[dim - 1];
	for (int i = 0; i < dim - 1; i++) feasible_point[i] = r.innerPoint[i];
    exitcode = qh_new_qhull_klevel(qh, dim, numpoints, halfspaces, ismalloc, qhull_cmd, feasible_point, NULL, NULL);

    //exitcode = qh_new_qhull(qh, dim, numpoints, halfspaces, ismalloc, qhull_cmd, NULL, NULL);

	r.vertices.clear();

	boolT zerodiv;
	FORALLfacets{
        vector<float> tmp_vertex;
		tmp_vertex.clear();
		for (int d = 0; d < qh->hull_dim; d++) {
			if (facet->offset < -qh->MINdenom) {
				tmp_vertex.push_back((facet->normal[d] / -facet->offset) + qh->feasible_point[d]);
			}
			else {
				tmp_vertex.push_back(qh_divzero(facet->normal[d], facet->offset, qh->MINdenom_1,
					&zerodiv) + qh->feasible_point[d]);
			}
		}
		r.vertices.emplace_back(tmp_vertex);
	}
	
#ifdef qh_NOmem
	qh_freeqhull(qh, qh_ALL);
#else
	qh_freeqhull(qh, !qh_ALL);
	qh_memfreeshort(qh, &curlong, &totlong);
	if (curlong || totlong)  /* could also check previous runs */
		fprintf(stderr, "qhull internal warning (user_eg, #3): did not free %d bytes of long memory (%d pieces)\n",
			totlong, curlong);
#endif
	//qh_freeqhull(qh, !qh_ALL);

    // memory free

    return;
}

/*
bool Halfspace::is_Feasible_cplex(Halfspace& r, int dimm) {
	// EPS
	float EPS_control = 0.000001; // discuss!!!
	int dim = dimm - 1;
	IloEnv env;
	IloModel model(env);
	IloNumVarArray vars(env);
	for (int i = 0; i < dim; i++) vars.add(IloNumVar(env, 0.0, 1.0)); // all dimension 0<=d_i<=1.0
	
	// d_0+d_1+...<=1.0
	IloExpr expr(env);
	for (int i = 0; i < dim; i++) {
		expr += vars[i];
	}
	model.add(expr <= 1.0);
	expr.end();
	
	
	for (auto it = r.HPs.begin(); it != r.HPs.end(); it++) {
		int HPid = it->first;
		//cout << HPid << endl;
		bool side = it->second;
		IloExpr expr(env);
		for (int d = 0; d < dim; d++) {
			expr += AllHP[HPid].w[d] * vars[d];
		}
		if (side = false) {
			model.add(expr <= AllHP[HPid].w[dim] - EPS_control);
		}
		else {
			model.add(expr >= AllHP[HPid].w[dim] + EPS_control);
		}
		expr.end();
	}

	IloObjective obj_max0 = IloMaximize(env, vars[0]);
	IloObjective obj_max1 = IloMaximize(env, vars[1]);
	IloObjective obj_max2 = IloMaximize(env, vars[2]);

	model.add(obj_max0);

	IloCplex cplex(model);
	//cplex.setParam(IloCplex::PreInd, false);
	//cplex.setParam(IloCplex::ScaInd, 1);
	
	//cplex.setParam(IloCplex::RootAlg, IloCplex::Barrier); // set optimizer used interior point method
	//cplex.setParam(IloCplex::Param::Parallel, CPX_PARALLEL_DETERMINISTIC); // multi thread
	cplex.setParam(IloCplex::RootAlg, IloCplex::Concurrent); // set optimizer used interior point method
	cplex.setParam(IloCplex::Param::Parallel, CPX_PARALLEL_OPPORTUNISTIC); // multi thread
	cplex.setParam(IloCplex::Param::Threads, 8);

	cplex.setOut(env.getNullStream());

	if (!cplex.solve()) {
		env.end();
		return false;
	}
	else {
		IloNumArray vals(env);
		cplex.getValues(vals, vars);
		for (int i = 0; i < dim; i++) r.innerPoint[i] = vals[i];

		if (dim > 1) {
			model.remove(obj_max0);
			model.add(obj_max1);
			cplex.solve();
			cplex.getValues(vals, vars);
			for (int i = 0; i < dim; i++) r.innerPoint[i] = (r.innerPoint[i] + vals[i]) / 2.0;
		}
		if (dim > 2) {
			model.remove(obj_max1);
			model.add(obj_max2);
			cplex.solve();
			cplex.getValues(vals, vars);
			for (int i = 0; i < dim; i++) r.innerPoint[i] = (r.innerPoint[i] + vals[i]) / 2.0;
		}
		env.end();
		return true;
	}
}
*/
