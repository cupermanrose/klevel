#include "object.h"
#include "header.h"

Obj::Obj() {
}

Obj::Obj(int i, int  d, float* a_w) {
	ID = i;
	if (a_w != 0)
		memcpy(w, a_w, sizeof(float) * d);
	else
		memset(w, 0, sizeof(float) * d);
}

Obj::Obj(int a_id, int  dim, vector<float>& object) {
    ID = a_id;
    for (int d=0;d<dim;d++)
        w[d]=object[d];
}

Obj::~Obj() {
}