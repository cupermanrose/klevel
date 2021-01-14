#ifndef _OBJECT_H
#define _OBJECT_H

#define Max_Dimension 8

#include <vector>
#include <iostream>

using namespace std;

class Obj {
public:
	int ID;
	float w[Max_Dimension];

public:
	Obj();
	Obj(int i, int d, float* w = 0);
	Obj(int i, int d, vector<float>& object);
	~Obj();

	/*static float GetIntersection(Point& A, Point& B);
	static bool A_larger_B_after_intersect(Point& A, Point& B);*/
};

#endif // !_OBJECT_H
