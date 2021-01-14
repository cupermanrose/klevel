/*<html><pre>  -<a                             href="../libqhull_r/qh-user_r.htm"
  >-------------------------------</a><a name="TOP">-</a>

  user_eg_r.c
  sample code for calling qhull() from an application.  Uses reentrant libqhull_r

  call with:

	 user_eg "cube/diamond options" "delaunay options" "halfspace options"

  for example:

	 user_eg                             # return summaries

	 user_eg "n" "o" "Fp"                # return normals, OFF, points

	 user_eg "n Qt" "o" "Fp"             # triangulated cube

	 user_eg "QR0 p" "QR0 v p" "QR0 Fp"  # rotate input and return points
										 # 'v' returns Voronoi
										 # transform is rotated for halfspaces

   main() makes three runs of qhull.

	 1) compute the convex hull of a cube

	 2a) compute the Delaunay triangulation of random points

	 2b) find the Delaunay triangle closest to a point.

	 3) compute the halfspace intersection of a diamond

 notes:

   For another example, see main() in unix_r.c and user_eg2_r.c.
   These examples, call qh_qhull() directly.  They allow
   tighter control on the code loaded with Qhull.

   For a C++ example, see user_eg3/user_eg3_r.cpp

   Summaries are sent to stderr if other output formats are used

   compiled by 'make bin/user_eg'

   see libqhull_r.h for data structures, macros, and user-callable functions.
*/

#include "qhull_ra.h"

/*-------------------------------------------------
-internal function prototypes
*/
void print_summary(qhT* qh);
void makecube(coordT* points, int numpoints, int dim);
void makeDelaunay(qhT* qh, coordT* points, int numpoints, int dim, int seed);
void findDelaunay(qhT* qh, int dim);
void makehalf(coordT* points, int numpoints, int dim);

/*-------------------------------------------------
-print_summary(qh)
*/
void print_summary(qhT* qh) {
	facetT* facet;
	int k;

	printf("\n%d vertices and %d facets with normals:\n",
		qh->num_vertices, qh->num_facets);
	FORALLfacets{
	  for (k = 0; k < qh->hull_dim; k++)
		printf("%6.2g ", facet->normal[k]);
	  printf("\n");
	}
}

/*--------------------------------------------------
-makecube- set points to vertices of cube
  points is numpoints X dim
*/
void makecube(coordT * points, int numpoints, int dim) {
	int j, k;
	coordT* point;

	for (j = 0; j < numpoints; j++) {
		point = points + j * dim;
		for (k = dim; k--; ) {
			if (j & (1 << k))
				point[k] = 1.0;
			else
				point[k] = -1.0;
		}
	}
} /*.makecube.*/

/*--------------------------------------------------
-makeDelaunay- set points for dim Delaunay triangulation of random points
  points is numpoints X dim.
notes:
  makeDelaunay() in user_eg2_r.c uses qh_setdelaunay() to project points in place.
*/
void makeDelaunay(qhT * qh, coordT * points, int numpoints, int dim, int seed) {
	int j, k;
	coordT* point, realr;

	printf("seed: %d\n", seed);
	qh_RANDOMseed_(qh, seed);
	for (j = 0; j < numpoints; j++) {
		point = points + j * dim;
		for (k = 0; k < dim; k++) {
			realr = qh_RANDOMint;
			point[k] = 2.0 * realr / (qh_RANDOMmax + 1) - 1.0;
		}
	}
} /*.makeDelaunay.*/

/*--------------------------------------------------
-findDelaunay- find the Delaunay triangle or adjacent triangle for [0.5,0.5,...]
  assumes dim < 100
notes:
  See <a href="../../html/qh-code.htm#findfacet">locate a facet with qh_findbestfacet()</a>
  calls qh_setdelaunay() to project the point to a parabaloid
warning:
  Errors if it finds a tricoplanar facet ('Qt').  The corresponding Delaunay triangle
  is in the set of tricoplanar facets or one of their neighbors.  This search
  is not implemented here.
*/
void findDelaunay(qhT * qh, int dim) {
	int k;
	coordT point[100];
	boolT isoutside;
	realT bestdist;
	facetT* facet;
	vertexT* vertex, ** vertexp;

	for (k = 0; k < dim; k++)
		point[k] = 0.5;
	qh_setdelaunay(qh, dim + 1, 1, point);
	facet = qh_findbestfacet(qh, point, qh_ALL, &bestdist, &isoutside);
	if (facet->tricoplanar) {
		fprintf(stderr, "findDelaunay: search not implemented for triangulated, non-simplicial Delaunay regions (tricoplanar facet, f%d).\n",
			facet->id);
		qh_errexit(qh, qh_ERRqhull, facet, NULL);
	}
	FOREACHvertex_(facet->vertices) {
		for (k = 0; k < dim; k++)
			printf("%5.2f ", vertex->point[k]);
		printf("\n");
	}
} /*.findDelaunay.*/

/*--------------------------------------------------
-makehalf- set points to halfspaces for a (dim)-dimensional diamond
  points is numpoints X dim+1

  each halfspace consists of dim coefficients followed by an offset
*/
void makehalf(coordT * points, int numpoints, int dim) {
	int j, k;
	coordT* point;

	for (j = 0; j < numpoints; j++) {
		point = points + j * (dim + 1);
		point[dim] = -1.0; /* offset */
		for (k = dim; k--; ) {
			if (j & (1 << k))
				point[k] = 1.0;
			else
				point[k] = -1.0;
		}
	}
} /*.makehalf.*/

#define DIM 3     /* dimension of points, must be < 31 for SIZEcube */
#define SIZEcube (1<<DIM)
#define SIZEdiamond (2*DIM)
#define TOTpoints (SIZEcube + SIZEdiamond)

/*--------------------------------------------------
-main- derived from Qhull-template in user_r.c

  see program header

  this contains three runs of Qhull for convex hull, Delaunay
  triangulation or Voronoi vertices, and halfspace intersection
  */
