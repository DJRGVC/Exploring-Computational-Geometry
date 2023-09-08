#ifndef __geom_h
#define __geom_h

#include <vector>
#include <queue>

using namespace std; 


/**
 * A point in 2D space.
 */
typedef struct _point2d {
  int x,y; 
} point2d;

/**
 * struct for storing a point, its distance from p0, and its radial angle
 */
typedef struct _radial_point2d {
  point2d p; 
  double distance; 
  double angle; 
  bool operator<(const _radial_point2d& p) const {
    if (angle > p.angle) return true;
    if (angle < p.angle) return false;
    return distance > p.distance;
  }
} radial_point2d;



/* returns 2 times the signed area of triangle abc. The area is
   positive if c is to the left of ab, 0 if a,b,c are collinear and
   negative if c is to the right of ab
 */
int signed_area2D(point2d a, point2d b, point2d c); 

/* finds the distance between two points */
int distance(point2d p1, point2d p2);

/* sorts the points in the vector by angle and distance from the
   origin. The points are sorted in counter-clockwise order. The
   origin is the first point in the vector. The points are sorted
   using a priority queue. The priority queue is a heap, so the
   complexity is O(n log n).
 */
void sortRadially(point2d p0, vector<point2d>& pts, priority_queue<radial_point2d>& sortedPoints);

/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2d p, point2d q, point2d r);

/* return 1 if c is  strictly left of ab; 0 otherwise */
int left_strictly (point2d a, point2d b, point2d c); 

/* return 1 if c is left of ab or on ab; 0 otherwise */
int left_on(point2d a, point2d b, point2d c); 

/* compute the convex hull of a set of points using the Graham
   scan. The points are sorted radially, and the convex hull is
   computed using a stack. The complexity is O(n log n).
 */
void graham_scan(vector<point2d>& pts, vector<point2d>& hull);

#endif
