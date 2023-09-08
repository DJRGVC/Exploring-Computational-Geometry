#ifndef __geom_h
#define __geom_h
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <vector>

using namespace std; 


typedef struct _point2d {
  double x,y; 
  bool isReflex;
  bool isLeft;
} point2d;

typedef struct _guard {
  point2d pos;
  GLfloat r;
  GLfloat g;
  GLfloat b;
  double dir;
} guard;

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

class HorizontalComparator {
public:
    bool operator()(point2d p1, point2d p2)
    {
        if (p1.x < p2.x) {
            return true;
        }
        if (p1.x > p2.x) {
          return false;
        }
        else if (p1.y > p2.y) {
          return true;
        }
        return false;
    }
};

class VerticalComparator {
public:
    bool operator()(point2d p1, point2d p2)
    {
        if (p1.y < p2.y) {
            return true;
        }
        if (p1.y > p2.y) {
          return false;
        }
        else if (p1.x > p2.x) {
          return true;
        }
        return false;
    }
};



/* returns 2 times the signed area of triangle abc. The area is
   positive if c is to the left of ab, 0 if a,b,c are collinear and
   negative if c is to the right of ab
 */
int signed_area2D(point2d a, point2d b, point2d c); 


/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2d p, point2d q, point2d r);


/* return 1 if c is  strictly left of ab; 0 otherwise */
int left_strictly (point2d a, point2d b, point2d c); 


/* return 1 if c is left of ab or on ab; 0 otherwise */
int left_on(point2d a, point2d b, point2d c); 



// compute the convex hull 
void graham_scan(vector<point2d>& pts, vector<point2d>& hull);
  

#endif
