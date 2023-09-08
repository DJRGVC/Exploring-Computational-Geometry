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
  // comparator
  bool operator==(const _point2d& p) const {
    return ((x == p.x && y == p.y) || (x == p.y && y == p.x));
  }
} point2d;

typedef struct _guard {
  point2d pos;
  GLfloat r;
  GLfloat g;
  GLfloat b;
  double dir;
} guard;

typedef struct _robot {
  point2d pos;
  GLfloat r;
  GLfloat g;
  GLfloat b;
} robot;

typedef struct _node {
  point2d* pos;
  double distance;
  bool onPolygon;

  // constructor
  _node(point2d* p, double d) {
    pos = p;
    distance = d;
    onPolygon = false;
  }
  _node(point2d* p, double d, bool onP) {
    pos = p;
    distance = d;
    onPolygon = onP;
  }

  // comparator
  bool operator==(const _node& n) const {
    return (pos == n.pos);
  }
} node;

typedef struct _goal {
  point2d pos;
  GLfloat r;
  GLfloat g;
  GLfloat b;
} goal;

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

class MyHashFunction {
public:
  // id is returned as hash function
  size_t operator()(const node& t) const {
    return 1000 * t.pos->x + t.pos->y;
    }
};

class PointerHashFunction {
public:
  // id is returned as hash function
  size_t operator()(const point2d* t) const {
    return (unsigned long long)t;
  }
};

// compare node objects, used for priority queue
class NodeComparator {
public:
  bool operator()(const node& n1, const node& n2) const {
    return n1.distance > n2.distance;
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
