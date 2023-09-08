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
  bool isJoint;
  // comparator

  // constructor
  _point2d() {
    isJoint = false;
  }

  _point2d(double x, double y) {
    this->x = x;
    this->y = y;
    isJoint = false;
  }

  bool operator==(const _point2d& p) const {
    return ((x == p.x && y == p.y) || (x == p.y && y == p.x));
  }

} point2d;

typedef struct _anglepoint {
  vector<double>* angles;
  point2d* pos;
  // constructor
  _anglepoint(vector<double> *a, point2d* p) {
    angles = a;
    pos = p;
  } 
  _anglepoint() {
    angles = new vector<double>();
    pos = new point2d();
  } 
} anglepoint;

typedef struct _rrtNode {
  anglepoint* orientation;
  struct _rrtNode* parent;
  bool partOfPath;

  // constructor
  _rrtNode(anglepoint* a, struct _rrtNode* par) {
    orientation = a;
    parent = par;
    partOfPath = false;
  }
  _rrtNode() {
    orientation = new anglepoint();
    parent = NULL;
    partOfPath = false;
  }
} rrtNode;


typedef struct _node {
  point2d* pos;
  double distance;

  // constructor
  _node(point2d* p, double d) {
    pos = p;
    distance = d;
  }
} node;

typedef struct _goal {
  point2d pos;
  GLfloat r;
  GLfloat g;
  GLfloat b;
} goal;

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

#endif
