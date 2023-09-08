#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>

using namespace std; 

/* **************************************** */
/* returns the signed area of triangle abc. The area is positive if c
   is to the left of ab, and negative if c is to the right of ab
 */
int signed_area2D(point2d a, point2d b, point2d c) {

  return 1; 
}



/* **************************************** */
/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2d p, point2d q, point2d r) {
  
  return 1; 
}



/* **************************************** */
/* return 1 if c is  strictly left of ab; 0 otherwise */
int left_strictly(point2d a, point2d b, point2d c) {
  
  return 1; 
}


/* return 1 if c is left of ab or on ab; 0 otherwise */
int left_on(point2d a, point2d b, point2d c) {

  return 1; 
}



// compute the convex hull 
void graham_scan(vector<point2d>& pts, vector<point2d>& hull ) {

  printf("hull2d (graham scan): start\n"); 
  hull.clear(); //should be empty, but cleara it to be safe

  //just for fun: at the moment the hull is the bounding box of pts!
  //erase this and insert your code instead
  int x1, x2, y1, y2;
  if (pts.size() > 0) {
    x1 = x2 = pts[0].x;
    y1 = y2 = pts[0].y;
    
    for (int i=1; i< pts.size(); i++) {
      if (pts[i].x < x1) x1 = pts[i].x;
      if (pts[i].x > x2) x2 = pts[i].x;
      if (pts[i].y < y1) y1 = pts[i].y;
      if (pts[i].y > y2) y2 = pts[i].y;
    }
    point2d p1 = {x1,y1}, p2 = {x2, y1}, p3 = {x2, y2}, p4 = {x1, y2}; 
    hull.push_back(p1);
    hull.push_back(p2);
    hull.push_back(p3);
    hull.push_back(p4);
  }
  
  printf("hull2d (graham scan): end\n"); 
  return; 
}

