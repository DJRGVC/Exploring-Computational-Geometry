#include "geom.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <queue>
#include <stack>
#include <math.h>     /* atan2 */

/**
 * @author:  Brian Liu, Daniel Grant
 * @date:    2023-02-11
 * @version: 1.0
 * @brief:   Implementation of the Graham Scan algorithm
 */

using namespace std; 

/* **************************************** */
/* returns the signed area of triangle abc. The area is positive if c
   is to the left of ab, and negative if c is to the right of ab
 */
int signed_area2D(point2d a, point2d b, point2d c) {
  return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x);
}

/* **************************************** */
/* return 1 if p,q,r collinear, and 0 otherwise */
int collinear(point2d p, point2d q, point2d r) {
  return (signed_area2D(p,q,r) == 0);
}


/* **************************************** */
/* return 1 if c is  strictly left of ab; 0 otherwise */
int left_strictly(point2d a, point2d b, point2d c) {
  return (signed_area2D(a,b,c) > 0);
}

/* **************************************** */
/* return 1 if c is left of ab or on ab; 0 otherwise */
int left_on(point2d a, point2d b, point2d c) {
  return (signed_area2D(a,b,c) >= 0);
}

/**
 * @brief:   Implementation of the Graham Scan algorithm
 * @param:   points - a vector of points
 * @param:   hull - a vector of points that will be filled with the convex hull
 * @return:  void
 */
void graham_scan(vector<point2d>& pts, vector<point2d>& hull ) {

  printf("hull2d (graham scan): start\n"); 

  // initialize priority queue of radial_point2d
  priority_queue<radial_point2d> sortedPoints;


  hull.clear(); // should be empty, but clear it to be safe

  // first find smallest y-coordinate point
  // this must be in the convex hull
  point2d p0 = pts[0];
  for (int i = 1; i < pts.size(); i++) {
    if (pts[i].y < p0.y) {
      p0 = pts[i];
    }
    if (pts[i].y == p0.y && pts[i].x < p0.x) {
      p0 = pts[i];
    }
  }

  // sort radially around the point with other points
  // if two points have same angle, then the one with smaller distance
  // put sorted point in a stack
  sortRadially(p0, pts, sortedPoints);
  
  // if there are less than 3 points, can not make a convex hull
  if (sortedPoints.size() < 3) {
    printf("hull2d (graham scan): error: less than 3 points\n");
    return;
  }

  point2d p1 = sortedPoints.top().p;
  sortedPoints.pop();

  // push the first three points onto the stack
  hull.push_back(p0);
  hull.push_back(p1);

  // loop through the rest of the points
  // if the next point is left of the line formed by the top two points
  // then push it onto the stack
  // otherwise, pop the top point off the stack and try again
  // until the next point is left of the line formed by the top two points
  // then push the next point onto the stack
  // this will ensure that the points on the stack form a convex hull
  while (!sortedPoints.empty()) {
    p0 = hull[hull.size() - 2];
    p1 = hull[hull.size() - 1];
    point2d p2 = sortedPoints.top().p;
    sortedPoints.pop();
    if (left_strictly(p0, p1, p2)) {
      hull.push_back(p2);
    } else {
      if (collinear(p0, p1, p2)) {
        hull.pop_back();
        hull.push_back(p2);
      } else {
        while (!left_strictly(p0, p1, p2)) {
          hull.pop_back();
          p0 = hull[hull.size() - 2];
          p1 = hull[hull.size() - 1];
        }
        hull.push_back(p2);
      }
    }
  }

  printf("hull2d (graham scan): end\n"); 
  return; 
}

/**
 * @brief:   Finds the distance between two points
 * @param:   p1 - the first point
 * @param:   p2 - the second point
 * @return:  the distance between the two points
 */
int distance(point2d p1, point2d p2) {
  return (p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y);
}

/**
 * @brief 
 * 
 * @param:   p0 - the point to sort around
 * @param:   pts - the vector of points to sort
 * @param:   sortedPoints - the priority queue to put the sorted points in
 * @return:  void
 */
void sortRadially(point2d p0, vector<point2d>& pts, priority_queue<radial_point2d>& sortedPoints) {
  // sort points radially into a priority queue
  for (int i = 0; i < pts.size(); i++) {
    if (pts[i].x == p0.x && pts[i].y == p0.y) {
      continue;
    }
    radial_point2d rp;
    rp.p = pts[i];
    rp.angle = atan2(pts[i].y - p0.y, pts[i].x - p0.x);
    rp.distance = distance(p0, pts[i]);
    sortedPoints.push(rp);
  }
}

