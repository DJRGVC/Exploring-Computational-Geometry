/*
 * closestPair.h
 *
 * Created on: 2023-1-26
 * Author: Brian Liu & Daniel Grant
*/

#include <map>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cstdlib>
#include <vector>

#ifndef CLOSESTPAIR_H

using namespace std;

// point struct
typedef struct {
    double x;
    double y;
} Point;

// closest pair struct
typedef struct {
    Point p1;
    Point p2;
    double distance;
} ClosestPair;

// constants for map size (can be changed)
const int MAP_LENGTH = 1000;
const int MAP_WIDTH = 1000;

// main function for calling naive and divide and conquer algorithms
int main(int argc, char** argv);

// find closest point in adjacent grids within distance of shortest pair
ClosestPair closest_around(int x, int y, double minDist, Point p);

// add points to the grid
void addToGrid(Point p);

// find the distance between two points
double distance(Point p1, Point p2);

// find the closest pair of points, brute force
ClosestPair closest_naive();

// find the closest pair of points, divide and conquer
ClosestPair closest_grid();

#endif /* CLOSESTPAIR_H_ */
