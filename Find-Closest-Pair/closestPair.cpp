/*
 * Project 1: Closest Pair
 * Authors: Brian Liu and Daniel Grant
 * Date: 2/5/2023
 * Description: This program finds the closest pair of points in a set of points
 *              using a naive approach and the gridding algorithm. */

#include "closestPair.h"
#include <string>
#include <map>
#include <ctime>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cstdlib>
#include <vector>

using namespace std; 

// for naive method
vector<Point> points;

// for gridding method
vector<Point> grid[MAP_LENGTH][MAP_WIDTH];

// number of points
int numPoints;

// number of grid divisions
int gridDivisions;

// random number
int rand();

/*
 * Function: main
 * Usage: main(argc, argv);
*/
int main(int argc, char** argv) {
  
  // make sure user inputs 1 or 2 arguments (second one is optional for grinding)
  if (argc < 2 || argc > 3) {
    printf("Usage: %s <Number of Dots> [Number of Grid Divisions]\n", argv[0]);
    exit(1);
  }


  // checking if the number of points is valid
  if (atoi(argv[1]) <= 0) {
    printf("Number of Dots must be a positive integer, not %d\n", atoi(argv[1]));
    exit(1);
  }

  // checking if the number of grid divisions is valid
  if (argc == 3 && atoi(argv[2]) <= 0) {
    printf("Number of Grid Divisions must be a positive integer, not %d\n", atoi(argv[2]));
    exit(1);
  }

  // set grid divisions
  if (argc == 2) {
    // if there is only one argument, set the number of grid divisions square root of the number of points (explained further in README)
    gridDivisions = (int) sqrt(atoi(argv[1]));
  } else {
    gridDivisions = atoi(argv[2]);
  }
  numPoints = atoi(argv[1]);

  // randomize the seed for the random number generator using hash function 
  int seed = time(NULL);
  seed = (int) ((1103515245 * ((unsigned int)seed) + 12345) & 0x7fffffffUL);
  srand(seed);

  // boolean to check if the point is already in the set
  bool isDupe;

  // generate random points
  for (int i = 0; i < numPoints; i++) {
    isDupe = false;
    Point p;
    p.x = (int)rand() % ((int)MAP_LENGTH);
    p.y = (int)rand() % ((int)MAP_WIDTH);

    // to make sure that the points are not the same
    for (int j = 0; j < i; j++) {
      if (p.x == points[j].x && p.y == points[j].y) {
        i--;
        isDupe = true;
      }
    }
    if (isDupe) {
      continue;
    }
    addToGrid(p);
    points.push_back(p);
  }

  // Naive approach
  //begin time measurement
  clock_t start = clock();

  printf("\nNaive Approach: ");
  ClosestPair closestPairNaive = closest_naive();

  clock_t end = clock();
  if (closestPairNaive.p1.x == -1) {
    printf("Time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
  } else {
    //end time measurement

    // print the closest pair, and the time it took to find it, and the 
    // distance between them
    printf("\nClosest Pair: (%f, %f) and (%f, %f)\n", closestPairNaive.p1.x, closestPairNaive.p1.y, closestPairNaive.p2.x, closestPairNaive.p2.y);
    printf("Distance: %f\n", closestPairNaive.distance);
    printf("Time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
  }


  // Gridding approach:
  //begin time measurement
  start = clock();

  printf("\nGridding Approach: ");
  ClosestPair closestPairGrid = closest_grid();
  end = clock();
  if (closestPairGrid.p1.x == -1) {
    printf("Time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
  } else {
    //end time measurement

    // print the closest pair, and the time it took to find it, and the 
    // distance between them
    printf("\nClosest Pair: (%f, %f) and (%f, %f)\n", closestPairGrid.p1.x, closestPairGrid.p1.y, closestPairGrid.p2.x, closestPairGrid.p2.y);
    printf("Distance: %f\n", closestPairGrid.distance);
    printf("Time: %f\n", (double)(end - start) / CLOCKS_PER_SEC);
    
 }
}

/*
 * Function: distance
 * Usage: distance(p1, p2);
*/

/*
 * Function: addToGrid
 * Usage: addToGrid(p);
*/
void addToGrid(Point p) {
  int x = p.x / (MAP_LENGTH / gridDivisions);
  int y = p.y / (MAP_WIDTH / gridDivisions);
  // print point and grid location
  grid[x][y].push_back(p);
}

double distance(Point p1, Point p2) {
  return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

/*
 * Function: naiveClosestPair
 * Description: This function finds the closest pair of points in a set of points
 *             using a naive approach.
*/
ClosestPair closest_naive() {
  // closestPair object initialized
  ClosestPair closestPairNaive;
  Point invalidPoint = {-1, -1};
  closestPairNaive.p1 = invalidPoint;
  closestPairNaive.p2 = invalidPoint;
  closestPairNaive.distance = distance(closestPairNaive.p1, closestPairNaive.p2);

  // base case if there are less than 2 points
  if (points.size() < 2) {
    printf("Not enough points to find a pair\n");
    return closestPairNaive;
  }

  closestPairNaive.p1 = points[0];
  closestPairNaive.p2 = points[1];
  closestPairNaive.distance = distance(points[0], points[1]);

  // now, loop through all the points, and find the closest pair
  for (int i = 0; i < numPoints; i++) {
    for (int j = i + 1; j < numPoints; j++) {
      // calculate the distance between the two points
      double dist = distance(points[i], points[j]);
      if (dist < closestPairNaive.distance) {
        closestPairNaive.distance = dist;
        closestPairNaive.p1 = points[i];
        closestPairNaive.p2 = points[j];
      }
    }
  }
  return closestPairNaive;
}

/*
 * Function: closest_grid
 * Description: This function finds the closest pair of points in a set of points
 *             using the divide and conquer algorithm.
*/
ClosestPair closest_grid() {

  // closestPair object initialized
  ClosestPair closestPairGrid;
  Point invalidPoint = {-1, -1};
  closestPairGrid.p1 = invalidPoint;
  closestPairGrid.p2 = invalidPoint;
  closestPairGrid.distance = distance(points[0], points[1]);

  // base case if there are less than 2 points
  if (points.size() < 2) {
    printf("Not enough points to find a pair\n");
    return closestPairGrid;
  }

  // closestPair object initialized
  closestPairGrid.p1 = points[0];
  closestPairGrid.p2 = points[1];
  closestPairGrid.distance = distance(points[0], points[1]);



  // loop through 2d array of points
  for (int i = 0; i < gridDivisions; i++) {
    for (int j = 0; j < gridDivisions; j++) {
      // if there are points in the grid
      if (grid[i][j].size() > 0) {
        // loop through the points in the grid
        for (int k = 0; k < grid[i][j].size(); k++) {
          // loop through the points in the grid
          for (int l = k + 1; l < grid[i][j].size(); l++) {
            // calculate the distance between the two points
            double dist = distance(grid[i][j][k], grid[i][j][l]);
            if (dist < closestPairGrid.distance) {
              closestPairGrid.distance = dist;
              closestPairGrid.p1 = grid[i][j][k];
              closestPairGrid.p2 = grid[i][j][l];
            }
            if (dist == 1) {
              return closestPairGrid;
            }
          }
          ClosestPair closestAdj = closest_around(i, j, closestPairGrid.distance, grid[i][j][k]);
          if (closestAdj.distance < closestPairGrid.distance) {
            closestPairGrid.distance = closestAdj.distance;
            closestPairGrid.p1 = closestAdj.p1;
            closestPairGrid.p2 = closestAdj.p2;
          }
        }
      }
    }
  }
  return closestPairGrid;
}

/*
 * Function: closest_around
 * Description: This function finds the closest pair of points in adjacent grids within the minumum distance away from the current grid
*/
ClosestPair closest_around(int x, int y, double minDist, Point p) {
// closestPair object initialized
 ClosestPair closestPairGrid;
 closestPairGrid.p1 = points[0];
 closestPairGrid.p2 = points[1];
 closestPairGrid.distance = minDist;

 // find which grids are within minDist away from point p
 int xMin = (p.x - minDist) / (MAP_LENGTH / gridDivisions);
 int xMax = (p.x + minDist) / (MAP_LENGTH / gridDivisions);
 int yMin = (p.y - minDist) / (MAP_WIDTH / gridDivisions);
 int yMax = (p.y + minDist) / (MAP_WIDTH / gridDivisions);

 // make sure the grid locations are within the bounds of the grid
 xMin = max(xMin, 0);
 xMax = min(xMax, gridDivisions - 1);
 yMin = max(yMin, 0);
 yMax = min(yMax, gridDivisions - 1);
 
 // loop through the grids within minDist away from point p
 for (int i = xMin; i <= xMax; i++) {
   for (int j = yMin; j <= yMax; j++) {
     // if there are points in the grid
     if ((i != x || j != y) && grid[i][j].size() > 0) {
       // loop through the points in the grid
       for (int k = 0; k < grid[i][j].size(); k++) {
         if (abs(grid[i][j][k].x - p.x) < minDist && abs(grid[i][j][k].y - p.y) < minDist) {
         // calculate the distance between the two points
           double dist = distance(p, grid[i][j][k]);
           if (dist < closestPairGrid.distance) {
             closestPairGrid.distance = dist;
             closestPairGrid.p1 = p;
             closestPairGrid.p2 = grid[i][j][k];
           }
         }
       }
     }
   }
 }
 return closestPairGrid;
}





