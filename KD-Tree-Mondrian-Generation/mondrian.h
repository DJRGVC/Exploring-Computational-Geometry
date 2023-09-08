#ifndef __mondrian_h
#define __mondrian_h

#include <iostream>
#include <string>

using namespace std;

// CONSTANTS
const int MAX_WIDTH = 1000000;
const int MAX_HEIGHT = 1000000;

/** 
 * point2d struct that sorts the x and y coordinates of a point
*/
typedef struct _point2D {
    double x;
    double y;
} point2D;

/**
 * Tree node type enum
*/
enum NodeType {
    HORIZONTAL = 1,
    VERTICAL = 2,
    LEAF = 3
};

/**
 * Tree node class
*/
class TreeNode {
  public: 
     point2D p; 
     NodeType type; 
     TreeNode *left, *right; 
     TreeNode();
     TreeNode(point2D point, NodeType nodeType); 
};

/**
 * kd tree class
*/
class Kdtree {
  public:
     TreeNode *root; 
     int count;  //number of leaves/points in the tree
     int height;     
     Kdtree();
     void initialize(vector<point2D> points);
};


/**
 * build the tree from the given points, recursively
*/
void buildTree(TreeNode *node, vector<point2D> left_right, vector<point2D> top_bottom);

/**
 * print points
*/
void printPoints(vector<point2D> points);

/**
 * print the tree
 */
void printTree(TreeNode *node);

/**
 * print node
 */
void printNode(TreeNode *node);

/**
 * make new tree
 */
Kdtree makeNewTree(int nbPoints);

#endif
