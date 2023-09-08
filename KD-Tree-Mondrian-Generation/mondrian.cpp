/**
 * Mondrian Class
 *
 * @author Brian Liu and Daniel Grant
 *
 * @date 22/2/2023
 * @version 1.0
*/

#include "mondrian.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <random>

using namespace std;

// number of points
int numPoints;


// vector of points
vector<point2D> points;

// Kdtree instance
Kdtree* tree = new Kdtree();

/**
 * Make a new tree
 */
Kdtree makeNewTree(int pointsNum) {
  // randomize the seed for the random number generator using hash function 
  // safe print hello
  int seed = time(NULL);
  seed = (int) ((1103515245 * ((unsigned int)seed) + 12345) & 0x7fffffffUL);
  srand(seed);

  vector<point2D> pointsForPaint;

  uniform_real_distribution<double> unif(0, MAX_HEIGHT);

  default_random_engine re;

  // generate random points
  for (int i = 0; i < pointsNum; i++) {
    bool isDupe;
    isDupe = false;
    point2D p;
    p.x = unif(re);
    p.y = unif(re);

    // to make sure that the points are not the same
    for (int j = 0; j < i; j++) {
      if (p.x == pointsForPaint[j].x && p.y == pointsForPaint[j].y) {
        i--;
        isDupe = true;
      }
    }
    if (isDupe) {
      continue;
    }
    pointsForPaint.push_back(p);
  }

  // create the tree
  tree->initialize(pointsForPaint);
  return *tree;
}

/**
 * Sort lexicographically by x coordinate
*/
bool bottom_top(const point2D &a, const point2D &b) {
  if (a.y == b.y) {
    return a.x <= b.x;
  }
  return a.y < b.y;
}

/**
 * Sort lexicographically by y coordinate
*/
bool left_right (const point2D &a, const point2D &b) {
  if (a.x == b.x) {
    return a.y <= b.y;
  }
  return a.x < b.x;
}

// constructor for Kdtree class
Kdtree::Kdtree() 
{
  height = 0;
  count = 0;
  root = NULL;
}

// print out points
void printPoints(vector<point2D> points) {
  for (int i = 0; i < numPoints; i++) {
    printf("Point #1: (%f, %f)\n", points[i].x, points[i].y);
  }
}

/**
 * Recursive method to build the tree
*/
void buildTree(TreeNode *node, vector<point2D> pByX, vector<point2D> pByY) 
{

  if (pByX.size() == 1) {
    node->p = pByX[0];
    node->left = NULL;
    node->right = NULL;
    node->type = LEAF;
    return;
  }
  
  vector<point2D> P1_x;
  vector<point2D> P1_y;

  vector<point2D> P2_x;
  vector<point2D> P2_y;

  point2D median;

  if (node->type == VERTICAL) {
    // find the median point by x coordinate
    if (pByX.size() % 2 == 0) {
      median = pByX[(pByX.size() / 2) - 1];
    } else {
      median = pByX[pByX.size() / 2];
    }

    // find new set of points for left subtree and right subtree
    for (int i = 0; i < pByX.size(); i++) {
      if (left_right(pByX[i], median)) {
        P1_x.push_back(pByX[i]);
      } else {
        P2_x.push_back(pByX[i]);
      }
      if (left_right(pByY[i], median)) {
        P1_y.push_back(pByY[i]);
      } else {
        P2_y.push_back(pByY[i]);
      }
    }

    // update the node
    node->p = median;
    TreeNode* left = new TreeNode();
    node->left = left;
    left->type = HORIZONTAL;
    TreeNode* right = new TreeNode();
    node->right = right;
    right->type = HORIZONTAL;

    // recursively build the left and right subtrees
    buildTree(node->left, P1_x, P1_y);
    buildTree(node->right, P2_x, P2_y);
  }

  if (node->type == HORIZONTAL) {
    // find the median point by y coordinate
    if (pByY.size() % 2 == 0) {
      median = pByY[(pByY.size() / 2) - 1];
    } else {
      median = pByY[pByY.size() / 2];
    }

    // find new set of points for left subtree and right subtree
    for (int i = 0; i < pByY.size(); i++) {
      if (bottom_top(pByY[i], median)) {
        P1_y.push_back(pByY[i]);
      } else {
        P2_y.push_back(pByY[i]);
      }
      if (bottom_top(pByX[i], median)) {
        P1_x.push_back(pByX[i]);
      } else {
        P2_x.push_back(pByX[i]);
      }
    }

    // update the node
    node->p = median;
    TreeNode* left = new TreeNode();
    node->left = left;
    left->type = VERTICAL;
    TreeNode* right = new TreeNode();
    node->right = right;
    right->type = VERTICAL;

    // recursively build the left and right subtrees
    buildTree(node->left, P1_x, P1_y);
    buildTree(node->right, P2_x, P2_y);
  }



}

/**
 * Initialize the tree
 */
void Kdtree::initialize(vector<point2D> points) 
{
  // sort the points
  // by y
  stable_sort(points.begin(), points.end(), bottom_top);
  vector<point2D> pByY = points;

  // by x
  stable_sort(points.begin(), points.end(), left_right);
  vector<point2D> pByX = points;



  // create the root node as the middle of the points sorted by x
  if (pByX.size() % 2 == 0) {
    root = new TreeNode(pByX[points.size() / 2 - 1], VERTICAL);
  } else {
    root = new TreeNode(pByX[points.size() / 2], VERTICAL);
  }

  // call the recursive method to create the tree
  buildTree(root, pByX, pByY);
}

/**
 * Tree node constructor from point and node type
 */
TreeNode::TreeNode(point2D point, NodeType nodeType) 
{
  p = point;
  type = nodeType;
  left = NULL;
  right = NULL;
}

/**
 * Tree node constructor (default)
 */
TreeNode::TreeNode()
{
  left = NULL;
  right = NULL;
  type = VERTICAL;
  p = point2D();
}

/**
 * function to print tree
 */
void printTree(TreeNode *node) {
  if (node == NULL) {
    return;
  }
  if (node->type == LEAF) {
    printf("Node type: LEAF");
  }
  if (node->type == VERTICAL) {
    printf("Node type: VERTICAL");
  }
  if (node->type == HORIZONTAL) {
    printf("Node type: HORIZONTAL");
  }
  printf(", Point: (%f, %f)", node->p.x, node->p.y);
  printf("\n");
  printTree(node->left);
  printTree(node->right);
}

/**
 * function to print node
 */
void printNode(TreeNode *node) {
  if (node == NULL) {
    return;
  }
  if (node->type == LEAF) {
    printf("Node type: LEAF");
  }
  if (node->type == VERTICAL) {
    printf("Node type: VERTICAL");
  }
  if (node->type == HORIZONTAL) {
    printf("Node type: HORIZONTAL");
  }
  printf(", Point: (%f, %f)", node->p.x, node->p.y);
  printf("\n");
}

