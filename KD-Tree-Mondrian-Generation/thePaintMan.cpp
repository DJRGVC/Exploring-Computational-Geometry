/**
 * @author: Brian Liu and Daniel Grant
 *
 * @date: 2/22/2023
 * @version: 1.0
 */
#define GL_SILENCE_DEPRECATION

#include "mondrian.h"

#include <random>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

using namespace std;

// Global variables

// window size 
const int WINDOWSIZE = 1000;

// edge width
const int EDGE_WIDTH = 4000;

// Colors for the painting
GLfloat yellow[3] = {(float)251/255, (float)202/255, (float)0/255}; 
GLfloat blue[3] = {(float)27/255, (float)79/255, (float)154/255}; 
GLfloat black[3] = {(float)2/255, (float)2/255, (float)2/255}; 
GLfloat white[3] = {(float)255/255, (float)255/255, (float)255/255}; 
GLfloat red[3] = {(float)188/255, (float)3/255, (float)25/255};

// add all colors to a vector
vector<GLfloat*> colors;

// Global variables
// The kd tree
Kdtree myTree;


// Function prototypes
void display(void);
void keypress(unsigned char key, int x, int y);
void drawMondrian(TreeNode *node, double bounds[4]);

// number of points (user inputted)
int NPOINTS;





/**
 * Main function
 */
int main(int argc, char** argv) {

  //read number of points from user
  if (argc!=2) {
    printf("usage: viewPoints <nbPoints>\n");
    exit(1); 
  }
  NPOINTS = atoi(argv[1]); 

  // create the tree
  myTree = makeNewTree(NPOINTS);

  fflush(stdout); 


 
  //start the rendering 
  /* initialize GLUT  */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  
  
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display); 
  glutKeyboardFunc(keypress);


  /* init GL */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);   


  /* give control to event handler */
  glutMainLoop();

  return 0;
}


/* ****************************** */
/* This is the main render function. We registered this function to be
   called by GL to render the window. 
 */
void display(void) {

  // Clear the screen
  glClear(GL_COLOR_BUFFER_BIT);

  //clear all modeling transformations
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity(); 



  /* The default GL window is [-1,1]x[-1,1]x[-1,1] with the origin in
     the center. The camera is at (0,0,0) looking down negative
     z-axis.  

     The points are in the range (0,0) to (WINSIZE,WINSIZE), so they
     need to be mapped to [-1,1]x [-1,1] */
  //then scale the points to [0,2]x[0,2]
  glScalef(2.0/MAX_WIDTH, 2.0/MAX_WIDTH, 1.0);  
  //first translate the points to [-WINDOWSIZE/2, WINDOWSIZE/2]
  glTranslatef(-MAX_WIDTH/2, -MAX_HEIGHT/2, 0.0);
  
  // print the location of bottom left of window
 
  //draw the kd-tree
  // create int array representing xmin, xmax, ymin, ymax
  double bounds[4] = {0, MAX_WIDTH, 0, MAX_HEIGHT};

  colors.push_back(red);
  colors.push_back(blue);
  colors.push_back(black);
  colors.push_back(white);
  colors.push_back(white);
  colors.push_back(white);
  colors.push_back(white);
  colors.push_back(yellow);

  drawMondrian(myTree.root, bounds);

  // Flush the pipeline, and swap the buffers
  glFlush();
}

/**
 * This function is called whenever a key is pressed
 */
void keypress(unsigned char key, int x, int y) {
  switch(key) {
  case 'q':
    exit(0);
    break;
  } //switch (key)
}

/**
 * This function draws the mondrian painting
 */
void drawMondrian(TreeNode* node, double bounds[4]) {
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  // if node is a leaf, draw a rectangle
  if (node->type == LEAF) {
    // draw interior
    // choose random color from vector
    int randColor = rand() % colors.size();
    glColor3fv(colors[randColor]);
    if (bounds[1] - bounds[0] > 6000 && bounds[3] - bounds[2] > 6000) {
      glBegin(GL_POLYGON);
        glVertex2f(bounds[0] + EDGE_WIDTH, bounds[2] + EDGE_WIDTH);
        glVertex2f(bounds[1] - EDGE_WIDTH, bounds[2] + EDGE_WIDTH);
        glVertex2f(bounds[1] - EDGE_WIDTH, bounds[3] - EDGE_WIDTH);
        glVertex2f(bounds[0] + EDGE_WIDTH, bounds[3] - EDGE_WIDTH);
      glEnd();
    }
    return;
  }

  // if vertical, split bounds into two rectangles (left and right)
  if (node->type == VERTICAL) {
    double leftBounds[4] = {bounds[0], node->p.x, bounds[2], bounds[3]};
    double rightBounds[4] = {node->p.x, bounds[1], bounds[2], bounds[3]};
    drawMondrian(node->left, leftBounds);
    drawMondrian(node->right, rightBounds);
  }

  // if horizontal, split bounds into two rectangles (top and bottom)
  if (node->type == HORIZONTAL) {
    double topBounds[4] = {bounds[0], bounds[1], node->p.y, bounds[3]};
    double bottomBounds[4] = {bounds[0], bounds[1], bounds[2], node->p.y};
    drawMondrian(node->right, topBounds);
    drawMondrian(node->left, bottomBounds);
  }
  
}
