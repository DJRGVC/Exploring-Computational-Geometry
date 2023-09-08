/**
 * @file vis.cpp
 * @author Brian Liu and Daniel Grant
 * 
 * 
 */
#include "geom.h" 

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif


#include <vector> 
#include <iostream>

#include <queue>
#include <unordered_map>
#include <unordered_set>

using namespace std; 



GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat cyan[3] = {0.0, 1.0, 1.0};



/* global variables */
const int WINDOWSIZE = 750; 

enum mode { RUN_DIJKSTRA, FINISH_POLYGON, NEW_POLYGON, CALC_VISIBLE_LINES, DRAW_OBSTACLES, RUN, DRAW_ROBOT, DRAW_GOAL, DRAW_VISIBLE_LINES };

mode curMode = NEW_POLYGON;


//the current polygon 
vector<point2d> poly;

//the current user polygon
vector<point2d> polygon;

// vector of polygon obstacles

vector<point2d> obstacle1;
vector<point2d> obstacle2;
vector<vector<point2d> > obstacles;

// the robot
robot rob;

// the current goal
goal gol;

// hash map of where key: point2d, value: vector of nodes
// use custom hashing function in geom.h
unordered_map<point2d*, vector<node>, PointerHashFunction> nodeMap;
unordered_set<point2d*, PointerHashFunction> processedPoints;



// for storing shortest path
vector<point2d*> shortestPath;
double shortestPathLength;


// unordered_map<point2d*, vector<node>, MyHashFunction> nodeMap;
// unordered_set<point2d*, MyHashFunction> visitedNodes;



//coordinates of the last mouse click
double mouse_x=-10, mouse_y=-10;  //initialized to a point outside the window
                                  //
double lastx = -1, lasty = -1;


/* forward declarations of functions */
void display(void);
void keypress(unsigned char key, int x, int y);
void mousepress(int button, int state, int x, int y);
void timerfunc(); 

void initialize_polygon(); 
void print_polygon(vector<point2d>& poly); 
void populate_obstacle(vector<point2d>& poly, int divisions, int sector, int numVertices);



void populate_obstacles_hardcoded();
void draw_obstacles();
void calc_visible_lines();
void draw_visible_lines();



/* ****************************** */
int main(int argc, char** argv) {

  cout << "Begin by creating a polygon.\n";
  polygon.clear();

  /* initialize GLUT  */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display); 
  glutKeyboardFunc(keypress);
  glutMouseFunc(mousepress); 
  //glutIdleFunc(timerfunc); //register this if you want it called aat every fraame

  /* init GL */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);
  
  /* give control to event handler */
  glutMainLoop();
  return 0;
}

void populate_obstacles_hardcoded() {

  point2d p1, p2, p3, p4, p5, p6, p7, p8;

  p1.x = 100; p1.y = 350;
  p2.x = 275; p2.y = 100;
  p3.x = 160; p3.y = 350;
  p4.x = 275; p4.y = 600;
  p5.x = 400; p5.y = 300;
  p6.x = 600; p6.y = 150;
  p7.x = 600; p7.y = 650;
  p8.x = 300; p8.y = 650;

  obstacle1.push_back(p1);
  obstacle1.push_back(p2);
  obstacle1.push_back(p3);
  obstacle1.push_back(p4);

  obstacle2.push_back(p5);
  obstacle2.push_back(p6);
  obstacle2.push_back(p7);
  obstacle2.push_back(p8);
}

void draw_obstacles() {
  glColor3fv(blue);
  for (int i = 0; i < obstacles.size(); i++) {
    // draw lines for each obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      if (j == obstacles[i].size() - 1) {
        glBegin(GL_LINES);
        glVertex2f(obstacles[i][j].x, obstacles[i][j].y); 
        glVertex2f(obstacles[i][0].x, obstacles[i][0].y); 
        glEnd();
      } else {
        glBegin(GL_LINES);
        glVertex2f(obstacles[i][j].x, obstacles[i][j].y); 
        glVertex2f(obstacles[i][j+1].x, obstacles[i][j+1].y); 
        glEnd();
      }
    }
  }
}


bool parallel_intersect(point2d p1, point2d p2, point2d p3, point2d p4) {
  // return true if coincident, otherwise false
  if (p1.x == p3.x) {
    if (p1.y == p3.y) {
      if (p2.x == p4.x) {
        if (p2.y == p4.y) {
          return true;
        }
      } 
    }
  }
  return false;
}


char seg_seg_intersect(point2d p1, point2d p2, point2d p3, point2d p4, point2d* p5) {
  double s, t; // parametric value of first and second line respectively
  double num, denom; // numerator and denominator of the equations
  char code = '?';

  // finding the determinant
  denom = p1.x * (p4.y - p3.y) + p2.x * (p3.y - p4.y) + p4.x * (p2.y - p1.y) + p3.x * (p1.y - p2.y);

  // determinant is zero, lines are parallel or coincident
  if (denom == 0) {
    return parallel_intersect(p1, p2, p3, p4);
  }

  // find determinant of two line segments and intersection point
  num = p1.x * (p4.y - p3.y) + p3.x * (p1.y - p4.y) + p4.x * (p3.y - p1.y);

  // coincident or collinear
  if ((num == 0) || (num == denom)) {
    code = 'v';
  }

  // find parametric value of intersection point
  s = num / denom;

  num = -1 * (p1.x * (p3.y - p2.y) + p2.x * (p1.y - p3.y) + p3.x * (p2.y - p1.y));

  // coincident or collinear
  if ((num == 0) || (num == denom)) {
    code = 'v';
  }

  t = num / denom;

// should be between 0 and 1 to indicate intersection at some point between the two line segments
  if ((0 < s) && (s < 1) && (0 < t) && (t < 1)) {
    code = '1';
  } else if (0 > s || s > 1 || 0 > t || t > 1) {
    code = '0';
  }

  p5->x = p1.x + s * (p2.x - p1.x);
  p5->y = p1.y + s * (p2.y - p1.y);

  return code;
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


bool is_simple_polygon() {
  point2d *p5 = new point2d ;
  p5->x = -1;
  p5->y = -1;
  for (int i = 0; i < polygon.size(); i++) {
    for (int j = i + 2; j < polygon.size(); j++) {
      if (seg_seg_intersect(polygon[i], polygon[(i + 1) % polygon.size()], polygon[j], polygon[(j + 1) % polygon.size()], p5) == '1') {
        return false;
      }
    }
  }
  return true;
}


/* function that computes if a point is within 10 distance away from a line segment */
bool calc_collision(point2d p0, point2d p1, point2d p2) {
  double A = p0.x - p1.x;
  double B = p0.y - p1.y;
  double C = p2.x - p1.x;
  double D = p2.y - p1.y;

  double dot = A * C + B * D;
  double len_sq = C * C + D * D;
  double param = -1;
  if (len_sq != 0) //in case of 0 length line
      param = dot / len_sq;

  double xx, yy;

  if (param < 0) {
    xx = p1.x;
    yy = p1.y;
  }
  else if (param > 1) {
    xx = p2.x;
    yy = p2.y;
  }
  else {
    xx = p1.x + param * C;
    yy = p1.y + param * D;
  }

  double dx = p0.x - xx;
  double dy = p0.y - yy;
  if (sqrt(dx * dx + dy * dy) < 5) {
    return true;
  }
  return false;
}

double getDoubleAngleNormal(double direction, point2d a, point2d b) {
  double lineAngle = atan2(b.y - a.y, b.x - a.x);
  double angleDiff = lineAngle - direction;

  // normalize the angle to be between 0 and 2pi
  while (angleDiff < 0) {
    angleDiff += 2 * M_PI;
  }
  while (angleDiff > 2 * M_PI) {
    angleDiff -= 2 * M_PI;
  }

  if (abs(angleDiff) < 1e-4) {
    // guard is traveling parallel to line segment, return angle opposite of line angle
    return 0;
  }

  double doubleAngleNormal = lineAngle + angleDiff;
  if (doubleAngleNormal < 0) {
    doubleAngleNormal += 2 * M_PI;
  }
  if (doubleAngleNormal > 2 * M_PI) {
    doubleAngleNormal -= 2 * M_PI;
  }

  return doubleAngleNormal;
}


/* ****************************** */
/* draw the polygon */
void draw_polygon(vector<point2d>& poly){
   if (poly.size() == 0) return; 

  //set color  and polygon mode 
  glColor3fv(yellow);   
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
 
  int i;
  for (i=0; i<poly.size()-1; i++) {
    //draw a line from point i to i+1
    glBegin(GL_LINES);
    glVertex2f(poly[i].x, poly[i].y); 
    glVertex2f(poly[i+1].x, poly[i+1].y);
    glEnd();
  }
  //draw a line from the last point to the first  
  int last=poly.size()-1; 
    glBegin(GL_LINES);
    glVertex2f(poly[last].x, poly[last].y); 
    glVertex2f(poly[0].x, poly[0].y);
    glEnd();
}


/* ****************************** */
/* Draw lines to make up museum polygon */
void add_lines_to_polygon(double x, double y) {
  if ((x==lastx && y==lasty) || x<0 || y<0) {
    return;
  }
  lastx = x;
  lasty = y;
  point2d p; 
  p.x = x; 
  p.y = y; 
  p.isReflex = false;
  p.isLeft = false;
  polygon.push_back(p); 
}

/* ****************************** */
/* Draw lines to make up museum polygon */
void draw_lines_to_polygon() {
  if (polygon.size() == 0) return; 

  //set color  and polygon mode 
  glColor3fv(yellow);   
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
 
  int i;
  for (i=0; i<polygon.size()-1; i++) {
    //draw a line from point i to i+1
    glBegin(GL_LINES);
    glVertex2f(polygon[i].x, polygon[i].y); 
    glVertex2f(polygon[i+1].x, polygon[i+1].y);
    glEnd();
  }
  //draw a line from the last point to the first  
  int last=polygon.size()-1; 
    glBegin(GL_LINES);
    glVertex2f(polygon[last].x, polygon[last].y); 
    glVertex2f(polygon[0].x, polygon[0].y);
    glEnd();
}

void add_robot(double x, double y) {
  if (x==lastx && y==lasty) {
    return;
  }
  lastx = x;
  lasty = y;

  rob.pos.x = x;
  rob.pos.y = y;

  // generate random GLfloat [3]:
  GLfloat r = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  GLfloat g = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  GLfloat b = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  rob.r = r;
  rob.g = g;
  rob.b = b;
}

void draw_robot(){
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  double r = 5;
  GLfloat color[3] = {rob.r, rob.g, rob.b};
  glColor3fv(color);   
  glBegin(GL_POLYGON);
  for (double theta = 0; theta < 2*M_PI; theta+=.3) {
    glVertex2f(rob.pos.x + r*cos(theta), rob.pos.y + r*sin(theta));
  }
  glEnd();
}

void add_goal(double x, double y) {
  if (x==lastx && y==lasty) {
    return;
  }
  lastx = x;
  lasty = y;

  gol.pos.x = x;
  gol.pos.y = y;

  // generate random GLfloat [3]:
  // change random color to be different from robot
  GLfloat r = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  GLfloat g = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  GLfloat b = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  gol.r = r;
  gol.g = g;
  gol.b = b;
}

void draw_goal(){
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  double r = 5;
  GLfloat color[3] = {gol.r, gol.g, gol.b};
  glColor3fv(color);   
  glBegin(GL_POLYGON);
  for (double theta = 0; theta < 2*M_PI; theta+= M_PI/4) {
    glVertex2f(gol.pos.x + r*cos(theta), gol.pos.y + r*sin(theta));
    if (r == 5) {
      r = 15;
    } else {
      r = 5;
    }
  }
  glEnd();
}



bool isPointToRightOfLine(point2d a, point2d b, point2d c) {
  // there is a line segment from firstPoint to middleVertex; if thirdPoint is to the right of that segment,
  // return true, otherwise return false
  // if the cross product of the vector from firstPoint to middleVertex and the vector from firstPoint to thirdPoint
  // is positive, then thirdPoint is to the right of the line segment
  // if the cross product is negative, then thirdPoint is to the left of the line segment
  // if the cross product is zero, then thirdPoint is on the line segment
  return ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) < 0.001;
}

/* our coordinate system is (0,0) x (WINDOWSIZE,WINDOWSIZE) with the
   origin in the lower left corner 
*/
//draw a circle with center = (x,y)
void draw_circle(double x, double y){

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glColor3fv(blue);   

  double r = 5;
  glBegin(GL_POLYGON);
  for(double theta = 0; theta < 2*M_PI; theta+=.3){
   glVertex2f(x + r*cos(theta), y + r*sin(theta));
  }
  glEnd();
}



/* ******************************** */
void print_polygon(vector<point2d>& poly) {
  printf("polygon:"); 
  for (unsigned int i=0; i<poly.size()-1; i++) {
    printf("\tedge %d: [(%f,%f), (%f,%f)]\n",
	   i, poly[i].x, poly[i].y, poly[i+1].x, poly[i+1].y);
  }
  //print last edge from last point to first point 
  int last = poly.size()-1; 
    printf("\tedge %d: [(%f,%f), (%f,%f)]\n",
	   last, poly[last].x, poly[last].y, poly[0].x, poly[0].y);

}

bool isVisible(point2d * p, point2d * d) {
  // want to check if there are any intersections between the line segment from p to d and any of the edges of the obstacles
  // loop through all obstacles
  point2d * newIntersect = new point2d;
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      // check if the line segment intersects with the edge from obstacles[i][j] to obstacles[i][(j+1)%obstacles[i].size()]
      char intersect = seg_seg_intersect(*p, *d, obstacles[i][j], obstacles[i][(j+1)%obstacles[i].size()], newIntersect);
      if (intersect == '1') {
        return false;
      }
    }
  }
  return true;
}

double find_distance(point2d p1, point2d p2) {
  return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

bool checkIfPointInObstacle(point2d p, int numObstacle) {
  int numberOfIntersects = 0;
  point2d outOfBounds;
  outOfBounds.x = 800;
  outOfBounds.y = 375;
  point2d *p5 = new point2d;
  p5->x = -1;
  p5->y = -1;
  for (int j = 0; j < obstacles[numObstacle].size(); j++) {
    char intersection = seg_seg_intersect(obstacles[numObstacle][j], obstacles[numObstacle][(j + 1) % obstacles[numObstacle].size()], p, outOfBounds, p5);
    if (intersection == '1' || intersection == 'v') {
        numberOfIntersects += 1;
    }
  }
  // cout << "number of intersects: " << numberOfIntersects << endl;
  if (numberOfIntersects % 2 == 0) {
    return false;
  }
  return true;
}


void find_lines(point2d * p, int numObstacle, int pointNum) {
  // loop through all points in the obstacles
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      // check if obstacles[i][j] has already been processed
      if (processedPoints.find(&obstacles[i][j]) != processedPoints.end() || (&obstacles[i][j] == p)) {
        continue;
      }

      if (i == numObstacle) {
        // if the points are neighbooring, then we can add the line
        if (abs(j - pointNum) == 1 || (j == 0 && pointNum == obstacles[i].size() - 1) || (j == obstacles[i].size() - 1 && pointNum == 0)) {
          nodeMap[p].push_back(node(&obstacles[i][j], find_distance(*p, obstacles[i][j]), true));
          nodeMap[&obstacles[i][j]].push_back(node(p, find_distance(*p, obstacles[i][j]), true));
          continue;
        }

        // now, since the points are not neighbooring, we need to check if the line segment intersects with any other line segments in the obstacle

        bool cont = true;

        // loop through all points in the obstacle
        for (int k = 0; k < obstacles[i].size(); k++) {
          // check if the line segment intersects with the edge from obstacles[i][j] to obstacles[i][(j+1)%obstacles[i].size()]
          if (k == j || k == pointNum) {
            continue;
          }
          point2d* temp = new point2d;
          char intersect = seg_seg_intersect(*p, obstacles[i][j], obstacles[i][k], obstacles[i][(k+1)%obstacles[i].size()], temp);
          if (intersect == '1') {
            cont = false;
          }
        }

        if (!cont) {
          continue;
        }

        // if there is no intersection, then the line is either entirely within the obstacle, or entirely outside the obstacle

        // so, want to find a point between the two points, and check if it is within the obstacle
        // first, find the midpoint between the two points
        point2d midpoint;
        midpoint.x = (p->x + obstacles[i][j].x) / 2;
        midpoint.y = (p->y + obstacles[i][j].y) / 2;

        // now, check if the midpoint is within the obstacle
        if (!checkIfPointInObstacle(midpoint, numObstacle)) {
          // check if this line segment intersects with any other obstacles (if nested non-convex polygons)
          for (int l = 0; l < obstacles.size(); l++) {
            // check if the line segment intersects with any of the edges of the obstacle
            for (int m = 0; m < obstacles[l].size(); m++) {
              // check if the line segment intersects with the edge from obstacles[i][j] to obstacles[i][(j+1)%obstacles[i].size()]
              if (l == i && (m == j || m == pointNum)) {
                continue;
              }
              point2d* temp = new point2d;
              char intersect = seg_seg_intersect(*p, obstacles[i][j], obstacles[l][m], obstacles[l][(m+1)%obstacles[l].size()], temp);
              if (intersect == '1') {
                cont = false;
              }
            }
          }

          if (cont) {
            nodeMap[p].push_back(node(&obstacles[i][j], find_distance(*p, obstacles[i][j])));
            nodeMap[&obstacles[i][j]].push_back(node(p, find_distance(*p, obstacles[i][j])));
          }
        }

        continue;
      }

      // if no intersection, add the line to the visible lines
      else if (isVisible(p, &obstacles[i][j])) {
        // if p is the robot, print "hello"
        // if (p == &rob.pos) {
        //   printf("hello");
        // }

        // add the line to the visible lines

        // check if p is not already in the hashtable
        if (nodeMap.find(p) == nodeMap.end()) {
          // add p to the hashtable, with a vector containing only point obstacles[i][j]
          nodeMap.insert(pair<point2d*, vector<node> >(p, vector<node>()));
        }
        if (nodeMap.find(&obstacles[i][j]) == nodeMap.end()) {
          // add p to the hashtable, with a vector containing only point obstacles[i][j]
          nodeMap.insert(pair<point2d*, vector<node> >(&obstacles[i][j], vector<node>()));
        }

        // add the line from p to obstacles[i][j]
        nodeMap[p].push_back(node(&obstacles[i][j], find_distance(*p, obstacles[i][j])));
        // also add the line from obstacles[i][j] to p
        nodeMap[&obstacles[i][j]].push_back(node(p, find_distance(*p, obstacles[i][j])));
      }
    }
  }
  // now, add point p to hashset for already processed points
  processedPoints.insert(p);
}

void calc_visible_lines() {
  // first, find all visible lines from the robot's position
  find_lines(&rob.pos, -1, -1);
  // also, want to find all visible lines from the goal's position
  find_lines(&gol.pos, -1, -1);
  // loop through all points in the obstacles
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      // call function to calculate visible lines from that vertex
      find_lines(&obstacles[i][j], i, j);
    }
  }
}



void draw_visible_lines() {
  // create a unordered_set to keep track of which points have already been processed
  unordered_set<point2d*> processedPoints;

  // loop through all points in the obstacles
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      // check if obstacles[i][j] has already been processed
      if (processedPoints.find(&obstacles[i][j]) != processedPoints.end()) {
        continue;
      }
      // check if the point is in the hashtable
      if (nodeMap.find(&obstacles[i][j]) != nodeMap.end()) {
        // loop through all the nodes in the vector
        for (int k = 0; k < nodeMap[&obstacles[i][j]].size(); k++) {
          // draw a line from obstacles[i][j] to nodeMap[&obstacles[i][j]][k].point
          glColor3f(.3, 0, .5);


          // want to see if pos is the point before &obstacles[i][j] in the ith obstacle, or the point after &obstacles[i][j] in the ith obstacle
          int indexBefore = (j - 1) % obstacles[i].size();
          int indexAfter = (j + 1) % obstacles[i].size();
          if (indexBefore < 0) {
            indexBefore += obstacles[i].size();
          }

          // check if the point is the point before &obstacles[i][j] in the ith obstacle
          if (nodeMap[&obstacles[i][j]][k].pos == &obstacles[i][indexBefore] || nodeMap[&obstacles[i][j]][k].pos == &obstacles[i][indexAfter]) {
              glColor3f(0.0, .2, .7);
          }

          // first, check to see if the edge is in the shortestPath; if so, then draw it in red
          if (shortestPath.size() != 0) {
            for (int l = 0; l < shortestPath.size() - 1; l++) {
              if ((shortestPath[l] == &obstacles[i][j] && shortestPath[l + 1] == nodeMap[&obstacles[i][j]][k].pos) || (shortestPath[l] == nodeMap[&obstacles[i][j]][k].pos && shortestPath[l + 1] == &obstacles[i][j])) {
                glColor3f(1.0, 0.0, 0.0);
                break;
              }
            }
          }


          glBegin(GL_LINES);
          glVertex2f(obstacles[i][j].x, obstacles[i][j].y);
          glVertex2f(nodeMap[&obstacles[i][j]][k].pos->x, nodeMap[&obstacles[i][j]][k].pos->y);
          glEnd();
        }
      }
      // add the point to the processedPoints hashset
      processedPoints.insert(&obstacles[i][j]);
    }
  }
}



// vector<point2d*> shortestPath;
// double shortestPathLength;

void dijktras_algorithm() {
  // first, create a priority queue
  priority_queue<node, vector<node>, NodeComparator> pq;

  // create a unordered set to keep track of which points have already been processed
  unordered_set<point2d*, PointerHashFunction> visitedPoints;

  // create a unordered map to keep track of the distance from the start to each point
  unordered_map<point2d*, double, PointerHashFunction> distanceMap;
  unordered_map<point2d*, point2d*, PointerHashFunction> fromMap;


  // initialize the distance from the start to each point to infinity
  // first, do so for the robot and goal
  distanceMap.insert(pair<point2d*, double>(&rob.pos, 0));
  distanceMap.insert(pair<point2d*, double>(&gol.pos, numeric_limits<double>::infinity()));

  //now, loop through all points in the obstacles, and do the same
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      // add the point to the distanceMap, with a distance of infinity
      distanceMap.insert(pair<point2d*, double>(&obstacles[i][j], numeric_limits<double>::infinity()));
    }
  }

  // initialize the fromMap, so that each point stores a null pointer
  // first, do so for the goal
  fromMap.insert(pair<point2d*, point2d*>(&gol.pos, nullptr));

  // now, loop through all points in the obstacles, and do the same
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      // add the point to the fromMap, with a null pointer
      fromMap.insert(pair<point2d*, point2d*>(&obstacles[i][j], nullptr));
    }
  }



  // now, we can start the algorithm
  // first, add the robot to the priority queue
  pq.push(node(&rob.pos, 0));



  // now, loop through the priority queue until it is empty
  while (!pq.empty()) {
    // get the top node from the priority queue
    node n = pq.top();
    pq.pop();

    // check if the node has already been visited
    if (visitedPoints.find(n.pos) != visitedPoints.end()) {
      // it has already been visited, so continue
      continue;
    }

    // add the node to the visitedNodes hashset
    visitedPoints.insert(n.pos);

    // check if the node is the goal
    if (n.pos == &gol.pos) {

      // it is the goal, so we are done
      // want to go through the nodes, adding each vertex that is part of the shortest
      // path to the ShortestPath, and the shortestPathLength value.
      shortestPathLength = distanceMap[&gol.pos];
      
      // loop through the fromMap, starting at the goal, and going back to the robot
      // add these vertices to the shortestPath vector
      point2d * currentPoint = &gol.pos;
      shortestPath.push_back(currentPoint);
      while (currentPoint != &rob.pos) {
        // get the next point
        currentPoint = fromMap[currentPoint];

        // add the current point to the shortestPath vector
        shortestPath.push_back(currentPoint);
      }

      // reverse the shortestPath vector
      reverse(shortestPath.begin(), shortestPath.end());
      break;
    }

    // now, we need to loop through all of the neighbors of the current node
    for (int i = 0; i < nodeMap[n.pos].size(); i++) {
      // get the current neighbor
      point2d * neighbor = nodeMap[n.pos][i].pos;

      // check if the neighbor has already been visited
      if (visitedPoints.find(neighbor) != visitedPoints.end()) {
        // it has already been visited, so continue
        continue;
      }

      // calculate the distance from the start to the neighbor
      double distance = distanceMap[n.pos] + nodeMap[n.pos][i].distance;

      // check if the distance is less than the current distance from the start to the neighbor
      if (distance < distanceMap[neighbor]) {
        // it is less, so update the distanceMap and fromMap
        distanceMap[neighbor] = distance;
        fromMap[neighbor] = n.pos;

        // add the neighbor to the priority queue
        pq.push(node(neighbor, distance));
      }
    }




  }




}



/* ****************************** */
/* This is the main render function. We registered this function to be
   called by GL to render the window. 
 */
void display(void) {

  glClear(GL_COLOR_BUFFER_BIT);
  //clear all modeling transformations 
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();


  /* The default GL window is [-1,1]x[-1,1] with the origin in the
     center.  The camera is at (0,0,0) looking down negative
     z-axis.  

     The points are in the range (0,0) to (WINSIZE,WINSIZE), so they
     need to be mapped to [-1,1]x [-1,1] */
  
  //First we scale down to [0,2] x [0,2] */ 
  glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);  
  /* Then we translate so the local origin goes in the middle of teh
     window to (-WINDOWSIZE/2, -WINDOWSIZE/2) */
  glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0); 
  
  //now we draw in our local coordinate system (0,0) to
  //(WINSIZE,WINSIZE), with the origin in the lower left corner.

  //draw a circle where the mouse was last clicked. Note that this
  //point is stored as a global variable and is modified by the mouse handler function 
  if (curMode == DRAW_OBSTACLES) {
    draw_obstacles();
  }

  if (curMode == RUN) {
    draw_lines_to_polygon();
    timerfunc();
  }

  if (curMode == NEW_POLYGON) {
    draw_circle(mouse_x, mouse_y);
    add_lines_to_polygon(mouse_x, mouse_y); 
    draw_lines_to_polygon();
    if (obstacles.size() > 0) {
      draw_obstacles();
    }
  }

  if (curMode == DRAW_ROBOT) {
    draw_obstacles();
    add_robot(mouse_x, mouse_y);
    draw_robot();
  }

  if (curMode == DRAW_GOAL) {
    draw_obstacles();
    draw_robot();
    add_goal(mouse_x, mouse_y);
    draw_goal();
  }

  if (curMode == CALC_VISIBLE_LINES) {
    calc_visible_lines();
    // // print robot's point and the number of nodes it has
    // cout << "Robot: " << rob.pos.x << ", " << rob.pos.y << " has " << nodeMap[&rob.pos].size() << " nodes" << endl;

    // // loop through points in the obstacles, and print the number of nodes each point has
    // for (int i = 0; i < obstacles.size(); i++) {
    //   for (int j = 0; j < obstacles[i].size(); j++) {
    //     cout << "Point: " << obstacles[i][j].x << ", " << obstacles[i][j].y << " has " << nodeMap[&obstacles[i][j]].size() << " nodes" << endl;
    //   }
    // }
    curMode = DRAW_VISIBLE_LINES;
  }

  if (curMode == DRAW_VISIBLE_LINES) {
    draw_obstacles();
    draw_goal();
    draw_robot();
    draw_visible_lines();
  }

  /* execute the drawing commands */
  glFlush();
  glutPostRedisplay();
}


/* ****************************** */
void keypress(unsigned char key, int x, int y) {


  // for checking if check if polygon is simple
  point2d *p5 = new point2d ;
  p5->x = -1;
  p5->y = -1;

  switch(key) {
  case 'q':	
    exit(0);
    break;
  case 'n':
    if (obstacles.size() > 0) {
      // check if the polygon is valid by seeing if it intersects with any of the other polygons
      for (int i = 0; i < obstacles.size(); i++) {
        for (int j = 0; j < obstacles[i].size(); j++) {
          // check if the line segment intersects with any of the edges of the obstacle
          for (int k = 0; k < polygon.size(); k++) {
            // check if the line segment intersects with the edge from obstacles[i][j] to obstacles[i][(j+1)%obstacles[i].size()]
            point2d* temp = new point2d;
            char intersect = seg_seg_intersect(obstacles[i][j], obstacles[i][(j+1)%obstacles[i].size()], polygon[k], polygon[(k+1)%polygon.size()], temp);
            if (intersect == '1') {
              // the polygon is not valid, so don't add it
              printf("polygon is not valid\n");
              polygon.clear();
              break;
            }
          }
        }
      }
    }

    for (int i = 0; i < polygon.size(); i++) {
      for (int j = i + 2; j < polygon.size(); j++) {
        if (seg_seg_intersect(polygon[i], polygon[(i + 1) % polygon.size()], polygon[j], polygon[(j + 1) % polygon.size()], p5) == '1') {
          printf("polygon is not simple\n");
          polygon.clear();
          break;
        }
      }
    }

    obstacles.push_back(polygon);
    polygon.clear();
    break;
  case 's':
    // printf("pushing back polygon\n");
    // print out all points in the polygon
    for (int i = 0; i < polygon.size(); i++) {
      // cout << polygon[i].x << ", " << polygon[i].y << endl;
    }
    obstacles.push_back(polygon);
    curMode = DRAW_ROBOT;
    break;
  case 'e':

    // check if the robot is in any of the obstacles
    for (int i = 0; i < obstacles.size(); i++) {
      if (checkIfPointInObstacle(rob.pos, i)) {
        printf("robot is in obstacle\n");
        return;
      }
    }

    curMode = DRAW_GOAL;
    break;
  case 'c':

    // check if the goal is in any of the obstacles
    for (int i = 0; i < obstacles.size(); i++) {
      if (checkIfPointInObstacle(gol.pos, i)) {
        printf("goal is in obstacle\n");
        return;
      }
    }

    curMode = CALC_VISIBLE_LINES;
    break;
  case 'd':
    dijktras_algorithm();
    break;
  }
}


/* 
void glutMouseFunc(void (*func)(int button, int state, int x, int y));

glutMouseFunc sets the mouse callback for the current window. When a
user presses and releases mouse buttons in the window, each press and
each release generates a mouse callback. The button parameter is one
of GLUT_LEFT_BUTTON, GLUT_MIDDLE_BUTTON, or GLUT_RIGHT_BUTTON. For
systems with only two mouse buttons, it may not be possible to
generate GLUT_MIDDLE_BUTTON callback. For systems with a single mouse
button, it may be possible to generate only a GLUT_LEFT_BUTTON
callback. The state parameter is either GLUT_UP or GLUT_DOWN
indicating whether the callback was due to a release or press
respectively. The x and y callback parameters indicate the window
relative coordinates when the mouse button state changed. If a
GLUT_DOWN callback for a specific button is triggered, the program can
assume a GLUT_UP callback for the same button will be generated
(assuming the window still has a mouse callback registered) when the
mouse button is released even if the mouse has moved outside the
window.
*/
void mousepress(int button, int state, int x, int y) {



  if(state == GLUT_DOWN) { //mouse click detected 

    //(x,y) are in window coordinates, where the origin is in the upper
    //left corner; our reference system has the origin in lower left
    //corner, this means we have to reflect y
    mouse_x = (double)x;
    mouse_y = (double)(WINDOWSIZE - y); 
    // printf("mouse pressed at (%.1f,%.1f)\n", mouse_x, mouse_y); 
  }
  
  glutPostRedisplay();
}





//this function is called every frame. Use for animations 
void timerfunc() {
  
  //if you want the window to be re-drawn, call this 
  // glutPostRedisplay(); 
  if (curMode == RUN) {
    glutPostRedisplay(); 
  }

}




/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
     
   // Set the viewport to cover the new window
   glViewport(0, 0, width, height);
 
   glMatrixMode(GL_PROJECTION);  // operate on the Projection matrix
   glLoadIdentity();             // reset
   gluOrtho2D(0.0, (GLdouble) width, 0.0, (GLdouble) height); 
}
