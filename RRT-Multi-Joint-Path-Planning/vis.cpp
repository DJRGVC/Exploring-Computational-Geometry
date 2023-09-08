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

#ifndef CALLBACK
#define CALLBACK
#endif


#include <vector> 
#include <iostream>
#include <iomanip>

#include <utility>
#include <queue>
#include <unordered_map>
#include <unordered_set>



using namespace std; 
using std::cout;
using std::cerr;
using std::endl;
using std::ends;
using std::stringstream;



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

const double epsilon = 6;

const double entropyThreshold = .17;

const double movementSpeed = .8;  // Adjust this value to control the movement speed

const int sizeRRT = 1000000;

int numSections;

enum mode { NOT_TREE, GO, COMPUTE_PATH, FINISH_POLYGON, NEW_POLYGON, CALC_VISIBLE_LINES, DRAW_OBSTACLES, RUN, DRAW_ROBOT, DRAW_GOAL, DRAW_VISIBLE_LINES, SHOW_ROT };

mode curMode = NEW_POLYGON;

//the current user polygon
vector<point2d> polygon;

// vector of polygon obstacles
vector<vector<point2d> > obstacles;

// the robot
vector< vector<point2d> > rob;
vector< vector<point2d> > origRob;

// the node that is the goal
rrtNode* goalNode;

// rotated and translated robots
vector< vector<vector<point2d> > > rotatedRobots;

// the current goal
goal gol;

// vector for storing RRT nodes
vector<rrtNode*> rrtNodes;

// vector for storing RRT nodes in final path
vector<rrtNode*> finalPath;

GLuint listId1, listId2, listId3;       // IDs of display lists
GLdouble vertices[64][6];               // arrary to store newly created vertices (x,y,z,r,g,b) by combine callback
int vertexIndex = 0;                    // array index for above array incremented inside combine callback




// hash map of where key: point2d, value: vector of nodes
// use custom hashing function in geom.h
// unordered_map<point2d*, vector<node>, PointerHashFunction> nodeMap;
// unordered_set<point2d*, PointerHashFunction> processedPoints;


void CALLBACK tessBeginCB(GLenum which);
void CALLBACK tessEndCB();
void CALLBACK tessErrorCB(GLenum errorCode);
void CALLBACK tessVertexCB(const GLvoid *data);
void CALLBACK tessVertexCB2(const GLvoid *data);
void CALLBACK tessCombineCB(const GLdouble newVertex[3], const GLdouble *neighborVertex[4],
                            const GLfloat neighborWeight[4], GLdouble **outData);
const char* getPrimitiveType(GLenum type);

vector<GLuint> obIds;
vector<GLuint> robIds;


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
void draw_obstacles();
char seg_seg_intersect(point2d p1, point2d p2, point2d p3, point2d p4, point2d* p5);
void draw_circle(double x, double y, bool user);
double find_distance(point2d p1, point2d p2);
double entropy(vector<double> *a1, vector<double> *a2);
bool pointInPolygon(point2d p, vector<point2d> poly);



/* ****************************** */
int main(int argc, char** argv) {

  cout << "press 's' to start entering the vertices of the polygon and 'e' to end.\n (Please enter the vertices counter-clockwise)\n";
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
  glutIdleFunc(timerfunc); //register this if you want it called aat every fraame


  glShadeModel(GL_SMOOTH);                    // shading mathod: GL_SMOOTH or GL_FLAT
  glPixelStorei(GL_UNPACK_ALIGNMENT, 4);      // 4-byte pixel alignment

  // enable /disable features
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_CULL_FACE);

    // track material ambient and diffuse from surface color, call it before glEnable(GL_COLOR_MATERIAL)
  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);

  glClearColor(0, 0, 0, 0);                   // background color
  glClearStencil(0);                          // clear stencil buffer
  glClearDepth(0.0f);                         // 0 is near, 1 is far
  glDepthFunc(GL_LEQUAL);

  /* init GL */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);
  
  /* give control to event handler */
  glutMainLoop();
  return 0;
}

void compute_path() {
  // clear the final path
  finalPath.clear();

  goalNode->partOfPath = true;

  // add the goal node to the final path
  finalPath.push_back(goalNode);

  // now add the parent of the goal node to the final path
  rrtNode* par = goalNode->parent;

  // while the parent of the last node in the final path is not NULL
  while (par->parent != NULL) {
    // add the parent of the last node in the final path to the final path
    par->partOfPath = true;
    finalPath.push_back(par);

    par = finalPath[finalPath.size() - 1]->parent;
  }

  finalPath.push_back(par);
  // reverse the final path
  reverse(finalPath.begin(), finalPath.end());

  printf("final path size: %lu\n", finalPath.size());
}

anglepoint* generate_random_point() {
  anglepoint* p = new anglepoint;

  // make angle between 0 and 2pi
  for (int i = 0; i < origRob.size(); i++) {
    p->angles->push_back(static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 2 * M_PI);
  }
  p->pos = new point2d;
  p->pos->x = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 750;
  p->pos->y = static_cast <double> (rand()) / static_cast <double> (RAND_MAX) * 750;
  return p;
}

pair<double, rrtNode*> min_dist_to_node(anglepoint* p) {
  double minDist = 1000000;
  rrtNode* closest_node = new rrtNode;
  for (int i = 0; i < rrtNodes.size(); i++) {
    double dist = sqrt(pow(p->pos->x - rrtNodes[i]->orientation->pos->x, 2) + pow(p->pos->y - rrtNodes[i]->orientation->pos->y, 2));
    if (dist < minDist) {
      minDist = dist;
      closest_node = rrtNodes[i];
    }
  }
  return make_pair(minDist, closest_node);
}

void find_epsilon_point(rrtNode* closestNode, anglepoint* p, rrtNode* newPoint) {
  double theta = atan2(p->pos->y - closestNode->orientation->pos->y, p->pos->x - closestNode->orientation->pos->x);
  newPoint->orientation->pos->x = closestNode->orientation->pos->x + epsilon * cos(theta);
  newPoint->orientation->pos->y = closestNode->orientation->pos->y + epsilon * sin(theta);
  newPoint->orientation->angles = p->angles;
  newPoint->parent = closestNode;

  // generate new random angles within entropy range of closestNode for each section
  for (int i = 0; i < origRob.size(); i++) {
    double minAngle = closestNode->orientation->angles->at(i) - entropyThreshold;
    double maxAngle = closestNode->orientation->angles->at(i) + entropyThreshold;
    double newAngle = minAngle + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (maxAngle - minAngle)));
    newPoint->orientation->angles->at(i) = newAngle;
  }
}

vector< vector<point2d> > rotateRobot(vector<double>* thetas) {
  // for each point in the robot, rotate it by theta
  vector< vector<point2d> > rotatedRobot;

  // fill rotatedRobot with the original robot
  for (int i = 0; i < origRob.size(); i++) {
    vector<point2d> rotatedSection;
    for (int j = 0; j < origRob[i].size(); j++) {
      rotatedSection.push_back(origRob[i][j]);
    }
    rotatedRobot.push_back(rotatedSection);
  }

  // cout << "numSections in rotate: " << numSections << endl;
  // cout << "thetas size in rotate: " << thetas->size() << endl;
  // cout << "origRob size in rotate: " << origRob.size() << endl;

  for (int i = 0; i < origRob.size(); i++) {
    // find pivot point
    point2d* pivot = &rotatedRobot[i][0];

    // now loop through all remaining sections and rotate them by theta
    for (int j = i; j < origRob.size(); j++) {
      for (int p = 0; p < rotatedRobot[j].size(); p++) {
        // if the point is the pivot point, don't rotate it
        if (p == 0 && j == i) {
          continue;
        }
        // rotate the point by theta
        double rotatedX = (rotatedRobot[j][p].x - pivot->x) * cos((*thetas)[i]) - (rotatedRobot[j][p].y - pivot->y) * sin((*thetas)[i]) + pivot->x;
        double rotatedY = (rotatedRobot[j][p].x - pivot->x) * sin((*thetas)[i]) + (rotatedRobot[j][p].y - pivot->y) * cos((*thetas)[i]) + pivot->y;

        rotatedRobot[j][p].x = rotatedX;
        rotatedRobot[j][p].y = rotatedY;
      }
    }
  }

  return rotatedRobot;
}

void translateRobot(vector< vector<point2d> >* rotatedRobot, double x, double y) {
  // find the difference between the first point in the robot and the point (x, y)
  double xDiff = x - (*rotatedRobot)[0][0].x;
  double yDiff = y - (*rotatedRobot)[0][0].y;


  // translate the robot by that difference
  for (int i = 0; i < (*rotatedRobot).size(); i++) {
    for (int j = 0; j < (*rotatedRobot)[i].size(); j++) {
      (*rotatedRobot)[i][j].x += xDiff;
      (*rotatedRobot)[i][j].y += yDiff;
    }
  }
}

bool doesRobotIntersect(vector< vector<point2d> > rotatedRobot) {
  // find if robot intersects itself
  for (int a = 0; a < rotatedRobot.size(); a++) {
    for (int b = 0; b < rotatedRobot[a].size(); b++) {
      for (int c = 0; c < rotatedRobot.size(); c++) {
        for (int d = 0; d < rotatedRobot[c].size(); d++) {
          if (a != c || b != d) {
            char intersect = seg_seg_intersect(rotatedRobot[a][b], rotatedRobot[a][(b+1)%rotatedRobot[a].size()], rotatedRobot[c][d], rotatedRobot[c][(d+1)%rotatedRobot[c].size()], NULL);
            if (intersect == '1') {
              return true;
            }
          }
        }
      }
    }
  }

  // find if robot inside itself (choose a random point from each segment, and check if it is inside prior segment)
  for (int i = 1; i < rotatedRobot.size(); i++) {
    for (int j = 0; j < rotatedRobot[i].size(); j++) {
      point2d randomPoint = rotatedRobot[i][rand() % rotatedRobot[i].size()];
      if (pointInPolygon(randomPoint, rotatedRobot[i-1])) {
        return true;
      }
    }
  }


  // check for intersections between the rotated and translated robot and the obstacles
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      for (int k = 0; k < rotatedRobot.size(); k++) {
        // check if the line segment intersects with the edge from obstacles[i][j] to obstacles[i][(j+1)%obstacles[i].size()]
        for (int q = 0; q < rotatedRobot[k].size(); q++) {
          char intersect = seg_seg_intersect(rotatedRobot[k][q], rotatedRobot[k][(q+1)%rotatedRobot[k].size()], obstacles[i][j], obstacles[i][(j+1)%obstacles[i].size()], NULL);
          if (intersect == '1' || intersect == 'v') {
            return true;
          }
        }
      }
    }
  }
  return false;
}


bool pointInPolygon(point2d p, vector<point2d> poly) {
  int numberOfIntersects = 0;
  point2d outOfBounds;
  outOfBounds.x = 800;
  outOfBounds.y = 375;
  for (int j = 0; j < poly.size(); j++) {
    char intersection = seg_seg_intersect(poly[j], poly[(j + 1) % poly.size()], p, outOfBounds, NULL);
    if (intersection == '1' || intersection == 'v') {
        numberOfIntersects += 1;
    }
  }
  if (numberOfIntersects % 2 == 0) {
    return false;
  }
  return true;
}

bool robotInObstacle(vector< vector<point2d> >* rotatedRobot) {
  // simply want to choose a random point in the robot and check if it is within the bounds of an obstacle
  // if it is, return true, otherwise return false
  point2d randomPoint = (*rotatedRobot)[0][rand() % (*rotatedRobot)[0].size()];

  // loop through all obstacles
  for (int i = 0; i < obstacles.size(); i++) {
    if (pointInPolygon(randomPoint, obstacles[i])) {
      return true;
    }
  }

  return false;
}

bool robotOutOfBound(vector< vector<point2d> >* rotatedRobot) {
  // loop through all points, and check that they are all between 0 and windowsize
  for (int i = 0; i < (*rotatedRobot).size(); i++) {
    for (int j = 0; j < (*rotatedRobot)[i].size(); j++) {
      if ((*rotatedRobot)[i][j].x < 0 || (*rotatedRobot)[i][j].x > WINDOWSIZE || (*rotatedRobot)[i][j].y < 0 || (*rotatedRobot)[i][j].y > WINDOWSIZE) {
        return true;
      }
    }
  }
  return false;
}

bool isFree(rrtNode* newNode) {

  // cout << "cinco" << endl;

  // print out angles in newNode
  // for (int i = 0; i < newNode->orientation->angles->size(); i++) {
  //   cout << (*newNode->orientation->angles)[i] << endl;
  // }

  // first, want to rotate the robot by the angle of the node
  vector< vector<point2d> > rotatedRobot = rotateRobot(newNode->orientation->angles);

  // cout << "seis" << endl;

  // then, want to translate the robot by the position of the node
  translateRobot(&rotatedRobot, newNode->orientation->pos->x, newNode->orientation->pos->y);

  // cout << "siete" << endl;

  // check for intersections between the rotated and translated robot and the obstacles
  if (doesRobotIntersect(rotatedRobot)) {
    return false;
  }

  // cout << "ocho" << endl;

  // now, want to check if the robot is within the bounds of an obstacle
  if (robotInObstacle(&rotatedRobot)) {
    return false;
  }

  // cout << "nueve" << endl;

  if (robotOutOfBound(&rotatedRobot)) {
    return false;
  }

  // cout << "diez" << endl;

  rotatedRobots.push_back(rotatedRobot);

  // if we get here, the robot is not intersecting any obstacles and is not within the bounds of an obstacle, so it is free
  return true;
}

void draw_rotated_robots() {
  for (int i = 0; i < rotatedRobots.size(); i++) {
    for (int j = 0; j < rotatedRobots[i].size(); j++) {
      for (int k = 0; k < rotatedRobots[i][j].size(); k++) {
        draw_circle(rotatedRobots[i][j][k].x, rotatedRobots[i][j][k].y, false);

        // connect lines between the points
        glColor3fv(green);
        if (k != rotatedRobots[i][j].size() - 1) {
          glBegin(GL_LINES);
          glVertex2f(rotatedRobots[i][j][k].x, rotatedRobots[i][j][k].y); 
          glVertex2f(rotatedRobots[i][j][k+1].x, rotatedRobots[i][j][k+1].y); 
          glEnd();
        } else {
          glBegin(GL_LINES);
          glVertex2f(rotatedRobots[i][j][k].x, rotatedRobots[i][j][k].y); 
          glVertex2f(rotatedRobots[i][j][0].x, rotatedRobots[i][j][0].y); 
          glEnd();
        }
      }
    }
  }
}

bool newLineIntersectObstacle(point2d* p1, point2d* p2) {
  // loop through all obstacles
  for (int i = 0; i < obstacles.size(); i++) {
    // loop through all points in the obstacle
    for (int j = 0; j < obstacles[i].size(); j++) {
      // check if the line segment intersects with the edge from obstacles[i][j] to obstacles[i][(j+1)%obstacles[i].size()]
      char intersect = seg_seg_intersect(*p1, *p2, obstacles[i][j], obstacles[i][(j+1)%obstacles[i].size()], NULL);
      if (intersect == '1' || intersect == 'v') {
        return true;
      }
    }
  }
  return false;
}

double entropy(vector<double> *a1, vector<double> *a2) {
  double entropy = 0;
  for (int i = 0; i < a1->size(); i++) {
    entropy += pow((*a1)[i] - (*a2)[i], 2);
  }
  return sqrt(entropy) / a1->size();
}


rrtNode* populate_RRT() {
  // cout << "uno" << endl;
  anglepoint* p = generate_random_point();
  pair<double, rrtNode*> result = min_dist_to_node(p);
  // cout << "dos" << endl;
  double minDistance = result.first;
  rrtNode* closestNode = result.second;

  while (minDistance <= epsilon) {
    p = generate_random_point();
    result = min_dist_to_node(p);
    minDistance = result.first;
    closestNode = result.second;
  }

  // cout << "tres" << endl;

  // create new node at epsilon distance from closestNode in direction of p
  rrtNode* newNode = new rrtNode;
  find_epsilon_point(closestNode, p, newNode);

  // cout << "cuatro" << endl;

  // find out if newPoint is a valid position and orientation for the robot
  // if it is, add it to the RRT
  if (isFree(newNode)) {
    // check if this connection would intersect any obstacles
    // if it does, don't add it to the RRT
    if (newLineIntersectObstacle(newNode->orientation->pos, newNode->parent->orientation->pos)) {
      return NULL;
      // pop last robot from rotatedRobots
      rotatedRobots.pop_back();
    }

    rrtNodes.push_back(newNode);

    // want to see if node is within epsilon distance of goal
    if (find_distance(*newNode->orientation->pos, gol.pos) <= epsilon) {
      // if it is, add it to the RRT and return true
      return newNode;
    }
  }

  // cout angles in newNode
  // for (int i = 0; i < newNode->orientation->angles->size(); i++) {
  //   cout << "angles: " << (*newNode->orientation->angles)[i] << endl;
  // }

  return NULL;
}

void fill_RRT() {
  int lastSize = rrtNodes.size();
  goalNode = NULL;
  while (goalNode == NULL && rrtNodes.size() < sizeRRT) {
    goalNode = populate_RRT();
    if (rrtNodes.size() != lastSize) {
      cout << "rrtNodes size: " << rrtNodes.size() << endl;
      lastSize = rrtNodes.size();
    }
    // draw position of last added node
    draw_circle(rrtNodes[rrtNodes.size() - 1]->orientation->pos->x, rrtNodes[rrtNodes.size() - 1]->orientation->pos->y, false);
  }
}

void draw_RRT() {
  for (int i = 0; i < rrtNodes.size(); i++) {
    // draw a circle at the node
    bool partOfPath = false;
    if (rrtNodes[i]->partOfPath) {
      partOfPath = true;
    } else {
      partOfPath = false;
    }
    draw_circle(rrtNodes[i]->orientation->pos->x, rrtNodes[i]->orientation->pos->y, partOfPath);
    // draw a line from the node to its parent

    GLfloat magenta[] = {.8, 0.0f, .8f};  // Assuming magenta color values

    glColor3f(magenta[0], magenta[1], magenta[2]);  
    if (i != 0) {
      glBegin(GL_LINES);
      glVertex2f(rrtNodes[i]->orientation->pos->x, rrtNodes[i]->orientation->pos->y); 
      glVertex2f(rrtNodes[i]->parent->orientation->pos->x, rrtNodes[i]->parent->orientation->pos->y); 
      glEnd();
    }
  }
}


void draw_obstacles() {
  // clear obids
  obIds.clear();
  // set glmode to polygon
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  // set color for tesselation filling
  glColor3f(.5, .5, .9);


  for (int i = 0; i < obstacles.size(); i++) {
    GLuint id = glGenLists(1);  // create a display list
    if(!id) cout << "ERROR: Cannot create a new list" << endl; 
    GLUtesselator *tess = gluNewTess(); // create a tessellator
    if (!tess) cout << "ERROR: Cannot create tessellation object" << endl;

    gluTessCallback(tess, GLU_TESS_BEGIN, (void (CALLBACK *)())tessBeginCB);
    gluTessCallback(tess, GLU_TESS_END, (void (CALLBACK *)())tessEndCB);
    gluTessCallback(tess, GLU_TESS_ERROR, (void (CALLBACK *)())tessErrorCB);
    gluTessCallback(tess, GLU_TESS_VERTEX, (void (CALLBACK *)())tessVertexCB);
    gluTessCallback(tess, GLU_TESS_COMBINE, (void (CALLBACK *)())tessCombineCB);

    gluTessProperty(tess, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO);
    glNewList(id, GL_COMPILE);
    gluTessBeginPolygon(tess, 0);                   // with NULL data
        gluTessBeginContour(tess);
          for (int j = 0; j < obstacles[i].size(); j++) {
            // convert to GLdouble *
            GLdouble *vert = new GLdouble[3];
            vert[0] = obstacles[i][j].x;
            vert[1] = obstacles[i][j].y;
            vert[2] = 0;
            gluTessVertex(tess, vert, vert);
          }
        gluTessEndContour(tess);
    gluTessEndPolygon(tess);
    glEndList();

    gluDeleteTess(tess);

    obIds.push_back(id);
  }

  // draw the obstacles
  for (int i = 0; i < obIds.size(); i++) {
    glCallList(obIds[i]);
  }

  glColor3f(.7, .7, 1);

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

  if (p5 != NULL) {
    p5->x = p1.x + s * (p2.x - p1.x);
    p5->y = p1.y + s * (p2.y - p1.y);
  }

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

bool is_simple_robot() {
  point2d *p5 = new point2d ;
  p5->x = -1;
  p5->y = -1;
  // if any of the segments intersect, return false, otherwise return true
  // loop through all segments in the robot, make sure none of them intersect, if they do, return false, otherwise return true
  for (int p = 0; p < rob.size(); p++) {
    for (int i = 0; i < rob[0].size(); i++) {
      for (int j = i + 2; j < rob[0].size(); j++) {
        if (seg_seg_intersect(rob[p][i], rob[p][(i + 1) % rob[p].size()], rob[p][j], rob[p][(j + 1) % rob[p].size()], p5) == '1') {
          return false;
        }
      }
    }
  }
  return true;
}

void checkIfGuardInPolygon() {
  // for (int i = 0; i < guards.size(); i++)  {
    int numberOfIntersects = 0;
    point2d outOfBounds;
    outOfBounds.x = 800;
    outOfBounds.y = 375;
    point2d *p5 = new point2d;
    p5->x = -1;
    p5->y = -1;
    for (int j = 0; j < polygon.size(); j++) {
      char intersection = seg_seg_intersect(polygon[j], polygon[(j + 1) % polygon.size()], gol.pos, outOfBounds, p5);
      if (intersection == '1' || intersection == 'v') {
          numberOfIntersects += 1;
      }
    }
    printf("# of intersects: %d", numberOfIntersects);
    if (numberOfIntersects % 2 == 0) {
      printf("guard outside bounds of polygon\n");
    }
    if (numberOfIntersects % 2 == 1) {
      printf("guard within bounds of polygon\n");
    }
  }
// }


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
  polygon.push_back(p); 
}

/* ******************************* */
void add_lines_to_robot(double x, double y) {
  if ((x==lastx && y==lasty) || x<0 || y<0) {
    return;
  }
  lastx = x;
  lasty = y;
  point2d p; 
  p.x = x; 
  p.y = y; 

  if (rob.size() <= numSections) {
    // create new vector of point2d
    // cout << "1" << endl;
    vector<point2d> section;

    // add point from previous section
    // cout << "2" << endl;
    if (numSections > 0) {
      section.push_back(rob[numSections - 1][rob[numSections - 1].size() - 1]);
      // make this point a isJoint
      section[0].isJoint = true;
    }

    // add p
    // cout << "2" << endl;
    section.push_back(p);


    // add section to rob
    // cout << "3" << endl;
    rob.push_back(section);

    // cout << "4" << endl;
    return;
  }

  rob[numSections].push_back(p);
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


void draw_robot(){
  if (rob.size() == 0) return; 

  robIds.clear();
  // set glmode to polygon
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  // set color for tesselation filling
  glColor3f(.3, 0, 0);


  for (int i = 0; i < rob.size(); i++) {
    GLuint id = glGenLists(1);  // create a display list
    if(!id) cout << "ERROR: Cannot create a new list" << endl; 
    GLUtesselator *tess = gluNewTess(); // create a tessellator
    if (!tess) cout << "ERROR: Cannot create tessellation object" << endl;

    gluTessCallback(tess, GLU_TESS_BEGIN, (void (CALLBACK *)())tessBeginCB);
    gluTessCallback(tess, GLU_TESS_END, (void (CALLBACK *)())tessEndCB);
    gluTessCallback(tess, GLU_TESS_ERROR, (void (CALLBACK *)())tessErrorCB);
    gluTessCallback(tess, GLU_TESS_VERTEX, (void (CALLBACK *)())tessVertexCB);
    gluTessCallback(tess, GLU_TESS_COMBINE, (void (CALLBACK *)())tessCombineCB);

    gluTessProperty(tess, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO);
    glNewList(id, GL_COMPILE);
    gluTessBeginPolygon(tess, 0);                   // with NULL data
        gluTessBeginContour(tess);
          for (int j = 0; j < rob[i].size(); j++) {
            // convert to GLdouble *
            GLdouble *vert = new GLdouble[3];
            vert[0] = rob[i][j].x;
            vert[1] = rob[i][j].y;
            vert[2] = 0;
            gluTessVertex(tess, vert, vert);
            // cout << "vert: " << vert[0] << ", " << vert[1] << endl;
          }
        gluTessEndContour(tess);
    gluTessEndPolygon(tess);
    glEndList();

    gluDeleteTess(tess);

   robIds.push_back(id);
  }

  // draw the obstacles
  for (int i = 0; i < robIds.size(); i++) {
    glCallList(robIds[i]);
  }

  //set color and polygon mode 
  glColor3fv(red);   
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  for (int j = 0; j < rob.size(); j++) {

    for (int i=0; i<rob[j].size()-1; i++) {
      if (rob[j][i].isJoint) {
        // draw a circle at the joint
        draw_circle(rob[j][i].x, rob[j][i].y, false);
      }

      glColor3fv(red);
      //draw a line from point i to i+1
      glBegin(GL_LINES);
      glVertex2f(rob[j][i].x, rob[j][i].y); 
      glVertex2f(rob[j][i+1].x, rob[j][i+1].y);
      glEnd();
    }
    //draw a line from the last point to the first  
    int last=rob[j].size()-1; 
    glBegin(GL_LINES);
    glVertex2f(rob[j][last].x, rob[j][last].y); 
    glVertex2f(rob[j][0].x, rob[j][0].y);
    glEnd();
  }

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
void draw_circle(double x, double y, bool user){

  double r = 5;
  if (user) {
    glColor3fv(blue);   
  } else {
    glColor3fv(white);   
    r = 2;
  }

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

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



double shortestAngularDistance(double a, double b) {
  double diff = fmod(( b - a + M_PI ),  2*M_PI) - M_PI;
  return diff < -M_PI ? diff + 2*M_PI : diff;
  // if (shortest)
  // double dist = fmod(b - a, 2 * M_PI);
  // if (abs(dist) > M_PI) {
  //   return dist - 2 * M_PI;
  // }
  // return dist;
  // if((a + d - M_PI) < 0 || (a + d - M_PI) > 2 * M_PI)
  //   return d;
  // else
  //   return d - 2 * M_PI;
  // 6 cases!

  if (a<b && ((a < M_PI && b < M_PI) || (a > M_PI && b > M_PI))) {
    return b - a;
  } else if (b<a && ((a<M_PI && b < M_PI) || (a > M_PI && b > M_PI))) {
    return b - a;
  } else if (a < b) {
    return -(2 * M_PI - (b - a));
  }
  return 2 * M_PI - (a - b);
}

// Interpolate the angle between a and b based on the given percent
double interpolateAngle(double a, double b, double percent) {
  double shortestDistance = shortestAngularDistance(a, b);
  return a + shortestDistance * percent;
}

int currentTargetIndex = 1;  // Global variable to track the index of the current target point
int stepsTaken = 1;
int stepsToTake = 0;

void moveRobotSmoothly() {
  // Check if there are more target positions to move to
  if (currentTargetIndex < finalPath.size()) {
    if (stepsTaken == 1) {
      // update number of steps to take
      stepsToTake = find_distance(*finalPath[currentTargetIndex - 1]->orientation->pos, *finalPath[currentTargetIndex]->orientation->pos) / movementSpeed;
    }

    // create a vector of angle differences between the last target orientation and the new target orientation

    vector<double> *angleDiffs = new vector<double>;
    for (int i = 0; i < finalPath[currentTargetIndex - 1]->orientation->angles->size(); i++) {
      angleDiffs->push_back(shortestAngularDistance(fmod((*finalPath[currentTargetIndex - 1]->orientation->angles)[i], 2 * M_PI), fmod((*finalPath[currentTargetIndex]->orientation->angles)[i], 2 * M_PI)));
    }

    // double lastAngle = fmod(finalPath[currentTargetIndex - 1]->orientation->angle, 2 * M_PI);
    // double newAngle = fmod(finalPath[currentTargetIndex]->orientation->angle + M_PI, 2 * M_PI);

    // first, find the angle between the current position and the target position
    double angle = atan2(finalPath[currentTargetIndex]->orientation->pos->y - finalPath[currentTargetIndex - 1]->orientation->pos->y, finalPath[currentTargetIndex]->orientation->pos->x - finalPath[currentTargetIndex - 1]->orientation->pos->x);

    // now, find the point that is steps taken * movementSpeed away from the current position in the direction of the target position
    point2d nextPosition;
    nextPosition.x = finalPath[currentTargetIndex - 1]->orientation->pos->x + stepsTaken * movementSpeed * cos(angle);
    nextPosition.y = finalPath[currentTargetIndex - 1]->orientation->pos->y + stepsTaken * movementSpeed * sin(angle);

    // cout << "whallah" << (stepsTaken / stepsToTake) << endl;

    // now, update angle diffs by multiplying each by stepsTaken / stepsToTake, and add to the last angle
    for (int i = 0; i < angleDiffs->size(); i++) {
      double temp = (*angleDiffs)[i] * ((double)stepsTaken / (double)stepsToTake) + (*finalPath[currentTargetIndex - 1]->orientation->angles)[i]; 
      (*angleDiffs)[i] = temp < 0 ? temp + 2 * M_PI : temp;
    }

    // print out angle diffs
    // for (int i = 0; i < angleDiffs->size(); i++) {
    //   cout << i << ": angle diff: " << (*angleDiffs)[i] << endl;
    // }

    // vector<point2d> newRob = rotateRobot(lastAngle + stepsTaken*(angleDiff / stepsToTake));
    vector< vector<point2d> > newRob = rotateRobot(angleDiffs);

    // now, want to translate the robot by the difference between the current position and the next position
    translateRobot(&newRob, nextPosition.x, nextPosition.y);

    rob = newRob;

    stepsTaken++;

    if (find_distance(*finalPath[currentTargetIndex]->orientation->pos, rob[0][0]) <= movementSpeed) {
      // if it is, increment the current target index
      currentTargetIndex++;
      stepsTaken = 1;
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

  if (curMode == RUN || curMode == COMPUTE_PATH || curMode == GO || curMode == SHOW_ROT || curMode == NOT_TREE) {
    draw_obstacles();
    draw_robot();
    draw_goal();
    if (curMode != NOT_TREE) {
      draw_RRT();
    }
  }

  if (curMode == NEW_POLYGON) {
    draw_circle(mouse_x, mouse_y, true);
    add_lines_to_polygon(mouse_x, mouse_y); 
    draw_lines_to_polygon();
    if (obstacles.size() > 0) {
      draw_obstacles();
    }
  }

  if (curMode == DRAW_ROBOT) {
    draw_circle(mouse_x, mouse_y, true);
    draw_obstacles();
    add_lines_to_robot(mouse_x, mouse_y);
    draw_lines_to_polygon();
    draw_robot();
  }

  if (curMode == DRAW_GOAL) {
    draw_obstacles();
    draw_robot();
    add_goal(mouse_x, mouse_y);
    draw_goal();
  }

  if (curMode == SHOW_ROT) {
    draw_rotated_robots();
  }

  /* execute the drawing commands */
  glFlush();
  glutPostRedisplay();
}

bool generateRobPos() {
  point2d* robInit = new point2d;
  robInit->x = rob[0][0].x;
  robInit->y = rob[0][0].y;
  // create vector of angles, all 0s, of length numSections
  vector<double>* angles = new vector<double>;
  for (int i = 0; i < rob.size(); i++) {
    angles->push_back(0);
  }
  anglepoint* robInitPos = new anglepoint(angles, robInit);
  rrtNode* initNode = new rrtNode(robInitPos, NULL);
  if (isFree(initNode)) {
    rrtNodes.push_back(initNode);
    return true;
  }
  return false;
}

/* ****************************** */
void keypress(unsigned char key, int x, int y) {
  // create vector of angles
  vector<double>* angles = new vector<double>;

  switch(key) {
  case 'q':	
    exit(0);
    break;
  case 'n':
    if (is_simple_polygon()) {
      cout << "polygon is simple\n";
    } else {
      cout << "polygon is not simple\n";
      polygon.clear();
      break;
    }
    obstacles.push_back(polygon);
    polygon.clear();
    break;
  case 's':
    if (curMode == DRAW_ROBOT) {
      // add 0s to the vector of angles
      for (int i = 0; i <= numSections; i++) {
        angles->push_back(0);
      }

      // check if the robot is free (make new anglepoint with vector of 0s of length numSections, check if it is free)

      // cout << "un" << endl;
      numSections++;

      origRob = rob;
      if (isFree(new rrtNode(new anglepoint(angles, new point2d(rob[0][0].x, rob[0][0].y)), NULL))) {
        cout << "rob is simple\n";
      } else {
        cout << "rob is not simple\n";
        // remove the last section of the robot

        // cout << "deux" << endl;
        rob.pop_back();
        origRob = rob;
        numSections--;
        // cout << "trois" << endl;
        break;
      }
    }
    else { 
      printf("pushing back polygon\n");
      // print out all points in the polygon
      // for (int i = 0; i < polygon.size(); i++) {
      //   cout << polygon[i].x << ", " << polygon[i].y << endl;
      // }
      obstacles.push_back(polygon);
      curMode = DRAW_ROBOT;

      numSections = 0;
    }

    origRob = rob;
    break;
  case 'e':
    if (generateRobPos()) {
      cout << "rob is free\n";
      origRob = rob;
    } else {
      cout << "rob is not free\n";
      rob.clear();
      origRob.clear();
      break;
    }
    curMode = DRAW_GOAL;
    break;
  case 'r':
    fill_RRT();
    curMode = RUN;
    break;
  case 'p':
    curMode = COMPUTE_PATH;
    compute_path();
    break;
  case 'g':
    curMode = GO;
    break;
  case 'c':
    if (curMode == SHOW_ROT) {
      curMode = RUN;
    }
    else {
      curMode = SHOW_ROT;
    }
    break;
  case 't':
    if (curMode == NOT_TREE) {
      curMode = GO;
    }
    else {
      curMode = NOT_TREE;
    }
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
    printf("mouse pressed at (%.1f,%.1f)\n", mouse_x, mouse_y); 
  }
  
  glutPostRedisplay();
}

//this function is called every frame. Use for animations 
void timerfunc() {
  
  //if you want the window to be re-drawn, call this 
  // glutPostRedisplay(); 
  if (curMode == GO || curMode == NOT_TREE) {
    draw_robot();
    moveRobotSmoothly();
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

///////////////////////////////////////////////////////////////////////////////
// GLU_TESS CALLBACKS
///////////////////////////////////////////////////////////////////////////////
void CALLBACK tessBeginCB(GLenum which)
{
    glBegin(which);

}



void CALLBACK tessEndCB()
{
    glEnd();

}



void CALLBACK tessVertexCB(const GLvoid *data)
{
    // cast back to double type
    const GLdouble *ptr = (const GLdouble*)data;

    glVertex3dv(ptr);

}



///////////////////////////////////////////////////////////////////////////////
// draw a vertex with color
///////////////////////////////////////////////////////////////////////////////
void CALLBACK tessVertexCB2(const GLvoid *data)
{
    // cast back to double type
    const GLdouble *ptr = (const GLdouble*)data;

    glColor3dv(ptr+3);
    glVertex3dv(ptr);

}



///////////////////////////////////////////////////////////////////////////////
// Combine callback is used to create a new vertex where edges intersect.
// In this function, copy the vertex data into local array and compute the
// color of the vertex. And send it back to tessellator, so tessellator pass it
// to vertex callback function.
//
// newVertex: the intersect point which tessellator creates for us
// neighborVertex[4]: 4 neighbor vertices to cause intersection (given from 3rd param of gluTessVertex()
// neighborWeight[4]: 4 interpolation weights of 4 neighbor vertices
// outData: the vertex data to return to tessellator
///////////////////////////////////////////////////////////////////////////////
void CALLBACK tessCombineCB(const GLdouble newVertex[3], const GLdouble *neighborVertex[4],
                            const GLfloat neighborWeight[4], GLdouble **outData)
{
    // copy new intersect vertex to local array
    // Because newVertex is temporal and cannot be hold by tessellator until next
    // vertex callback called, it must be copied to the safe place in the app.
    // Once gluTessEndPolygon() called, then you can safly deallocate the array.
    vertices[vertexIndex][0] = newVertex[0];
    vertices[vertexIndex][1] = newVertex[1];
    vertices[vertexIndex][2] = newVertex[2];

    // compute vertex color with given weights and colors of 4 neighbors
    // the neighborVertex[4] must hold required info, in this case, color.
    // neighborVertex was actually the third param of gluTessVertex() and is
    // passed into here to compute the color of the intersect vertex.
    vertices[vertexIndex][3] = neighborWeight[0] * neighborVertex[0][3] +   // red
                               neighborWeight[1] * neighborVertex[1][3] +
                               neighborWeight[2] * neighborVertex[2][3] +
                               neighborWeight[3] * neighborVertex[3][3];
    vertices[vertexIndex][4] = neighborWeight[0] * neighborVertex[0][4] +   // green
                               neighborWeight[1] * neighborVertex[1][4] +
                               neighborWeight[2] * neighborVertex[2][4] +
                               neighborWeight[3] * neighborVertex[3][4];
    vertices[vertexIndex][5] = neighborWeight[0] * neighborVertex[0][5] +   // blue
                               neighborWeight[1] * neighborVertex[1][5] +
                               neighborWeight[2] * neighborVertex[2][5] +
                               neighborWeight[3] * neighborVertex[3][5];


    // return output data (vertex coords and others)
    *outData = vertices[vertexIndex];   // assign the address of new intersect vertex

    ++vertexIndex;  // increase index for next vertex
}



void CALLBACK tessErrorCB(GLenum errorCode)
{
    const GLubyte *errorStr;

    errorStr = gluErrorString(errorCode);
    cerr << "[ERROR]: " << errorStr << endl;
}


const char* getPrimitiveType(GLenum type)
{
    switch(type)
    {
    case 0x0000:
        return "GL_POINTS";
        break;
    case 0x0001:
        return "GL_LINES";
        break;
    case 0x0002:
        return "GL_LINE_LOOP";
        break;
    case 0x0003:
        return "GL_LINE_STRIP";
        break;
    case 0x0004:
        return "GL_TRIANGLES";
        break;
    case 0x0005:
        return "GL_TRIANGLE_STRIP";
        break;
    case 0x0006:
        return "GL_TRIANGLE_FAN";
        break;
    case 0x0007:
        return "GL_QUADS";
        break;
    case 0x0008:
        return "GL_QUAD_STRIP";
        break;
    case 0x0009:
        return "GL_POLYGON";
        break;
    }
    return NULL;
}

