/* mouse1.c 
Laura Toma

Example using the mouse in OpenGL.  First the mouse is registered via
a callback function. Once registered, this function will be called on
any mouse event in the window.  The user can use this function to
respond to mouse events. This code will print the coordinates of the
point where the mouse is clicked, and will draw a small blue disk at
the point where the mouse is pressed.

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

enum mode { BEFORE, DRAW_LINES, DRAW_GAURDS, RUN};

mode curMode = BEFORE;

//the current polygon 
vector<point2d>  poly;

//the current user polygon
vector<point2d>  polygon;

//the current vector of gaurds
vector<guard>  guards;

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






/* ****************************** */
int main(int argc, char** argv) {

  cout << "press 's' to start entering the vertices of the polygon and 'e' to end.\n (Please enter the vertices counter-clockwise)\n";

  // print_polygon(poly);

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



void checkIfGuardInPolygon() {
  for (int i = 0; i < guards.size(); i++)  {
    int numberOfIntersects = 0;
    point2d outOfBounds;
    outOfBounds.x = 800;
    outOfBounds.y = 375;
    point2d *p5 = new point2d;
    p5->x = -1;
    p5->y = -1;
    for (int j = 0; j < polygon.size(); j++) {
      char intersection = seg_seg_intersect(polygon[j], polygon[(j + 1) % polygon.size()], guards[i].pos, outOfBounds, p5);
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


void move_guards() {
  for (int i = 0; i < guards.size(); i++) {
    double newAngle = 10;
    // check for collision with polygon:
    for (int j = 0; j < polygon.size(); j++) {
      point2d endpointOne = polygon[j];
      point2d endpointTwo = polygon[(j + 1) % polygon.size()];
      if (calc_collision(guards[i].pos, endpointOne, endpointTwo)) {
        newAngle = getDoubleAngleNormal(guards[i].dir, endpointOne, endpointTwo);
        break;
      }
    }

    if (newAngle != 10) {
      guards[i].dir = newAngle;
    }

    guards[i].pos.x += .7 * cos(guards[i].dir);
    guards[i].pos.y += .7 * sin(guards[i].dir);
  }
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
  if (x==lastx && y==lasty) {
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


void add_guard(double x, double y){
  if (x==lastx && y==lasty) {
    return;
  }
  lastx = x;
  lasty = y;

  point2d p; 
  p.x = x; 
  p.y = y; 
  // generate random GLfloat [3]:
  GLfloat r = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  GLfloat g = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  GLfloat b = static_cast <GLfloat> (rand()) / static_cast <GLfloat> (RAND_MAX);
  guard guar;
  guar.r = r;
  guar.g = g;
  guar.b = b;
  p.isReflex = false;
  p.isLeft = false;
  guar.pos = p;
  // choose a random direction, between 0 and 2pi
  guar.dir = static_cast <double> (rand()) / (static_cast <double> (RAND_MAX/(2*M_PI)));
  printf("guard direction: %f", guar.dir);
  guards.push_back(guar);

  cout << "added guard\n";
}

void draw_gaurds(){
  if (guards.size() == 0) return; 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  double r = 5;
  for (int i = 0; i < guards.size(); i++) {
    guard temp = guards[i];
    GLfloat tempColor[3] = {temp.r, temp.g, temp.b};
    glColor3fv(tempColor);   
    glBegin(GL_POLYGON);
    for(double theta = 0; theta < 2*M_PI; theta+=.3){
     glVertex2f(temp.pos.x + r*cos(theta), temp.pos.y + r*sin(theta));
    }
    glEnd();
  }
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
  // check if first point is a reflex point
  
  for (int i = 0; i < pts.size(); i++) {
    if (pts[i].x == p0.x && pts[i].y == p0.y) {
      continue;
    }
    radial_point2d rp;
    rp.p = pts[i];
    rp.angle = atan2(pts[i].y - p0.y, pts[i].x - p0.x);
    // truncate double rp.angle to 5 decimal points
    rp.angle = (float) rp.angle;
    if (pts[i].isLeft) {
      rp.distance = 0;
    } else {
      rp.distance = distance(p0, pts[i]);
    }
    sortedPoints.push(rp);
  }
}

bool isPointToRightOfLine(point2d a, point2d b, point2d c) {
  // there is a line segment from firstPoint to middleVertex; if thirdPoint is to the right of that segment,
  // return true, otherwise return false
  // if the cross product of the vector from firstPoint to middleVertex and the vector from firstPoint to thirdPoint
  // is positive, then thirdPoint is to the right of the line segment
  // if the cross product is negative, then thirdPoint is to the left of the line segment
  // if the cross product is zero, then thirdPoint is on the line segment
  return ((b.x - a.x)*(c.y - a.y) - (b.y - a.y)*(c.x - a.x)) < 0;
}

void calculate_visible_polygons() {
  // for each guard
  for (int i = 0; i < guards.size(); i++) {
    // get the guard's position
    // guard tempGuard = guards[i];
    // point2d guardPoint = tempGuard.pos;
    vector<point2d> visiblePoints;

    // first find visible points:
    for (int j = 0; j < polygon.size(); j++) {
      // first find closest visible point along that line segment
      point2d endPoint;

      double changeInY = (polygon[j].y - guards[i].pos.y);
      double changeInX = (polygon[j].x - guards[i].pos.x);

      endPoint.y = guards[i].pos.y + 20*changeInY;
      endPoint.x = guards[i].pos.x + 20*changeInX;
      endPoint.isReflex = false;


      bool passThroughReflexVertex = false;
      bool arePointsLeft = false;
      point2d leftPoint = j == 0 ? polygon[polygon.size() - 1] : polygon[j-1];
      point2d rightPoint = j == polygon.size() - 1 ? polygon[0] : polygon[j+1];
      endPoint.isLeft = false;
      passThroughReflexVertex = (isPointToRightOfLine(guards[i].pos, endPoint, leftPoint) != isPointToRightOfLine(guards[i].pos, endPoint, rightPoint));
      arePointsLeft = !isPointToRightOfLine(guards[i].pos, endPoint, leftPoint);

      // now, check line segment (gaurdpoint -> endPoint), and see if any intersection points occur
      point2d intersectionPoint;
      intersectionPoint.x = endPoint.x;
      intersectionPoint.y = endPoint.y;
      intersectionPoint.isReflex = false;


      if (!polygon[j].isReflex || passThroughReflexVertex) {
        intersectionPoint.x = polygon[j].x;
        intersectionPoint.y = polygon[j].y;
        intersectionPoint.isLeft = false;
        for (int k = 0; k < polygon.size(); k++) {
          point2d *newIntersect = new point2d;
          newIntersect->x = polygon[j].x;
          newIntersect->y = polygon[j].y;
          char intSect = seg_seg_intersect(guards[i].pos, endPoint, polygon[k], polygon[(k+1)%polygon.size()], newIntersect);
          if (intSect == '1' || intSect == 'v') {
            // if the new intersection point is closer than the current one, replace it
            if (distance(guards[i].pos, *newIntersect) < distance(guards[i].pos, intersectionPoint)) {
              intersectionPoint.x = newIntersect->x;
              intersectionPoint.y = newIntersect->y;
            }
          }
        }
        visiblePoints.push_back(intersectionPoint);
      }

      // is reflex vertex AND passes tangentially through the point, staying in the polygon (bouncin)
      else {
        point2d *newIntersect = new point2d;
        intersectionPoint.isLeft = false;
        for (int k = 0; k < polygon.size(); k++) {
          newIntersect->x = endPoint.x;
          newIntersect->y = endPoint.y;
          if (k == j || (k + 1)%polygon.size() == j) {
            continue;
          }
          char intSect = seg_seg_intersect(guards[i].pos, endPoint, polygon[k], polygon[(k+1)%polygon.size()], newIntersect);
          if (intSect == '1' || intSect == 'v') {
            // if the new intersection point is closer than the current one, replace it
            if (distance(guards[i].pos, *newIntersect) < distance(guards[i].pos, intersectionPoint)) {
              intersectionPoint.x = newIntersect->x;
              intersectionPoint.y = newIntersect->y;
              
            }
          }
        }
        if ((polygon[j].x == intersectionPoint.x && polygon[j].y == intersectionPoint.y) || (intersectionPoint.x == endPoint.x && intersectionPoint.y == endPoint.y )) {
          visiblePoints.push_back(polygon[j]);
        }
        else if (distance(guards[i].pos, intersectionPoint) < distance(guards[i].pos, polygon[j])) {
          visiblePoints.push_back(intersectionPoint);
        }
        else if (arePointsLeft) {
          intersectionPoint.isLeft = true;
          visiblePoints.push_back(intersectionPoint);
          intersectionPoint.isLeft = false;
          visiblePoints.push_back(polygon[j]);
        }
        else {
          visiblePoints.push_back(polygon[j]);
          visiblePoints.push_back(intersectionPoint);
        }
      }
      // now, we have the closest intersection point, so we can add it to the visible points

    }

    // now, sort the visible points radially
    priority_queue<radial_point2d> sortedPoints;
    sortRadially(guards[i].pos, visiblePoints, sortedPoints);

    // now, we can draw the visible polygon
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


    // want to make a filled in polygon; therefore, can draw triangles with the guard point, and two points from radially sorted points
    // will do this with gl_triangle_fan

    radial_point2d p1;
    radial_point2d p2 = sortedPoints.top();
    radial_point2d firstPoint = p2;
    sortedPoints.pop();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(guards[i].r, guards[i].g, guards[i].b, 0.5f);

    while (sortedPoints.size() != 0){
      p1 = p2;
      p2 = sortedPoints.top();
      sortedPoints.pop();
      glBegin(GL_POLYGON);
       glVertex2f(guards[i].pos.x, guards[i].pos.y);
       glVertex2f(p1.p.x, p1.p.y);
       glVertex2f(p2.p.x, p2.p.y);
      glEnd();
    }
    glBegin(GL_POLYGON);
       glVertex2f(guards[i].pos.x, guards[i].pos.y);
       glVertex2f(firstPoint.p.x, firstPoint.p.y);
       glVertex2f(p2.p.x, p2.p.y);
      glEnd();
    glDisable(GL_BLEND);

  }
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


  // draw_polygon(poly); 

  //draw a circle where the mouse was last clicked. Note that this
  //point is stored as a global variable and is modified by the mouse handler function 
  if (curMode == DRAW_LINES) {
    draw_circle(mouse_x, mouse_y);
    add_lines_to_polygon(mouse_x, mouse_y); 
    draw_lines_to_polygon();
  }

  if (curMode == DRAW_GAURDS) {
    draw_lines_to_polygon();
    add_guard(mouse_x, mouse_y);
    draw_gaurds();
  }

  if (curMode == RUN) {
    draw_lines_to_polygon();
    timerfunc();
  }

  /* execute the drawing commands */
  glFlush();
}


void findReflexVertices() {
  // loop through all points in polgon, finding out if each vertex is a reflex vertex or not;
  // since we know that the points are all inputting counter-clockwise, just need to see if 
  // a third point is to the right of two given points; if so, then it is a reflex vertex
  for (unsigned int i = 0; i < polygon.size(); i++) {
    // get the two points that are adjacent to the current point
    point2d p1;
    point2d p2;
    if (i == 0) {
      p1 = polygon[polygon.size() - 1];
      p2 = polygon[i + 1];
    } else if (i == polygon.size() - 1) {
      p1 = polygon[i - 1];
      p2 = polygon[0];
    } else {
      p1 = polygon[i - 1];
      p2 = polygon[i + 1];
    }

    // now, we have the two points that are adjacent to the current point
    // we need to see if the third point is to the right of the line formed by the first two points
    // if so, then it is a reflex vertex
    if (isPointToRightOfLine(p1, polygon[i], p2)) {
      // it is a reflex vertex
      polygon[i].isReflex = true;
    }
  }
}

/* ****************************** */
void keypress(unsigned char key, int x, int y) {
  switch(key) {
  case 'q':	
    exit(0);
    break;
  case 's':
    polygon.clear();
    curMode = DRAW_LINES;
    break;
  case 'e':
    if (is_simple_polygon()) {
      cout << "polygon is simple\n";
    } else {
      cout << "polygon is not simple\n";
      exit(1);
    }
    findReflexVertices();
    curMode = DRAW_GAURDS;
    break;
  case 'g':
    checkIfGuardInPolygon();
    curMode = RUN;
    glutPostRedisplay();
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
  if (curMode == RUN) {
    move_guards();
    draw_gaurds();
    calculate_visible_polygons();
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
