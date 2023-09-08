/* viewpoints.cpp, Laura Toma
   
   What it does: Draws a set of points in the default 2D projection.
   Includes a tentative function for printing and drawing a list of
   points (assumed to be a convex hull).
   
   This code is provided as a startup for your 2d hull.  Change it as
   needed to work with your project.
*/

#include "geom.h"
#include "rtimer.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <strings.h>

//to compile on both apple and unix platform
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <vector>
using namespace std; 




//pre-defined colors for convenience 
GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat cyan[3] = {0.0, 1.0, 1.0};
/* from https://www.opengl.org/discussion_boards/showthread.php/132502-Color-tables  */
GLfloat brown[3] = { 0.647059, 0.164706, 0.164706}; 
GLfloat DarkBrown[3] = { 0.36, 0.25, 0.20}; 
GLfloat DarkTan[3] = { 0.59, 0.41, 0.31};
GLfloat Maroon[3]= { 0.556863, 0.137255, 0.419608}; 
GLfloat DarkWood[3] = { 0.52, 0.37, 0.26}; 
GLfloat  Copper[3] = { 0.72,  0.45,  0.20};
GLfloat green1[3] = {.5, 1, 0.5};
GLfloat green2[3] = {0.0, .8, 0.0};
GLfloat green3[3] = {0.0, .5, 0.0};
GLfloat ForestGreen[3] = { 0.137255, 0.556863, 0.137255};
GLfloat MediumForestGreen[3] = { 0.419608 , 0.556863 , 0.137255}; 
GLfloat LimeGreen[3] ={ 0.196078,  0.8 , 0.196078}; 
GLfloat Orange[3] = { 1, .5, 0}; 
GLfloat Silver[3] = { 0.90, 0.91, 0.98};
GLfloat Wheat[3] = { 0.847059 , 0.847059, 0.74902}; 





/* global variables */

//desired number of points, entered by the user on the command line
int NPOINTS;

//the vector of points
//note: needs to be global in order to be rendered
vector<point2d>  points;

//the convex hull 
//needs to be global in order to be rendered
vector<point2d>  hull; 



//window size for the graphics window
const int WINDOWSIZE = 500; 

/* currently there are 4 different ways to initialize points.  The
   user can cycle through them by pressing 'i'. Check out the display()
   function to see how that's implemented.
*/
int NB_INIT_CHOICES = 13; 
int  POINT_INIT_MODE = 0; //the first inititalizer





/********************************************************************/
/* forward declarations of functions */

//print label, then the vector 
void print_vector(const char* label, vector<point2d> p); 



/* render the array of points stored in global variable points.
   Each point is drawn as a small square.  */
void draw_points(vector<point2d> pts); 

/* Render the hull; the points on the hull are expected to be in
   boundary order (either ccw or cw), otherwise it will look
   zig-zaagged.  
*/
  void draw_hull(vector<point2d> hull); 


void display(void);
void keypress(unsigned char key, int x, int y);

// initializer function
void initialize_points_circle(vector<point2d>& pts, int n); 
void initialize_points_horizontal_line(vector<point2d>&pts, int n);
void initialize_points_random(vector<point2d>&pts, int n) ;
void initialize_points_cross(vector<point2d>&pts, int n) ;
void initialize_points_spiral(vector<point2d>&pts, int n) ;
void initialize_points_six_point_star(vector<point2d>&pts, int n) ;
void initialize_points_four_clusters(vector<point2d> &pts, int n) ;
void initialize_points_spaceship(vector<point2d>& pts, int n) ;
void initialize_points_polar(vector<point2d>& pts, int n) ;
void initialize_points_balloon(vector<point2d>& pts, int n) ;
void initialize_points_heart(vector<point2d>& pts, int n) ;
void initialize_points_multiple_lowest(vector<point2d> &pts, int n) ;
void initialize_intermediate_collinear(vector<point2d> &pts, int n) ;





//you'll add more 


/********************************************************************/


/* ****************************** */
/* Initializes pts with an intermediate collinear point on the hull
   and duplicate points
*/
void initialize_intermediate_collinear(vector<point2d> &pts, int n) {
  printf("\ninitialize points intermediate collinear\n");
  pts.clear();

  point2d far;
  far.x = 365;
  far.y = 251;
  pts.push_back(far);

  point2d q;
  q.x = 147;
  q.y = 400;
  pts.push_back(q);

  point2d r;
  r.x = 147;
  r.y = 400;
  pts.push_back(r);

  point2d s;
  s.x = 145 + 22;
  s.y = 101 + 15;
  pts.push_back(s);

  point2d t;
  t.x = 145;
  t.y = 101;
  pts.push_back(t);

  point2d x;
  x.x = 258;
  x.y = 261;
  pts.push_back(x);

  point2d z;
  z.x = 287;
  z.y = 400;
  pts.push_back(z);

  point2d k;
  k.x = 88;
  k.y = 201;
  pts.push_back(k);

  // Initializes n random points inside of the hull
  for (int i = 0; i < n; i++) {
    point2d pt;
    pt.y = (random() % 129) + 251;
    pt.x = (random() % 140) + 147;
    pts.push_back(pt);
  }
}

/* ****************************** */
/* Initializes pts with multiple points of the same lowest y value.
  Initializes the middle low point first to check that collinear points
  are also handled correctly.
*/
void initialize_points_multiple_lowest(vector<point2d> &pts, int n) {
  printf("\ninitialize points multiple lowest\n");
  pts.clear();

  // Create three lowest points with the same y value
  point2d first_low;
  first_low.x = 100;
  first_low.y = 24;
  pts.push_back(first_low);

  point2d second_low;
  second_low.x = 110;
  second_low.y = 24;
  pts.push_back(second_low);

  point2d third_low;
  third_low.x = 80;
  third_low.y = 24;
  pts.push_back(third_low);

  // Create random points above the lowest points
  for (int i = 0; i < n; i++) {
    point2d p;
    p.x = (int)(.3 * WINDOWSIZE) / 2 + random() % ((int)(.7 * WINDOWSIZE) + 24);
    p.y = (int)(.3 * WINDOWSIZE) / 2 + random() % ((int)(.7 * WINDOWSIZE));
    pts.push_back(p);
  }
}


void initialize_points_heart(vector<point2d>& pts, int n) {
  printf("\ninitialize points heart with stem and leaf\n");
  pts.clear();
  double interval = M_PI / n;
  double t = 0;
  point2d p;
  for (int i = 0; i < n; i++) {
    double x = sin(t) * i * WINDOWSIZE / (2 * n) + (WINDOWSIZE / 2);
    double y = cos(t) * i * WINDOWSIZE / (2 * n) + (WINDOWSIZE / 2);
    p.x = x;
    p.y = y+100;
    pts.push_back(p);
    t += interval;
  }
  t = M_PI;
  for (int i = 0; i < n; i++) {
    double x = sin(t) * i * WINDOWSIZE / (2 * n) + (WINDOWSIZE / 2);
    double y = -cos(t) * i * WINDOWSIZE / (2 * n) + (WINDOWSIZE / 2);
    p.x = x;
    p.y = y+100;
    pts.push_back(p);
    t += interval;
  }
  pts.push_back(p);
}

void initialize_points_balloon(vector<point2d>& pts, int n) {
  printf("\ninitialize points balloon\n");
  //clear the vector just to be safe
  pts.clear();
  point2d p;
  float a = (float)WINDOWSIZE / 4;
  float b = (float)WINDOWSIZE / 8;
  //allocate a percentage of the n points to the balloon and to the string
  int balloon_n = n*0.9;
  int string_n = n - balloon_n;
  for (int i=0; i<balloon_n; i++) {
    float angle = 2 * M_PI * i / balloon_n;
    p.x = b * cos(angle) + (WINDOWSIZE / 2);
    p.y = a * sin(angle) + (WINDOWSIZE / 2) + 100;
    pts.push_back(p);
  }
  // Find the point closest to the center of the balloon
  int closest_idx = 0;
  float closest_distance = INT_MAX;
  for (int i = 0; i < balloon_n; i++) {
    float distance = sqrt(pow(pts[i].x - (WINDOWSIZE / 2), 2) + pow(pts[i].y - (WINDOWSIZE / 2) + 100, 2));
    if (distance < closest_distance) {
      closest_distance = distance;
      closest_idx = i;
    }
  }
  int bottom_x = pts[closest_idx].x;
  int bottom_y = pts[closest_idx].y;
  int string_spacing = ((WINDOWSIZE-30)/2)/string_n;
  for (int i=0; i<string_n; i++) {
    p.x = bottom_x;
    p.y = bottom_y - string_spacing*i;
    pts.push_back(p);
  }
}


/* ****************************** */
/* Initializes pts in the shape of a space ship. 
The points are in the range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_spaceship(vector<point2d>& pts, int n) {
  printf("\ninitialize points spaceship\n"); 
  //clear the vector just to be safe 
  pts.clear(); 
  n = n/9; //we'll generate a space ship
  double  step = 2* M_PI/n; 
  int radius = 10; 
  point2d p; 
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ radius*cos(i*step) - 8*radius; 
    p.y = WINDOWSIZE/2+ radius*sin(i*step); 
    pts.push_back(p); 
  }
 
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ radius*cos(i*step) -4*radius; 
    p.y = WINDOWSIZE/2+ radius*sin(i*step) - 1.2*(radius); 
    pts.push_back(p); 
  }
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ radius*cos(i*step); 
    p.y = WINDOWSIZE/2+ radius*sin(i*step) - 5*(radius)/3; 
    pts.push_back(p); 
  }
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ radius*cos(i*step) + 4*radius; 
    p.y = WINDOWSIZE/2+ radius*sin(i*step) - 1.2*(radius); 
    pts.push_back(p); 
  }
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ radius*cos(i*step) + 8*radius; 
    p.y = WINDOWSIZE/2+ radius*sin(i*step); 
    pts.push_back(p); 
  }
  
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ 22*radius*cos(i*step)*cos(i*step) - 11*radius;
    p.y = WINDOWSIZE/2+ 7*radius*cos(i*step)*sin(i*step);
    pts.push_back(p); 
  }
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ 5*radius*cos(i*step/2);
    p.y = WINDOWSIZE/2+ 4*radius*sin(i*step/2)+2*radius;
    pts.push_back(p); 
  }
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ 3*radius*cos(i*step)+2*radius;
    p.y = WINDOWSIZE/2 -1.1*radius*sin(i*step/2)+2*radius;
    pts.push_back(p); 
  }
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2 -3*radius*cos(i*step)-2*radius;
    p.y = WINDOWSIZE/2 -1.1*radius*sin(i*step/2)+2*radius;
    pts.push_back(p); 
  }
}
/* ****************************** */
/* Initializes pts in cool polar shapes. 
The points are in the range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_polar(vector<point2d>& pts, int n) {
  printf("\ninitialize points polar\n"); 
  //clear the vector just to be safe 
  pts.clear(); 
  n = n/10; //we'll generate a space ship
  double  step = 2 * M_PI/n; 
  
  point2d p;
  for (int i = 0; i < 5*n; i++) {
    double r = 3*(sin(2.4 * i * step)*sin(2.4 * i * step))+3*(cos(2.4*i*step)*cos(2.4*i*step)*cos(2.4*i*step)*cos(2.4*i*step));
    p.x = WINDOWSIZE/2 + 50*r*cos(i*step);
    p.y = WINDOWSIZE/2 + 50*r*sin(i*step);
    pts.push_back(p);
  }
  for (int i = 0; i < 5*n; i++) {
    double r = 3*(sin(1.2 * i * step)*sin(1.2 * i * step))+3*(cos(6*i*step)*cos(6*i*step)*cos(6*i*step));
    p.x = WINDOWSIZE/2 + 33*r*cos(i*step);
    p.y = WINDOWSIZE/2 + 33*r*sin(i*step);
    pts.push_back(p);
  }
}


/*
Initializes points in 4 clusters in the corners, with two lines intersecting
in the middle to make an 'X'

Author: Drew Hofer & Deven Kanwal
*/
void initialize_points_four_clusters(vector<point2d> &pts, int n) {
  pts.clear();
  for (int i = 0; i < n; i++) {
    int section;
    point2d p;
    section = random() % 6;
    if (section == 0) {
      p.x = 100 + random() % 40;
      p.y = 100 + random() % 40;
      pts.push_back(p);
    }
    if (section == 1) {
      p.x = 400 + random() % 40;
      p.y = 100 + random() % 40;
      pts.push_back(p);
    }
    if (section == 2) {
      p.x = 100 + random() % 40;
      p.y = 400 + random() % 40;
      pts.push_back(p);
    }
    if (section == 3) {
      p.x = 400 + random() % 40;
      p.y = 400 + random() % 40;
      pts.push_back(p);
    }
    if (section == 4) {
      p.x = (rand() % (301)) + 100;
      p.y = p.x;
      pts.push_back(p);
    }
    if (section == 5) {
      p.x = (rand() % (301)) + 100;
      p.y = (p.x * -1) + 530;
      pts.push_back(p);
    }
  }
}





/* ****************************** */
/* Initializes pts with n points on two circles.  The points are in the
   range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_circle(vector<point2d>& pts, int n) {

  printf("\ninitialize points circle\n"); 
  //clear the vector just to be safe 
  pts.clear(); 

  n = n/2; //we'll generaate two circles, n/2 points each
  double  step = 2* M_PI/n; 
  int radius = 100; 

  point2d p; 
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ radius*cos(i*step); 
    p.y = WINDOWSIZE/2+ radius*sin(i*step); 
    pts.push_back(p); 
  }

  radius /= 2; 
  for (int i=0; i<n; i++) {
    p.x = WINDOWSIZE/2+ radius*cos(i*step); 
    p.y = WINDOWSIZE/2+ radius*sin(i*step); 
    pts.push_back(p); 
  }
  
}





/* ****************************** */
/* Initializes pts with n points on a line.  The points are in the
   range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_horizontal_line(vector<point2d>& pts, int n) {

  printf("\ninitialize points line\n"); 
  //clear the vector just to be safe 
  pts.clear(); 
  
  point2d p; 
  for (int i=0; i<n; i++) {
    p.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
    p.y =  WINDOWSIZE/2; 
    pts.push_back(p); 
  }
}




/* ****************************** */
/* Initializes pts with n random points.  The points are in the
   range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_random(vector<point2d>& pts, int n) {

   printf("\ninitialize points random\n"); 
  //clear the vector just to be safe 
  pts.clear(); 

  point2d p; 
  for (int i=0; i<n; i++) {
    p.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
    p.y =  (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
    pts.push_back(p); 
  }
}




/* ****************************** */
/* Initializes pts with n points on a cross-like shape.  The points are
   in the range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_cross(vector<point2d>& pts, int n) {
  
  printf("\ninitialize points cross\n"); 
  //clear the vector just to be safe 
  pts.clear(); 
  
  point2d p; 
  for (int i=0; i<n; i++) {
    if (i%2 == 0) {
      
      p.x = (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE)); 
      p.y =  random() % ((int)(.7*WINDOWSIZE))  / 5;
      p.y += (int)((1-.7/5)*WINDOWSIZE/2);
    };
    if (i%2 == 1)  {
      
      p.x = random() % ((int)(.7*WINDOWSIZE)) / 5; 
      p.x +=  (int)((1-.7/5)*WINDOWSIZE/2);
      p.y =  (int)(.3*WINDOWSIZE)/2 + random() % ((int)(.7*WINDOWSIZE));
    }
   
    pts.push_back(p); 
    
  }//for i

}


/* ****************************** */
/* Initializes pts with n points in a spiral pattern.  The points are
   in the range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_spiral(vector<point2d>& pts, int n) {
  
  printf("\ninitialize points spiral\n"); 
  //clear the vector just to be safe 
  pts.clear(); 

  // random value for creating some randomness in the spiral
  int r = random() % 10;

  // radius of the spiral, to increase as points are added
  point2d p;
  for (int i=0; i<n; i++) {
    // increase the radius
    double radius = double(i)/(1.2*n) * WINDOWSIZE/2;
    // calculate the angle of the point
    double angle = (double)3*i/n * 2 * M_PI;
    // calculate the x and y coordinates of the point
    r = random() % 30 - 15;
    p.x = WINDOWSIZE/2 + radius * cos(angle) + r;
    r = random() % 30 - 15;
    p.y = WINDOWSIZE/2 + radius * sin(angle) + r;
    pts.push_back(p);
  }
  
  
}


/* ****************************** */
/* Initializes pts with n points in a six point star pattern.  The points are
   in the range (0,0) to (WINSIZE,WINSIZE).
*/ 
void initialize_points_six_point_star(vector<point2d>& pts, int n) {
  
  printf("\ninitialize points six point star\n"); 
  //clear the vector just to be safe 
  pts.clear(); 

  // random value for creating some randomness in the spiral
  int r = random() % 10;

  // current radius
  double radius = WINDOWSIZE/8;

  // radius of the spiral, to increase as points are added
  point2d p;
  for (int i=0; i<n; i++) {
    if (i % (n/6) < n/12) {
      // increase the radius
      radius += (3*(WINDOWSIZE/(8))/(n/12));
    }
    else {
      // decrease the radius
      radius -= (3*(WINDOWSIZE/(8))/(n/12));
    }

    // calculate the angle of the point
    double angle = (double)i/n * 2 * M_PI;
    // calculate the x and y coordinates of the point
    r = random() % 20 - 10;
    p.x = WINDOWSIZE/2 + radius * cos(angle) + r;
    r = random() % 20 - 10;
    p.y = WINDOWSIZE/2 + radius * sin(angle) + r;
    pts.push_back(p);
  }
}



/* ****************************** */
/* print the vector of points */
void print_vector(const char* label, vector<point2d> points) {
  
  printf("%s ", label);
  for (int i=0; i< points.size(); i++) {
    printf("[%3d,%3d] ", points[i].x, points[i].y);
  }
  printf("\n");
}





/* ****************************** */
int main(int argc, char** argv) {

  //read number of points from user
  if (argc!=2) {
    printf("usage: viewPoints <nbPoints>\n");
    exit(1); 
  }
  NPOINTS = atoi(argv[1]); 
  printf("you entered n=%d\n", NPOINTS);
  assert(NPOINTS >0); 

  //populate the points 
  initialize_points_random(points, NPOINTS);
  //print_vector("points:", points);

  //compute the convex hull 
  Rtimer rt1; 
  rt_start(rt1); 
  graham_scan(points, hull); 
  rt_stop(rt1); 
  print_vector("hull:", hull);
  
  //print the timing 
  char buf [1024]; 
  rt_sprint(buf,rt1);
  printf("hull time:  %s\n\n", buf);
  fflush(stdout); 

 
  //start the rendering 
  /* initialize GLUT  */
  // cout argc and argv
  printf("argc: %d\n", argc);
  for (int i=0; i<argc; i++) {
    printf("argv[%d]: %s\n", i, argv[i]);
  }
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
  glScalef(2.0/WINDOWSIZE, 2.0/WINDOWSIZE, 1.0);  
  //first translate the points to [-WINDOWSIZE/2, WINDOWSIZE/2]
  glTranslatef(-WINDOWSIZE/2, -WINDOWSIZE/2, 0); 
 
  draw_points(points);
  draw_hull(hull); 

  /* execute the drawing commands */
  glFlush();
}




/* ****************************** */
/* draw the points. each point is drawn as a small square
*/
void draw_points(vector<point2d> points){

  const int R= 1;
  //draw polygon filled or line 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  //set drawing color 
  glColor3fv(yellow);   
  
  for (int i=0; i< points.size(); i++) {
    //draw a small square centered at (points[i].x, points[i].y)
    glBegin(GL_POLYGON);
    glVertex2f(points[i].x -R,points[i].y-R);
    glVertex2f(points[i].x +R,points[i].y-R);
    glVertex2f(points[i].x +R,points[i].y+R);
    glVertex2f(points[i].x -R,points[i].y+R);
    glEnd();
  }
} //draw_points 





/* ****************************** */
/* Draaaw the hull; the points on the hull are expected to be in
   boundary order (either ccw or cw) or else it will look
   zig-zaagged. To render the hull we'll draw lines between
   consecutive points */
void draw_hull(vector<point2d> hull){

  //set color 
  glColor3fv(red);   //this should be a constant
  
  if (hull.size() >0) {
    int i; 
    for (i=0; i< hull.size()-1; i++) {
      
      //draw a line from  i to i+1
      glBegin(GL_LINES);
      glVertex2f(hull[i].x, hull[i].y); 
      glVertex2f(hull[i+1].x, hull[i+1].y); 
      glEnd();
    }
    
    //draw a line from the last point to the first point
    i =  hull.size()-1; 
    glBegin(GL_LINES);
    glVertex2f(hull[i].x, hull[i].y); 
    glVertex2f(hull[0].x, hull[0].y); 
    glEnd();
  }//if (hull not empty)
}



/* ****************************** */
/* Handler for key presses. called whenever a key is spressed */
void keypress(unsigned char key, int x, int y) {
  switch(key) {
  case 'q':
    exit(0);
    break;

  case 'i':
    //when the user presses 'i', we want to re-initialize the points and
    //recompute the hull
    POINT_INIT_MODE = (POINT_INIT_MODE+1) % NB_INIT_CHOICES; 
    switch (POINT_INIT_MODE) {
    case 0: 
      initialize_points_circle(points, NPOINTS); 
      break; 
    case 1: 
      initialize_points_cross(points, NPOINTS); 
      break; 
    case 2: 
      initialize_points_horizontal_line(points, NPOINTS); 
      break; 
    case 3: 
      initialize_points_random(points, NPOINTS); 
      break; 
    case 4: 
      initialize_points_spiral(points, NPOINTS);
      break;
    case 5:
      initialize_points_six_point_star(points, NPOINTS);
      break;
    case 6:
      initialize_points_four_clusters(points, NPOINTS);
      break;
    case 7:
      initialize_points_spaceship(points, NPOINTS);
      break;
    case 8:
      initialize_points_polar(points, NPOINTS);
      break;
    case 9:
      initialize_points_balloon(points, NPOINTS);
      break;
    case 10:
      initialize_points_heart(points, NPOINTS);
      break;
    case 11:
      initialize_points_multiple_lowest(points, NPOINTS);
      break;
    case 12:
      initialize_intermediate_collinear(points, NPOINTS);
      break;
    } //switch 
    //we changed the points, so we need to recompute the hull
    graham_scan(points, hull); 

    //we changed stuff, so we need to tell GL to redraw
    glutPostRedisplay();

  } //switch (key)

}//keypress



/* Handler for window re-size event. Called back when the window first appears and
   whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
     
   // Set the viewport to cover the new window
   glViewport(0, 0, width, height);
 
   glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
   glLoadIdentity();             // Reset
   gluOrtho2D(0.0, (GLdouble) width, 0.0, (GLdouble) height); 
}


