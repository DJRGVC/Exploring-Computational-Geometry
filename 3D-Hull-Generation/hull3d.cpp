
/* hull3d.cpp

   This code is provided as a startup for your 3d hull.  Change it as
   needed to work with your project. 
   OpenGL 1.x
   Laura Toma
*/
#include "geom.h"

#include <OpenGL/OpenGL.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

//this allows this code to compile both on apple and linux platforms
#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif


#include <vector>

//to read the mesh from a file 
#include <fstream>
#include <iostream>
#include <sstream>
#include <random>

#include<chrono>
#include<thread>

using namespace std; 

// points have coordinates in range -RANGE..RANGE
int  RANGE = 1000; 



//pre-defined colors for convenience 
GLfloat red[3] = {1.0, 0.0, 0.0};
GLfloat green[3] = {0.0, 1.0, 0.0};
GLfloat blue[3] = {0.0, 0.0, 1.0};
GLfloat black[3] = {0.0, 0.0, 0.0};
GLfloat white[3] = {1.0, 1.0, 1.0};
GLfloat yellow[3] = {1.0, 1.0, 0.0};
GLfloat gray[3] = {0.5, 0.5, 0.5};
GLfloat lightgrey[3] = { 0.8, 0.8, 0.8};
GLfloat darkgrey[3] = { 0.6, 0.6, 0.6};
GLfloat magenta[3] = {1.0, 0.0, 1.0};
GLfloat darkmagenta[3] = { 0.5, 0.0, 0.5};
GLfloat cyan[3] = {0.0, 1.0, 1.0};
GLfloat darkcyan[3] = { 0.0, 0.4, 0.4};
GLfloat lightcyan[3] = { 0.0, 0.6, 0.6};




/********************************************************************/
/* global variables */
/********************************************************************/



int n;  //desired number of points

//the vector  of n points; note: needs to be global in order to be rendered
vector<point3d>  points;

//the convex hull. note: needs to be global in order to be  rendered
vector<triangle3d>  hull; 


///window size for the graphics window
const int WINDOWSIZE = 500; 

//global translation matrix. updated when the user translates the scene 
GLfloat pos[3] = {0,0,0};

//global rotation matrix. updated when the user translates the scene 
GLfloat theta[3] = {0,0,0};

// draw polygons line or filled.  
bool fill_mode = 0; 

//render the hull one face at a time, continuously 
bool animate_mode = false; 

//render the hull one face at a time, continuously 
bool dot_mode = true; 

bool background_mode = false;

//used by animate 
int hull_index_render = 0;










/********************************************************************/
/* forward declarations of functions */
/********************************************************************/

void display(void);
void keypress(unsigned char key, int x, int y);
void idlefunc();

// generate n random points with int coords in -RANGE..RANGE
void initialize_points_random_RANGE(vector<point3d>& points, int n);
// generate n  points on a sphere  with int coords in -RANGE..RANGE
void initialize_points_sphere_RANGE(vector<point3d>& points, int n, double rad);
// initialize n points from a mesh 
int initialize_points_from_mesh(vector<point3d>& pts, char* fpath); 

				
//print label, then the vector 
void print_points(const char* label, vector<point3d> p);

// render the points; each point is drawn as a small square. 
void draw_points(vector<point3d> pts);

// render the hull 
void draw_hull(vector<triangle3d>& hull);


//various helper functions 
void draw_xy_rect(GLfloat z, GLfloat* col); 
void draw_xz_rect(GLfloat y, GLfloat* col); 
void draw_yz_rect(GLfloat x, GLfloat* col); 
void cube(GLfloat side); 
void filledcube(GLfloat side); 
void draw_axes(); 




/* ************************************************************ */
int main(int argc, char** argv) {

  printf("hello!");
  //read number of points from user
  if (argc!=2) {
    printf("usage: hull3d <nbPoints>\n");
    exit(1); 
  }
  n = atoi(argv[1]); 
  printf("you entered n=%d\n", n);
  assert(n>0); 

  //populate the points 
  initialize_points_random_RANGE(points, n);
  //print_points("points:", points);

  //compute the hull 
  naive_hull(points, hull);
  //gift_wrapping_hull(points, hull); 
  //print_hull(hull);

  
  /* OPEN GL STUFF */
    /* initialize GLUT and open a window */
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(WINDOWSIZE, WINDOWSIZE);
  glutInitWindowPosition(100,100);
  glutCreateWindow(argv[0]);

  /* register callback functions */
  glutDisplayFunc(display); 
  glutKeyboardFunc(keypress);
  glutIdleFunc(idlefunc);
   
  /* GL init */
  /* set background color black*/
  glClearColor(0, 0, 0, 0);  
  //when depth test is enabled, GL determines which objects are in
  //front/behind and renders them correctly
  glEnable(GL_DEPTH_TEST); 

  // setup the projection  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60, 1 /* aspect */, 1, 10.0);
  /* the frustrum is from z=-1 to z=-10 */
  /* camera is at (0,0,0) looking along negative y axis */
  
  //initialize the translation to bring the points in the view
  //frustrum which is [-1, -10]
  pos[2] = -3;


  /* give control to the event handler */
  glutMainLoop();

  return 0;
}





/* ************************************************************ */
/* This is the function that renders the window. We registered this
   function as the "displayFunc". It will be called by GL everytime
   the window needs to be rendered.
 */
void display(void) {
  
  //clear the screen
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* The default GL window is x=[-1,1], y= [-1,1] with the origin in
     the center.  The view frustrum was set up from z=-1 to z=-10. The
     camera is at (0,0,0) looking along negative z axis.
  */
  
  //clear all modeling transformations
  glMatrixMode(GL_MODELVIEW); 
  glLoadIdentity();

  /*  First we translate and rotate our local reference system with the
      user transformation. pos[] represents the cumulative translation
      entered by the user, and theta[] the cumulative rotation entered
      by the user */
  glTranslatef(pos[0], pos[1], pos[2]);  
  glRotatef(theta[0], 1,0,0); //rotate theta[0] around x-axis, etc 
  glRotatef(theta[1], 0,1,0);
  glRotatef(theta[2], 0,0,1);
   
  //draw a cube, nice for perspective 
  cube(1); 

  //our points are in the range [-RANGE, RANGE], scale them to screen coords [-1, 1]
  float s = 1/(float)RANGE; 
  glScalef(s, s, s);

  draw_points(points);
  draw_hull(hull); 
  
  // execute the drawing commands
  glFlush();
}


/* this is the function we registered as idleFunc. It's called every
   frame or so.  */
void idlefunc() {
  if (animate_mode) {
    hull_index_render++;
    glutPostRedisplay();
  }
}


/* ************************************************************ */
/* This is the function that handles key presses, We registered this
   function as the keypressFunc.  It will be called by GL whenever a
   key is pressed */
void keypress(unsigned char key, int x, int y) {
  
  switch(key) {
    
  case 'i': 
    //re-initialize  RANDOM
    initialize_points_random_RANGE(points, n); 
    //re-compute the hull
    // naive_hull(points, hull);
    hull_index_render = 0;
    giftwrapping_hull(points, hull); 
    glutPostRedisplay(); 

    break; 
    
  case 's': 
    //re-initialize SPHERE 
    initialize_points_sphere_RANGE(points, n, .8); 
    //re-compute the hull
    // naive_hull(points, hull);
    hull_index_render = 0;
    giftwrapping_hull(points, hull); 
    glutPostRedisplay(); 
    break; 
    
  case 'm':
    //re-initialize from mesh 
    if (initialize_points_from_mesh(points, "./meshes/spot/spot_control_mesh.obj")) {
      //re-compute the hull
      // naive_hull(points, hull);
      giftwrapping_hull(points, hull); 
      glutPostRedisplay(); 
    }
    break;
    
  case 'c': 
    //fillmode
    fill_mode = !fill_mode; 
    glutPostRedisplay();
    break;
    
  case 'a': //animate 
    animate_mode = !animate_mode;
    glutPostRedisplay();
    break;

  case 'o': //dot mode
    dot_mode = !dot_mode;
    glutPostRedisplay();
    break;

  case 'p': //dot mode
    background_mode = !background_mode;
    glutPostRedisplay();
    break;
    
    //ROTATIONS
  case 'x':
    theta[0] += 5.0; 
    glutPostRedisplay();
    break;
  case 'y':
    theta[1] += 5.0;
    glutPostRedisplay();
    break;
  case 'z':
    theta[2] += 5.0;
    glutPostRedisplay();
    break;
  case 'X':
    theta[0] -= 5.0; 
    glutPostRedisplay();
    break;
  case 'Y':
    theta[1] -= 5.0; 
    glutPostRedisplay();
    break;
  case 'Z':
    theta[2] -= 5.0; 
    glutPostRedisplay();
    break;

    //TRANSLATIONS 
case 'b': //backward (zoom out)
    pos[2] -= 0.1; 
    glutPostRedisplay();
    break;
    //forward (zoom in)
  case 'f': //forward (zoom in) 
    pos[2] += 0.1; 
    glutPostRedisplay();
    break;
  case 'd':      //down 
    pos[1] -= 0.1; 
    glutPostRedisplay();
    break;
    //up
  case 'u':  //up
    pos[1] += 0.1; 
    glutPostRedisplay();
    break;
  case 'l': //left 
    pos[0] -= 0.1; 
    glutPostRedisplay();
    break;
  case 'r': //right 
    pos[0] += 0.1; 
    glutPostRedisplay();
    break;

  case '0':
    //reset all transformations (back to original position)
    pos[0] = pos[1] = 0; pos[2] = -3;
    theta[0] = theta[1] = theta[2] = 0;
    glutPostRedisplay();
    break;
    
  case 'q':
    exit(0);
    break;
  } 

}//keypress






/* ************************************************************** */
/* populate the vector with n random points with int coords in the
   range [-RANGE, RANGE]
 */
void initialize_points_random_RANGE(vector<point3d>& points, int n) {
  printf("initialize points random\n"); 
  //clear the vector just to be safe 
  points.clear(); 

  srand((unsigned)time(0)); 
      
  point3d p;
  for (int i=0; i<n; i++) {
    p.x = (int) rand()% RANGE  - (RANGE/2); 
    p.y = (int) rand()% RANGE  - (RANGE/2); 
    p.z = (int) rand()% RANGE  - (RANGE/2); 
    //cout << "point: " << p.x << ","<< p.y << "," <<  p.z << endl; 
    points.push_back(p); 
  }
}




/* ************************************************************** */
/* populate the vector with n points on a sphere with int coords in
   the range [-RANGE, RANGE]
 */
void initialize_points_sphere_RANGE(vector<point3d>& points, int n, double rad) {

  printf("initialize points sphere\n"); 
  //clear the vector just to be safe 
  points.clear();

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dist(0, 2*M_PI);
  
  double phi, theta;
  point3d p;
  for (int i=0; i<n/2; i++) {
    phi = dist(gen); 
    theta = dist(gen);
    p.x = RANGE * rad *sin(phi) * cos(theta);
    p.y = RANGE*rad * sin(phi) * sin(theta);
    p.z = RANGE*rad * cos(phi); 
    //cout << p.x << ","<< p.y << "," <<  p.z << endl; 
    points.push_back(p); 
  }

  //let's add a second smaller sphere inside first one
  rad = rad/2; 
  for (int i=0; i<n/2; i++) {
    phi = dist(gen); 
    theta = dist(gen);
    p.x = RANGE * rad *sin(phi) * cos(theta);
    p.y = RANGE*rad * sin(phi) * sin(theta);
    p.z = RANGE*rad * cos(phi); 
    //cout << p.x << ","<< p.y << "," <<  p.z << endl; 
    points.push_back(p); 
  }
}







/* ************************************************************ */
/* credit: Juan Atehortua and Lily Smith, fall 2021 */
int initialize_points_from_mesh(vector<point3d>& pts, char* fpath) {
  printf("initialize points mesh\n"); 
  pts.clear();
  
  string line;
  ifstream mesh;
  //mesh.open("./meshes/spot/spot_control_mesh.obj", ios::in);
  mesh.open(fpath, ios::in);
  point3d point;
  if (mesh.is_open()) {
    while (getline(mesh, line)) {
      if (line.substr(0, 2) == ("v ")){
        istringstream v(line.substr(2));
        v >> point.x;
        v >> point.y;
        v >> point.z;
	point.x *= RANGE;
	point.y *= RANGE;
	point.z *= RANGE;
        pts.push_back(point);
      }//if
    }//while 
  } else {
    cerr << "Could not open file.\n";
    return 0;
  }
  mesh.close();
  return 1;
}

  
  
/* ******************************************************** */
/* Draw the array of points, Each point is drawn as a small cube.
   The points are in the range x, y, z in [-1,1],
*/
  void draw_points(vector<point3d> pts){
    if (dot_mode) {
    //filled 
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //set color 
    glColor3fv(yellow);   
    
    for (int i=0; i < points.size(); i++) {
      
      //draw small filled cube at (points[i].x, points[i].y, points[i].z)
      
      //first save local coordinate system 
      glPushMatrix(); 
      //translate our local coordinate system to the point that we want to draw
      glTranslatef(points[i].x, points[i].y, points[i].z);
      //draw the cube 
      filledcube(5); 
      //go  back to where we were
      glPopMatrix(); 
    } //for 

    }

    
  }//draw_points
  
  
  



/* ********************************************************* */
// draw the faces of the hull.
void draw_hull(vector<triangle3d>& hull){
  
  //if there is no hull, skip rendering 
  if (hull.size() == 0) return;
  
  //set the fill mode 
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  
  //N represents how many faces to render 
  if (hull_index_render >= hull.size()) hull_index_render = 0;
  int N = (animate_mode) ? hull_index_render : hull.size();

  for (int i=0; i<N; i++) {

    //set color of this face
    GLfloat color[3]; 
    for (int x=0; x<3; x++) color[x] = (GLfloat)hull[i].color[x];
    glColor3fv(color); 

    //draw the triangle
    glBegin(GL_TRIANGLES);
    glVertex3f(hull[i].a->x, hull[i].a->y,hull[i].a->z);
    glVertex3f(hull[i].b->x, hull[i].b->y,hull[i].b->z);
    glVertex3f(hull[i].c->x, hull[i].c->y,hull[i].c->z);
    glEnd();

  }//for 
}//draw-hull 


//draw a square x=[-side,side] x y=[-side,side] at depth z
void draw_xy_rect(GLfloat z, GLfloat side, GLfloat* color) {

  if (background_mode) {
    glColor3fv(color);
    glBegin(GL_POLYGON);
    glVertex3f(-side,-side, z);
    glVertex3f(-side,side, z);
    glVertex3f(side,side, z);
    glVertex3f(side,-side, z);
    glEnd();
  }

}

//draw a square y=[-side,side] x z=[-side,side] at given x
void draw_yz_rect(GLfloat x, GLfloat side, GLfloat* color) {
  
  if (background_mode) {
    glColor3fv(color);
    glBegin(GL_POLYGON);
    glVertex3f(x,-side, side);
    glVertex3f(x,side, side);
    glVertex3f(x,side, -side);
    glVertex3f(x,-side, -side);
    glEnd();

  }
}

//draw a square x=[-side,side] x z=[-side,side] at given y
void draw_xz_rect(GLfloat y, GLfloat side, GLfloat* color) {

  if (background_mode) {
    glColor3fv(color);
    glBegin(GL_POLYGON);
    glVertex3f(-side,y, side);
    glVertex3f(-side,y, -side);
    glVertex3f(side,y, -side);
    glVertex3f(side,y, side);
    glEnd();
  }
}

//draw a cube with middle planes 
void cube(GLfloat side) {
  GLfloat f = side, b = -side;

  if (fill_mode) glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  else  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  /* back face*/
  draw_xy_rect(b,side, lightcyan);
  /* side faces */
  draw_yz_rect(b, side, darkcyan);
  draw_yz_rect(f, side, darkcyan);
  //front, up, down faces missing to be able to see inside 

  /* middle z=0 face */
  draw_xy_rect(0, side, lightgrey);
  /* middle x=0 face */
  draw_yz_rect(0,side, darkgrey);
  /* middle y=0 face  */
  draw_xz_rect(0, side, darkmagenta);
}

//draw a filled cube [-side,side]^3
void filledcube(GLfloat side) {
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  
  /* back, front faces */
  draw_xy_rect(-side,side, yellow);
  draw_xy_rect(side,side, yellow);
  
  /* left, right faces*/
  draw_yz_rect(-side, side, yellow);
  draw_yz_rect(side, side, yellow);
  
  /* up, down  faces  */
  draw_xz_rect(side,side, yellow);
  draw_xz_rect(-side,side, yellow);
}
