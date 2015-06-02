
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vec3.hpp"

static double alpha_ = 0;
static double beta_ = 0;
static double angleXAxis_=20;
static double angleYAxis_=45;
static double angleZAxis_=0;
static double window_width_ = 1000;
static double window_height_ =1000;
double speed=1;


void DrawPyramid(const Vec3& P,double a, double b,double h){
  Vec3 normal,p1,p2,p3,p4,p5;

  glBegin(GL_TRIANGLES);            // Start Drawing A Triangle
  glNormal3f( 0.0,  0.0, 1.0);			// Set Top Point Of Triangle To Red
  glVertex3f( 0.0,  2.0, 0.0);      // First Point Of The Triangle
  glVertex3f(-2.0, -2.0, 0.0);      // Second Point Of The Triangle
  glVertex3f( 2.0, -2.0, 0.0);      // Third Point Of The Triangle
  glEnd();
}

// draw a sphere composed of triangles
void DrawSphere(const Vec3& ctr, double r){
  int     i, j,
          n1 = 20, n2 = 60;
  Vec3    normal, v1;
  double  a1, a1d = M_PI / n1,
          a2, a2d = M_PI / n2,
          s1, s2,
          c1, c2;

  glShadeModel(GL_SMOOTH);
  for(i = 0; i < n1; i++){
    a1 = i * a1d;

    glBegin(GL_TRIANGLE_STRIP);
    for(j = 0; j <= n2; j++){
      a2 = (j + .5 * (i % 2)) * 2 * a2d;

      s1 = sin(a1);
      c1 = cos(a1);
      s2 = sin(a2);
      c2 = cos(a2);
      normal = c1 * XVec3 + s1 * (c2 * YVec3 + s2 * ZVec3);
      v1 = ctr + r * normal;
      glNormal3dv(normal.p);
      glVertex3dv(v1.p);

      s1 = sin(a1 + a1d);
      c1 = cos(a1 + a1d);
      s2 = sin(a2 + a2d);
      c2 = cos(a2 + a2d);
      normal = c1 * XVec3 + s1 * (c2 * YVec3 + s2 * ZVec3);
      v1 = ctr + r * normal;
      glNormal3dv(normal.p);
      glVertex3dv(v1.p);
    }
    glEnd();
  }
}

void SetMaterialColor(int side, double r, double g, double b) {
  float	amb[4], dif[4], spe[4];
  int mat;

  dif[0] = r;
  dif[1] = g;
  dif[2] = b;

  for(int i = 0; i < 3; i++) {
    amb[i] = .1 * dif[i];
    spe[i] = .5;
  }
  amb[3] = dif[3] = spe[3] = 1.0;

  switch(side){
    case 1:	mat = GL_FRONT;
      break;
    case 2:	mat = GL_BACK;
      break;
    default: mat = GL_FRONT_AND_BACK;
  }

  glMaterialfv(mat, GL_AMBIENT, amb);
  glMaterialfv(mat, GL_DIFFUSE, dif);
  glMaterialfv(mat, GL_SPECULAR, spe);
  glMaterialf( mat, GL_SHININESS, 20);
}



// set viewport transformations and draw objects
void InitLighting() {
  GLfloat lp1[4]  = { 10,  20,  10,  0};
  GLfloat lp2[4]  = { -5,  5, -10,  0};
  GLfloat red[49]  = {1.0, .8,  .8,  1};
  GLfloat blue[4] = { .8, .8, 1.0,  1};

  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);
  glClearColor(1, 1, 1, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glShadeModel(GL_SMOOTH);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
  glEnable(GL_LIGHTING);

  glLightfv(GL_LIGHT1, GL_POSITION, lp1);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  red);
  glLightfv(GL_LIGHT1, GL_SPECULAR, red);
  glEnable(GL_LIGHT1);

  glLightfv(GL_LIGHT2, GL_POSITION, lp2);
  glLightfv(GL_LIGHT2, GL_DIFFUSE,  blue);
  glLightfv(GL_LIGHT2, GL_SPECULAR, blue);
  glEnable(GL_LIGHT2);

  glClearColor(1, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // init viewport to canvassize
  glViewport(0, 0, window_width_, window_height_);

  // init coordinate system
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-15, 15, -10, 10, -20, 20);

//  double n=20;
//  double f=10;
//  double a=1;
//  double fov=180;
//  double l=-1*n*tan((a*fov)/2);
//  double r=f*tan((a*fov)/2);
//  double b=-1*n*tan(fov/2);
//  double t=f*tan(fov/2);
//
//  glFrustum(l, r, b, t, n, f);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

void DrawQuad(Vec3 a,Vec3 b,Vec3 c,Vec3 d){
	  glBegin(GL_QUADS);
	  glNormal3f( 0.0,  0.0, 1.0);
	  glVertex3dv(a.p);
	  glVertex3dv(b.p);
	  glVertex3dv(c.p);
	  glVertex3dv(d.p);
	  glEnd();
}
void DrawTriangle(Vec3 a,Vec3 b,Vec3 c){
	  glBegin(GL_TRIANGLES);
	  glNormal3f( 0.0,  0.0, 1.0);
	  glVertex3dv(a.p);
	  glVertex3dv(b.p);
	  glVertex3dv(c.p);
	  glEnd();
}
void DrawPyramidOrigami(float seitenlaenge, float angle){
	float hoehe=sqrt(seitenlaenge*seitenlaenge-(seitenlaenge/2)*(seitenlaenge/2));
	Vec3 P1(0,0,0);
	Vec3 P2(seitenlaenge,0,0);
	Vec3 P3(seitenlaenge/2,hoehe,0);

	DrawTriangle(P1,P2,P3);

	glPushMatrix();
	glTranslatef(0,0,0);
	glRotatef(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(seitenlaenge,0,0);
	glRotated(180,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();

}
void DrawIkosaederOrigami(float seitenlaenge, float angle){
	float hoehe=sqrt(seitenlaenge*seitenlaenge-(seitenlaenge/2)*(seitenlaenge/2));
	Vec3 P1(0,0,0);
	Vec3 P2(seitenlaenge,0,0);
	Vec3 P3(seitenlaenge/2,hoehe,0);

	//glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(0,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	//glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();

	//glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	//glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();

	//glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	//glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();

	//glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	//glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();

	//glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	//glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//glPopMatrix();
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
}

void DrawTTPOrigami(float seitenlaenge, float angle){//fl�chen richtig aber es muss verschiedene winkel gefaltet werden

	float hoehe=sqrt(seitenlaenge*seitenlaenge-(seitenlaenge/2)*(seitenlaenge/2));
	Vec3 P1(0,0,0);
	Vec3 P2(seitenlaenge,0,0);
	Vec3 P3(seitenlaenge/2,hoehe,0);
	//Root
	//draw tri 1
	DrawTriangle(P1,P2,P3);
	//draw group 2 3 4
	// draw 2
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	// draw 3
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	// draw 4
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	// draw 5
	glPushMatrix();
	glTranslatef(seitenlaenge,0,0);
	glRotated(180,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	//draw group 6 7 8 9 10 11 12 13 14
	//draw 6

	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);

	//draw group 7 8 9 10
	//draw 7
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);

	//draw 10
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	//draw group 8 9
	//draw 8

	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	//draw 9

	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();
	glPopMatrix();
	SetMaterialColor(3,1,1,1);
	//draw group 11 12 13 14
	//draw 11
	glPushMatrix();
	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);

	//draw 12

	glTranslatef(seitenlaenge/2,hoehe,0);
	glRotated(-60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);

	//draw 13

	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);

	//draw 14

	glTranslatef(0,0,0);
	glRotated(60,0,0,1);
	glRotated(angle,1,0,0);
	DrawTriangle(P1,P2,P3);
	glPopMatrix();



}
void DrawCubeOrigami(float seitenlaenge, float angle){
	Vec3 P1(0,0,0);
	Vec3 P2(0,seitenlaenge,0);
	Vec3 P3(seitenlaenge,seitenlaenge,0);
	Vec3 P4(seitenlaenge,0,0);

	DrawQuad(P1,P2,P3,P4);
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(90,0,0,1);
	glRotated(angle,1,0,0);
	DrawQuad(P1,P2,P3,P4);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0,seitenlaenge,0);
	glRotated(0,0,0,1);
	glRotated(angle,1,0,0);
	DrawQuad(P1,P2,P3,P4);

	glTranslatef(0,seitenlaenge,0);
	glRotated(0,0,0,1);
	glRotated(angle,1,0,0);
	DrawQuad(P1,P2,P3,P4);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(seitenlaenge,seitenlaenge,0);
	glRotated(-90,0,0,1);
	glRotated(angle,1,0,0);
	DrawQuad(P1,P2,P3,P4);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(seitenlaenge,0,0);
	glRotated(180,0,0,1);
	glRotated(angle,1,0,0);
	DrawQuad(P1,P2,P3,P4);
	glPopMatrix();

}
void DrawCube(double seitenlaenge){
	//Eckpunkte
	Vec3 P1(-seitenlaenge/2,seitenlaenge/2,-seitenlaenge/2);
	Vec3 P2(seitenlaenge/2,seitenlaenge/2,-seitenlaenge/2);
	Vec3 P3(seitenlaenge/2,seitenlaenge/2,seitenlaenge/2);
	Vec3 P4(-seitenlaenge/2,seitenlaenge/2,seitenlaenge/2);
	Vec3 P5(-seitenlaenge/2,-seitenlaenge/2,-seitenlaenge/2);
	Vec3 P6(seitenlaenge/2,-seitenlaenge/2,-seitenlaenge/2);
	Vec3 P7(seitenlaenge/2,-seitenlaenge/2,seitenlaenge/2);
	Vec3 P8(-seitenlaenge/2,-seitenlaenge/2,seitenlaenge/2);
	//Seiten
	DrawQuad(P1,P2,P3,P4);
	DrawQuad(P3,P2,P6,P7);
	DrawQuad(P4,P3,P7,P8);
	DrawQuad(P1,P4,P8,P5);
	DrawQuad(P2,P1,P5,P6);
	DrawQuad(P5,P6,P7,P8);
}

void DrawBox(double seitenlaenge){
	//Eckpunkte
	Vec3 P1(-seitenlaenge/2,seitenlaenge/2,-seitenlaenge/2);
	Vec3 P2(seitenlaenge/2,seitenlaenge/2,-seitenlaenge/2);
	Vec3 P3(seitenlaenge/2,seitenlaenge/2,seitenlaenge/2);
	Vec3 P4(-seitenlaenge/2,seitenlaenge/2,seitenlaenge/2);
	Vec3 P5(-seitenlaenge/2,-seitenlaenge/2,-seitenlaenge/2);
	Vec3 P6(seitenlaenge/2,-seitenlaenge/2,-seitenlaenge/2);
	Vec3 P7(seitenlaenge/2,-seitenlaenge/2,seitenlaenge/2);
	Vec3 P8(-seitenlaenge/2,-seitenlaenge/2,seitenlaenge/2);
	//Die 6 seiten zeichnen
	glPushMatrix();
	glTranslated(-seitenlaenge/2,seitenlaenge/2,-seitenlaenge/2);
	glRotated(beta_, 0, 0, 1);
	DrawQuad(Vec3 (0,0,0),Vec3 (seitenlaenge/2,0,0),Vec3 (seitenlaenge/2,0,seitenlaenge),Vec3 (0,0,seitenlaenge));
	glPushMatrix();
	glTranslated(seitenlaenge/2,0,0);
	glRotated(-2*beta_, 0, 0, 1);
	DrawQuad(Vec3 (0,0,0),Vec3 (seitenlaenge/2,0,0),Vec3 (seitenlaenge/2,0,seitenlaenge),Vec3 (0,0,seitenlaenge));
	glPopMatrix();
	glPopMatrix();
	//SetMaterialColor(2, 0, 0, 0);
	DrawQuad(P3,P2,P6,P7);
	//SetMaterialColor(2, 0, 0, 1);
	DrawQuad(P4,P3,P7,P8);
	//SetMaterialColor(2, 1, 0, 1);
	DrawQuad(P1,P4,P8,P5);
	//SetMaterialColor(2, 1, 1, 1);
	DrawQuad(P2,P1,P5,P6);
	//SetMaterialColor(2, 0, 1, 1);
	DrawQuad(P5,P8,P7,P6);
}
// draw the entire scene
void Preview() {
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();						    // Reset The Current Modelview Matrix
  glTranslated(0, 0, -10.0);      // Move 10 units backwards in z,
                                  // since camera is at origin

  //glPushMatrix();
  glRotated(angleXAxis_,0,1,0);
  glRotated(angleYAxis_,1,0,0);
  SetMaterialColor(1, 1, 0, 0);
  SetMaterialColor(2, 0, 0, 1);
  DrawIkosaederOrigami(10,beta_);
  //glPopMatrix();

}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	  if (key == GLFW_KEY_UP)angleYAxis_+=5;
	  if (key == GLFW_KEY_DOWN)angleYAxis_-=5;
	  if (key == GLFW_KEY_LEFT)angleXAxis_+=5;
	  if (key == GLFW_KEY_RIGHT)angleXAxis_-=5;
	  if (key == GLFW_KEY_O){
		  if(beta_<90)beta_ +=1;
	  }
	  if (key == GLFW_KEY_C){
		  if(beta_>0)beta_ -=5;
	  }
	  if (key == GLFW_KEY_R){
		  beta_ =0;
		  angleXAxis_=0;
		  angleYAxis_=0;
	  }
	  if(key == GLFW_KEY_ESCAPE)exit(0);
}


int main() {
  GLFWwindow* window = NULL;

  printf("Here we go!\n");

  if(!glfwInit()){
    return -1;
  }

  window = glfwCreateWindow(window_width_, window_height_,
                            "Simple 3D Animation", NULL, NULL);
  if(!window) {
    glfwTerminate();
    return -1;
  }

  glfwMakeContextCurrent(window);

  while(!glfwWindowShouldClose(window)) {
    // switch on lighting (or you don't see anything)
    InitLighting();

    // set background color
    glClearColor(0.8, 0.8, 0.8, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //KeyFunctions
    glfwSetKeyCallback(window, key_callback);
    // draw the scene
    Preview();

    // make it appear (before this, it's hidden in the rear buffer)
    glfwSwapBuffers(window);

    glfwPollEvents();
  }

  glfwTerminate();

  printf("Goodbye!\n");

  return 0;
}
