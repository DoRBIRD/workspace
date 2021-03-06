
#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "vec3.hpp"
static double alpha_ = 0;
static double beta_ = 0;
static double gamma_= 0;
static double delta_ = 0;
static double epsilon_ = 0;
static double zetta_ = 0;
static double dt_ = 1;
static double angleXAxis_=20;
static double angleYAxis_=45;
static double angleZAxis_=0;
static double posX_=0;
static double posY_=0;
static double posZ_=0;

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
  GLfloat lp2[4]  = { -50,  50, -10,  0};
  GLfloat lp3[4]  = { -5,  -50, -10,  0};
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

  glLightfv(GL_LIGHT2, GL_POSITION, lp3);
  glLightfv(GL_LIGHT2, GL_DIFFUSE,  blue);
  glLightfv(GL_LIGHT2, GL_SPECULAR, blue);
  glEnable(GL_LIGHT3);

  glClearColor(1, 0, 0, 1);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  // init viewport to canvassize
  glViewport(0, 0, window_width_, window_height_);

  // init coordinate system
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(-30, 30, -30, 30, -50, 50);

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
	Vec3 norm(b%a);
	//norm.Normalize();
	norm=norm/(sqrt(pow(norm.p[0],2)+pow(norm.p[1],2)+pow(norm.p[2],2)));
	glBegin(GL_QUADS);
	glNormal3dv(norm.p);
	//glNormal3f( 0.0,  0.0, 1.0);
	glVertex3dv(a.p);
	glVertex3dv(b.p);
	glVertex3dv(c.p);
	glVertex3dv(d.p);
	glEnd();
}
void DrawTriangle(Vec3 a,Vec3 b,Vec3 c){
	Vec3 norm(b%a);
	norm.Normalize();
	glBegin(GL_TRIANGLES);
	//glNormal3dv(norm.p);
	glNormal3f( 0.0,  0.0, 1.0);
	glVertex3dv(a.p);
	glVertex3dv(b.p);
	glVertex3dv(c.p);
	glEnd();
}
void DrawFancyQuad(Vec3 a,Vec3 b,Vec3 c,Vec3 d,double pro,double hoehe){
	Vec3 h(0,0,hoehe);
	Vec3 a2=((a-b)*0.1-(a-d)*0.1)+h;
	Vec3 b2=((b-c)*0.1-(b-a)*0.1)+h;
	Vec3 c2=((c-d)*0.1-(c-b)*0.1)+h;
	Vec3 d2=((d-a)*0.1-(d-c)*0.1)+h;

	DrawQuad(a,d,c,b);
	DrawQuad(a2,b2,c2,d2);
	DrawQuad(c2,b2,b,c);
	DrawQuad(d2,c2,c,d);
	DrawQuad(a2,d2,d,a);
	DrawQuad(b2,a2,a,b);

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
	glTranslatef(0,-3*seitenlaenge,0);
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
	DrawQuad(P1,P4,P3,P4);
	glPopMatrix();

}
void DrawFancyCubeOrigami(float seitenlaenge, float angle){
	Vec3 P1(0,0,0);
	Vec3 P2(0,seitenlaenge,0);
	Vec3 P3(seitenlaenge,seitenlaenge,0);
	Vec3 P4(seitenlaenge,0,0);

	DrawFancyQuad(P1,P2,P3,P4,0.1,1);
	glPushMatrix();
	glTranslatef(0,0,0);
	glRotated(90,0,0,1);
	glRotated(angle,1,0,0);
	DrawFancyQuad(P1,P2,P3,P4,0.1,1);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(0,seitenlaenge,0);
	glRotated(0,0,0,1);
	glRotated(angle,1,0,0);
	DrawFancyQuad(P1,P2,P3,P4,0.1,1);

	glTranslatef(0,seitenlaenge,0);
	glRotated(0,0,0,1);
	glRotated(angle,1,0,0);
	DrawFancyQuad(P1,P2,P3,P4,0.1,1);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(seitenlaenge,seitenlaenge,0);
	glRotated(-90,0,0,1);
	glRotated(angle,1,0,0);
	DrawFancyQuad(P1,P2,P3,P4,0.1,1);
	glPopMatrix();

	glPushMatrix();
	glTranslatef(seitenlaenge,0,0);
	glRotated(180,0,0,1);
	glRotated(angle,1,0,0);
	DrawFancyQuad(P1,P2,P3,P4,0.1,1);
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
	DrawQuad(Vec3 (0,0,0),Vec3 (0,0,seitenlaenge),Vec3 (seitenlaenge/2,0,seitenlaenge),Vec3 (seitenlaenge/2,0,0));
	glPushMatrix();
	glTranslated(seitenlaenge/2,0,0);
	glRotated(-2*beta_, 0, 0, 1);
	DrawQuad(Vec3 (0,0,0),Vec3 (0,0,seitenlaenge),Vec3 (seitenlaenge/2,0,seitenlaenge),Vec3 (seitenlaenge/2,0,0));
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
void DrawZylinder(double radius,double hoehe,double res){
        int i;
        double a1=2* M_PI /res;
        Vec3 M1(0,0,0);
        Vec3 M2(0,0,hoehe);
        for(i=0; i<=res; i++){
                Vec3 P0(radius*sin(a1*(i+0)),radius*cos(a1*(i+0)),hoehe);
                Vec3 P1(radius*sin(a1*(i+0.5)),radius*cos(a1*(i+0.5)),0);
                Vec3 P2(radius*sin(a1*(i+1)),radius*cos(a1*(i+1)),hoehe);
                Vec3 P3(radius*sin(a1*(i+1.5)),radius*cos(a1*(i+1.5)),0);
                DrawTriangle(P0,P2,P1);
                DrawTriangle(P1,P2,P3);
                DrawTriangle(P2,P0,M2);
                DrawTriangle(P1,P3,M1);
        }
}
void DrawKrane(Vec3 traget){
        double baseSize=2;
        double radius=1;
        double firstArmLenght=5;
        double secondArmLenght=5;
        double thirdArmLength=5;
        double greifer=3;
        glPushMatrix();
        glRotated(alpha_,0,0,1);
        DrawCube(baseSize);
        glTranslated(0,0,baseSize/2);
        glRotated(beta_,0,1,0);
        DrawSphere(Vec3(0,0,0),radius*1.2);
        DrawZylinder(radius,firstArmLenght,100);
        glTranslated(0,0,firstArmLenght);
        glRotated(gamma_,0,1,0);
        DrawSphere(Vec3(0,0,0),radius*1.2);
        DrawZylinder(radius,secondArmLenght,100);
        glTranslated(0,0,secondArmLenght);
        glRotated(delta_,0,1,0);
        DrawSphere(Vec3(0,0,0),radius*1.2);
        DrawZylinder(radius,thirdArmLength,100);
        glTranslated(-greifer/2,-0.3*sqrt(greifer*greifer-(greifer/2)*(greifer/2)),secondArmLenght);
        DrawPyramidOrigami(greifer,110);
        glPopMatrix();
}

bool cFC2D(double seitenLaenge, double radius, double px, double py){
	//test ob im x schlauch
	if(py<seitenLaenge/2&&py>-seitenLaenge/2){
		if(abs(px-seitenLaenge/2)-radius<0)return true;
		if(abs(-px-seitenLaenge/2)-radius<0)return true;
	}
	//test ob im y schlauch
	if(px<seitenLaenge/2&&px>-seitenLaenge/2){
		if(abs(py-seitenLaenge/2)-radius<0)return true;
		if(abs(-py-seitenLaenge/2)-radius<0)return true;
	}
	//testen nach Ecken
	if(px>seitenLaenge/2&&py>seitenLaenge/2&&sqrt(pow(px-seitenLaenge/2,2)+pow(py-seitenLaenge/2,2))-radius<0)return true;
	if(px>seitenLaenge/2&&py<-seitenLaenge/2&&sqrt(pow(px-seitenLaenge/2,2)+pow(-py-seitenLaenge/2,2))-radius<0)return true;
	if(px<-seitenLaenge/2&&py<-seitenLaenge/2&&sqrt(pow(-px-seitenLaenge/2,2)+pow(-py-seitenLaenge/2,2))-radius<0)return true;
	if(px<-seitenLaenge/2&&py>seitenLaenge/2&&sqrt(pow(-px-seitenLaenge/2,2)+pow(py-seitenLaenge/2,2))-radius<0)return true;
	return false;
}
bool cFC3D(double seitenLaenge, double radius, double posX, double posY,double posZ){
	int ccounter=0;
	if(cFC2D(seitenLaenge,radius,posX, posY))ccounter++;
	if(cFC2D(seitenLaenge,radius,posZ, posY))ccounter++;
	if(cFC2D(seitenLaenge,radius,posZ, posX))ccounter++;
	if(ccounter>=2)return true;
	return false;
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
  SetMaterialColor(3, 1, 0.5, 0);
  double seitenLaenge=10;
  double radius=3;
  DrawCube(seitenLaenge);
  if(cFC3D(seitenLaenge,radius,posX_,posY_,posZ_))SetMaterialColor(3, 1, 0, 0.5);
  DrawSphere(Vec3(posX_,posY_,posZ_),radius);
  //DrawIkosaederOrigami(5,beta_);
  //DrawKrane(Vec3(0,0,0));

  glPopMatrix();

}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	  if (key == GLFW_KEY_N)dt_+=0.1;
	  if (key == GLFW_KEY_M)dt_-=0.1;
	  if (key == GLFW_KEY_UP)angleYAxis_+=dt_;
	  if (key == GLFW_KEY_DOWN)angleYAxis_-=dt_;
	  if (key == GLFW_KEY_LEFT)angleXAxis_+=dt_;
	  if (key == GLFW_KEY_RIGHT)angleXAxis_-=dt_;
	  if (key == GLFW_KEY_W)posY_+=dt_;
	  if (key == GLFW_KEY_S)posY_-=dt_;
	  if (key == GLFW_KEY_A)posX_-=dt_;
	  if (key == GLFW_KEY_D)posX_+=dt_;
	  if (key == GLFW_KEY_Q)posZ_+=dt_;
	  if (key == GLFW_KEY_E)posZ_-=dt_;
      if (key == GLFW_KEY_1)alpha_+=dt_;
      if (key == GLFW_KEY_2)alpha_-=dt_;
      if (key == GLFW_KEY_3)beta_+=dt_;
      if (key == GLFW_KEY_4)beta_-=dt_;
      if (key == GLFW_KEY_5)gamma_+=dt_;
      if (key == GLFW_KEY_6)gamma_-=dt_;
      if (key == GLFW_KEY_7)delta_+=dt_;
      if (key == GLFW_KEY_8)delta_-=dt_;
      if (key == GLFW_KEY_9)epsilon_+=dt_;
      if (key == GLFW_KEY_0)epsilon_-=dt_;

	  if (key == GLFW_KEY_O){
//		  if(beta_<90)beta_ +=1;
		  beta_ +=1;
	  }
	  if (key == GLFW_KEY_C){
//		  if(beta_>0)beta_ -=5;
		  beta_ -=1;
	  }
	  if (key == GLFW_KEY_R){
		  beta_ =0;
		  angleXAxis_=0;
		  angleYAxis_=0;
		  posX_=0;
		  posY_=0;
		  posZ_=0;
		  dt_=1;
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
