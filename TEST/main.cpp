#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <vector>
//#include <glm\glm.hpp>

using namespace std;

#include "vec3.hpp"
static double alpha_ = 0;
static double beta_ = 0;
static double gamma_= 0;
static double delta_ = 0;
static double epsilon_ = 0;
static double zetta_ = 0;
static double bla_ = 0;
static double dt_ = 1;
static double angleXAxis_=20;
static double angleYAxis_=45;
static double angleZAxis_=0;
static double posX_=0;
static double posY_=0;
static double posZ_=0;
static double posX2_=0;
static double posY2_=0;
static double posZ2_=0;

static Vec3 e1(1,0,0);
static Vec3 e2(0,1,0);
static Vec3 e3(0,0,1);
static Vec3 KAX_(1,0,0);
static Vec3 KAY_(0,1,0);
static Vec3 KAZ_(0,0,1);


vector<Vec3>Kugel;
vector<Vec3>KugelSpeeds;
vector<double>KugelRadiuses;
vector<double>KugelMasses;


static Vec3 kugel1(10, 0, 0);
static Vec3 kugel2(-10, 0, 3);
static Vec3 speedKugel1(-0.05, 0, 0);
static Vec3 speedKugel2(0.05, 0, 0);
static double radiusKugel1 = 3;
static double radiusKugel2 = 2;
static double massKugel1 = 0.75 * M_PI*pow(radiusKugel1, 3);
static double massKugel2 = 0.75 * M_PI*pow(radiusKugel2, 3);




static Vec3 mouseSpeed(0,0,0);

static double window_width_ = 1000;
static double window_height_ =1000;
double speed=1;

double getDistanceToLine(Vec3 normale, Vec3 punkt, double verschiebung){
	return abs((normale/normale.Length())*punkt-verschiebung);
}

double getAngle2D(double Ax, double  Ay, double Bx, double By){
	return atan2(Ay,Ax)-atan2(By,Bx)/M_PI*180;
}
Vec3 rotationUmVector(Vec3 N,double alpha){
	Vec3 result(0,0,0);



	return result;

}
Vec3 spiegelungUrsprungsgerade(double alpha, Vec3 V){
	Vec3 result(0,0,0);
	result.p[0]=V.p[0]*cos(2*alpha)+V.p[2]*sin(2*alpha);
	result.p[2]=V.p[0]*sin(2*alpha)+V.p[2]*(-cos(2*alpha));
	return result;
}
Vec3 spiegelungUmNormale(Vec3 V,Vec3 A){
	Vec3 result = V-2*((A*V)/(A*A))*A;
	return result;
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
//    if(i<=n1*0.5)SetMaterialColor(3,0,0,0);
//    if(i>n1*0.5)SetMaterialColor(3,1,1,1);

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

//  double n=10;
//  double f=40;
//  double a=0.5;
//  double fov=90;
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
	Vec3 norm((a-b)%(a-c));
	norm=norm/norm.Length();
	glBegin(GL_QUADS);
	glNormal3dv(norm.p);
	glVertex3dv(a.p);
	glVertex3dv(b.p);
	glVertex3dv(c.p);
	glVertex3dv(d.p);
	glEnd();
}
void DrawTriangle(Vec3 a,Vec3 b,Vec3 c){
	Vec3 norm((a-b)%(a-c));
	norm=norm/norm.Length();
	glBegin(GL_TRIANGLES);
	glNormal3dv(norm.p);
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
	DrawQuad(P3,P2,P6,P7);
	DrawQuad(P4,P3,P7,P8);
	DrawQuad(P1,P4,P8,P5);
	DrawQuad(P2,P1,P5,P6);
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
                DrawTriangle(P0,P1,P2);
                DrawTriangle(P1,P2,P3);
                DrawTriangle(P0,M2,P2);
                DrawTriangle(P3,P1,M1);
        }
}
void DrawKrane(){
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
bool cFCC(double seite,double punkt, double radius){
	return abs(punkt-seite/2)-radius<0;
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

void wandCollision(double x, double z){
	for(int i=0; i < Kugel.size(); i++){
		if(abs(Kugel[i].p[0]-x)-KugelRadiuses[i]<0){
			KugelSpeeds[i].p[0]=-KugelSpeeds[i].p[0];
		}
		if(abs(-Kugel[i].p[0]-x)-KugelRadiuses[i]<0){
			KugelSpeeds[i].p[0]=-KugelSpeeds[i].p[0];
		}
		if(abs(Kugel[i].p[2]-z)-KugelRadiuses[i]<0){
			KugelSpeeds[i].p[2]=-KugelSpeeds[i].p[2];
		}
		if(abs(-Kugel[i].p[2]-z)-KugelRadiuses[i]<0){
			KugelSpeeds[i].p[2]=-KugelSpeeds[i].p[2];
		}
	}
}
void wandCollision2(double x, double z){
    if(abs(kugel1.p[0]-x)-radiusKugel1<0){
    	speedKugel1.p[0]=-speedKugel1.p[0];
    }
    if(abs(-kugel1.p[0]-x)-radiusKugel1<0){
    	speedKugel1.p[0]=-speedKugel1.p[0];
    }
    if(abs(kugel1.p[2]-z)-radiusKugel1<0){
    	speedKugel1.p[2]=-speedKugel1.p[2];
    }
    if(abs(-kugel1.p[2]-z)-radiusKugel1<0){
    	speedKugel1.p[2]=-speedKugel1.p[2];
    }

    if(abs(kugel2.p[0]-x)-radiusKugel2<0){
    	speedKugel2.p[0]=-speedKugel2.p[0];
    }
    if(abs(-kugel2.p[0]-x)-radiusKugel1<0){
    	speedKugel2.p[0]=-speedKugel2.p[0];
    }
    if(abs(kugel2.p[2]-z)-radiusKugel2<0){
    	speedKugel2.p[2]=-speedKugel2.p[2];
    }
    if(abs(-kugel2.p[2]-z)-radiusKugel2<0){
    	speedKugel2.p[2]=-speedKugel2.p[2];
    }
}

void kugelnBewegen(){
	for(int i=0; i < Kugel.size(); i++){
		Kugel[i]+=KugelSpeeds[i];
	}
}
void kugelnZeichnen(){
	for(int i=0; i < Kugel.size(); i++){
		DrawSphere(Kugel[i],KugelRadiuses[i]);
	}
}
void kugelCollision(){
	for(int i=0; i < Kugel.size(); i++){
		for(int j=0; j < Kugel.size(); j++){
			if(i != j){
				if (abs((Kugel[i] - Kugel[j]).Length()) - (KugelRadiuses[i] + KugelRadiuses[j]) <= 0){
					Vec3 neuSpeedKugel1 = KugelSpeeds[i] - ((2 * KugelMasses[j]) / (KugelMasses[i] + KugelMasses[j])) * (((KugelSpeeds[i] - KugelSpeeds[j]).Dot(Kugel[i] - Kugel[j])) / ((Kugel[i] - Kugel[j]).Length2())) * (Kugel[i] - Kugel[j]);
					Vec3 neuSpeedKugel2 = KugelSpeeds[j] - ((2 * KugelMasses[i]) / (KugelMasses[i] + KugelMasses[j])) * (((KugelSpeeds[j] - KugelSpeeds[i]).Dot(Kugel[j] - Kugel[i])) / ((Kugel[j] - Kugel[i]).Length2())) * (Kugel[j] - Kugel[i]);
					KugelSpeeds[i] = neuSpeedKugel1;
					KugelSpeeds[j] = neuSpeedKugel2;
				}
			}

		}
	}
}


double kugelrotation;
void DrawBillardtisch(){
	double x=20+beta_;
	double z=15+alpha_;
	double h=2+delta_;
	double rb=1+gamma_;
	double kugelradius=1+epsilon_;
	double scaling=0.0001;
	double friction=1-0.0000;

	Vec3 hoehe(0,h,0);
	Vec3 P1(-x,0,z);
	Vec3 P2(x,0,z);
	Vec3 P3(x,0,-z);
	Vec3 P4(-x,0,-z);

	//Boden
	SetMaterialColor(3, 0, 0.5, 0);
	DrawQuad(P1,P2,P3,P4);
	//seiten
	SetMaterialColor(3, 0, 0.6, 0);
	DrawQuad(P1+e2*h,P2+e2*h,P2,P1);
	DrawQuad(P2+e2*h,P3+e2*h,P3,P2);
	DrawQuad(P3+e2*h,P4+e2*h,P4,P3);
	DrawQuad(P4+e2*h,P1+e2*h,P1,P4);
	//rahmen
	SetMaterialColor(3, 0.6, 0.6, 0);
	DrawQuad(P1+e2*h-e1*rb+e3*rb,  P2+e2*h+e1*rb+e3*rb,  P2+e2*h,  P1+e2*h);
	DrawQuad(P2+e2*h+e1*rb+e3*rb,  P3+e2*h+e1*rb-e3*rb,  P3+e2*h,  P2+e2*h);
	DrawQuad(P3+e2*h+e1*rb-e3*rb,  P4+e2*h-e1*rb-e3*rb,  P4+e2*h,  P3+e2*h);
	DrawQuad(P4+e2*h-e1*rb-e3*rb,  P1+e2*h-e1*rb+e3*rb,  P1+e2*h,  P4+e2*h);

	DrawQuad(P2+e2*h+e1*rb+e3*rb,  P1+e2*h-e1*rb+e3*rb,  P1-e1*rb+e3*rb,  P2+e1*rb+e3*rb);
	DrawQuad(P3+e2*h+e1*rb-e3*rb,  P2+e2*h+e1*rb+e3*rb,  P2+e1*rb+e3*rb,  P3+e1*rb-e3*rb);
	DrawQuad(P4+e2*h-e1*rb-e3*rb,  P3+e2*h+e1*rb-e3*rb,  P3+e1*rb-e3*rb,  P4-e1*rb-e3*rb);
	DrawQuad(P1+e2*h-e1*rb+e3*rb,  P4+e2*h-e1*rb-e3*rb,  P4-e1*rb-e3*rb,  P1-e1*rb+e3*rb);

	DrawQuad(P1-e1*rb+e3*rb,  P2+e1*rb+e3*rb,  P2,  P1);
	DrawQuad(P2+e1*rb+e3*rb,  P3+e1*rb-e3*rb,  P3,  P2);
	DrawQuad(P3+e1*rb-e3*rb,  P4-e1*rb-e3*rb,  P4,  P3);
	DrawQuad(P4-e1*rb-e3*rb,  P1-e1*rb+e3*rb,  P1,  P4);


	wandCollision(x,z);
	wandCollision2(x,z);
	kugelCollision();
	kugelnBewegen();
	kugelnZeichnen();



	if ((kugel1-kugel2).Length() - (radiusKugel1 + radiusKugel2) <= 0){
		Vec3 neuSpeedKugel1 = speedKugel1 - ((2 * massKugel2) / (massKugel1 + massKugel2)) * (((speedKugel1 - speedKugel2).Dot(kugel1 - kugel2)) / ((kugel1 - kugel2).Length2())) * (kugel1 - kugel2);
		Vec3 neuSpeedKugel2 = speedKugel2 - ((2 * massKugel1) / (massKugel1 + massKugel2)) * (((speedKugel2 - speedKugel1).Dot(kugel2 - kugel1)) / ((kugel2 - kugel1).Length2())) * (kugel2 - kugel1);
		speedKugel1 = neuSpeedKugel1;
		speedKugel2 = neuSpeedKugel2;
	}


	kugel1+=speedKugel1;
	kugel2+=speedKugel2;

	DrawSphere(kugel1,radiusKugel1);
	DrawSphere(kugel2,radiusKugel2);



}
void DrawCubeKugel(double seitenLaenge,double kugelradius){
	  DrawCube(seitenLaenge);
	  if(cFC3D(seitenLaenge, kugelradius, posX_, posY_,posZ_))SetMaterialColor(3, 1, 0, 0.5);
	  DrawSphere(Vec3(posX_,posY_,posZ_),kugelradius);
}

// draw the entire scene
void Preview() {
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();						    // Reset The Current Modelview Matrix
  glTranslated(0, 0, -30.0);      // Move 10 units backwards in z,
                                  // since camera is at origin

  //glPushMatrix();
  glRotated(angleXAxis_,0,1,0);
  glRotated(angleYAxis_,1,0,0);
  SetMaterialColor(3, 1, 0.5, 0);
  //DrawKrane();
  DrawBillardtisch();
  //DrawZylinder(10,10,100);
  //DrawCubeKugel(5,2);


  glPopMatrix();

}
void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS){

 	   double xpos, ypos;
 	   glfwGetCursorPos(window, &xpos, &ypos);
 	   mouseSpeed.p[0] = xpos;
 	   mouseSpeed.p[2] = ypos;

    }

    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE){

 	   double xpos, ypos;
 	   glfwGetCursorPos(window, &xpos, &ypos);
 	   mouseSpeed.p[0] = xpos - mouseSpeed.p[0];
 	   mouseSpeed.p[2] = ypos - mouseSpeed.p[2];

// 	   cout << mouseSpeed.p[0] << endl << mouseSpeed.p[2] << endl;
//
// 	  SpherePos.p[0] = mouseSpeed.p[0];
// 	  SpherePos.p[2] = mouseSpeed.p[2];
    }


}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	  if (key == GLFW_KEY_N)dt_+=0.1;
	  if (key == GLFW_KEY_M)dt_-=0.1;
	  if (key == GLFW_KEY_UP)angleYAxis_+=dt_;
	  if (key == GLFW_KEY_DOWN)angleYAxis_-=dt_;
	  if (key == GLFW_KEY_LEFT)angleXAxis_+=dt_;
	  if (key == GLFW_KEY_RIGHT)angleXAxis_-=dt_;
	  if (key == GLFW_KEY_W)posZ2_+=dt_;
	  if (key == GLFW_KEY_S)posZ2_-=dt_;
	  if (key == GLFW_KEY_A)posX2_-=dt_;
	  if (key == GLFW_KEY_D)posX2_+=dt_;
	  if (key == GLFW_KEY_Q)posY2_+=dt_;
	  if (key == GLFW_KEY_E)posY2_-=dt_;
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
      if (key == GLFW_KEY_B)zetta_+=0.10;
      if (key == GLFW_KEY_V)zetta_-=0.10;
      if (key == GLFW_KEY_F)bla_+=0.1;
      if (key == GLFW_KEY_G)bla_-=0.1;

	  if (key == GLFW_KEY_O){
//		  if(beta_<90)beta_ +=1;
		  beta_ +=1;
	  }
	  if (key == GLFW_KEY_C){
//		  if(beta_>0)beta_ -=5;
		  beta_ -=1;
	  }
	  if (key == GLFW_KEY_R){
		  kugel1 = Vec3(10, 2, 0);
		  kugel2 = Vec3(0, 0, 10);
		  speedKugel1 = Vec3(-0.05, 0, 0);
		  speedKugel2 = Vec3(0, 0, -0.05);
		  radiusKugel1 = 2;
		  radiusKugel2 = 1;
		  massKugel1 = 0.75 * M_PI*pow(radiusKugel1, 3);
		  massKugel2 = 0.75 * M_PI*pow(radiusKugel2, 3);

	  }
	  if(key == GLFW_KEY_ESCAPE)exit(0);
}

void deleteBall(int pointer){
	Kugel.erase(Kugel.begin()+pointer);
	KugelSpeeds.erase(KugelSpeeds.begin()+pointer);;
	KugelRadiuses.erase(KugelRadiuses.begin()+pointer);;
	KugelMasses.erase(KugelMasses.begin()+pointer);
}

void createBall(Vec3 pos, Vec3 speed, double radius){
	Kugel.push_back(pos);
	KugelSpeeds.push_back(speed);
	KugelRadiuses.push_back(radius);
	KugelMasses.push_back(0.75 * M_PI*pow(radius, 3));
}

void initGame(){
	createBall(Vec3(2, 0, 0),Vec3(0, 0, 0),4);
	createBall(Vec3(0, 0, 10),Vec3(0, 0, -0.1),3);
	createBall(Vec3(-10, 0, 10),Vec3(0.1, 0, -0.1),2);
	createBall(Vec3(10, 0, -10),Vec3(0.1, 0, -0.1),1);
}


int main() {
  GLFWwindow* window = NULL;

  printf("Here we go!\n");

  if(!glfwInit()){
    return -1;
  }
  //initGame();
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
    glfwSetMouseButtonCallback(window, mouse_button_callback);
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
