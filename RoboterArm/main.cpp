/*
	Ben�tigt matrix.h erh�ltlich hier http://matrix.drque.net/
	Tastenbelegung:
		W,S,A,D,Q,E - Steuerung der Kugel/Ziel-punkt
		X,Y,Z - Kameraausrichtung auf XY,YZ und ZX Ebenen(Front-, Seiten- und Draufsicht)
		O,P - Punkte menge durchlaufen
		1,2,3,4 Arml�nge �ndern
		i - start/stop animate
		u - richtung der animation

*/

#define GLFW_DLL
#include <GLFW/glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

using namespace std;

#include "vec3.hpp"
#include "matrix.h"

static double alpha_ = 0;
static double beta_ = 0;
static double gamma_= 0;
static double delta_ = 0;
static double epsilon_ = 0;
static double zetta_ = 0;
static double bla_ = 0;
static double dt_ = 1;
static double angleXAxis_=0;
static double angleYAxis_=-90;
static double angleZAxis_=0;

static double a[2] = {15,15};
Vec3 p( 5, 0, 5);
Vec3 q(0,0,0);
Vec3 theta(M_PI/2, M_PI/2, M_PI/2);
//static double q[3] = {0, 0, 0};
int counter=0;
bool animate;
bool richtung;

double J[3][3];



static Vec3 mouseSpeed(0,0,0);

static double window_width_ = 1000;
static double window_height_ =1000;
static double zoomIn = 1;

vector<Vec3> getTorus(double radius1,double radius2,int resolution1,int resolution2){
	vector<Vec3>punkte;
	for (int i=0;i<=resolution1;i++){
		double alpha=i*2*M_PI/resolution1;
		for (int j=0;j<=resolution2;j++){
			double beta=j*2*M_PI/resolution2;
			punkte.push_back(Vec3((radius1+radius2*cos(beta))*cos(alpha),(radius1+radius2*cos(beta))*sin(alpha),radius2*sin(beta)));
		}
	}
	return punkte;
}


vector<Vec3> getBall(double size,int slices,double resolution){
	vector<Vec3>punkte;
	for (int i=-slices;i<=slices;i++){
		double aktslice=i*size/slices;
		double R=sqrt(pow(size,2)-pow(aktslice,2));
		for (int j=0;j<=resolution;j++){
			double aktPunkt=j*2*M_PI/resolution;
			punkte.push_back(Vec3(R*cos(aktPunkt),R*sin(aktPunkt),aktslice));
		}
	}
	return punkte;
}

void SetMaterialColor(int side, double r, double g, double b) {
  float amb[4], dif[4], spe[4];
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
    case 1:     mat = GL_FRONT;
      break;
    case 2:     mat = GL_BACK;
      break;
    default: mat = GL_FRONT_AND_BACK;
  }

  glMaterialfv(mat, GL_AMBIENT, amb);
  glMaterialfv(mat, GL_DIFFUSE, dif);
  glMaterialfv(mat, GL_SPECULAR, spe);
  glMaterialf( mat, GL_SHININESS, 20);
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
	glOrtho(-30*zoomIn, 30*zoomIn, -30*zoomIn, 30*zoomIn, -50*zoomIn, 50*zoomIn);

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



void DrawKrane(){
	double radius=1;
	glPushMatrix();
	glRotated(180+(180/M_PI) * theta.p[2],0,0,1);

	glRotated(-90+(180/M_PI) * theta.p[0],0,1,0);
	DrawSphere(Vec3(0,0,0),radius);
	DrawZylinder(radius,a[0],50);
	glTranslated(0,0,a[0]);
	glRotated((180/M_PI) * theta.p[1],0,1,0);
	DrawSphere(Vec3(0,0,0),radius);
	DrawSphere(Vec3(0,0,a[1]),radius);
	DrawZylinder(radius,a[1],50);

	glPopMatrix();
}

vector<Vec3>Punkte=getTorus(10,5,100,100);
void IK(){
	for(int i = 0; i <= 5; i++){
		Vec3 r=a[0]*Vec3(cos(theta.p[0]) * cos(theta.p[2]),
							cos(theta.p[0]) * sin(theta.p[2]),
							sin(theta.p[0]));

		Vec3 s=a[1]*Vec3(cos(theta.p[0] + theta.p[1]) * cos(theta.p[2]),
							cos(theta.p[0] + theta.p[1]) * sin(theta.p[2]),
							sin(theta.p[0] + theta.p[1]));

		Vec3 f = q + r + s;

		Vec3 F_theta1(-a[0] * sin(theta.p[0]) * cos(theta.p[2]) - a[1] * sin(theta.p[0] + theta.p[1]) * cos(theta.p[2]),
						-a[0] * sin(theta.p[0]) * sin(theta.p[2]) - a[1] * sin(theta.p[0] + theta.p[1]) * sin(theta.p[2]),
						a[0] * cos(theta.p[0]) + a[1] * cos(theta.p[0] + theta.p[1]));

		Vec3 F_theta2(-a[1] * sin(theta.p[0] + theta.p[1]) * cos(theta.p[2]),
						-a[1] * sin(theta.p[0] + theta.p[1]) * sin(theta.p[2]),
						a[1] * cos(theta.p[0] + theta.p[1]));

		Vec3 F_theta3(-a[0]*cos(theta.p[0]) *sin(theta.p[2]) - a[1]*cos(theta.p[0] + theta.p[1] ) * sin(theta.p[2]),
						a[0]*cos(theta.p[0]) *cos(theta.p[2]) + a[1]*cos(theta.p[0] + theta.p[1] ) * cos(theta.p[2]),
						0);

		for(int i=0;i<=2;i++){
			J[i][0] = F_theta1.p[i];
			J[i][1] = F_theta2.p[i];
			J[i][2] = F_theta3.p[i];
		}


		Matrix<double> Jacobi(3,3);
		for(int k=0;k<=2;k++){
			for(int l=0;l<=2;l++){
				Jacobi.put(k,l,J[k][l]);
			}
		}

		Matrix<double> JacobiRes(3,3);
		JacobiRes = (Jacobi.getTranspose() * Jacobi).getInverse() * Jacobi.getTranspose();

		Vec3 dtheta(0,0,0);
		for(int i = 0; i < 3; i++){
				for(int j = 0; j < 3; j++){
						dtheta.p[i]+=JacobiRes.get(i,j)*(p-f).p[j];
				}
		}
		theta = theta + dtheta;
	}
}


// draw the entire scene
void Preview() {
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();                                                 // Reset The Current Modelview Matrix
	glTranslated(0, 0, -30.0);      // Move 10 units backwards in z,
								  // since camera is at origin

	//glPushMatrix();
	glRotated(angleXAxis_,0,1,0);
	glRotated(angleYAxis_,1,0,0);
	glRotated(angleZAxis_,0,0,1);
	SetMaterialColor(3, 1, 0.5, 0);
	for(int g=0;g<counter;g++){
	//for(int g=0;g<(int)Punkte.size();g++){
		Vec3 p1=Punkte.at(g);
		glPushMatrix();
		glTranslated(p1.p[0],p1.p[1],p1.p[2]);
		DrawCube(0.2);
		glPopMatrix();
	}
	  if(animate){
		  if (counter<Punkte.size()-1 && !richtung){
			  p=Punkte.at(counter);
			  counter++;}
		  if(counter>0 && richtung){
			  p=Punkte.at(counter);
			  counter--;
		  }
	  }

	IK();
	DrawKrane();
	SetMaterialColor(3, 1, 0, 0);
	DrawSphere(p, 1.1);
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
	}


}
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods){
	if (key == GLFW_KEY_N)dt_+=0.1;
	if (key == GLFW_KEY_M)dt_-=0.1;
	if (key == GLFW_KEY_UP)angleYAxis_+=dt_;
	if (key == GLFW_KEY_DOWN)angleYAxis_-=dt_;
	if (key == GLFW_KEY_LEFT)angleXAxis_+=dt_;
	if (key == GLFW_KEY_RIGHT)angleXAxis_-=dt_;

	if (key == GLFW_KEY_J){
		zoomIn += 0.01;
	}

	if (key == GLFW_KEY_K){
		zoomIn -= 0.01;
	}


	if (key == GLFW_KEY_O){
		if(counter<Punkte.size()-1){
		  	  p=Punkte.at(counter);
		  	  counter++;
		}
	}
	if (key == GLFW_KEY_P){
		  if(counter>0){
			  p=Punkte.at(counter);
			  counter--;
		  }
	}
	if (key == GLFW_KEY_I && action == GLFW_PRESS){
		animate = !animate;
	}
	if (key == GLFW_KEY_U && action == GLFW_PRESS){
		richtung = !richtung;
	}



	if (key == GLFW_KEY_Y){
		angleXAxis_=0;
		angleYAxis_=-90;
		angleZAxis_=0;
	}
	if (key == GLFW_KEY_Z){
		angleXAxis_=0;
		angleYAxis_=0;
		angleZAxis_=0;
	}
	if (key == GLFW_KEY_X){
		angleXAxis_=0;
		angleYAxis_=-90;
		angleZAxis_=-90;
	}

	if (key == GLFW_KEY_W)p.p[2] +=dt_;
	if (key == GLFW_KEY_S)p.p[2] -=dt_;
	if (key == GLFW_KEY_A)p.p[0] -=dt_;
	if (key == GLFW_KEY_D)p.p[0] +=dt_;
	if (key == GLFW_KEY_Q)p.p[1] +=dt_;
	if (key == GLFW_KEY_E)p.p[1] -=dt_;

	if (key == GLFW_KEY_1)a[0]+=dt_;
	if (key == GLFW_KEY_2)a[0]-=dt_;
	if (key == GLFW_KEY_3)a[1]+=dt_;
	if (key == GLFW_KEY_4)a[1]-=dt_;


//      if (key == GLFW_KEY_1)alpha_+=dt_;
//      if (key == GLFW_KEY_2)alpha_-=dt_;
//      if (key == GLFW_KEY_3)beta_+=dt_;
//      if (key == GLFW_KEY_4)beta_-=dt_;
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

	if (key == GLFW_KEY_R){
		alpha_=0;
		beta_ =0;
		gamma_=0;
		delta_=0;
		epsilon_=0;
		zetta_=0;
		angleXAxis_=0;
		angleYAxis_=-90;
		angleZAxis_=0;
		dt_=1;
		alpha_=0;
		a[0]=15;
		a[1]=15;
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
