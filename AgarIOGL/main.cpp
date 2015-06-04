
#include "main.h"
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

std::vector<Blob> blopList;

static Vec3 mouseSpeed(0,0,0);

static double window_width_ = 1000;
static double window_height_ = 1000;
double speed=1;



int randomInt(int min, int max){
	//srand(time(NULL)); // Seed the time
	int finalNum = rand()%(max-min)+min; // Generate the number, assign to variable.
	return finalNum;
}
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
Vec3 spiegelungUmNomrale(Vec3 V,Vec3 A){
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
  glOrtho(-50, 50, -50, 50, -50, 50);

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

//class blop{
//	double r=0;
//	double g=0;
//	double b=0;
//	double mass=0;
//	static Vec3 KAY_(0,1,0);
//
//    blob(double r, double g, double b, double mass, Vec3 postion){
//
//    };
//	public:
//    	void set_values (int,int);
//    	int area() {return width*height;}
//};

Blob::Blob(Vec3 farbe, double mass, Vec3 position){
	this->farbe=farbe;
	this->mass=mass;
	this->position=position;
}
//Blob::Blob(){
//	this->farbe(randomInt(0,255)/255,randomInt(0,255)/255,randomInt(0,255)/255);
//	this->mass=randomInt(1,2);
//	this->position(randomInt(0,255)/255,randomInt(0,255)/255,randomInt(0,255)/255);
//}
double Blob::calcRadius(){
	return pow((3*this->mass)/(4*M_PI),1/3);
}
bool Blob::checkForCollision(Blob andererBlob){
	return ((this->position-andererBlob.position).Length())<=this->calcRadius()+andererBlob.calcRadius();
}
void Blob::drawBlop(){
	DrawSphere(this->position,calcRadius());
}
void Blob::eatBlob(Blob andererBlob){
	this->mass+=andererBlob.mass;
}

void checkForBlobCollisions(){
	vector<int> intlist=new vector<int>();
	for(int i : intlist){

	}
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



double kugelrotation;
void DrawBillardtisch(){
	double x=20+beta_;
	double z=15+alpha_;
	double h=2+delta_;
	double rb=1+gamma_;
	double kugelradius=1+epsilon_;
	double scaling=0.0001;
	double friction=1-0.0000;

	Vec3 hoehe(0,0,h);
	Vec3 P1(-x,z,0);
	Vec3 P2(x,z,0);
	Vec3 P3(x,-z,0);
	Vec3 P4(-x,-z,0);

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

	Vec3 W1(-x,z/2,0);
	Vec3 W2(x,z/2,0);
	DrawQuad(W1,W1+e2*h,W2+e2*h,W2);

	Vec3 W3(bla_,-z,0);
	Vec3 W4(-x,zetta_,0);
	DrawQuad(W3,W3+e2*h,W4+e2*h,W4);

	Vec3 W5(x,z,0);
	Vec3 W6(-x,0,0);
	DrawQuad(W5,W5+e2*h,W6+e2*h,W6);

	Vec3 Pfeiler(x/2+posX2_,0,-z/2+posZ2_);
	double pfeilerRadius=2;
	glPushMatrix();
	glTranslated(Pfeiler.p[0],Pfeiler.p[1],Pfeiler.p[2]);
	glRotated(-90,1,0,0);
	DrawZylinder(pfeilerRadius,h,100);
	glPopMatrix();

	mouseSpeed.p[0]=mouseSpeed.p[0]*friction;
	mouseSpeed.p[2]=mouseSpeed.p[2]*friction;

	posX_+=mouseSpeed.p[0]*scaling;
	posZ_+=mouseSpeed.p[2]*scaling;

	if(abs(posX_-x)-kugelradius<0){
		mouseSpeed.p[0]=-mouseSpeed.p[0];
	}
	if(abs(-posX_-x)-kugelradius<0){
		mouseSpeed.p[0]=-mouseSpeed.p[0];
	}
	if(abs(posZ_-z)-kugelradius<0){
		mouseSpeed.p[2]=-mouseSpeed.p[2];
	}
	if(abs(-posZ_-z)-kugelradius<0){
		mouseSpeed.p[2]=-mouseSpeed.p[2];
	}

	Vec3 W1W2ORTHO((W1-W2).p[2],0,-(W1-W2).p[0]);
	if(getDistanceToLine(W1W2ORTHO,Vec3(posX_,0,posZ_),W1W2ORTHO/W1W2ORTHO.Length()*W1)<=0){
				//mouseSpeed=spiegelungUrsprungsgerade(getAngle2D((W1-W2).p[0],(W1-W2).p[2],1,0),mouseSpeed);
				mouseSpeed=spiegelungUmNomrale(mouseSpeed,W1W2ORTHO);
		}

	Vec3 W3W4ORTHO((W3-W4).p[2],0,-(W3-W4).p[0]);
	if(getDistanceToLine(W3W4ORTHO,Vec3(posX_,0,posZ_),W3W4ORTHO/W3W4ORTHO.Length()*W3)<=0){
				//mouseSpeed=spiegelungUrsprungsgerade(getAngle2D((W3-W4).p[0],(W3-W4).p[2],1,0),mouseSpeed);
				mouseSpeed=spiegelungUmNomrale(mouseSpeed,W3W4ORTHO);
		}
	Vec3 W5W6ORTHO((W5-W6).p[2],0,-(W5-W6).p[0]);
	if(getDistanceToLine(W5W6ORTHO,Vec3(posX_,0,posZ_),W5W6ORTHO/W5W6ORTHO.Length()*W5)<=0){
				//mouseSpeed=spiegelungUrsprungsgerade(getAngle2D((W5-W6).p[0],(W5-W6).p[2],1,0),mouseSpeed);
				mouseSpeed=spiegelungUmNomrale(mouseSpeed,W5W6ORTHO);
		}

	if((Pfeiler-Vec3(posX_,0,posZ_)).Length()-kugelradius-pfeilerRadius<=0){
		mouseSpeed=spiegelungUrsprungsgerade(getAngle2D(-(Pfeiler-Vec3(posX_,0,posZ_)).p[2],(Pfeiler-Vec3(posX_,0,posZ_)).p[0],1,0),mouseSpeed);

	}



	glPushMatrix();
	kugelrotation+=(sqrt(pow(mouseSpeed.p[0],2)+pow(mouseSpeed.p[2],2))/50)/(2*M_PI*kugelradius);

	glTranslated(posX_,kugelradius,posZ_);
	//glRotated(getAngleToXAxis(mouseSpeed.p[0],mouseSpeed.p[2]),0,1,0);
	glRotated(kugelrotation,mouseSpeed.p[2],0,-mouseSpeed.p[0]);
	DrawSphere(Vec3(0,0,0),kugelradius);
	glPopMatrix();
	//getAngleToXAxis(1,0);

}
// draw the entire scene
void Preview() {
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();						    // Reset The Current Modelview Matrix
  glTranslated(0, 0, -30.0);      // Move 10 units backwards in z,
                                  // since camera is at origin

  //glPushMatrix();
//  glRotated(angleXAxis_,0,1,0);
//  glRotated(angleYAxis_,1,0,0);
  SetMaterialColor(3, 1, 0.5, 0);

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
		  alpha_=0;
		  beta_ =0;
		  gamma_=0;
		  delta_=0;
		  epsilon_=0;
		  zetta_=0;
		  angleXAxis_=20;
		  angleYAxis_=45;
		  posX_=0;
		  posY_=0;
		  posZ_=0;
		  dt_=1;
		  alpha_=0;

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
