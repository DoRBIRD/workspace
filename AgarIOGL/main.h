/*
 * main.h
 *
 *  Created on: 18.03.2015
 *      Author: Jonas
 */

#ifndef MAIN_H_
#define MAIN_H_
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

class Blob{
	private:
		Vec3 farbe;
		double mass;
		Vec3 position;
	public:
		Blob(Vec3 farbe, double mass, Vec3 position);
		Blob();
		double calcRadius();
		bool checkForCollision(Blob andererBlob);
		void drawBlop();
		void eatBlob(Blob andererBlob);
};


#endif /* MAIN_H_ */
