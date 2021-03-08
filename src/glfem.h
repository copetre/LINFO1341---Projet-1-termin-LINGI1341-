/*
 *  glfem.h - BOV version
 *  Library for EPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 *  GLFW  http://www.glfw.org/ (version utilisée 3.3.2)
 *  BOV   https://git.immc.ucl.ac.be/hextreme/NGP/-/tree/master/deps/BOV
 *
 */
 
 

#ifndef _GLFEM_H_
#define _GLFEM_H_

#include <stdlib.h>
#include "BOV.h"
#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>
#include "fem.h"


static float GLFEM_BLACK[4] = {0.0,0.0,0.0,1.0};
static float GLFEM_BLUE[4]  = {0.0,0.0,1.0,1.0};
static float GLFEM_RED[4]   = {1.0,0.0,0.0,1.0};
static float GLFEM_WHITE[4]   = {1.0,1.0,1.0,1.0};

void          	glfemSetColor(float color[4]);
void          	glfemSetTextColor(float color[4]);
void 		glfemSetLineWidth(float width);

static void   	glfemKeyCallback(GLFWwindow* self,int key,int scancode,int action,int mods);

void          	glfemWindowCreate(const char *windowName,int w,int h,int n,double *x,double *y);
void          	glfemWindowUpdate();
void          	glfemWindowFree();
int           	glfemWindowShouldClose();

void          	glfemDrawMessage(char *message, double pos[2]);

void 		glfemDrawElement(double *x, double *y, int n);
void          	glfemPlotMesh(femMesh *theMesh);
double        	glfemScale(double minimum, double maximum, double value);
void 	     	glfemPlotSolution(femMesh* theMesh, double *u);

#endif


