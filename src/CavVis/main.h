//---------------------------------------------------------------------------------------
//
//	CAVVIS - A Field-of-View Geometric Algorithm for Protein Cavity Detection
//
//  Copyright (C) 2019 Instituto de Telecomunicações (www.it.pt)
//  Copyright (C) 2019 Universidade da Beira Interior (www.ubi.pt)
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Contacts:
//      Tiago Simões <tiago.simoes@it.ubi.pt>
//      Abel Gomes <agomes@di.ubi.pt>
//---------------------------------------------------------------------------------------
#ifndef CavVis_main_h
#define CavVis_main_h
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#ifdef __APPLE__ // Mac OS X
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else // Linux
#include <GL/glut.h>
#endif
#include "zpr.h"
#include "DataStructures.h"
#include "Parsers.h"
#include "Utils.h"
#include "CavVis.h"
#include "Graphics.h"
using namespace std;
using namespace CavVis;

string version                      = "1.0";  // Please check Changelog.txt for details.

// ---------------------------------------
// User defined parameters
// ---------------------------------------
float gridSpacing                   = 0.6f;         // 0.6 is the optimal grid spacing value for accurate results,
                                                    // however, 0.8 may be used for better time performance results.
int top                             = 10;           // Default top
bool dareas                         = false;        // Computes cavity areas
float tdistance                     = 4.0;          // Default distance, 1.4 can also be used
bool ffilter                        = true;         // Filter filling spheres
string input                        = "";           // Input protein filename
string output                       = "";           // Output directory
bool opengl_viewer                  = false;        // Use embedded opengl viewer
float v_gridSpacing                 = 0.3f;         // Viewer default grid spacing (also used to generate cavity surfaces)


// ---------------------------------------
// Embedded opengl viewer parameters
// ---------------------------------------
int  window_width                   = 900;
int  window_height                  = 800;
static GLfloat light_ambient[]  = { 0.0, 0.0, 0.0, 1.0 };
static GLfloat light_diffuse[]  = { 1.0, 1.0, 1.0, 1.0 };
static GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
static GLfloat light_position0[]= { 1.0, 1.0, 1.0, 0.0 };
static GLfloat mat_ambient[]    = { 0.7, 0.7, 0.7, 1.0 };
static GLfloat mat_diffuse[]    = { 0.8, 0.8, 0.8, 1.0 };
static GLfloat mat_specular[]   = { 1.0, 1.0, 1.0, 1.0 };
static GLfloat high_shininess[] = { 100.0 };
// --------------------------------------

Protein protein;                    // Supplied Protein
string PDB_id;                      // id of the protein
vector<mTriangle> * mblob;          // Protein surface
vector<Point *> points;             // main point vector of point objects
vector<vector<mTriangle>> surfaces; // Surface of each cavity

/* Generate Cavity Surfaces */
void cavitySurfaces();
/* OpenGL core functions */
void glut(int argc, char * argv[], string title, int width, int height);
/* OpenGL Draw objects */
void drawObjects();
/* OpenGL core functions */
void displayFunc(void);
/* OpenGL core functions */
void init();
/* Handle terminal related stuff */
void printHelp();
/* Handle terminal related stuff */
bool processArgs(int argc, char** argv);

#endif
