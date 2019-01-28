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
#ifndef __CavVis__Graphics__
#define __CavVis__Graphics__
#include <stdio.h>
#ifdef __APPLE__ // Mac OS X
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else // Linux
#include <GL/glut.h>
#endif
#include "DataStructures.h"
using namespace std;
using namespace CavVis;

namespace Viewer{
    class Graphics{
    public:
        /* Draw sphere - Point object version */
        static void drawSphere(Point center, float size, float R, float G, float B, bool wiresphere);
        /* Draw sphere - Vec3 object version */
        static void drawSphere(Vector4 center, float size, float R, float G, float B, bool wiresphere);
        /* Draw a sphere - Atom object version */
        static void drawSphere(Atom * atom, float size, float R, float G, float B, bool wiresphere);
        /* Draw protein surface */
        static void drawSurface(vector<mTriangle> * blob, Color color, bool wireframe, bool debug);
        Graphics(){};
        ~Graphics(){};
    private:
        //
    };
}

#endif /* defined(__CavVis__Graphics__) */
