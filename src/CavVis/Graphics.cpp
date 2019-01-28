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
#include "Graphics.h"

namespace Viewer{
    
    /* Draw sphere */
    void Graphics::drawSphere(Vector4 center, float size, float R, float G, float B, bool wiresphere){
        glPushMatrix();
        glTranslatef(center.x, center.y, center.z);
        glColor3f(R, G, B);
        if (wiresphere) glutWireSphere(size, 30.0, 30.0);
        else glutSolidSphere(size, 30.0, 30.0);
        glPopMatrix();
    }
    
    /* Draw a sphere - Point object version */
    void Graphics::drawSphere(Point center, float size, float R, float G, float B, bool wiresphere){
        glPushMatrix();
        glTranslatef(center.coord.x, center.coord.y, center.coord.z);
        glColor3f(R, G, B);
        if (wiresphere) glutWireSphere(size, 30.0, 30.0);
        else glutSolidSphere(size, 30.0, 30.0);
        glPopMatrix();
    }
    
    /* Draw a sphere - Atom object version */
    void Graphics::drawSphere(Atom * atom, float size, float R, float G, float B, bool wiresphere){
        glPushMatrix();
        glTranslatef(atom->coord.x, atom->coord.y, atom->coord.z);
        glColor3f(R, G, B);
        if (wiresphere) glutWireSphere(size, 30.0, 30.0);
        else glutSolidSphere(size, 30.0, 30.0);
        glPopMatrix();
    }
    
    /* Draw protein surface */
    void Graphics::drawSurface(vector<mTriangle> * blob, Color color, bool wireframe, bool debug){
        if (wireframe) glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
        int s = (int) blob->size();
        if (debug) printf("Blob.size()=%d:\n", s);
        glBegin(GL_TRIANGLES);
        
        // Iterate triangles of the blob
        for(int i=0; i < blob->size(); i++){
            mTriangle * triangle = &blob->at(i);
            for(int j=0; j < 3; j++){
                mVertex * vertex = triangle->vertices.at(j);
                glColor3f(color.R, color.G, color.B);
                glNormal3f(vertex->normal.x, vertex->normal.y, vertex->normal.z);
                glVertex3f(vertex->coord.x, vertex->coord.y, vertex->coord.z);
                //printf("(%.1f, %.1f, %.1f)\n", vertex->coord.x, vertex->coord.y, vertex->coord.z);
            }
        }
        glEnd();
        if (wireframe) glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    }
    
}
