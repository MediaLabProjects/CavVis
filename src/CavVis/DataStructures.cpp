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
#include "DataStructures.h"
#include "Surface.h"
namespace CavVis{
    
    Vector4::Vector4 (float _x, float _y, float _z){
        x = _x;
        y = _y;
        z = _z;
    }
    
    Vector4::Vector4 (float _x, float _y, float _z, int _i, int _j, int _k){
        x = _x;
        y = _y;
        z = _z;
        i = _i;
        j = _j;
        k = _k;
    }
    
    /* Create vector using two points */
    Vector4 Vector4::Vector(Vector4 * b){
        Vector4 res;
        // [a] is the object point
        res.x = b->x-x;
        res.y = b->y-y;
        res.z = b->z-z;
        return(res);
    }
    
    /* Compute distance between two points */
    float Vector4::Distance(Vector4 * a){
        float xd = x - a->x;
        float yd = y - a->y;
        float zd = z - a->z;
        return sqrt(xd*xd + yd*yd  + zd*zd);
    }
    
    /* Normalize current Vector4 object */
    void Vector4::Normalize(){
        if (length == 0.0) Length();
        x = x/length;
        y = y/length;
        z = z/length;
    }
    
    /* Compute length */
    void Vector4::Length () { length = sqrt(x * x + y * y + z * z); }
    
    Atom::Atom(int _id, string _type, float _x, float _y, float _z, float _radius, Color _color){
        id      = _id;
        type    = _type;
        coord.x       = _x;
        coord.y       = _y;
        coord.z       = _z;
        radius  = _radius;
        color   = _color;
    }
    
    /* Check if point is inside the current atom/sphere */
    bool Atom::pointInside(Point * p){
        if (((coord.x - p->coord.x)*(coord.x - p->coord.x)) + ((coord.y - p->coord.y)*(coord.y - p->coord.y)) + ((coord.z - p->coord.z)*(coord.z - p->coord.z)) <= radius*2){
            return true;
        }else
            return false;
    }
    
    /* Interpolates atom center to a parsed grid */
    Point Atom::interpolate(Grid * grid){
        Point interpolated;
        // This follows the same implementation detail on the marchingCubes algorithm
        // (i.e. When the grid is initialized on the startGrid() method)
        interpolated.coord.i = roundf( (coord.x - (grid->minx-grid->padding)) / grid->gridSpacing);
        interpolated.coord.j = roundf( ((grid->maxy+grid->padding) - coord.y) / grid->gridSpacing);
        interpolated.coord.k = roundf( ((grid->maxz+grid->padding) - coord.z) / grid->gridSpacing);
        interpolated.coord.x = grid->X[interpolated.coord.i];
        interpolated.coord.y = grid->Y[interpolated.coord.j];
        interpolated.coord.z = grid->Z[interpolated.coord.k];
        return(interpolated);
    }

    /* Check if point is inside the current atom/sphere */
    bool Atom::pointInside(Vector4 * p){
        if (((coord.x - p->x)*(coord.x - p->x)) + ((coord.y - p->y)*(coord.y - p->y)) + ((coord.z - p->z)*(coord.z - p->z)) <= radius*2){
            return true;
        }else
            return false;
    }

    Color::Color(float _R, float _G, float _B){
        R = _R;
        G = _G;
        B = _B;
    }
    
    /* Generates a random color to the current color object */
    void Color::random(){
        R = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        G = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        B = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    }
    
    Point::Point(float _x, float _y, float _z){
        coord.x = _x;
        coord.y = _y;
        coord.z = _z;
    }
    
    /* Interpolates point to the parsed grid */
    Point Point::interpolate(Grid * grid){
        Point interpolated;
        // This follows the same implementation detail on the marchingCubes algorithm
        // (i.e. When the grid is initialized on the startGrid() method)
        interpolated.coord.i = roundf( (coord.x - (grid->minx-grid->padding)) / grid->gridSpacing);
        interpolated.coord.j = roundf( ((grid->maxy+grid->padding) - coord.y) / grid->gridSpacing);
        interpolated.coord.k = roundf( ((grid->maxz+grid->padding) - coord.z) / grid->gridSpacing);
        interpolated.coord.x = grid->X[interpolated.coord.i];
        interpolated.coord.y = grid->Y[interpolated.coord.j];
        interpolated.coord.z = grid->Z[interpolated.coord.k];
        return(interpolated);
    }
    
    /* Creates point object with (x,y,z) equal to the geometric center from a set of objects */
    void Point::computeGeometricCenter(vector<Atom *> points){
        double sumX=0.0, sumY=0.0, sumZ=0.0;
        // Sum
        for(int i=0; i <points.size(); i++){
            sumX += points[i]->coord.x;
            sumY += points[i]->coord.y;
            sumZ += points[i]->coord.z;
        }
        // Divide by number of available atoms
        coord.x = sumX / points.size();
        coord.y = sumY / points.size();
        coord.z = sumZ / points.size();
    }
    
    /* Creates point object with (x,y,z) equal to the geometric center from a set of objects */
    void Point::computeGeometricCenter(vector<Atom> * points){
        double sumX=0.0, sumY=0.0, sumZ=0.0;
        // Sum
        for(int i=0; i <points->size(); i++){
            sumX += points->at(i).coord.x;
            sumY += points->at(i).coord.y;
            sumZ += points->at(i).coord.z;
        }
        // Divide by number of available atoms
        coord.x = sumX / points->size();
        coord.y = sumY / points->size();
        coord.z = sumZ / points->size();
    }
    
    /* Creates point object with (x,y,z) equal to the geometric center from a set of objects */
    void Point::computeGeometricCenter(vector<Point> * points){
        double sumX=0.0, sumY=0.0, sumZ=0.0;
        // Sum
        for(int i=0; i <points->size(); i++){
            Point p = (*points).at(i);
            sumX += p.coord.x;
            sumY += p.coord.y;
            sumZ += p.coord.z;
        }
        // Divide by number of available atoms
        coord.x = sumX / points->size();
        coord.y = sumY / points->size();
        coord.z = sumZ / points->size();
    }
    
}
