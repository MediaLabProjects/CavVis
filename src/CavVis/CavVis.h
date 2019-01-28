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
#ifndef __CavVis__CavVis__
#define __CavVis__CavVis__
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <algorithm>
#include <set>
#include <float.h>
#include <math.h>
#include "DataStructures.h"
#include "Utils.h"
#include "Surface.h"
#include "tree.h"
#include "neighbors.h"
#include "dbscan.h"
#define PI 3.14159265
namespace CavVis{
    
    class Algorithm{
    public:
        bool variable_fov=true;
        Algorithm(vector<Point *> * points, Protein * protein, float gridSpacing, int top, float distance, bool areas, bool ffilter, bool viewer, bool output);
        /* Filter points. Check if the normal of the point is pointing towards the surface (valid)
         or to the outside of the surface (invalid). */
        static bool validPoint(Grid * grid, Point * point, float threshold);
        /* Voxel Ray Tracer implementation, returns true if the first intersected voxel is a surface one */
        static bool voxelRayTracer(Point * pt, Grid * grid, float &distance);
        /* Voxel Ray Tracer implementation, return the first intersected voxel */
        Point voxelRayTracer_alt(Point * pt, Grid * grid);
    private:
        /* Compute cavity areas */
        void computeAreas(Protein * protein, float cgridspacing);
        /* Auxiliar function of the visibility method, run tests (adopted solution 3)*/
        bool Visibility_tests(Protein * protein, Grid * grid, Point * a, bool debug);
        /* Main method */
        void Visibility(Protein * protein, Grid * grid, vector<Point *> * points);
        /* Test if point [b] is over the plane of point [a] */
        bool test1(Point * a, Point * b, bool debug);
        /* Test if point [a] is facing point [b] by using their normal vectors [na, nb] */
        static bool test2(Vector4 * na, Vector4 * nb);
        /* Check if point [b] is visible to [a] */
        static bool test3(Point * a, Point * b);
        /* Recursive method that finds points visible from an initial point [a] and its visibles [a->final] */
        void find(vector<Point *> * R, Point * a, Point * original, int& inc, bool debug, int& gid, set<int> * joint);
        /* Union all set of visible point that are visible to each other */
        void Union(vector<Point *> * points, vector<Group> * groups, bool debug, bool output);
        /* Auxiliar method for formCavities */
        bool intersect(vector<int> * A, vector<int> * B);
        /* Form cavities by joining and filtering groups */
        void formCavities(float gridspacing, vector<Group> * visibles, vector<Cavity> * cavities, bool debug, bool output);
        /* Identify the set of atoms that fills each cavity */
        void IdentifyFillingAtoms(Grid * grid, vector<Cavity> * cavities, float gridSpacing, bool filter_spheres, bool debug, bool output);
        /* Cluster final list of cavities using a algorithm based on N-ary Trees - Legacy support */
        void cluster(Protein * protein, vector<Cavity> * cavities, int topc, float distance, bool debug, bool output);
    };
}
#endif /* defined(__CavVis__CavVis__) */
