//---------------------------------------------------------------------------------------
//
//	CAVVIS - A Field-of-View Geometric Algorithm for Protein Cavity Detection
//
//  Custom Neighbors Library
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
#ifndef __CavVis__neighbors__
#define __CavVis__neighbors__
#include <stdio.h>
#include "DataStructures.h"
#include "Surface.h"
#include "CavVis.h"
using namespace CavVis;
namespace Neighbors{
    struct Voxel{
    public:
        int id;                                         // id of the voxel (can be misleading)
        int i,j,k;                                      // Position of the current voxel in the grid
        int inside=0;                                   // Number of points inside the voxel
        Point centerPoint;                              // Center point of the voxel
        vector<Point> vertices;                         // Vertices of the current voxel (to avoid re-computation)
        vector<Voxel> neighbors;                        // Voxel neighbors
        float max_x, max_y, max_z, min_x, min_y, min_z; // Max and min values when the voxel is a bounding box
    private:
        //
    };
    
    /* Standard main method of the Neighbors Library */
    vector<Point> getNeighborsOfVerticeAt(Grid * grid, Point currentVertice, int at);
    
}
#endif /* defined(__CavVis__neighbors__) */
