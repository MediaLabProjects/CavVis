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
#include "Surface.h"
namespace CavVis{
        
    /* Main MC method */
    void Surface::MarchingCubes(bool viewer, Protein * _protein, Grid * _grid, bool output){
        grid    = _grid;
        protein = _protein;
        // Dynamically compute parameters (ioc and seg)
        grid->ioc=(protein->max_radius/grid->gridSpacing)*2;
        grid->seg=(int) protein->max_radius/grid->gridSpacing;
        
        startGrid(output);
        vertexIntensity();
        
        int vt; int i=0; int j=0; int k=0; int m=0;
        float edge[12][9];
        // The first 3 positions are coordinates (x,y,z) of the interpolated triangle vertex
        // The last 6  positions are coordinates (i,j,k) of points of the edge used to interpolate the triangle vertex
        if (output) printf("  Starting to Map Triangle Vertices\n");
        for(i=1;i<grid->size_x-1;i++){
            for(j=1;j<grid->size_y-1;j++){
                for(k=1;k<grid->size_z-1;k++){
                    
                    // Associates each cell with an edge and then calculates the interpolation in the vertex
                    // Avoid the first indice, it was init with 0
                    if(indexs[i][j][k]>0){
                        if (edges_table[cube[i][j][k]] & 1){
                            edge[0][0] =vertex_interpolation_x(i,j+1,k+1, i+1,j+1,k+1);
                            edge[0][1] =vertex_interpolation_y(i,j+1,k+1, i+1,j+1,k+1);
                            edge[0][2] =vertex_interpolation_z(i,j+1,k+1, i+1,j+1,k+1);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[0][3] = i;
                            edge[0][4] = j+1;
                            edge[0][5] = k+1;
                            //
                            edge[0][6] = i+1;
                            edge[0][7] = j+1;
                            edge[0][8] = k+1;
                            
                            //printf("Storing (%f, %f, %f)<->(%f, %f, %f)\n", edge[0][3], edge[0][4], edge[0][5], edge[0][6], edge[0][7], edge[0][8]);
                            
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point(); // triangle vertex
                            tVertex->coord.x = edge[0][0];
                            tVertex->coord.y = edge[0][1];
                            tVertex->coord.z = edge[0][2];
                            // One position can be mapped to more than one triangle vertex
                            grid->map[i][j+1][k+1].tPoints->insert(tVertex);
                            grid->map[i][j+1][k+1].surface      = true;
                            //
                            grid->map[i+1][j+1][k+1].tPoints->insert(tVertex);
                            grid->map[i+1][j+1][k+1].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 2){
                            edge[1][0] =vertex_interpolation_x(i+1,j+1,k+1, i+1,j+1,k);
                            edge[1][1] =vertex_interpolation_y(i+1,j+1,k+1, i+1,j+1,k);
                            edge[1][2] =vertex_interpolation_z(i+1,j+1,k+1, i+1,j+1,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[1][3] = i+1;
                            edge[1][4] = j+1;
                            edge[1][5] = k+1;
                            //
                            edge[1][6] = i+1;
                            edge[1][7] = j+1;
                            edge[1][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[1][0];
                            tVertex->coord.y = edge[1][1];
                            tVertex->coord.z = edge[1][2];
                            grid->map[i+1][j+1][k+1].tPoints->insert(tVertex);
                            grid->map[i+1][j+1][k+1].surface  = true;
                            //
                            grid->map[i+1][j+1][k].tPoints->insert(tVertex);
                            grid->map[i+1][j+1][k].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 4){
                            edge[2][0] =vertex_interpolation_x(i,j+1,k, i+1,j+1,k);
                            edge[2][1] =vertex_interpolation_y(i,j+1,k, i+1,j+1,k);
                            edge[2][2] =vertex_interpolation_z(i,j+1,k, i+1,j+1,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[2][3] = i;
                            edge[2][4] = j+1;
                            edge[2][5] = k;
                            //
                            edge[2][6] = i+1;
                            edge[2][7] = j+1;
                            edge[2][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[2][0];
                            tVertex->coord.y = edge[2][1];
                            tVertex->coord.z = edge[2][2];
                            grid->map[i][j+1][k].tPoints->insert(tVertex);
                            grid->map[i][j+1][k].surface  = true;
                            //
                            grid->map[i+1][j+1][k].tPoints->insert(tVertex);
                            grid->map[i+1][j+1][k].surface= true;
                        }
                        if (edges_table[cube[i][j][k]] & 8){
                            edge[3][0] =vertex_interpolation_x(i,j+1,k+1, i,j+1,k);
                            edge[3][1] =vertex_interpolation_y(i,j+1,k+1, i,j+1,k);
                            edge[3][2] =vertex_interpolation_z(i,j+1,k+1, i,j+1,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[3][3] = i;
                            edge[3][4] = j+1;
                            edge[3][5] = k+1;
                            //
                            edge[3][6] = i;
                            edge[3][7] = j+1;
                            edge[3][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[3][0];
                            tVertex->coord.y = edge[3][1];
                            tVertex->coord.z = edge[3][2];
                            grid->map[i][j+1][k+1].tPoints->insert(tVertex);
                            grid->map[i][j+1][k+1].surface  = true;
                            //
                            grid->map[i][j+1][k].tPoints->insert(tVertex);
                            grid->map[i][j+1][k].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 16){
                            edge[4][0] =vertex_interpolation_x(i,j,k+1, i+1,j,k+1);
                            edge[4][1] =vertex_interpolation_y(i,j,k+1, i+1,j,k+1);
                            edge[4][2] =vertex_interpolation_z(i,j,k+1, i+1,j,k+1);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[4][3] = i;
                            edge[4][4] = j;
                            edge[4][5] = k+1;
                            //
                            edge[4][6] = i+1;
                            edge[4][7] = j;
                            edge[4][8] = k+1;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[4][0];
                            tVertex->coord.y = edge[4][1];
                            tVertex->coord.z = edge[4][2];
                            grid->map[i][j][k+1].tPoints->insert(tVertex);
                            grid->map[i][j][k+1].surface   = true;
                            //
                            grid->map[i+1][j][k+1].tPoints->insert(tVertex);
                            grid->map[i+1][j][k+1].surface = true;
                        }
                        if (edges_table[cube[i][j][k]] & 32){
                            edge[5][0] =vertex_interpolation_x(i+1,j,k+1, i+1,j,k);
                            edge[5][1] =vertex_interpolation_y(i+1,j,k+1, i+1,j,k);
                            edge[5][2] =vertex_interpolation_z(i+1,j,k+1, i+1,j,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[5][3] = i+1;
                            edge[5][4] = j;
                            edge[5][5] = k+1;
                            //
                            edge[5][6] = i+1;
                            edge[5][7] = j;
                            edge[5][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[5][0];
                            tVertex->coord.y = edge[5][1];
                            tVertex->coord.z = edge[5][2];
                            grid->map[i+1][j][k+1].tPoints->insert(tVertex);
                            grid->map[i+1][j][k+1].surface  = true;
                            //
                            grid->map[i+1][j][k].tPoints->insert(tVertex);
                            grid->map[i+1][j][k].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 64){
                            edge[6][0] =vertex_interpolation_x(i,j,k, i+1,j,k);
                            edge[6][1] =vertex_interpolation_y(i,j,k, i+1,j,k);
                            edge[6][2] =vertex_interpolation_z(i,j,k, i+1,j,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[6][3] = i;
                            edge[6][4] = j;
                            edge[6][5] = k;
                            //
                            edge[6][6] = i+1;
                            edge[6][7] = j;
                            edge[6][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[6][0];
                            tVertex->coord.y = edge[6][1];
                            tVertex->coord.z = edge[6][2];
                            grid->map[i][j][k].tPoints->insert(tVertex);
                            grid->map[i][j][k].surface  = true;
                            //
                            grid->map[i+1][j][k].tPoints->insert(tVertex);
                            grid->map[i+1][j][k].surface= true;
                        }
                        if (edges_table[cube[i][j][k]] & 128){
                            edge[7][0] =vertex_interpolation_x(i,j,k+1, i,j,k);
                            edge[7][1] =vertex_interpolation_y(i,j,k+1, i,j,k);
                            edge[7][2] =vertex_interpolation_z(i,j,k+1, i,j,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[7][3] = i;
                            edge[7][4] = j;
                            edge[7][5] = k+1;
                            //
                            edge[7][6] = i;
                            edge[7][7] = j;
                            edge[7][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[7][0];
                            tVertex->coord.y = edge[7][1];
                            tVertex->coord.z = edge[7][2];
                            grid->map[i][j][k+1].tPoints->insert(tVertex);
                            grid->map[i][j][k+1].surface  = true;
                            //
                            grid->map[i][j][k].tPoints->insert(tVertex);
                            grid->map[i][j][k].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 256){
                            edge[8][0] =vertex_interpolation_x(i,j+1,k+1, i,j,k+1);
                            edge[8][1] =vertex_interpolation_y(i,j+1,k+1, i,j,k+1);
                            edge[8][2] =vertex_interpolation_z(i,j+1,k+1, i,j,k+1);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[8][3] = i;
                            edge[8][4] = j+1;
                            edge[8][5] = k+1;
                            //
                            edge[8][6] = i;
                            edge[8][7] = j;
                            edge[8][8] = k+1;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[8][0];
                            tVertex->coord.y = edge[8][1];
                            tVertex->coord.z = edge[8][2];
                            grid->map[i][j+1][k+1].tPoints->insert(tVertex);
                            grid->map[i][j+1][k+1].surface  = true;
                            //
                            grid->map[i][j][k+1].tPoints->insert(tVertex);
                            grid->map[i][j][k+1].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 512){
                            edge[9][0] =vertex_interpolation_x(i+1,j+1,k+1, i+1,j,k+1);
                            edge[9][1] =vertex_interpolation_y(i+1,j+1,k+1, i+1,j,k+1);
                            edge[9][2] =vertex_interpolation_z(i+1,j+1,k+1, i+1,j,k+1);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[9][3] = i+1;
                            edge[9][4] = j+1;
                            edge[9][5] = k+1;
                            //
                            edge[9][6] = i+1;
                            edge[9][7] = j;
                            edge[9][8] = k+1;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[9][0];
                            tVertex->coord.y = edge[9][1];
                            tVertex->coord.z = edge[9][2];
                            grid->map[i+1][j+1][k+1].tPoints->insert(tVertex);
                            grid->map[i+1][j+1][k+1].surface  = true;
                            //
                            grid->map[i+1][j][k+1].tPoints->insert(tVertex);
                            grid->map[i+1][j][k+1].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 1024){
                            edge[10][0] =vertex_interpolation_x(i+1,j+1,k, i+1,j,k);
                            edge[10][1] =vertex_interpolation_y(i+1,j+1,k, i+1,j,k);
                            edge[10][2] =vertex_interpolation_z(i+1,j+1,k, i+1,j,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[10][3] = i+1;
                            edge[10][4] = j+1;
                            edge[10][5] = k;
                            //
                            edge[10][6] = i+1;
                            edge[10][7] = j;
                            edge[10][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[10][0];
                            tVertex->coord.y = edge[10][1];
                            tVertex->coord.z = edge[10][2];
                            grid->map[i+1][j+1][k].tPoints->insert(tVertex);
                            grid->map[i+1][j+1][k].surface  = true;
                            //
                            grid->map[i+1][j][k].tPoints->insert(tVertex);
                            grid->map[i+1][j][k].surface    = true;
                        }
                        if (edges_table[cube[i][j][k]] & 2048){
                            edge[11][0] =vertex_interpolation_x(i,j+1,k, i,j,k);
                            edge[11][1] =vertex_interpolation_y(i,j+1,k, i,j,k);
                            edge[11][2] =vertex_interpolation_z(i,j+1,k, i,j,k);
                            // Store edge information (node positions) used to interpolate current triangle vertex
                            edge[11][3] = i;
                            edge[11][4] = j+1;
                            edge[11][5] = k;
                            //
                            edge[11][6] = i;
                            edge[11][7] = j;
                            edge[11][8] = k;
                            //
                            // Store interpolated point on the corresponding position of the grid
                            Point * tVertex = new Point();; // triangle vertex
                            tVertex->coord.x = edge[11][0];
                            tVertex->coord.y = edge[11][1];
                            tVertex->coord.z = edge[11][2];
                            grid->map[i][j+1][k].tPoints->insert(tVertex);
                            grid->map[i][j+1][k].surface  = true;
                            //
                            grid->map[i][j][k].tPoints->insert(tVertex);
                            grid->map[i][j][k].surface    = true;
                        }
                        
                        // Map triangle vertices
                        // The main idea is that equal vertex that are mapped as the same vertex to each triangle.
                        // For example, triangle A was mapped with vertex 1 and triangle B was also mapped with the same
                        // vertex - In this case the vertex object is the same for both triangles (pointer).
                        for (m=0;case_table[cube[i][j][k]][m]!=-1;m+=3) {
                            vt=indexs[i][j][k]+m;
                            
                            // --------------------------------------
                            // ### Vertex [1] of the current triangle
                            // --------------------------------------
                            // Check if the current vertex was allocated to other triangle
                            map<vid, mVertex *>::iterator found = mapVertices.find(tie(edge[case_table[cube[i][j][k]][m]][0],
                                                                                       edge[case_table[cube[i][j][k]][m]][1],
                                                                                       edge[case_table[cube[i][j][k]][m]][2]));
                            if (found != mapVertices.end()){
                                found->second->repeated++; // Increment the number of triangles that are sharing the current vertex
                                // IF found, set the current vertex of the triangle pointing to that one
                                fvertices[vt] = found->second;
                            }else{
                                // IF NOT found add the new vertex to the triangle and to the map of vertices
                                fvertices[vt]->coord.x = edge[case_table[cube[i][j][k]][m]][0];
                                fvertices[vt]->coord.y = edge[case_table[cube[i][j][k]][m]][1];
                                fvertices[vt]->coord.z = edge[case_table[cube[i][j][k]][m]][2];
                                // Store edges used to interpolate the current vertex
                                // Edge 0
                                fvertices[vt]->edge_0.i = edge[case_table[cube[i][j][k]][m]][3];
                                fvertices[vt]->edge_0.x = grid->X[fvertices[vt]->edge_0.i];
                                fvertices[vt]->edge_0.j = edge[case_table[cube[i][j][k]][m]][4];
                                fvertices[vt]->edge_0.y = grid->Y[fvertices[vt]->edge_0.j];
                                fvertices[vt]->edge_0.k = edge[case_table[cube[i][j][k]][m]][5];
                                fvertices[vt]->edge_0.z = grid->Z[fvertices[vt]->edge_0.k];
                                // Edge 1
                                fvertices[vt]->edge_1.i = edge[case_table[cube[i][j][k]][m]][6];
                                fvertices[vt]->edge_1.x = grid->X[fvertices[vt]->edge_1.i];
                                fvertices[vt]->edge_1.j = edge[case_table[cube[i][j][k]][m]][7];
                                fvertices[vt]->edge_1.y = grid->Y[fvertices[vt]->edge_1.j];
                                fvertices[vt]->edge_1.k = edge[case_table[cube[i][j][k]][m]][8];
                                fvertices[vt]->edge_1.z = grid->Z[fvertices[vt]->edge_1.k];
                                
                                // Map it
                                mapVertices.insert(make_pair(tie(fvertices[vt]->coord.x,
                                                                 fvertices[vt]->coord.y,
                                                                 fvertices[vt]->coord.z), fvertices[vt]));
                            }
                            
                            // --------------------------------------
                            // ### Vertex [2] of the current triangle
                            // --------------------------------------
                            vt=indexs[i][j][k]+m+1;
                            // Check if the current vertex was allocated to other triangle
                            found = mapVertices.find(tie(edge[case_table[cube[i][j][k]][m+1]][0],
                                                                                       edge[case_table[cube[i][j][k]][m+1]][1],
                                                                                       edge[case_table[cube[i][j][k]][m+1]][2]));
                            if (found != mapVertices.end()){
                                found->second->repeated++;
                                // IF found, set the current vertex of the triangle pointing to that one
                                fvertices[vt] = found->second;
                            }else{
                                // IF NOT found add the new vertex to the triangle and to the map of vertices
                                fvertices[vt]->coord.x = edge[case_table[cube[i][j][k]][m+1]][0];
                                fvertices[vt]->coord.y = edge[case_table[cube[i][j][k]][m+1]][1];
                                fvertices[vt]->coord.z = edge[case_table[cube[i][j][k]][m+1]][2];
                                // Store edges used to interpolate the current vertex
                                // Edge 0
                                fvertices[vt]->edge_0.i = edge[case_table[cube[i][j][k]][m+1]][3];
                                fvertices[vt]->edge_0.x = grid->X[fvertices[vt]->edge_0.i];
                                fvertices[vt]->edge_0.j = edge[case_table[cube[i][j][k]][m+1]][4];
                                fvertices[vt]->edge_0.y = grid->Y[fvertices[vt]->edge_0.j];
                                fvertices[vt]->edge_0.k = edge[case_table[cube[i][j][k]][m+1]][5];
                                fvertices[vt]->edge_0.z = grid->Z[fvertices[vt]->edge_0.k];
                                // Edge 1
                                fvertices[vt]->edge_1.i = edge[case_table[cube[i][j][k]][m+1]][6];
                                fvertices[vt]->edge_1.x = grid->X[fvertices[vt]->edge_1.i];
                                fvertices[vt]->edge_1.j = edge[case_table[cube[i][j][k]][m+1]][7];
                                fvertices[vt]->edge_1.y = grid->Y[fvertices[vt]->edge_1.j];
                                fvertices[vt]->edge_1.k = edge[case_table[cube[i][j][k]][m+1]][8];
                                fvertices[vt]->edge_1.z = grid->Z[fvertices[vt]->edge_1.k];
                                
                                mapVertices.insert(make_pair(tie(fvertices[vt]->coord.x,
                                                                 fvertices[vt]->coord.y,
                                                                 fvertices[vt]->coord.z), fvertices[vt]));
                            }
                            
                            // --------------------------------------
                            // ### Vertex [3] of the current triangle
                            // --------------------------------------
                            vt=indexs[i][j][k]+m+2;
                            // Check if the current vertex was allocated to other triangle
                            found = mapVertices.find(tie(edge[case_table[cube[i][j][k]][m+2]][0],
                                                         edge[case_table[cube[i][j][k]][m+2]][1],
                                                         edge[case_table[cube[i][j][k]][m+2]][2]));
                            if (found != mapVertices.end()){
                                found->second->repeated++;
                                // IF found, set the current vertex of the triangle pointing to that one
                                fvertices[vt] = found->second;
                            }else{
                                // IF NOT found add the new vertex to the triangle and to the map of vertices
                                fvertices[vt]->coord.x = edge[case_table[cube[i][j][k]][m+2]][0];
                                fvertices[vt]->coord.y = edge[case_table[cube[i][j][k]][m+2]][1];
                                fvertices[vt]->coord.z = edge[case_table[cube[i][j][k]][m+2]][2];
                                // Store edges used to interpolate the current vertex
                                // Edge 0
                                fvertices[vt]->edge_0.i = edge[case_table[cube[i][j][k]][m+2]][3];
                                fvertices[vt]->edge_0.x = grid->X[fvertices[vt]->edge_0.i];
                                fvertices[vt]->edge_0.j = edge[case_table[cube[i][j][k]][m+2]][4];
                                fvertices[vt]->edge_0.y = grid->Y[fvertices[vt]->edge_0.j];
                                fvertices[vt]->edge_0.k = edge[case_table[cube[i][j][k]][m+2]][5];
                                fvertices[vt]->edge_0.z = grid->Z[fvertices[vt]->edge_0.k];
                                // Edge 1
                                fvertices[vt]->edge_1.i = edge[case_table[cube[i][j][k]][m+2]][6];
                                fvertices[vt]->edge_1.x = grid->X[fvertices[vt]->edge_1.i];
                                fvertices[vt]->edge_1.j = edge[case_table[cube[i][j][k]][m+2]][7];
                                fvertices[vt]->edge_1.y = grid->Y[fvertices[vt]->edge_1.j];
                                fvertices[vt]->edge_1.k = edge[case_table[cube[i][j][k]][m+2]][8];
                                fvertices[vt]->edge_1.z = grid->Z[fvertices[vt]->edge_1.k];
                                
                                mapVertices.insert(make_pair(tie(fvertices[vt]->coord.x,
                                                                 fvertices[vt]->coord.y,
                                                                 fvertices[vt]->coord.z), fvertices[vt]));
                            }
                        }
                    }
                }
            }
        }
        if (output) printf("  End Mapping of Triangle Vertices\n");
        
        // Generate triangles of the blob
        genTriangles();
        // Compute normals
        normalization(viewer, false);
        
        // debug only
        /*
        for(int i=0; i < mblob.size(); i++){
            mTriangle * triangle = &mblob.at(i);
            // Get vertices of the current triangle object
            for(int j=0; j < 3; j++){
                mVertex * vertex = triangle->vertices[j];
                // Check how many times the current [vertex] is shared among the set of triangles
                cout << vertex << " Shared = " << vertex->repeated << "x times\n";
            }
        }
        */
        
        if (output){
            printf("  Number of Triangles:%d\n",nvertices*3);
            printf("  Number of Points:%d\n",nvertices);
        }
    
        // free internal resources
        free(cube);
        free(indexs);
        free(fvertices);
        mapVertices.clear();
    }
    
    /* Populate the a grid object */
    void Surface::startGrid(bool output){
        // Obtain the minimum and maximum values of the proteintein (i.e. in terms of atoms)
        grid->maxx = protein->max_x;
        grid->minx = protein->min_x;
        grid->maxy = protein->max_y;
        grid->miny = protein->min_y;
        grid->maxz = protein->max_z;
        grid->minz = protein->min_z;
        
        // Compute the max number of cells
        grid->size_x = ((grid->maxx-grid->minx)/grid->gridSpacing)+2*grid->seg;
        grid->size_y = ((grid->maxy-grid->miny)/grid->gridSpacing)+2*grid->seg;
        grid->size_z = ((grid->maxz-grid->minz)/grid->gridSpacing)+2*grid->seg;
        
        if (output){
            printf("  Grid spacing=%.2f, seg=%d, ioc=%d, decay=%.2f, iso=%.2f\n", grid->gridSpacing, grid->seg, grid->ioc, grid->d, grid->iso);
            printf("  Grid size x=%d, y=%d, z=%d \n", grid->size_x, grid->size_y, grid->size_z);
        }
        
        /* Attention, by using a malloc on c++ class and not the keyword New you don't call the constructor, 
           so there is no object. However, we only want to allocated some space for those variables inside a 
           class. That's why we use new std::vector<Point>() later because by using malloc the vector was never 
           properly initialized. */
        
        // Allocate space for grid values
        grid->X=(float *)malloc((grid->size_x+1)*sizeof(float));
        grid->Y=(float *)malloc((grid->size_y+1)*sizeof(float));
        grid->Z=(float *)malloc((grid->size_z+1)*sizeof(float));
        // Allocate auxiliar
        grid->map=(mPoint ***)malloc((grid->size_x+1)*sizeof(mPoint**));
        // Allocate auxiliar
        cube=(int ***)malloc((grid->size_x+1)*sizeof(int**));
        indexs=(int ***)malloc((grid->size_x+1)*sizeof(int**));
        // Continuing allocation for remaining...
        for(int i=0; i<grid->size_x+1;i++){
            grid->map[i]=(mPoint **)malloc((grid->size_y+1)*sizeof(mPoint*));
            cube[i]=(int **)malloc((grid->size_y+1)*sizeof(int*));
            indexs[i]=(int **)malloc((grid->size_y+1)*sizeof(int*));
            for(int j=0; j<grid->size_y+1; j++){
                grid->map[i][j]=(mPoint *)malloc((grid->size_z+1)*sizeof(mPoint));
                cube[i][j]=(int *)malloc((grid->size_z+1)*sizeof(int));
                indexs[i][j]=(int *)malloc((grid->size_z+1)*sizeof(int));
                for(int k=0; k<grid->size_z+1; k++){
                    grid->map[i][j][k].tPoints   = new std::set<Point *>();
                    grid->map[i][j][k].atoms     = new std::vector<Atom *>();
                    grid->map[i][j][k].intensity = 0.0;
                    grid->map[i][j][k].surface   = false;
                    cube[i][j][k]=0;
                    indexs[i][j][k]=0;
                }
            }
        }
        
        // Init grid values
        grid->X[0]=grid->minx-grid->padding;
        grid->Y[0]=grid->maxy+grid->padding;
        grid->Z[0]=grid->maxz+grid->padding;
        
        // Iterate each cell of the grid, add assign value (x coordinate)
        for(int i=1;i<=grid->size_x;i++) grid->X[i]=grid->X[i-1]+grid->gridSpacing;
        
        // Iterate each cell of the grid, add assign value (y coordinate)
        for(int i=1;i<=grid->size_y;i++) grid->Y[i]= grid->Y[i-1]-grid->gridSpacing;
        
        // Iterate each cell of the grid, add assign value (Z coordinate)
        for(int i=1;i<=grid->size_z;i++) grid->Z[i]= grid->Z[i-1]-grid->gridSpacing;

    }
    
    /* Gaussian Function (see documentation/paper for more details) */
    float Surface::Function(int i, float x, float y, float z, float radius){
        float R =(x-protein->atoms[i].coord.x)*(x-protein->atoms[i].coord.x) + (y-protein->atoms[i].coord.y)*(y-protein->atoms[i].coord.y) + (z-protein->atoms[i].coord.z)*(z-protein->atoms[i].coord.z);
        return expf(-grid->d*((R/(radius*radius))-1));
    }
    
    /* Gradient absolute sum value of the first Derivatives of the Gaussian Function defined in Function */
    Vector4 Surface::Derivate(float x, float y, float z){
        Vector4 dr;
        dr.x=0.0, dr.y=0.0, dr.z=0.0;
        float d=grid->d;
        for(int a=0; a < protein->atoms.size(); a++){
            dr.x+= -(((2*d*(x-protein->atoms[a].coord.x))*Function(a, x, y, z, protein->atoms[a].radius))/(protein->atoms[a].radius*protein->atoms[a].radius));
            
            dr.y+= -(((2*d*(y-protein->atoms[a].coord.y))*Function(a, x, y, z, protein->atoms[a].radius))/(protein->atoms[a].radius*protein->atoms[a].radius));
            
            dr.z+= -(((2*d*(z-protein->atoms[a].coord.z))*Function(a, x, y, z, protein->atoms[a].radius))/(protein->atoms[a].radius*protein->atoms[a].radius));
        }
        return(dr);
    }
    
    /* Vertex interpolation related methods */
    float Surface::vertex_interpolation_x(int x1, int y1, int z1, int x2, int y2, int z2){
        float iso   = grid->iso;
        float mu    = (iso - grid->map[x1][y1][z1].intensity) / (grid->map[x2][y2][z2].intensity - grid->map[x1][y1][z1].intensity);
        float point = grid->X[x1] + mu * (grid->X[x2]-grid->X[x1]);
        return point;
    }
    float Surface::vertex_interpolation_y(int x1, int y1, int z1, int x2, int y2, int z2){
        float iso   = grid->iso;
        float mu    = (iso - grid->map[x1][y1][z1].intensity) / (grid->map[x2][y2][z2].intensity - grid->map[x1][y1][z1].intensity);
        float point = grid->Y[y1] + mu * (grid->Y[y2]-grid->Y[y1]);
        return point;
    }
    float Surface::vertex_interpolation_z(int x1, int y1, int z1, int x2, int y2, int z2){
        float iso   = grid->iso;
        float mu    = (iso - grid->map[x1][y1][z1].intensity) / (grid->map[x2][y2][z2].intensity - grid->map[x1][y1][z1].intensity);
        float point = grid->Z[z1] + mu * (grid->Z[z2]-grid->Z[z1]);
        return point;
    }
    
    /* Compute the intensity on each point/vertex of the grid taking into account the set of atoms of the protein */
    void Surface::vertexIntensity(){
        int p,q,r;
        int xl,xu,yl,yu,zl,zu;
        int n     = (int) protein->atoms.size();
        float iso = grid->iso;
        float ioc = grid->ioc;
        // Iterate each atom of the protein
        for(int i=0;i<n;i++){
            p = ((int)(protein->atoms[i].coord.x-grid->X[0])/grid->gridSpacing)+1;
            q = ((int)(grid->Y[0]-protein->atoms[i].coord.y)/grid->gridSpacing)+1;
            r = ((int)(grid->Z[0]-protein->atoms[i].coord.z)/grid->gridSpacing)+1;
            if( p > ioc ) xl=p-ioc;
            else xl=0;
            if( p+ioc < grid->size_x+1 ) xu = p+ioc;
            else xu= grid->size_x+1;
            if( q > ioc ) yl=q-ioc;
            else yl=0;
            if( q+ioc < grid->size_y+1 ) yu = q+ioc;
            else yu= grid->size_y+1;
            if( r > ioc ) zl=r-ioc;
            else zl=0;
            if( r+ioc < grid->size_z+1 ) zu = r+ioc;
            else zu= grid->size_z+1;
            
            for(int x=xl;x<xu;x++){
                for(int y=yl;y<yu;y++){
                    for(int z=zl;z<zu;z++){
                        grid->map[x][y][z].intensity+=Function(i,grid->X[x],grid->Y[y],grid->Z[z], protein->atoms[i].radius);
                        Vector4 tmp(grid->X[x], grid->Y[y], grid->Z[z], x, y, z);
                        if (protein->atoms[i].pointInside(&tmp)) grid->map[x][y][z].atoms->push_back(&protein->atoms[i]);
                        
                    }
                }
            }
            
        }
        
        // Compute intensities in the vertices of the cubes
        for(int i=0;i<grid->size_x-1;i++){
            for(int j=0;j<grid->size_y-1;j++){
                for(int k=0;k<grid->size_z-1;k++){
                    // 8 vertices of the cube
                    if(grid->map[i][j+1][k+1].intensity  <iso)	{cube[i][j][k] |=1;}//(1<<0);
                    if(grid->map[i+1][j+1][k+1].intensity<iso)  {cube[i][j][k] |=2;}//(1<<1);
                    if(grid->map[i+1][j+1][k].intensity  <iso)	{cube[i][j][k] |=4;}//(1<<2);
                    if(grid->map[i][j+1][k].intensity    <iso)  {cube[i][j][k] |=8;}//(1<<3);
                    if(grid->map[i][j][k+1].intensity    <iso)  {cube[i][j][k] |=16;}//(1<<4);
                    if(grid->map[i+1][j][k+1].intensity  <iso)	{cube[i][j][k] |=32;}//(1<<5);
                    if(grid->map[i+1][j][k].intensity    <iso)  {cube[i][j][k] |=64;}//(1<<6);
                    if(grid->map[i][j][k].intensity      <iso) 	{cube[i][j][k] |=128;}//(1<<7);
                }
            }
        }
        
        //Number of vertices that exist on the surface
        int i,j,k;
        for(i=0;i<grid->size_x;i++){
            for(j=0;j<grid->size_y;j++){
                for(k=0;k<grid->size_z;k++){
                    indexs[i][j][k]=0;
                    if(vertex_table[cube[i][j][k]]>0){
                        indexs[i][j][k]=nvertices;
                        nvertices=nvertices+vertex_table[cube[i][j][k]];
                    }
                }
            }
        }
        // Creates the array of vertices by using the size of nvertices
        fvertices=(mVertex **)malloc(nvertices*sizeof(mVertex));
        for(int j=0; j<nvertices; j++) fvertices[j] = new mVertex();
    }
    
    /* Compute Vector direction */
    Vector4 Surface::getvDirection(mVertex * vertex, Vector4 * edge_a, Vector4 * edge_b){
        bool debug=false;
        // ## Get intensity of each edge node, and find out the max and min node
        float edge_ai = grid->map[edge_a->i][edge_a->j][edge_a->k].intensity;
        float edge_bi = grid->map[edge_b->i][edge_b->j][edge_b->k].intensity;
        Vector4 * max, * min;
        if (edge_ai > edge_bi){ max = edge_a; min = edge_b; }else{ max = edge_b; min = edge_a; }
        
        if (debug){
            printf("Vertex                        = (%.1f, %.1f, %.1f)\n", vertex->coord.x, vertex->coord.y, vertex->coord.z);
            printf("Edge 0 (%d, %d, %d) intensity = %f\n", edge_a->i, edge_a->j, edge_a->k, edge_ai);
            printf("Edge 1 (%d, %d, %d) intensity = %f\n", edge_b->i, edge_b->j, edge_b->k, edge_bi);
            printf("Max                           = (%d, %d, %d)(%.1f, %.1f, %.1f)\n", max->i, max->j, max->k, max->x, max->y, max->z);
            printf("Min                           = (%d, %d, %d)(%.1f, %.1f, %.1f)\n", min->i, min->j, min->k, min->x, min->y, min->z);
        }
        
        // ## Create vector between max and min, why not vertex and min ?
        Vector4 direction = min->Vector(max);
        direction.Length();
        if (debug){
            printf("Vector direction              = (%1.f, %1.f, %1.f)\n", direction.x, direction.y, direction.z);
            printf("Vector direction magnitude    =  %f\n", direction.length);
            printf("\n\n");
        }
        
        return(direction);
    }
    
    /* Compute normal of triangle vertex using the cross product */
    Vector4 Surface::computeNormal(Vector4 * direction, mVertex * _A, mVertex * _B, mVertex * _C){
        bool debug=false;
        // ## Compute cross product with random order of vertices
        Vector4 * A = &_A->coord;
        Vector4 * B = &_B->coord;
        Vector4 * C = &_C->coord;
        
        Vector4 AB = *A - *B;
        Vector4 CB = *C - *B;
        Vector4 normal = CB % AB;
        normal.Length();
        
        // ## Check if the angle between vector [direction] and [normal] is greater than 90º, if yes then invert [normal]
        float angle = ((*direction) * normal) / (direction->length * normal.length);
        if (debug) printf("Angle = %f\n", angle);
        if (angle >= -0.01745240643){ // cos 91 = 0.017
            // Inverted
            //printf("inverted\n");
            normal.x = -normal.x;
            normal.y = -normal.y;
            normal.z = -normal.z;
        }
        return(normal);
    }
    
    /* Add triangles genereated by the marching cubes algorithm into the blob */
    void Surface::genTriangles(){
        // ----------------------------------------
        // ## Create triangles and add it to the marching cubes blob
        // ----------------------------------------
        for(int i=0; i < nvertices; i+=3){
            if (i==0) continue; // Avoid the first dummy inicialized value in vector vertices that is set to V=[0,0,0]
            
            // Get mapped vertices
            mVertex * A = fvertices[i];
            mVertex * B = fvertices[i+1];
            mVertex * C = fvertices[i+2];
            
            // Avoid
            if (A->coord.x == 0 && A->coord.y == 0 && A->coord.z == 0){continue;}
            if (B->coord.x == 0 && B->coord.y == 0 && B->coord.z == 0){continue;}
            if (C->coord.x == 0 && C->coord.y == 0 && C->coord.z == 0){continue;}
            
            /* debug only
            printf("(%.1f, %.1f, %.1f) interpolated from [%d][%d][%d]=(%.1f, %.1f, %.1f) AND [%d][%d][%d]=(%.1f, %.1f, %.1f)\n",
                   A->coord.x, A->coord.y, A->coord.z,
                   A->edge_0.i, A->edge_0.j, A->edge_0.k,
                   A->edge_0.x, A->edge_0.y, A->edge_0.z,
                   A->edge_1.i, A->edge_1.j, A->edge_1.k,
                   A->edge_1.x, A->edge_1.y, A->edge_1.z);
            */
            mTriangle triangle;
            // Add vertices to triangle
            triangle.vertices.push_back(A);
            triangle.vertices.push_back(B);
            triangle.vertices.push_back(C);
            
            // Add triangle to the marching cubes blob
            mblob.push_back(triangle);
        }
        
        
    }
    
    /* Compute the normal vector associated with each vertex of the blob */
    void Surface::normalization(bool viewer, bool debug){
        float threshold=0.001;
        bool useDerivatives=false; // Compute normal vector of each vertex of the surface using derivates, This is much slower!
        
        // ----------------------------------------
        // ## Compute normals using the cross product (smooth shading)
        // ----------------------------------------
        for(int i=0; i < mblob.size() && !useDerivatives; i++){
            mTriangle * triangle = &mblob.at(i);
            // Get vertices of the current triangle object
            mVertex * A = triangle->vertices.at(0);
            mVertex * B = triangle->vertices.at(1);
            mVertex * C = triangle->vertices.at(2);
            
            // Get edges of vertices A, B, and C
            Vector4 * edge0_A = &A->edge_0;
            Vector4 * edge1_A = &A->edge_1;
            
            Vector4 * edge0_B = &B->edge_0;
            Vector4 * edge1_B = &B->edge_1;
            
            Vector4 * edge0_C = &C->edge_0;
            Vector4 * edge1_C = &C->edge_1;
            
            /* debug only
             printf("(%.1f, %.1f, %.1f) interpolated from [%d][%d][%d]=(%.1f, %.1f, %.1f) AND [%d][%d][%d]=(%.1f, %.1f, %.1f)\n",
             A->coord.x, A->coord.y, A->coord.z,
             A->edge_0.i, A->edge_0.j, A->edge_0.k,
             A->edge_0.x, A->edge_0.y, A->edge_0.z,
             A->edge_1.i, A->edge_1.j, A->edge_1.k,
             A->edge_1.x, A->edge_1.y, A->edge_1.z);
             */
            
            // Compute Vector Direction, later need to validate the normal of vertex A, B, and C.
            Vector4 A_direction = getvDirection(A, edge0_A, edge1_A);
            Vector4 B_direction = getvDirection(B, edge0_B, edge1_B);
            Vector4 C_direction = getvDirection(C, edge0_C, edge1_C);
            
            
            // Compute vertex [A, B, C] normals and validated it against the vector direction found for vertex [A, B, C]
            // - Invert or not to invert that is the question
            Vector4 normal = computeNormal(&A_direction, C, B, A);
            A->normal.x += normal.x;
            A->normal.y += normal.y;
            A->normal.z += normal.z;
            
            normal = computeNormal(&B_direction, C, B, A);
            B->normal.x += normal.x;
            B->normal.y += normal.y;
            B->normal.z += normal.z;
            
            normal = computeNormal(&C_direction, C, B, A);
            C->normal.x += normal.x;
            C->normal.y += normal.y;
            C->normal.z += normal.z;
            
        }
        
        // Finally, divide the normal vector by the number of times the same vertex is mapped to a triangle
        for(int i=0; i < mblob.size() && !useDerivatives; i++){
            mTriangle * triangle = &mblob.at(i);
            for(int j=0; j < 3; j++){
                mVertex * vertex = triangle->vertices.at(j);
                vertex->normal.x = vertex->normal.x / vertex->repeated;
                vertex->normal.y = vertex->normal.y / vertex->repeated;
                vertex->normal.z = vertex->normal.z / vertex->repeated;
            }
        }
        
        
        // ----------------------------------------
        // ## Compute normals using derivatives (alternative)
        // ----------------------------------------
        for(int i=0; i < mblob.size() && useDerivatives; i++){
            mTriangle * triangle = &mblob.at(i);
            for(int j=0; j < 3; j++){
                mVertex * vertex = triangle->vertices.at(j);
                Vector4 dr = Derivate(vertex->coord.x, vertex->coord.y, vertex->coord.z);
                vertex->normal.x = -dr.x;
                vertex->normal.y = -dr.y;
                vertex->normal.z = -dr.z;
            }
        }

        // Populate CavVis structures
        int id=0;
        for(int i=0; i < grid->size_x+1; i++){
            for(int j=0; j < grid->size_y+1; j++){
                for(int k=0; k < grid->size_z+1; k++){
                    if (!grid->map[i][j][k].surface) continue;
                    for(set<Point *>::iterator it = grid->map[i][j][k].tPoints->begin(); it != grid->map[i][j][k].tPoints->end(); it++){
                        Point * vertex = (*it);
                        // Find this vertex on the map structure
                        map<vid, mVertex *>::iterator found = mapVertices.find(tie(vertex->coord.x, vertex->coord.y, vertex->coord.z));
                        if (found != mapVertices.end()){
                            //printf("FOUND IT\n");
                            vertex->id         = id;
                            vertex->normal.id  = id; // the [id] is the same of the vertex
                            vertex->normal.x   = found->second->normal.x;
                            vertex->normal.y   = found->second->normal.y;
                            vertex->normal.z   = found->second->normal.z;
                            vertex->normal.length = sqrt( powf(vertex->normal.x, 2.0) +
                                                          powf(vertex->normal.y, 2.0) +
                                                          powf(vertex->normal.z, 2.0));
                            //printf("--->(%.1f, %.1f, %.1f)-(%.1f, %.1f, %.1f)\n", vertex->coord.x, vertex->coord.y, vertex->coord.z, vertex->normal.x, vertex->normal.y, vertex->normal.z);
                       
                            // # Filter domain
                            // Check if the normal of the point is pointing towards the surface (valid)
                            // or to the outside of the surface (invalid)
                            // ! The main optimization idea is to reduce the number of points that will be processed by the algorithm
                            if (!Algorithm::validPoint(grid, vertex, threshold)) vertex->invalid=true;
                            id++; // increment point ids
                        }
                    }
                }
            }
        }
        
    }
    

}
