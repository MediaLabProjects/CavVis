//---------------------------------------------------------------------------------------
//
//	CAVVIS - A Field-of-View Geometric Algorithm for Protein Cavity Detection
//
//  Custom Neighbors Library
//
//  Copyright (C) 2018 Instituto de Telecomunicações (www.it.pt)
//  Copyright (C) 2018 Universidade da Beira Interior (www.ubi.pt)
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
#include "neighbors.h"
namespace Neighbors{
    
    /* Check if [nodes] belongs to the grid */
    bool isOutOfBounds(Grid * grid, Point * node){
        // Check if [v] is a valid vertice on the parse grid
        if (node->coord.i <= grid->size_x && node->coord.j <= grid->size_y && node->coord.k <= grid->size_z) return false;
        return true;
    }
    
    /* Check if [node] was used, by the marching cubes algorithm, to interpolate triangle vertices */
    bool isSurface(Grid * grid, Point * node){
        if (grid->map[node->coord.i][node->coord.j][node->coord.k].surface) return true;
        return false;
    }
    
    /* Check if [node] is valid through a set of tests */
    bool isValid(Point * a, Point * node, Grid * grid){
        // Teste 1
        // Check if [v] is a valid vertice on the parse grid
        if (isOutOfBounds(grid, node)) return false;
        
        // Teste 2
        // Check if its a surface vertice
        if (!isSurface(grid, node)) return false;
        
        return (true);
    }
    
    /* Check if the vertice to add to the list of neighbors already exists on the vector and is not out of bounds of the grid */
    bool existsValidVertice(vector<Point> * vertices, Point * v, Grid * grid){
        // Check if the vertice was already stored
        for(int i=0; i < vertices->size(); i++){
            Point * p = &vertices->at(i);
            if (p->coord.x == v->coord.x && p->coord.y == v->coord.y && p->coord.z == v->coord.z) return true;
        }
        // Check if [v] is a valid vertice on the parse grid
        if (v->coord.i <= grid->size_x && v->coord.j <= grid->size_y && v->coord.k <= grid->size_z){
            if (v->coord.i > 0 && v->coord.j > 0 && v->coord.k > 0){
                return false;
            }else return true;
        }else
            return true;
    }
    
    /* Auxiliar Method */
    void getNeighborsOfVertice2dAt(vector<Point> * neighbors, Grid * grid, Point currentVertice, int at, bool change_k){
        
        // using this for we can save/show all the voxels until [at]
        //for(int c=1; c <= at; c++){
        for(;true;){
            
            int c=at;
            Voxel voxel;
            Point initialVertice;
            if (change_k){
                
                // Voxel 0
                // - starting/initial voxel z+ z-
                initialVertice.coord.i = currentVertice.coord.i;
                initialVertice.coord.j = currentVertice.coord.j;
                initialVertice.coord.k = currentVertice.coord.k;
                initialVertice.coord.x = grid->X[currentVertice.coord.i];
                initialVertice.coord.y = grid->Y[currentVertice.coord.j];
                initialVertice.coord.z = grid->Z[currentVertice.coord.k];
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                
                // Remaining voxels at !=1 near voxel 0
                for(int e=1; e <= at && c!=1 && change_k; e++){
                    initialVertice.coord.i = currentVertice.coord.i-(c-e);
                    initialVertice.coord.x = grid->X[initialVertice.coord.i];
                    // Add current voxel to the pool of voxels
                    if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                }
                
                for(int e=1; e <= at && c!=1 && change_k; e++){
                    initialVertice.coord.i = currentVertice.coord.i+(c-e);
                    initialVertice.coord.x = grid->X[initialVertice.coord.i];
                    // Add current voxel to the pool of voxels
                    if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                }
                
                
            }
            
            
            // Voxel 1
            // save index position of grid
            initialVertice.coord.i = currentVertice.coord.i-c;
            initialVertice.coord.j = currentVertice.coord.j+c;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            // Remaining voxels at !=1 near voxel 1
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i-(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            }
            
            // Remaining voxels at !=1 near voxel 1
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i-c;
                initialVertice.coord.j = currentVertice.coord.j+(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                initialVertice.coord.y = grid->Y[initialVertice.coord.j];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                
                for(int e=1; e <= at && c!=1 && change_k; e++){
                    initialVertice.coord.i = currentVertice.coord.i-(c-e);
                    initialVertice.coord.x = grid->X[initialVertice.coord.i];
                    // Add current voxel to the pool of voxels
                    if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                }
                
            }
            
            // Voxel 2
            initialVertice.coord.i = currentVertice.coord.i-c;
            initialVertice.coord.j = currentVertice.coord.j;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            // Remaining voxels at !=1 near voxel 2
            for(int e=1; e <= at && c!=1 && change_k; e++){
                initialVertice.coord.i = currentVertice.coord.i-(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            }
            
            
            // Voxel 3
            initialVertice.coord.i = currentVertice.coord.i-c;
            initialVertice.coord.j = currentVertice.coord.j-c;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            // Remaining voxels at !=1 near voxel 3
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i-(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            }
            
            // Remaining voxels at !=1 near voxel 3
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i-c;
                initialVertice.coord.j = currentVertice.coord.j-(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                initialVertice.coord.y = grid->Y[initialVertice.coord.j];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                
                for(int e=1; e <= at && c!=1 && change_k; e++){
                    initialVertice.coord.i = currentVertice.coord.i-(c-e);
                    initialVertice.coord.x = grid->X[initialVertice.coord.i];
                    // Add current voxel to the pool of voxels
                    if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                }
                
            }
            
            // Voxel 4
            initialVertice.coord.i = currentVertice.coord.i;
            initialVertice.coord.j = currentVertice.coord.j-c;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            
            // Voxel 5
            initialVertice.coord.i = currentVertice.coord.i+c;
            initialVertice.coord.j = currentVertice.coord.j-c;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            // Remaining voxels at !=1 near voxel 5
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i+(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            }
            // Remaining voxels at !=1 near voxel 5
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i+c;
                initialVertice.coord.j = currentVertice.coord.j-(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                initialVertice.coord.y = grid->Y[initialVertice.coord.j];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                
                for(int e=1; e <= at && c!=1 && change_k; e++){
                    initialVertice.coord.i = currentVertice.coord.i+(c-e);
                    initialVertice.coord.x = grid->X[initialVertice.coord.i];
                    // Add current voxel to the pool of voxels
                    if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                }
                
            }
            
            // Voxel 6
            initialVertice.coord.i = currentVertice.coord.i+c;
            initialVertice.coord.j = currentVertice.coord.j;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            // Remaining voxels at !=1 near voxel 6
            for(int e=1; e <= at && c!=1 && change_k; e++){
                initialVertice.coord.i = currentVertice.coord.i+(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                
            }
            
            // Voxel 7
            initialVertice.coord.i = currentVertice.coord.i+c;
            initialVertice.coord.j = currentVertice.coord.j+c;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            // Remaining voxels at !=1 near voxel 6
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i+(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            }
            
            for(int e=1; e <= at && c!=1; e++){
                initialVertice.coord.i = currentVertice.coord.i+c;
                initialVertice.coord.j = currentVertice.coord.j+(c-e);
                initialVertice.coord.x = grid->X[initialVertice.coord.i];
                initialVertice.coord.y = grid->Y[initialVertice.coord.j];
                // Add current voxel to the pool of voxels
                if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                
                
                for(int e=1; e <= at && c!=1 && change_k; e++){
                    initialVertice.coord.i = currentVertice.coord.i+(c-e);
                    initialVertice.coord.x = grid->X[initialVertice.coord.i];
                    // Add current voxel to the pool of voxels
                    if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
                }
                
            }
            
            // Voxel 8
            initialVertice.coord.i = currentVertice.coord.i;
            initialVertice.coord.j = currentVertice.coord.j+c;
            initialVertice.coord.k = currentVertice.coord.k;
            initialVertice.coord.x = grid->X[initialVertice.coord.i];
            initialVertice.coord.y = grid->Y[initialVertice.coord.j];
            initialVertice.coord.z = grid->Z[initialVertice.coord.k];
            // Add current voxel to the pool of voxels
            if (!existsValidVertice(neighbors, &initialVertice, grid)) neighbors->push_back(initialVertice);
            
            break;
        }
        
        //return(voxels);
    }
    
    /* Standard main method of the Neighbors Library */
    vector<Point> getNeighborsOfVerticeAt(Grid * grid, Point currentVertice, int at){
        bool debug=false;
        int initial_k=currentVertice.coord.k;
        int k=currentVertice.coord.k-at;
        int km=currentVertice.coord.k+at;
        vector<Point> neighbors;
        if (debug) cout << "INITIAL currentVoxel.k=" << currentVertice.coord.k << " and k=" << k << "\n";
        
        // current
        getNeighborsOfVertice2dAt(&neighbors, grid, currentVertice, at, false);
        
        // front
        for(int i=k; i != initial_k;i++){
            currentVertice.coord.k = i;
            if (debug)cout << "currentVoxel2.k=" << currentVertice.coord.k << "\n";
            bool flag=false;
            if (k==i) flag=true;
            getNeighborsOfVertice2dAt(&neighbors, grid, currentVertice, at, flag);
            
        }
        
        // behind
        for(int i=km; i != initial_k; i--){
            currentVertice.coord.k = i;
            if (debug)cout << "currentVoxel3.k=" << currentVertice.coord.k << "\n";
            bool flag=false;
            if(i==km) flag=true;
            getNeighborsOfVertice2dAt(&neighbors, grid, currentVertice, at, flag);
            
        }
        
        return (neighbors);
        
    }
    
}