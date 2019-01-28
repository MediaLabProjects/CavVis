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
#include "CavVis.h"
#define QUICKHULL_IMPLEMENTATION
#include "quickhull.h"
#define VOXELIZER_IMPLEMENTATION
#include "voxelizer.h"
namespace CavVis{
    
    Algorithm::Algorithm(vector<Point *> * points, Protein * protein, float gridSpacing, int top, float distance, bool areas, bool ffilter, bool viewer, bool output){
        // If the are of each cavity is requested beware that at the current stage, more computational time is needed to perform this step.
    
        printf("- Marching Cubes Data Computation...\n");
        // Assign parameters
        Grid grid; // Main protein grid
        grid.gridSpacing = gridSpacing;
        grid.iso = 1.0;
        grid.C   = 0.33; // Smooth parameter
        grid.d   = 2.35; // Decay rate parameter of the gaussian surface
        Surface MC;
        MC.MarchingCubes(viewer, protein, &grid, true);
        
        // #### Cavity Detection
        printf("\n- Visibility tests...\n");
        Visibility(protein, &grid, points);
        
        // #### Cavity Formation
        printf("\n- Creating groups of visible points...\n");
        vector<Group> visibles;
        Union(points, &visibles, false, false);
        printf("\n- Forming cavities from detected groups...\n");
        vector<Cavity> cavities;
        formCavities(grid.gridSpacing, &visibles, &cavities, false, output);
        printf("\n- Generating dummy atoms that fills each top-%d cavity...\n", top);
        IdentifyFillingAtoms(&grid, &cavities, gridSpacing, ffilter, false, output);
        
        // Create Generic content for the result file
        for(int i=0; i < cavities.size(); i++){
            Cavity * cav = &cavities[i];
            if (!cav->top) continue;
            Cavity nc;
            nc.id = i;
            for(int j=0; j < cav->chs.size(); j++){
                ConvexHull * ch = &cav->chs[j];
                nc.volume = ch->volume;
                nc.area   = ch->area;
                for(int k=0; k < ch->spheres.size(); k++){
                    Atom dummyAtom;
                    dummyAtom.id = (int) nc.atoms.size();
                    dummyAtom.coord.x = ch->spheres[k].coord.x;
                    dummyAtom.coord.y = ch->spheres[k].coord.y;
                    dummyAtom.coord.z = ch->spheres[k].coord.z;
                    // avoid rep.
                    bool exists=false;
                    for(int e=0; e < nc.atoms.size(); e++){
                        if (nc.atoms[e].coord.x == dummyAtom.coord.x &&
                            nc.atoms[e].coord.y == dummyAtom.coord.y &&
                            nc.atoms[e].coord.z == dummyAtom.coord.z){ exists=true; break;}
                    }
                    if (!exists) nc.atoms.push_back(dummyAtom);
                }
            }
            
            nc.gc.computeGeometricCenter(&nc.atoms);
            if (nc.atoms.size() != 0) protein->temp.push_back(nc);
        }
        
        cout << "\n- Clustering final list of cavities (d=" << distance << ")...\n\n";
        // Final cluster of cavities and output only those [top] cavities basead on volume
        cluster(protein, &protein->temp, top, distance, false, false);
        
        // ## If requested compute area of each cavity
        if (areas) computeAreas(protein, grid.gridSpacing);
        
        cout << "- Number of putative binding sites detected: " << protein->cavities.size() << "\n\n";
        for(int i=0; i < protein->cavities.size(); i++){
            if (areas) printf("  Protein Cavity [%d] volume = %f and area = %f\n",i,protein->cavities[i].volume, protein->cavities[i].area);
            else printf("  Protein Cavity [%d] volume = %f\n", i, protein->cavities[i].volume);
        }
        
        printf("\n");
        
        // Free resources
        free(grid.X); free(grid.Y); free(grid.Z); free(grid.map);
    }
    
    /* Compute cavity areas */
    void Algorithm::computeAreas(Protein * protein, float cgridspacing){
        // Compute area of each cavity, using the marching cubes to construct the surface triangles of the fillings atoms of the cavity
        // -------------------------------------
        bool debug_area     = false;
        bool output_mc      = false;    // Output marching cubes information
        float radius        = 1.0;      // Filling/dummy atom radius
        vector<vector<mTriangle>> areas; // Each cavity areas
        vector<Protein> acavities;

        float minx=FLT_MAX,miny=FLT_MAX,minz=FLT_MAX;
        float maxx=FLT_MIN,maxy=FLT_MIN,maxz=FLT_MIN;

        for(int i=0; i < protein->cavities.size(); i++){
            Protein acavity; // store cavity data as a protein object to the use the marching cubes
            acavity.pcavities.push_back(&protein->cavities[i]); // store pointer to the current cavity
            
            for(int j=0; j < protein->cavities[i].atoms.size(); j++){
                Atom a = protein->cavities[i].atoms[j];
                a.id = j;
                a.radius = radius;
                minx=a.coord.x<minx?a.coord.x:minx;
                maxx=a.coord.x>maxx?a.coord.x:maxx;
                miny=a.coord.y<miny?a.coord.y:miny;
                maxy=a.coord.y>maxy?a.coord.y:maxy;
                minz=a.coord.z<minz?a.coord.z:minz;
                maxz=a.coord.z>maxz?a.coord.z:maxz;
                acavity.atoms.push_back(a);
                if (debug_area) printf("H   %f  %f  %f\n", a.coord.x, a.coord.y, a.coord.z);
            }
            // Store min, max values
            acavity.min_x = minx;
            acavity.max_x = maxx;
            acavity.min_y = miny;
            acavity.max_y = maxy;
            acavity.min_z = minz;
            acavity.max_z = maxz;
            acavity.max_radius = 2.0;
            acavity.numberOfatoms = (int) acavity.atoms.size();
            // Store cavity
            acavities.push_back(acavity);
            // reset
            minx=FLT_MAX;miny=FLT_MAX;minz=FLT_MAX;
            maxx=FLT_MIN;maxy=FLT_MIN;maxz=FLT_MIN;
            if (debug_area){
                printf(" - Number of atoms = %d\n", (int) acavity.numberOfatoms);
                printf(" - min_x=%f, max_x=%f\n", acavity.min_x, acavity.max_x);
                printf(" - min_y=%f, max_y=%f\n", acavity.min_y, acavity.max_y);
                printf(" - min_z=%f, max_z=%f\n", acavity.min_z, acavity.max_z);
                printf("\n");
            }
        }

        // Generate surface for each cavity previously populated and compute its area by its triangles
        for(int i=0; i < acavities.size(); i++){
            Protein * acavity = &acavities[i];
            Cavity * cavity = acavity->pcavities[0];
            vector<mTriangle> cblob;
            Grid cgrid; // Main protein grid
            cgrid.gridSpacing = cgridspacing;
            cgrid.iso = 1.0;
            cgrid.C   = 0.33; // Smooth parameter
            cgrid.d   = 2.35; // Decay rate parameter of the gaussian surface
            Surface MC;
            MC.MarchingCubes(true, acavity, &cgrid, output_mc);
            MC.mblob[0].color = acavity->color;
            
            // Iterate triangles of the created surface and compute AREA of the current cavity surface
            vector<mTriangle> triangles = MC.mblob;
            for(int j=0; j < triangles.size(); j++){
                mTriangle * triangle = &triangles.at(j);
                mVertex * v1 = triangle->vertices.at(0);
                mVertex * v2 = triangle->vertices.at(1);
                mVertex * v3 = triangle->vertices.at(2);
                // https://math.stackexchange.com/questions/128991/how-to-calculate-area-of-3d-triangle
                float aux = 0.0; // Area of the current triangle at coordinates (x,y,z) for each vertex
                aux = powf(((v2->coord.x * v1->coord.y) - (v3->coord.x * v1->coord.y) - (v1->coord.x * v2->coord.y) +
                            (v3->coord.x * v2->coord.y) + (v1->coord.x * v3->coord.y) - (v2->coord.x * v3->coord.y)), 2.0);
                aux += powf(((v2->coord.x * v1->coord.z) - (v3->coord.x * v1->coord.z) - (v1->coord.x * v2->coord.z) +
                             (v3->coord.x * v2->coord.z) + (v1->coord.x * v3->coord.z) - (v2->coord.x * v3->coord.z)), 2.0);
                aux += powf(((v2->coord.y * v1->coord.z) - (v3->coord.y * v1->coord.z) - (v1->coord.y * v2->coord.z) +
                             (v3->coord.y * v2->coord.z) + (v1->coord.y * v3->coord.z) - (v2->coord.y * v3->coord.z)),2.0);
                cavity->area = cavity->area + (0.5*sqrtf(aux));
            }
        }
    }
    
    /* Auxiliar function of the visibility method */
    bool Algorithm::Visibility_tests(Protein * protein, Grid * grid, Point * a, bool debug){
        debug=false;
        
        if (debug) printf(" Running tests for point a[%d]:\n", a->id);
        
        // -------------------------------------------------------------
        // 1 - Get 1st intersected voxel by [a] and variable FoV per surface vertex
        // -------------------------------------------------------------
        Point ivoxel = voxelRayTracer_alt(a, grid); // First intersected voxel vertice by [a]
        
        // -------------------------------------------------------------
        // 2 - Computing the variable FoV using the water molecule radius
        // -------------------------------------------------------------
        if (variable_fov){
            Vector4 PR = a->coord.Vector(&ivoxel.coord);
            PR.Length();
            a->fov_angle = atan( 1.4 / PR.length);
            // printf("a->fov_angle = %f\n", a->fov_angle);
            // Although we have a variable FOV, we need to guarantee a minimum FOV of 20
            // in order to avoid problems in the detection of cavity regions
            // That is if fov angle < that 20 degrees then set it to 20 degrees
            if (a->fov_angle < 0.93969262) a->fov_angle = 0.93969262;
        }else a->fov_angle = 0.93969262;
        
        //printf("a->fov_angle = %f\n", a->fov_angle);
        
        // -------------------------------------------------------------
        // 3 - Get the group of neighbor voxels, of the first intersected voxel, that inscribes the water molecule
        // -------------------------------------------------------------
        if (ivoxel.coord.i == -1) return false;         // None voxel was intersected
        int nsize = ((int) (1.4/grid->gridSpacing))-1;  // Get the number of neighbor-levels necessary to inscribe the water molecule
                                                        // Set to one if better time performance is needed
        //printf("nsize=%d\n", nsize);
        vector<Point> voxel = Neighbors::getNeighborsOfVerticeAt(grid, ivoxel, nsize);
        
        // -------------------------------------------------------------
        // 4 - Get Vertices Inside FoV and back-face culling valid
        // -------------------------------------------------------------
        for(int j=0; j < voxel.size(); j++){
            Point * tmp = &voxel.at(j);
            if (!grid->map[tmp->coord.i][tmp->coord.j][tmp->coord.k].surface) continue;
    
            // Get those triangle vertex interpolated using the current voxel node/point [tmp]
            for(set<Point *>::iterator it = grid->map[tmp->coord.i][tmp->coord.j][tmp->coord.k].tPoints->begin(); it != grid->map[tmp->coord.i][tmp->coord.j][tmp->coord.k].tPoints->end(); it++){
                Point * b = (*it);
                if (b->invalid) continue;
                
                // [FOV Test] -> Is [b] in the field of view of [a] ?
                if (!test3(a, b)) continue;
                
                // [Back-Face Culling Test] -> Is point [a] roughly facing point [b] ?
                if (Algorithm::test2(&(a->normal), &(b->normal))){
                    // Check if already exists (need to change this to set, in order to avoid the search below)
                    bool exists=false;
                    for(int e=0; e < a->final.size(); e++){
                        if (a->final[e]->coord.x == b->coord.x &&
                            a->final[e]->coord.y == b->coord.y &&
                            a->final[e]->coord.z == b->coord.z){exists=true; break;}
                    }
                    // if not exists, added it to the final list of visible points of [a]
                    if (!exists) a->final.push_back(b);
                }
            }
        }
        
        if (a->final.size() != 0) return true; else return false;
    }
    
    /* Main method */
    void Algorithm::Visibility(Protein * protein, Grid * grid, vector<Point *> * points){
        set<Point *> tpoints;
        for(int i=0; i < grid->size_x+1; i++){
            for(int j=0; j < grid->size_y; j++){
                for(int k=0; k < grid->size_z; k++){
                    // Avoid those voxels that were not marked as surface one (i.e. with triangles)
                    if (!grid->map[i][j][k].surface) continue;
                    
                    for(set<Point *>::iterator it = grid->map[i][j][k].tPoints->begin(); it != grid->map[i][j][k].tPoints->end(); it++){
                        Point * tVertex = (*it);
                       
                        // Avoid those voxels that are not valid
                        if (tVertex->invalid) continue;
                        // Avoid anomalous situations
                        if (tVertex->coord.x == 0 && tVertex->coord.y == 0 && tVertex->coord.z == 0) continue;
                        
                        // Check if it was successfully added to the temporary set (i.e. not repeated)
                        if (tpoints.insert(tVertex).second){
                            // Run tests and store only those with relevant results
                            if (Visibility_tests(protein, grid, tVertex, false)) points->push_back(tVertex);
                        }
                    }
                }
            }
        }
        printf("  Number of points processed: %d\n", (int) points->size());
        tpoints.clear();
    }

    /* Test if point [b] is over the plane of point [a] - Legacy Support */
    bool Algorithm::test1(Point * a, Point * b, bool debug){
        // --> Ax + By + Cz + D
        Vector4 * N = &(a->normal);  // Get the normal vector [N] of point [a]
        Point * T = a;            // Get point of the plan where the normal vector [N] is
        Point * toTest = b;       // Check if point [b] is on the plane defined by [N]
        
        // Compute [D]
        // D=-(Ax+Bx+Cz+D)
        float D=-1*( N->x*T->coord.x + N->y*T->coord.y + N->z*T->coord.z);
        // Compute Ax+Bx+Cz+D
        float res= N->x*toTest->coord.x + N->y*toTest->coord.y + N->z*toTest->coord.z + D;
        if (debug) printf("%f*x+%f*y+%f*z+%f=%f\n", N->x, N->y, N->z, D, res);
        // If [res]>0 then point [b] is over the plane, thus visible to point [a]
        if (res > 0) return true; else return false;
    }
    
    /* Test if point [a] is facing point [b] by using their normal vectors [na, nb] */
    bool Algorithm::test2(Vector4 * na, Vector4 * nb){
        // --> a.b = |a| * |b| * cos()
        
        // [a.b] -> Compute dot product between normal vector [na] and [nb]
        float dotProduct = (na->x * nb->x) + (na->y * nb->y) + (na->z * nb->z);
        // cos() -> Compute cos using the dot product
        float cos_ = dotProduct / (na->length * nb->length);
        
        // We only consider those points visible using an angle of 170 degress
        // which represents a cos value of approx -0.985
        if (cos_ < -0.985) return true; // Point [b] is visible to [a]
        else return false; // Point [b] is not visible to [a]
    }
    
    /* Check if point [b] is visible to [a] */
    bool Algorithm::test3(Point * a, Point * b){
        // Create vector between point [a] and [b]
        Vector4 res = a->coord.Vector(&b->coord);
        Vector4 * normalA = &(a->normal);
        float dotProduct2   = (normalA->x * res.x) + (normalA->y * res.y) + (normalA->z * res.z);
        float magnitudeRes  = sqrtf(powf(res.x, 2.0)+powf(res.y, 2.0)+powf(res.z, 2.0));
        float magnitudeA    = sqrtf(powf(normalA->x, 2.0)+powf(normalA->y, 2.0)+powf(normalA->z, 2.0));
        float cAngle = dotProduct2 / (magnitudeA*magnitudeRes);
        if (cAngle >= a->fov_angle) return true; else return false;
    }
    
    /* Voxel Ray Tracer implementation, returns true if the first intersected voxel is a surface one */
    bool Algorithm::voxelRayTracer(Point * pt, Grid * grid, float &distance){
        bool debug=false;
        // Origin point coordinates of the ray
        Point A(pt->coord.x, pt->coord.y, pt->coord.z);
        
        // Ray vector components (i.e. normal vector of surface point [pt])
        Point n(pt->normal.x, pt->normal.y, pt->normal.z);
        
        // Voxel where point [pt] falls in
        Point sv = pt->interpolate(grid); // Start voxel
        
        // Starting position of the voxel of surface point [pt]
        int i = sv.coord.i;
        int j = sv.coord.j;
        int k = sv.coord.k;
        
        // In which direction the voxel ids are incremented. Needed to later update indices i, j, k
        int di = (n.coord.x >= 0) ? 1:-1;
        int dj = (n.coord.y >= 0) ? 1:-1;
        int dk = (n.coord.z >= 0) ? 1:-1;
        
        if (debug){
            printf("Surface Point              = (%f, %f, %f)\n", A.coord.x, A.coord.y, A.coord.z);
            printf("Normal Vector              = (%f, %f, %f)\n", n.coord.x, n.coord.y, n.coord.z);
            printf("Direction                  = (%d, %d, %d)\n", di, dj, dk);
        }
        
        // Compute next Voxel. That is, update bounding locations of X, Y, and Z
        float nextX = grid->X[i+di];
        float nextY = grid->Y[j-dj];
        float nextZ = grid->Z[k-dk];
        
        // Store the computed next voxel
        Point voxel;
        voxel.coord.x = nextX;
        voxel.coord.y = nextY;
        voxel.coord.z = nextZ;
        voxel.coord.i = i+di;
        voxel.coord.j = j-dj;
        voxel.coord.k = k-dk;
        //pt->path.push_back(voxel);
        
        if (debug){
            printf("Start Voxel                = [%d][%d][%d] -> (%f, %f, %f)\n", sv.coord.i, sv.coord.j, sv.coord.k, sv.coord.x, sv.coord.y, sv.coord.z);
            printf(" > Next voxel = [%d][%d][%d] -> (%f, %f, %f)\n", voxel.coord.i, voxel.coord.j, voxel.coord.k, voxel.coord.x, voxel.coord.y, voxel.coord.z);
        }
        
        // Compute tx, ty, tz that is the distance until next intersection with voxel-border;
        // the value of the parameter t at which the ray crosses the first vertical voxel boundary;
        // the expression of the ray is p + t . n = q, where p is the origin, n the nromal, and q i any ray point beyond p
        float tx = (n.coord.x!=0) ? (nextX - A.coord.x)/n.coord.x : FLT_MAX;
        float ty = (n.coord.y!=0) ? (nextY - A.coord.y)/n.coord.y : FLT_MAX;
        float tz = (n.coord.z!=0) ? (nextZ - A.coord.z)/n.coord.z : FLT_MAX;
        
        // Compute dx, dy, dz that's how far along the ray we must move for the horizontal
        // component to equal the width of a voxel the direction in which we traverse the grid
        // can only be FLT_MAX if we never go in that direction
        float dx = (n.coord.x!=0) ? grid->gridSpacing/n.coord.x*di : FLT_MAX;
        float dy = (n.coord.y!=0) ? grid->gridSpacing/n.coord.y*dj : FLT_MAX;
        float dz = (n.coord.z!=0) ? grid->gridSpacing/n.coord.z*dk : FLT_MAX;
        
        // Starting looping and storing remaining intersection voxels
        if (tx < ty) {
            if (tx < tz){
                i += di;
                tx += dx;
            }else{
                k -= dk;
                tz += dz;
            }
        } else {
            if (ty < tz){
                j -= dj;
                ty += dy;
            }else{
                k -= dk;
                tz += dz;
            }
        }
        int cont=0;
        while (i >= 0 && i <= grid->size_x && j>=0 && j<=grid->size_y && k>=0 && k<=grid->size_z){
            Point voxel;
            voxel.coord.x = grid->X[i];
            voxel.coord.y = grid->Y[j];
            voxel.coord.z = grid->Z[k];
            voxel.coord.i = i;
            voxel.coord.j = j;
            voxel.coord.k = k;
            //pt->path.push_back(voxel);
            // The current voxel of the path is a surface one, return true
            if (cont != 0) if (grid->map[voxel.coord.i][voxel.coord.j][voxel.coord.k].surface){
                distance = voxel.coord.Distance(&pt->coord);
                return true;
            }
            if (debug)
                printf(" > Next voxel = [%d][%d][%d] -> (%f, %f, %f)\n", voxel.coord.i, voxel.coord.j, voxel.coord.k, voxel.coord.x, voxel.coord.y, voxel.coord.z);
            if (tx < ty) {
                if (tx < tz){
                    i += di;
                    tx += dx;
                }else{
                    k -= dk;
                    tz += dz;
                }
            } else {
                if (ty < tz){
                    j -= dj;
                    ty += dy;
                }else{
                    k -= dk;
                    tz += dz;
                }
            }
            cont++;
        }
        return(false);
    }
    
    /* Voxel Ray Tracer implementation, return the first intersected voxel */
    Point Algorithm::voxelRayTracer_alt(Point * pt, Grid * grid){
        bool debug=false;
        // Origin point coordinates of the ray
        Point A(pt->coord.x, pt->coord.y, pt->coord.z);
        
        // Ray vector components (i.e. normal vector of surface point [pt])
        Point n(pt->normal.x, pt->normal.y, pt->normal.z);
        
        // Voxel where point [pt] falls in
        Point sv = pt->interpolate(grid); // Start voxel
        
        // Starting position of the voxel of surface point [pt]
        int i = sv.coord.i;
        int j = sv.coord.j;
        int k = sv.coord.k;
        
        // In which direction the voxel ids are incremented. Needed to later update indices i, j, k
        int di = (n.coord.x >= 0) ? 1:-1;
        int dj = (n.coord.y >= 0) ? 1:-1;
        int dk = (n.coord.z >= 0) ? 1:-1;
        
        if (debug){
            printf("Surface Point              = (%f, %f, %f)\n", A.coord.x, A.coord.y, A.coord.z);
            printf("Normal Vector              = (%f, %f, %f)\n", n.coord.x, n.coord.y, n.coord.z);
            printf("Direction                  = (%d, %d, %d)\n", di, dj, dk);
        }
        
        // Compute next Voxel. That is, update bounding locations of X, Y, and Z
        float nextX = grid->X[i+di];
        float nextY = grid->Y[j-dj];
        float nextZ = grid->Z[k-dk];
        
        // Store the computed next voxel
        Point voxel;
        voxel.coord.x = nextX;
        voxel.coord.y = nextY;
        voxel.coord.z = nextZ;
        voxel.coord.i = i+di;
        voxel.coord.j = j-dj;
        voxel.coord.k = k-dk;
        //pt->path.push_back(voxel);
        
        if (debug){
            printf("Start Voxel                = [%d][%d][%d] -> (%f, %f, %f)\n", sv.coord.i, sv.coord.j, sv.coord.k, sv.coord.x, sv.coord.y, sv.coord.z);
            printf(" > Next voxel = [%d][%d][%d] -> (%f, %f, %f)\n", voxel.coord.i, voxel.coord.j, voxel.coord.k, voxel.coord.x, voxel.coord.y, voxel.coord.z);
        }
        
        // Compute tx, ty, tz that is the distance until next intersection with voxel-border;
        // the value of the parameter t at which the ray crosses the first vertical voxel boundary;
        // the expression of the ray is p + t . n = q, where p is the origin, n the nromal, and q i any ray point beyond p
        float tx = (n.coord.x!=0) ? (nextX - A.coord.x)/n.coord.x : FLT_MAX;
        float ty = (n.coord.y!=0) ? (nextY - A.coord.y)/n.coord.y : FLT_MAX;
        float tz = (n.coord.z!=0) ? (nextZ - A.coord.z)/n.coord.z : FLT_MAX;
        
        // Compute dx, dy, dz that's how far along the ray we must move for the horizontal
        // component to equal the width of a voxel the direction in which we traverse the grid
        // can only be FLT_MAX if we never go in that direction
        float dx = (n.coord.x!=0) ? grid->gridSpacing/n.coord.x*di : FLT_MAX;
        float dy = (n.coord.y!=0) ? grid->gridSpacing/n.coord.y*dj : FLT_MAX;
        float dz = (n.coord.z!=0) ? grid->gridSpacing/n.coord.z*dk : FLT_MAX;
        
        // Starting looping and storing remaining intersection voxels
        if (tx < ty) {
            if (tx < tz){
                i += di;
                tx += dx;
            }else{
                k -= dk;
                tz += dz;
            }
        } else {
            if (ty < tz){
                j -= dj;
                ty += dy;
            }else{
                k -= dk;
                tz += dz;
            }
        }
        int cont=0;
        while (i >= 0 && i <= grid->size_x && j>=0 && j<=grid->size_y && k>=0 && k<=grid->size_z){
            Point voxel;
            voxel.coord.x = grid->X[i];
            voxel.coord.y = grid->Y[j];
            voxel.coord.z = grid->Z[k];
            voxel.coord.i = i;
            voxel.coord.j = j;
            voxel.coord.k = k;
            //pt->path.push_back(voxel);
            // The current voxel of the path is a surface one, return true
            if (cont != 0) if (grid->map[voxel.coord.i][voxel.coord.j][voxel.coord.k].surface) return voxel;
            
            if (debug)
                printf(" > Next voxel = [%d][%d][%d] -> (%f, %f, %f)\n", voxel.coord.i, voxel.coord.j, voxel.coord.k, voxel.coord.x, voxel.coord.y, voxel.coord.z);
            if (tx < ty) {
                if (tx < tz){
                    i += di;
                    tx += dx;
                }else{
                    k -= dk;
                    tz += dz;
                }
            } else {
                if (ty < tz){
                    j -= dj;
                    ty += dy;
                }else{
                    k -= dk;
                    tz += dz;
                }
            }
            cont++;
        }
        Point empty;
        empty.coord.i = -1;
        empty.coord.j = -1;
        empty.coord.k = -1;
        return(empty);
    }
    
    /* Filter points. Check if the normal of the point is pointing towards the surface (valid) 
       or to the outside of the surface (invalid). */
    bool Algorithm::validPoint(Grid * grid, Point * point, float threshold){
        bool debug=false;
        
        // Check if there is at least one surface voxel
        // If yes, then the [point] is valid (true) - It means that the point is pointing towards the surface
        // Otherwise is invalid (false) - It means that the point is pointing to the outside of the surface
        float distance; // distance to the intersected voxel from [point]
        bool existsVoxelOfSurface = voxelRayTracer(point, grid, distance);
        
        // There is at least one voxel, on the path of voxels from the starting point, that is a surface voxel
        // -> I.e. Surface [point] is pointing towards the surface and not to the outside of the surface
        if (existsVoxelOfSurface){
            point->dist = distance;
            // There are situations where the distance between a surface point and the first surface voxel, of the path, is not
            // large enough to be relevant to the algorithm. That is, is not large enough to fit a water molecule.
            // -> In order to measure this distance we use the test below in conjunction with a threshold
            Point B(point->coord.x + point->normal.x, point->coord.y + point->normal.y, point->coord.z + point->normal.z);
            Point BB=B.interpolate(grid);
            if (BB.coord.i >= 0 && BB.coord.i <= grid->size_x && BB.coord.j>=0 && BB.coord.j<=grid->size_y && BB.coord.k>=0 && BB.coord.k<=grid->size_z){
                float intensity = grid->map[BB.coord.i][BB.coord.j][BB.coord.k].intensity;
                if (intensity > threshold && intensity <= 1.0){ // Only store those with intensity greater than [threshold]
                    if (debug)
                        printf(" [BB] is pointing towards the surface [valid] (intensity[%d][%d][%d]=%f).\n", BB.coord.i, BB.coord.j, BB.coord.k, intensity);
                    return(true);
                }else{
                    if (debug)
                        printf(" [BB] is not pointing towards the surface [invalid] (intensity[%d][%d][%d]=%f).\n", BB.coord.i, BB.coord.j, BB.coord.k, intensity);
                    return(false);
                }
            }
        }else{
            // There isn't any voxel, on the path of voxels from the starting point, that is a surface voxel
            // -> Surface [point] is pointing to the outside of the surface
            return false;
        }
        // point invalid
        return false;
    }

    /* Recursive method that finds points visible from an initial point [a] and its visibles [a->final] */
    void Algorithm::find(vector<Point *> * R, Point * a, Point * original, int& inc, bool debug, int& gid, set<int> * joint){
        if (debug)
            printf("  Processing point [%d](%f, %f, %f) - inc=%d (%d)\n", a->id, a->coord.x, a->coord.y, a->coord.z, inc, (int) a->final.size());
        // Iterate visible points of point [a]
        for(vector<Point *>::iterator it = a->final.begin(); it != a->final.end(); it++){
            Point * add = *it;
            
            // Check if the current point was previously added to other group or the same one
            // If yes, skip it
            if (add->group_id != -1){ // If the group_id is != -1 means that the point was already added to the same or other group
                if (debug) printf("   -> Skipping point [%d] was already added to group [%d] \n", add->id, add->group_id);
                if (gid != add->group_id) joint->insert(add->group_id);
                continue;
            }
            if (debug) printf("   -> Add point [%d] to group [%d] \n", add->id, gid);
            
            // Union without repetitions
            bool exists=false;
            for(int c=0; c < R->size(); c++) if (R->at(c)->id == add->id){exists=true; break;}
            
            if (!exists){
                add->group_id = gid;  // Set the id of the group that the point is going to be added to
                R->push_back(add);
            }
            
            inc++;
            if (debug){
                printf("   -> R[%d]={", gid);
                for(int e=0; e < R->size(); e++) printf("%d, ", R->at(e)->id);
                printf("}\n");
            }
            find(R, add, original, inc, debug, gid, joint);
        }
        
    }
    
    /* Union all set of visible point that are visible to each others */
    void Algorithm::Union(vector<Point *> * points, vector<Group> * groups, bool debug, bool output){
        int inc=0;
        vector<Group> tmp_groups;
        if (debug) printf(" List of not-empty visible points of each point:\n");
        for(int i=0; i < points->size() && debug; i++){
            Point * a = (points->at(i));
            if (a->final.size()== 0) continue;
            printf("  Original set of visibles points of a[%d]={", a->id);
            for(vector<Point *>::iterator it = a->final.begin(); it != a->final.end(); it++)
                printf("%d, ", (*it)->id);
            printf("}\n");
            
        }
        
        if (debug) printf("\n\n");
        
        for(int i=0; i < points->size(); i++){
            vector<Point *> R; // vector of new visible points
            Group group;  // Current group
            group.id = (int) tmp_groups.size(); // id of the current group to be created
            Point * a = (points->at(i));
            
            // Recursively find the set of visible points (groups) of each point in [a->final]
            // ----------------------------------------------------------------------
            // For example:
            // a={b, c, d}, where b, c, and d are points visible to point [a]
            // b={e, f, g}, where e, f, and g are points visible to point [b]
            //
            // By Recursively finding visible points (group) of point [a] we get:
            // a={b, e, f, g, b, c, d}
            //
            // Due to the fact that the set of points of point [b] were included in the group of
            // [a], the group of points [b] will be empty:
            // b={}
            if (debug) printf("PROCESSING GROUP [%d] of a[%d]:\n", group.id, a->id);
            find(&R, a, a, inc, debug, group.id, &group.joint);
            
            if (debug){
            printf("RESULTS\n");
            printf("-------------------------------------------------------\n");
            }
            
            // Only store groups with visible points
            if (R.size() > 1){
                
                if (debug){
                if (group.joint.size() != 0)
                    printf(" Group [%d] should be joint to:", group.id);
                for (set<int>::iterator  it=group.joint.begin(); it!=group.joint.end(); ++it) printf("[%d]", *it);
                printf("\n");
                }
                
                group.a = a;
                group.color.random();
                group.points = R;
                tmp_groups.push_back(group);
                if (debug){
                    printf(" Recursive group of visible points of a[%d]=\n{", a->id);
                    for(int j=0; j < R.size(); j++) printf("%d,", R[j]->id);
                    printf("}\n");
                }
            }
            
            if (debug) printf("-------------------------------------------------------\n");
            R.clear(); //failsafe
        }
        
        if (debug){
        printf("\nList of Points for each Group:\n");
        for(int i=0; i < tmp_groups.size(); i++){
            printf("> Group [%d]={", tmp_groups[i].id);
            for(int j=0; j < tmp_groups[i].points.size(); j++){
                printf("%d,", tmp_groups[i].points[j]->id);
            }
            printf("\n");
        }
        printf("\n");
        }
        
        // Procedure to join groups
        // When a group [X], of surface point [a], is being created and there is at least one
        // point visible of [a] that has been added / belongs to another group [Y], the current
        // group [X] must be joined with group [Y].
        for(int i=0; i < tmp_groups.size(); i++){
            Group * A = &tmp_groups[i];
            Group ng;
            ng.pid      = A->id;
            ng.id       = (int) groups->size();
            ng.a        = A->a;
            ng.points   = A->points;
            ng.color.random();
            
            bool joint=false; // Flag to check if the group was joined with other or not
            if (debug) printf("Group [%d] will be joint to:", A->id);
            for (set<int>::iterator it=A->joint.begin(); it!=A->joint.end(); ++it){
                Group * B = &tmp_groups[*it]; // group id to joint to
                if (debug) printf("[%d]{", *it);
                
                // Search for group in the main vector of groups, to add points
                // we can't use those in [tmp_groups]
                for(int u=0; u < groups->size(); u++){
                    if (groups->at(u).pid == B->id){
                        // Iterate points of group
                        // - All points of the current group [A] will be moved to [B]
                        // - e.g. Group [1] will be joint to:[0] --> This means that all
                        //   points of group [1] will be copied to group [0]
                        for(int j=0; j < A->points.size(); j++){
                            groups->at(u).points.push_back(A->points[j]);
                            if (debug) printf("%d,", A->points[j]->id);
                        }
                        if (debug) printf("}");
                        joint=true;
                    }
                }
                
            }
            if (debug) printf("\n");
            
            if (!joint){
                // Add group to main vector of groups
                groups->push_back(ng);
                if (debug) printf(" > Group, with new id=[%d], added to the main vector of groups.\n", ng.id);
            }
        }
        
        if (output){
            printf("  Final list of Points for each Group:\n");
            for(int i=0; i < groups->size(); i++){
                printf("   > Group [%d]={", groups->at(i).id);
                for(int j=0; j < groups->at(i).points.size(); j++){
                    printf("%d,", groups->at(i).points[j]->id);
                }
                printf("}\n");
            }
        }
    }
    
    /* Auxiliar method for formCavities */
    bool Algorithm::intersect(vector<int> * A, vector<int> * B){
        int contB=0;
        int contA=(int) A->size();
        for(int i=0; i < A->size(); i++){
            for(int j=0; j < B->size(); j++){
                if (A->at(i) == B->at(j)) contB++;
            }
        }
        if (contB == contA) return true; else return false;
    }
    
    /* Form cavities by joining and filtering groups:
     - First , a group [X] is joined with group [Y] if [Y] is the closest to [X];
     - Second, cavities with isolated points are identified. Those points are removed from each identifed cavity.
     - Third , cavities containing points which do not point to any of them are removed
     - Fourth, Only top-n cavities are considered, being n the number of points of each cavity */
    void Algorithm::formCavities(float gridspacing,vector<Group> * visibles, vector<Cavity> * cavities, bool debug, bool output){
        // Default top to form cavities, the top is different from that provided by the user.
        // Instead of filtering "top" volumes, in this case, we are concerned to form cavities
        // with the most number of points.
        int top=10;
        // Joint groups
        struct Min{ Group * from; Group * to; float distance=-1.0; };
        struct Joint{ vector<int> ids; bool flag=false; };
        
        
        // -----------------------------------------
        // 1 - A group [X] is joined with group [Y] if [Y] is the closest to [X];
        // -----------------------------------------
        vector<Min> mins;
        for(int i=0; i < visibles->size(); i++){
            Group * resultA = &(visibles->at(i));
            
            if (debug){
                printf(" Group [%d]: {", resultA->id);
                for(int i=0; i < resultA->points.size(); i++) printf("%d,", resultA->points[i]->id);
                printf("}\n");
            }
            
            Min min;
            for(int j=0; j < resultA->points.size(); j++){
                Point * a = resultA->points.at(j);
                
                for(int ii=0; ii < visibles->size(); ii++){
                    if (i==ii) continue; // avoid the same
                    
                    Group * resultB = &(visibles->at(ii));
                    for(int jj=0; jj < resultB->points.size(); jj++){
                        //
                        Point * b = resultB->points.at(jj);
                        // Compute the distance between point [a] and [b] of each group
                        //float dist=Utils::Distances::distanceBetweenTwoPoints(a, b);
                        float dist = a->coord.Distance(&b->coord);
                        
                        // After the reduction of individual groups
                        // sometimes there are situations where the main group is
                        // flag to be joint to a much further group, due to the fact that individual groups
                        // are no longer in the vicinity.
                        // -> This is fixed by only joining groups where it's points are no more than the distance of the point to
                        //    to the closest voxel intersected.
                        
                        //if (dist > b->dist+gridspacing) continue;
                        if (dist > b->dist) continue;
                        
                        // Uncomment bellow if you need more debug data
                        // if (debug) printf("  Distance of A[%d] to B[%d] of group [%d]=%f\n", a->id, b->id, ii, dist);
                        if (min.distance == -1){
                            min.from = resultA;
                            min.to = resultB;
                            min.distance = dist;
                        }else{
                            if (min.distance > dist){
                                min.from = resultA;
                                min.to = resultB;
                                min.distance = dist;
                                
                            }
                        }
                    }
                    
                }
                
            }
            
            if (min.distance != -1){
                min.to->flag=true;
                mins.push_back(min);
            }
            
            if (debug)
                printf(" -> Min distance between two points was = %f. From group [%d] to group [%d]\n\n", min.distance, min.from->id, min.to->id);
            
            
        }
        
        // Uncomment bellow if you need more debug data
        /*
         for(int i=0; i < mins.size() && debug; i++){
         printf(" Group [%d] was joint with ={", mins[i].from->id);
         printf("%d", mins[i].to->id);
         printf("}\n");
         }
         */
        
        // Performs unions between groups
        vector<Joint> js;
        for(int i=0; i < mins.size(); i++){
            Min a = mins[i];
            // Uncomment bellow if you need more debug data
            // if (debug) printf("#->> Join [%d] to [%d]\n", a.from->id, a.to->id);
            bool found=false;
            for(int e=0; e < js.size(); e++){
                for(int j=0; j < js[e].ids.size(); j++){
                    if (js[e].ids[j] == a.to->id){
                        // Uncomment bellow if you need more debug data
                        // if (debug) printf("   > Found [%d] to add\n", a.from->id);
                        // avoid repetitions
                        bool rep=false;
                        for(int t=0; t < js[e].ids.size(); t++) if (js[e].ids[t] == a.from->id) rep=true;
                        if (!rep) js[e].ids.push_back(a.from->id);
                        found=true;
                    }
                }
            }
            if (!found){
                Joint j;
                j.ids.push_back(a.from->id);
                j.ids.push_back(a.to->id);
                js.push_back(j);
            }
        }
        // Uncomment bellow if you need more debug data
        /*
         if (debug) printf(" Initial before intersection check:\n");
         for(int i=0; i < js.size() && debug; i++){
         
         printf("[%d]={",i);
         for(int ii=0; ii < js[i].ids.size(); ii++){
         printf("%d,", js[i].ids[ii]);
         }
         printf("}\n\n");
         }
         */
        
        for(int i=0; i < js.size(); i++){
            for(int j=0; j < js.size(); j++){
                if (i==j) continue;
                if (intersect(&js[i].ids, &js[j].ids)){
                    js[i].flag=true;
                }
            }
            
        }
        
        if (output) printf(" Final list of joint groups (Clusters):\n");
        for(int i=0; i < js.size() && output; i++){
            if (js[i].flag) continue;
            printf(" -> Cluster [%d] is formed by groups={",i);
            for(int ii=0; ii < js[i].ids.size(); ii++){
                printf("%d,", js[i].ids[ii]);
            }
            printf("}\n");
        }
        if (output) printf("\n");
        
        // -----------------------------------------
        // 2 - Apply filters to form cavities
        // -----------------------------------------
        for(int i=0; i  < js.size(); i++){
            if (js[i].flag) continue; // avoid duplicates
            Cavity cav;
            cav.id = i;
            cav.color.random();
            
            // 2.1 -> Add the set of points of the group to the cavity
            for(int ii=0; ii < js[i].ids.size(); ii++){
                int id = js[i].ids[ii];
                Group * result = &visibles->at(id);
                for(int j=0; j < result->points.size(); j++)
                    cav.temp.push_back(result->points[j]);
                // Store the information related to set of groups that forms the cavities
                // (i.e. groups that were joint together to form the cavity)
                cav.groups.push_back(*result);
            }
            
            
            // 2.2 -> Cavities with isolated points are identified. Those points are removed from each identifed cavity.
            //    An isolated point is considered to be a point that does not intersects another.
            for(int j=0; j < cav.temp.size(); j++){
                Point * a = cav.temp[j];
                bool intersects=false;
                for(int jj=0; jj < cav.temp.size(); jj++){
                    if (j==jj) continue; // avoid the same
                    Point * b = cav.temp[jj];
                    //if (Utils::Distances::distanceBetweenTwoPoints(a, b) <= (0.8*2)){intersects=true;break;}
                    if (a->coord.Distance(&b->coord) <= (0.8*2)){intersects=true;break;}
                }
                if (intersects){
                    cav.points.push_back(a);
                }
            }
            
            // 2.3 -> Remove cavities containing points which do not point to any of them
            //        This is accomplished by using visibility test1, test2, and test3 (see documentation for more details)
            int nothingCount=0;
            for(int j=0; j < cav.points.size(); j++){
                Point * a = cav.points[j];
                bool pointsTo=false;
                for(int jj=0; jj < cav.points.size(); jj++){
                    if (j==jj) continue; // avoid the same
                    Point * b = cav.points[jj];
                    if (test1(a, b, false)){
                        if (test2(&a->normal, &b->normal)){
                            if (test3(a, b)){
                                pointsTo=true; break;
                            }
                        }
                    }
                }
                if (!pointsTo) nothingCount++;
            }
            
            // Cavities where ALL the set of points are not pointing to any point of the cavity are removed
            if (nothingCount != cav.points.size()){
                cavities->push_back(cav);
            }
            
        }
        
        
        // -----------------------------------------
        // 3 - Only top-n cavities are considered, being n the number of points of each cavity
        // See first comment of this method for more details.
        // -----------------------------------------
        if (top == -1){
            printf ("  Not using tops\n");
            for(int i=0; i < cavities->size(); i++) cavities->at(i).top = true;
            return;
        }
        
        if (output) printf(" Orded cavities by the number of points:\n");
        
        // Order by number of points in each cavity
        std::sort(cavities->begin(), cavities->end(), [](Cavity a, Cavity b) { return a.points.size() > b.points.size(); });
        
        // Identified individual sizes (i.e. number of points) of each cavity
        vector<int> sizes; // list of sizes
        for(int i=0; i < cavities->size(); i++){
            int size=(int) cavities->at(i).points.size();
            bool exists=false;
            vector<int>::iterator it = find_if(sizes.begin(), sizes.end(),  [size](int a){ return (a==size);});
            if (it != sizes.end())exists=true;
            if (!exists) sizes.push_back(size);
        }
        
        for(int i=0; i < sizes.size(); i++){
            int size=sizes[i];
            //printf("size=%d\n", size);
            if (i < top){
                for(int j=0; j < cavities->size(); j++){
                    if (cavities->at(j).points.size() == size){
                        if (output) printf(" -> Cavity [%d], size=%d (number of points of the cavity)\n", cavities->at(j).id, (int) cavities->at(j).points.size());
                        cavities->at(j).top=true;
                    }
                }
            }
        }
        
        if (output) printf("\n Final List of cavities:\n");
        for(int i=0; i < cavities->size() && output; i++){
            if (cavities->at(i).top)
                printf(" -> Cavity [%d] with %d points (TOP)\n", cavities->at(i).id, (int) cavities->at(i).points.size());
            else  printf(" -> Cavity [%d] with %d points\n", cavities->at(i).id, (int) cavities->at(i).points.size());
        }
        
    }
    
    /* Auxiliar method for IdentifyFillingAtoms */
    void combinePoints(qh_vertex_t * pp, Group * a, Group * b){
        int k=0;
        for(int i=0; i < a->points.size(); i++){
            pp[k].x = a->points[i]->coord.x;
            pp[k].y = a->points[i]->coord.y;
            pp[k].z = a->points[i]->coord.z;
            k++;
        }
        for(int i=0; i < b->points.size(); i++){
            pp[k].x = b->points[i]->coord.x;
            pp[k].y = b->points[i]->coord.y;
            pp[k].z = b->points[i]->coord.z;
            k++;
        }
        
    }
    
    /* Identify the set of atoms that fills each cavity */
    void Algorithm::IdentifyFillingAtoms(Grid * grid, vector<Cavity> * cavities, float gridSpacing, bool filter_spheres, bool debug, bool output){
        // Beware that at the current stage, more computational time is needed to perform this step.
        if (filter_spheres) printf("  Using filling spheres filtering (this step can take more time)...\n");
        
        // Vector of sets to map combined groups in order to avoid the recomputation of a pair of
        // groups that were previously processed.
        // e.g. cav->map[0]={1, 2} -> In position 0, means group 0, was combined with group 1 and 2
        // -
        // Take this opportunity also to reset group ids.
        for(int i=0; i < cavities->size(); i++){
            Cavity * cav = &(cavities->at(i));
            for(int j=0; j < cavities->at(i).groups.size(); j++){
                cavities->at(i).groups[j].pid = cavities->at(i).groups[j].id; // store previous id
                cavities->at(i).groups[j].id  = j; // add new id
                set<int> cb;
                cav->map.push_back(cb);
            }
        }
        
        // Use the set of points in every group to construct a convex polygon (i.e. global convex hull)
        for(int i=0; i < cavities->size(); i++){
            Cavity * cav = &(cavities->at(i));
            ConvexHull ch; // Convex Hull data;
            if (!cav->top) continue; // Process only the top cavities in terms of points
            if (debug) printf(" Processing Cavity [%d] of %d with %d groups\n", cav->id, (int) cavities->size(), (int) cav->groups.size());
            
            // Set of points of the global convex hull
            set<Point *> allPoints;
            
            // Get all group points and store it into one vector
            for(int j=0; j < cav->groups.size(); j++){
                Group * group = &(cav->groups[j]);
                for(int k=0; k < group->points.size(); k++) allPoints.insert(group->points[k]);
            }
            
            // Populate necessary data structures for the Convex Hull Library
            int n2=(int) allPoints.size();
            qh_vertex_t pp[n2];
            int k=0;
            // Store points into internal Convex Hull structure
            for(set<Point *>::iterator i = allPoints.begin(); i != allPoints.end(); i++){
                Point * point = (*i);
                pp[k].x = point->coord.x;
                pp[k].y = point->coord.y;
                pp[k].z = point->coord.z;
                k++;
            }
            qh_mesh_t mesh_ch = qh_quickhull3d(pp, n2);
            
            // The set of points was invalid to compute the convex hull
            if (mesh_ch.error) continue;
            
            // Add convex hull data into CavVis' structures
            for(int k=0; k < mesh_ch.nvertices; k++){
                Point point;
                point.coord.x = mesh_ch.vertices[k].x;
                point.coord.y = mesh_ch.vertices[k].y;
                point.coord.z = mesh_ch.vertices[k].z;
                // Add points of the computed convex hull
                ch.points.push_back(point);
            }
            
            // ###########################
            // VOXELIZE Convex Hull
            // ###########################
            vx_mesh_t* mesh_v;
            vx_point_cloud_t* result;
            mesh_v = vx_mesh_alloc(mesh_ch.nvertices, mesh_ch.nindices);
            
            mesh_v->indices = mesh_ch.indices;
            for (int i = 0; i < mesh_ch.nvertices; ++i) {
                qh_vertex_t v = mesh_ch.vertices[i];
                mesh_v->vertices[i].x = v.x;
                mesh_v->vertices[i].y = v.y;
                mesh_v->vertices[i].z = v.z;
            }
            
            // Add vertices and indices from the original mesh you want to voxelize
            // [...]
            // Later, Filling spheres are placed in the center of each
            // voxel, that in turn are computed by the voxelization library
            float m_gridSpacing=gridSpacing;
            // Precision factor to reduce "holes" artifact
            float precision = m_gridSpacing / 10.0; // Typically, precision = voxelsize / 10.
            
            // Run voxelization
            // If the cubes are not inside the convex hull it means we need a smaller resolution (grid spacing) that was set below
            result = vx_voxelize_pc(mesh_v, m_gridSpacing, m_gridSpacing, m_gridSpacing, precision);
            
            ch.spheres.clear(); // failsafe
            vector<Point> spheres_tmp;
            for(int i=0; i < result->nvertices; i++){
                Point centerOfVoxel;
                centerOfVoxel.coord.x = result->vertices[i].x;
                centerOfVoxel.coord.y = result->vertices[i].y;
                centerOfVoxel.coord.z = result->vertices[i].z;
                
                // Save spheres centers, that is = to center of the voxel
                if (filter_spheres){
                    // Only store spheres that are not inside the surface
                    Point BB=centerOfVoxel.interpolate(grid);
                    if (BB.coord.i >= 0 && BB.coord.i <= grid->size_x && BB.coord.j>=0 && BB.coord.j<=grid->size_y && BB.coord.k>=0 && BB.coord.k<=grid->size_z){
                        float intensity = grid->map[BB.coord.i][BB.coord.j][BB.coord.k].intensity;
                        // ## Filter 1 - Avoid those filling spheres that are intersecting the surface
                        if (intensity < 1.0) spheres_tmp.push_back(centerOfVoxel);
                        
                    }
                }else ch.spheres.push_back(centerOfVoxel);
            }
            if (filter_spheres){
                
                // ## Filter 2 - Use the DBScan algorithm in order to validate the filling spheres of the cavity.
                // Is important to have a continuous flow of filling spheres that delineates the cavity (i.e., only one cluster).
                // In the case where 2 clusters are found, the cluster with the most number filling spheres is the one considered.
                vector<Cluster> clusters = dbscan(1.4, 3, (int) spheres_tmp.size(), &spheres_tmp);
                // printf("clusters = %d\n", (int) clusters.size());
                int max_cluster_position=0;
                if (clusters.size() >= 2){
                    // Get the largest cluster position
                    float max=-1;
                    for(int jj=0; jj < clusters.size(); jj++){
                        if (clusters[jj].points.size() > max){
                            max = (int) clusters[jj].points.size();
                            max_cluster_position = jj;
                        }
                    }
                }
                // Store the final list of filling spheres
                for(int jj=0; jj < clusters[max_cluster_position].points.size(); jj++)
                    ch.spheres.push_back(*clusters[max_cluster_position].points[jj]);
            }
            
            
            // Compute the volume of the cavity based on the voxel of each filling spheres
            for(int jj=0; jj < ch.spheres.size(); jj++) ch.volume += m_gridSpacing*m_gridSpacing*m_gridSpacing;
            
            //printf("volume = %f with grid spacing = %f\n", ch.volume,m_gridSpacing);
            
            // Store global convex hull in the vector of convex hulls
            cav->chs.push_back(ch);
            
            // Free resources
            qh_free_mesh(mesh_ch);
            vx_mesh_free(mesh_v, false);
            vx_point_cloud_free(result);
            
        }
    }
    
    /* Cluster final list of cavities using a algorithm based on N-ary Trees - Legacy support */
    void Algorithm::cluster(Protein * protein, vector<Cavity> * cavities, int topc, float distance, bool debug, bool output){
        float threshold=distance;
        vector<tree<int>> trees; // Trees content
        vector<int> tree_roots;  // Trees ids
        vector<Cavity> fcavities; // not top cavities
        // Create trees, where the root/father is the current cavity and their childrens are those set of cavities
        // at a distance less than [threshold] of the father
        for(int i=0; i < cavities->size(); i++){
            Cavity * a = &cavities->at(i);
            if (a->flag) continue;
            Point gc_a = a->gc;
            tree<int> tr;
            tree<int>::iterator top, root;
            top=tr.begin();
            // Root/father
            root = tr.insert(top, a->id);
            
            for(int j=0; j < cavities->size(); j++){
                if (i==j) continue;
                Cavity * b = &cavities->at(j);
                if (b->flag) continue; // already in a tree
                Point gc_b = b->gc;
                //float res = Utils::Distances::distanceBetweenTwoPoints(gc_a, gc_b);
                float res = gc_a.coord.Distance(&gc_b.coord);
                if (debug) printf("  Distance between cav [%d] and [%d] = %.2f\n", a->id, b->id, res);
                if (res <= threshold){
                    // Add current cavity as a children of the current father (i.e. cavity [a] is near cavity [b]
                    if (debug) printf("  -> cav [%d] added as a children of [%d] (%.2f)\n", b->id, a->id, res);
                    tr.append_child(root, b->id);
                    b->flag=true;
                }
                
            }
            trees.push_back(tr);
            tree_roots.push_back(a->id);
            if (debug) printf("---------------------------------------------\n");
        }
        if (debug) printf("\n\n");
        
        for(int i=0; i < tree_roots.size() && output; i++){
            tree<int> tr = trees[i];
            tree<int>::iterator loc=::find(tr.begin(), tr.end(), tree_roots[i]);
            if(loc!=tr.end()) {
                tree<int>::sibling_iterator sib=tr.begin(loc);
                printf(" Cavity Parent [%d]={", *loc);
                while(sib!=tr.end(loc)) {
                    cout << (*sib) << ",";
                    ++sib;
                }
                cout << "}" << endl;
            }
        }
        
        // Iterate trees and join fathers and sons as one cavity
        for(int i=0; i < tree_roots.size(); i++){
            Cavity cav;
            cav.id = i;
            cav.color.random();
            tree<int> tr = trees[i];
            tree<int>::iterator loc=::find(tr.begin(), tr.end(), tree_roots[i]);
            if(loc!=tr.end()) {
                Cavity * result_initial;
                try{
                    result_initial = &cavities->at(*loc);
                }catch (const std::out_of_range& oor){ continue;}
                cav.volume += result_initial->volume;
                cav.area   += result_initial->area;
                
                for(int j=0; j < result_initial->atoms.size(); j++)
                    cav.atoms.push_back(result_initial->atoms[j]);
                
                tree<int>::sibling_iterator sib=tr.begin(loc);
                while(sib!=tr.end(loc)) {
                    //printf("%d", *sib);
                    Cavity * result = &cavities->at(*sib);
                    cav.volume += result->volume;
                    cav.area += result->area;
                    //printf("[%d](size=%d),", result->id, result->atoms.size());
                    for(int j=0; j < result->atoms.size(); j++)
                        cav.atoms.push_back(result->atoms[j]);
                    ++sib;
                }
                //cout << "}" << endl;
            }
            
            //cav.gc.computeGeometricCenter(&cav.atoms);
            fcavities.push_back(cav);
        }
        
        // # Create final cavity tops - sort by volume
        std::sort(fcavities.begin(), fcavities.end(), [](Cavity a, Cavity b) { return a.volume > b.volume; });
        for(int i=0; i < fcavities.size(); i++){
            if (i < topc) protein->cavities.push_back(fcavities.at(i));
        }
        
    }
    
    
}
