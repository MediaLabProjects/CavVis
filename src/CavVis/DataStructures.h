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
#ifndef __CavVis__DataStructures__
#define __CavVis__DataStructures__
#include <stdio.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <set>
#include <cmath>
using namespace std;
namespace CavVis{
    /* Atom Types */
    const float hydrogen_H      = 1.2;
    const float carbon_C        = 1.7;
    const float nitrogen_N      = 1.55;
    const float oxygen_O        = 1.52;
    const float fluorine_F      = 1.47;
    const float phosphorous_P   = 1.8;
    const float sulfur_S        = 1.8;
    const float chlorine_Ci     = 1.75;
    const float copper_Cu       = 1.4;
    class Atom;
    class Protein;
    class Grid;
    class Point;
    
    struct Vector4{
    public:
        int id=-1;
        float x, y, z;
        int i=-1, j=-1, k=-1;
        float length=0.0;
        //
        Vector4 (): x(0), y(0), z(0), i(-1), j(-1), k(-1) {}
        Vector4 (float _x, float _y, float _z);
        Vector4 (float _x, float _y, float _z, int _i, int _j, int _k);
        /* Create vector using two points */
        Vector4 Vector(Vector4 * b);
        /* Compute length */
        void Length ();
        /* Compute distance between two points */
        float Distance(Vector4 * a);
        /* Normalize current Vector4 object */
        void Normalize();
    };
    
#pragma region Operators
    /* Sum of 2 vectors */
    inline Vector4 operator+(const Vector4& lhs, const Vector4& rhs){
        Vector4 result;
        result.x = lhs.x + rhs.x;
        result.y = lhs.y + rhs.y;
        result.z = lhs.z + rhs.z;
        return result;
    }
    
    /* Substraction of 2 vectors */
    inline Vector4 operator-(const Vector4& lhs, const Vector4& rhs){
        Vector4 result;
        result.x = lhs.x - rhs.x;
        result.y = lhs.y - rhs.y;
        result.z = lhs.z - rhs.z;
        return result;
    }
    
    /* Multiplication of vector and number */
    inline Vector4 operator*(const Vector4& lhs, const float& rhs){
        Vector4 result;
        result.x = lhs.x * rhs;
        result.y = lhs.y * rhs;
        result.z = lhs.z * rhs;
        return result;
    }
    
    /* multiplication of number and vector */
    inline Vector4 operator*(const float& lhs, const Vector4& rhs){
        return rhs*lhs;
    }
    
    /* Division of vector and number */
    inline Vector4 operator/(const Vector4& lhs, const float& rhs){
        Vector4 result;
        result.x = lhs.x / rhs;
        result.y = lhs.y / rhs;
        result.z = lhs.z / rhs;
        return result;
    }
    
    /* Dot Product */
    inline float operator*(const Vector4& lhs, const Vector4& rhs){
        return (lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z);
    }
    
    /* 3D Cross Product */
    inline Vector4 operator%(const Vector4& lhs, const Vector4& rhs){
        Vector4 result;
        result.x = lhs.y * rhs.z - lhs.z * rhs.y;
        result.y = lhs.z * rhs.x - lhs.x * rhs.z;
        result.z = lhs.x * rhs.y - lhs.y * rhs.x;
        return result;
    }
    
#pragma endregion Operators
    
    class Color{
    public:
        float R, G, B;
        //
        Color(){}
        ~Color(){}
        Color(float _R, float _G, float _B);
        /* Generates a random color to the current color object */
        void random();
    private:
    };
    
    class Point{
    public:
        int id;
        Vector4 coord;                              // Coordinates position.
        Vector4 normal;                             // Normal vector of the point.
        float dist      = 0.0;
        int group_id    = -1;                       // Identification of the current group that point belongs to
        bool invalid    = false;                    // Filter domain flag, is a valid point to be processed ?
        float fov_angle = -1.0;                     // Variable FoV angle
        vector<Point *> final;                      // Final list of points visible to the current point (in the FOV and valid back-face).
        //
        Point(){}
        ~Point(){}
        Point(float _x, float _y, float _z);
        /* For std::set.insert() */
        bool operator<(Point const & rhs) const{ return std::tie(coord.x, coord.y, coord.z) < std::tie(rhs.coord.x, rhs.coord.y, rhs.coord.z);}
        /* Interpolates point to the parsed grid */
        Point interpolate(Grid * grid);
        /* Creates point object with (x,y,z) equal to the geometric center of a vector of Atom * */
        void computeGeometricCenter(vector<Atom *> points);
        void computeGeometricCenter(vector<Atom> * points);
        void computeGeometricCenter(vector<Point> * points);
    private:
        //
    };
    
    /* Group of points */
    class Group{
    public:
        int id  =-1;                                // Current id.
        int pid =-1;                                // Previous id.
        Vector4 gc;                                 // Geometric center of the group.
        Point * a;                                  // Starting point of the group. Each group is formed by the set of visible points of [a]
        // and their visibles, and so forth...
        vector<Point *> points;                     // Points beloning to the current group object.
        set<int> joint;                             // List of group ids that the current group should be joint to.
        bool flag=false;                            // Generic flag.
        Color color;                                // Color of the group.
        //
        Group(){};
        ~Group(){};
    private:
        //
    };
    
    class ConvexHull{
    public:
        int id;
        vector<Point> spheres;                      // Spheres that fills the convex hull.
        float volume=0.0;
        float area=0.0;
        vector<Point> points;
        //
        ConvexHull(){}
        ~ConvexHull(){}
    private:
        //
    };
    
    class Cavity{
    public:
        int id=0;                                   // Id of the cavity.
        float volume=0.0;                           // Volume of the cavity
        float area=0.0;                             // Area of the cavity (only if user-requested)
        Point gc;                                   // Geometric center of the cavity, calculated using filling atoms.
        vector<Atom> atoms;                         // Filling atoms, dummy atoms that fill the entire cavity region.
                                                    // These atoms are not presented on the protein.
        Color color;                                // Atom colors.
        Color * regionColor;                        // Color of the cavity region.
        vector<Group> groups;                       // Each cavity is formed by a set of groups.
        vector<Point *> temp;
        vector<Point *> points;                     // Final list of points of the cavity.
        vector<ConvexHull> chs;                     // Group of convex hulls formed by the set of groups.
        vector<set<int>> map;                       // Vector of sets to map combined groups in order to avoid the recomputation of convex hulls.
        bool top=false;                             // Default top to form cavities, this top is different from that provided by the user.
                                                    // See formCavities() for more details.
        bool flag=false;
        //
        Cavity(){};
        ~Cavity(){};
    private:
        //
    };
    
    // Marchinb Cubes Map Point
    class mPoint{
    public:
        float intensity;                            // Intensity value on the current [i][j][k] position.
        bool surface=false;                         // Flags that indicates if the current [i][j][k] position of
                                                    // the grid used to interpolate a triangle point.
        set<Point *> * tPoints;                     // Interpolated triangle points on the current [i][j][k] position.
        vector<Atom *> * atoms;                     // Set of atoms that intersected the current [i][j][k] position.
        //
        mPoint(){}
        ~mPoint(){}
    private:
        //
    };
    
    class Grid{
    public:
        /* Generic grid parameters */
        float gridSpacing=0.6;                      // Grid spacing value.
        float padding=2.0;                          // The padding of the grid is typically the maximum size of an atom (i.e. ~ 2.0).
        int size_x;                                 // Number of cells in each axis (x).
        int size_y;                                 // Number of cells in each axis (y).
        int size_z;                                 // Number of cells in each axis (z).
        float minx,miny,minz,maxx,maxy,maxz;        // Max, min values.
        int mini,minj,mink,maxi,maxj,maxk;          // Max, min index.
        float *X,*Y,*Z;                             // The triplet represents each voxel in x, y, z
        // (e.g. X[0], Y[0], Z[0] is the first voxel of the grid).
        mPoint ***map;                              // Maps the Intensity on each vertixe x,y,z (using the gaussian function).
        // and the triangle interpolated point at a [i][j][j] position of the grid.
        
        /* Marching cubes (MC) parameters */
        float iso = 1.0;                            //
        float C   = 0.33;                           // Smooth parameter.
        float d   = 2.35;                           // Decay rate parameter of the gaussian surface.
        int seg   = 80;                             // Charge influence (Dynamically computed in the MC algorithm).
        int ioc   = 40;                             // Dynamically computed in the MC algorithm.
        //
        Grid(){};
        ~Grid(){};
    private:
        //
    };
    
    class Protein{
    public:
        enum type {APO, HOLO, UNSPECIFIED} type;
        string id;                                  // Identification/id of the protein (PDB id).
        long numberOfatoms=0;                       // Number of atoms of the current protein.
        long numberOfHeteroAtoms=0;                 // Number of hetero atoms of the current protein.
        float max_radius=2.0;                       // The max radius of an atom of the protein.
        vector<Atom> atoms;                         // Atoms of the current id protein.
        vector<Cavity> temp;                        // Temp cavities.
        vector<Cavity> cavities;                    // Final set of cavities detected.
        vector<Cavity *> pcavities;                 // Temporary pointers
        Color * color;                              // Generic color.
        float min_x,min_y,min_z;
        float max_x,max_y,max_z;
        //
        Protein(){};
        ~Protein(){};
    private:
        //
    };
    
    class Atom{
    public:
        int id;
        string type;                                // Atom type.
        Vector4 coord;                              // Atom coordinates.
        float radius;                               // Atom radius.
        Color color;                                // Atom color
        //
        Atom(){};
        ~Atom(){};
        Atom(int _id, string _type, float _x, float _y, float _z, float _radius, Color _color);
        /* Check if Vector4 is inside the current atom/sphere */
        bool pointInside(Vector4 * p);
        /* Check if point is inside the current atom/sphere */
        bool pointInside(Point * p);
        /* Interpolates atom center to a parsed grid */
        Point interpolate(Grid * grid);
        /* For std::set.insert() */
        bool operator<(Atom const & rhs) const{ return id < rhs.id;}
    private:
        //
    };
    
    // Marching Cubes vertices
    struct mVertex{
        Vector4 coord;
        Vector4 normal;                            // Normal of [vertex].
        Vector4 edge_0;                            // The first point used to interpolate [vertex].
        Vector4 edge_1;                            // The second point used to interpolate [vertex].
        int repeated=1;                            // How many times this mVertex is shared with other triangles.
    };
    
    // Marching Cubes Triangles
    struct mTriangle{
        Color * color;                            // Triangle color.
        // Each triangle is composed by Vector4 vertex objects
        vector<mVertex *> vertices;
    };
}

#endif /* defined(__CavVis__DataStructures__) */

