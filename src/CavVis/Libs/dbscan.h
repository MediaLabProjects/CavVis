#ifndef __PLib__dbscan__
#define __PLib__dbscan__
/* Copyright 2015 Gagarine Yaikhom (MIT License) and modified by Tiago Simões */
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "Utils.h"
#include "DataStructures.h"
using namespace std;
using namespace CavVis;
/* Cluster structure (dbscan) */
struct Cluster{
    int id=-1;              // ID of the cluster
    int dbscan_id;          // ID provided by the dbscan algorithm
    Color color;            // Cluster color
    bool top=false;         
    Vector4 gc;                // Geometric Center of the cluster
    float distance=-1.0;    // Distance from this cluster to something
    vector<Point *> points; // References to points belonging to this cluster
    bool flag=false;        // Generic flag
};

vector<Cluster> dbscan(double epsilon_, int minpts_, int numberpoints_, vector<Point> * parsed_points);

#endif /* defined(__CavVis__dbscan__) */
