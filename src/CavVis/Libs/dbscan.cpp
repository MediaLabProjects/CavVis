#include "dbscan.h"
bool debug_dbscan=false;
/* Copyright 2015 Gagarine Yaikhom (MIT License) and modified by Tiago Simões */

#define UNCLASSIFIED -1
#define NOISE -2

#define CORE_POINT 1
#define NOT_CORE_POINT 0

#define SUCCESS 0
#define FAILURE -3

typedef struct point_s point_t;
struct point_s {
    int id;
    double x, y, z;
    int cluster_id;
    Vector4 normal;
    Point * ref;
};

typedef struct node_s node_t;
struct node_s {
    unsigned int index;
    node_t *next;
};

typedef struct epsilon_neighbours_s epsilon_neighbours_t;
struct epsilon_neighbours_s {
    unsigned int num_members;
    node_t *head, *tail;
};

node_t *create_node(unsigned int index);
int append_at_end(
                  unsigned int index,
                  epsilon_neighbours_t *en);
epsilon_neighbours_t *get_epsilon_neighbours(
                                             unsigned int index,
                                             point_t *points,
                                             unsigned int num_points,
                                             double epsilon,
                                             double (*dist)(point_t *a, point_t *b));
void print_epsilon_neighbours(
                              point_t *points,
                              epsilon_neighbours_t *en);
void destroy_epsilon_neighbours(epsilon_neighbours_t *en);
void dbscan(
            point_t *points,
            unsigned int num_points,
            double epsilon,
            unsigned int minpts,
            double (*dist)(point_t *a, point_t *b));
int expand(
           unsigned int index,
           unsigned int cluster_id,
           point_t *points,
           unsigned int num_points,
           double epsilon,
           unsigned int minpts,
           double (*dist)(point_t *a, point_t *b));
int spread(
           unsigned int index,
           epsilon_neighbours_t *seeds,
           unsigned int cluster_id,
           point_t *points,
           unsigned int num_points,
           double epsilon,
           unsigned int minpts,
           double (*dist)(point_t *a, point_t *b));
double euclidean_dist(point_t *a, point_t *b);
double adjacent_intensity_dist(point_t *a, point_t *b);
unsigned int parse_input(
                         FILE *file,
                         point_t **points,
                         double *epsilon,
                         unsigned int *minpts, int np);
void print_points(
                  point_t *points,
                  unsigned int num_points, vector<Cluster> * clusters);

node_t *create_node(unsigned int index)
{
    node_t *n = (node_t *) calloc(1, sizeof(node_t));
    if (n == NULL)
        perror("Failed to allocate node.");
    else {
        n->index = index;
        n->next = NULL;
    }
    return n;
}

int append_at_end(
                  unsigned int index,
                  epsilon_neighbours_t *en)
{
    node_t *n = create_node(index);
    if (n == NULL) {
        free(en);
        return FAILURE;
    }
    if (en->head == NULL) {
        en->head = n;
        en->tail = n;
    } else {
        en->tail->next = n;
        en->tail = n;
    }
    ++(en->num_members);
    return SUCCESS;
}

epsilon_neighbours_t *get_epsilon_neighbours(
                                             unsigned int index,
                                             point_t *points,
                                             unsigned int num_points,
                                             double epsilon,
                                             double (*dist)(point_t *a, point_t *b))
{
    epsilon_neighbours_t *en = (epsilon_neighbours_t *)
    calloc(1, sizeof(epsilon_neighbours_t));
    if (en == NULL) {
        perror("Failed to allocate epsilon neighbours.");
        return en;
    }
    for (int i = 0; i < num_points; ++i) {
        if (i == index)
            continue;
        if (dist(&points[index], &points[i]) > epsilon)
            continue;
        else {
            if (append_at_end(i, en) == FAILURE) {
                destroy_epsilon_neighbours(en);
                en = NULL;
                break;
            }
        }
    }
    return en;
}

void print_epsilon_neighbours(
                              point_t *points,
                              epsilon_neighbours_t *en)
{
    if (en) {
        node_t *h = en->head;
        while (h) {
            printf("(%lfm, %lf, %lf)\n",
                   points[h->index].x,
                   points[h->index].y,
                   points[h->index].z);
            h = h->next;
        }
    }
}

void destroy_epsilon_neighbours(epsilon_neighbours_t *en)
{
    if (en) {
        node_t *t, *h = en->head;
        while (h) {
            t = h->next;
            free(h);
            h = t;
        }
        free(en);
    }
}

void dbscan(
            point_t *points,
            unsigned int num_points,
            double epsilon,
            unsigned int minpts,
            double (*dist)(point_t *a, point_t *b))
{
    unsigned int i, cluster_id = 0;
    for (i = 0; i < num_points; ++i) {
        if (points[i].cluster_id == UNCLASSIFIED) {
            if (expand(i, cluster_id, points,
                       num_points, epsilon, minpts,
                       dist) == CORE_POINT)
                ++cluster_id;
        }
    }
}

int expand(
           unsigned int index,
           unsigned int cluster_id,
           point_t *points,
           unsigned int num_points,
           double epsilon,
           unsigned int minpts,
           double (*dist)(point_t *a, point_t *b))
{
    int return_value = NOT_CORE_POINT;
    epsilon_neighbours_t *seeds =
    get_epsilon_neighbours(index, points,
                           num_points, epsilon,
                           dist);
    if (seeds == NULL)
        return FAILURE;
    
    if (seeds->num_members < minpts)
        points[index].cluster_id = NOISE;
    else {
        points[index].cluster_id = cluster_id;
        node_t *h = seeds->head;
        while (h) {
            points[h->index].cluster_id = cluster_id;
            h = h->next;
        }
        
        h = seeds->head;
        while (h) {
            spread(h->index, seeds, cluster_id, points,
                   num_points, epsilon, minpts, dist);
            h = h->next;
        }
        
        return_value = CORE_POINT;
    }
    destroy_epsilon_neighbours(seeds);
    return return_value;
}

int spread(
           unsigned int index,
           epsilon_neighbours_t *seeds,
           unsigned int cluster_id,
           point_t *points,
           unsigned int num_points,
           double epsilon,
           unsigned int minpts,
           double (*dist)(point_t *a, point_t *b))
{
    epsilon_neighbours_t *spread =
    get_epsilon_neighbours(index, points,
                           num_points, epsilon,
                           dist);
    if (spread == NULL)
        return FAILURE;
    if (spread->num_members >= minpts) {
        node_t *n = spread->head;
        point_t *d;
        while (n) {
            d = &points[n->index];
            if (d->cluster_id == NOISE ||
                d->cluster_id == UNCLASSIFIED) {
                if (d->cluster_id == UNCLASSIFIED) {
                    if (append_at_end(n->index, seeds)
                        == FAILURE) {
                        destroy_epsilon_neighbours(spread);
                        return FAILURE;
                    }
                }
                d->cluster_id = cluster_id;
            }
            n = n->next;
        }
    }
    
    destroy_epsilon_neighbours(spread);
    return SUCCESS;
}

double euclidean_dist(point_t *a, point_t *b)
{
    return sqrt(pow(a->x - b->x, 2) +
                pow(a->y - b->y, 2) +
                pow(a->z - b->z, 2));
}

void addToCluster(int id, vector<Cluster> * clusters, Point * np){
    for(int i=0; i < clusters->size(); i++){
        Cluster * cluster = &(clusters)->at(i);
        if (cluster->dbscan_id == id)
            cluster->points.push_back(np);
    }
}

void print_points(
                  point_t *points,
                  unsigned int num_points, vector<Cluster> * clusters)
{
    unsigned int i = 0;

    
    if (debug_dbscan) printf("Number of points: %u\n"
           " x     y     z     cluster_id\n"
           "-----------------------------\n"
           , num_points);
    
    // Create id of the clusters
    while (i < num_points) {
        Cluster tmp;
        int cluster_id = points[i].cluster_id;
        // check if already exits
        bool exists=false;
        if (cluster_id != -2){
            for(int j=0; j < clusters->size(); j++)
                if (clusters->at(j).dbscan_id==cluster_id){ exists=true; break;}
            if (!exists){
                Cluster tmp;
                tmp.dbscan_id = cluster_id;
                clusters->push_back(tmp);
            }
        }
        i++;
    }
    i=0;
    
    while (i < num_points) {
        if (debug_dbscan) printf("%5.2lf %5.2lf %5.2lf: %d\n",
               points[i].x,
               points[i].y, points[i].z,
               points[i].cluster_id);
        
        /*
        Point np;
        np.id = points[i].id;
        np.x = points[i].x;
        np.y = points[i].y;
        np.z = points[i].z;
        np.normal = points[i].normal;
        */
        Point * np;
        np = points[i].ref;
        if (points[i].cluster_id == -2){
            Cluster tmp;
            tmp.dbscan_id = -2;
            tmp.points.push_back(np);
            clusters->push_back(tmp);
        }else
            addToCluster(points[i].cluster_id, clusters, np);
        ++i;
    }
    
    // separate clusters with id=-2;
    // 1 - get max id in clusters
    int max=clusters->at(0).dbscan_id;
    for(int j=1; j < clusters->size(); j++){
        if (clusters->at(j).dbscan_id > max) max=clusters->at(j).dbscan_id;
    }
    if (debug_dbscan) printf("max=%d\n", max);
    
    // 2 -
}

unsigned int parse_input(
                         FILE *file,
                         point_t **points,
                         double *epsilon,
                         unsigned int *minpts, int np, vector<Point> * parsed_points)
{
    unsigned int num_points=np, i = 0;
    
    point_t *p = (point_t *)
    calloc(num_points, sizeof(point_t));
    if (p == NULL) {
        perror("Failed to allocate points array");
        return 0;
    }
    while (i < num_points) {
        //fscanf(file, "%lf %lf %lf\n",
               
          //     &(p[i].x), &(p[i].y), &(p[i].z));
        p[i].id = parsed_points->at(i).id;
        p[i].x = parsed_points->at(i).coord.x;
        p[i].y = parsed_points->at(i).coord.y;
        p[i].z = parsed_points->at(i).coord.z;
        p[i].ref = &parsed_points->at(i);
        p[i].normal.x = parsed_points->at(i).normal.x;
        p[i].normal.y = parsed_points->at(i).normal.y;
        p[i].normal.z = parsed_points->at(i).normal.z;
        
        p[i].cluster_id = UNCLASSIFIED;
        ++i;
    }
    *points = p;
    return num_points;
}


vector<Cluster> dbscan(double epsilon_, int minpts_,int numberpoints_, vector<Point> * parsed_points) {
    vector<Cluster> clusters;
    
    point_t *points;
    double epsilon=epsilon_;
    unsigned int minpts=minpts_;
    unsigned int num_points =
    parse_input(stdin, &points, &epsilon, &minpts, numberpoints_, parsed_points);
    if (num_points) {
        dbscan(points, num_points, epsilon,
               minpts, euclidean_dist);
        if (debug_dbscan) printf("Epsilon: %lf\n", epsilon);
        if (debug_dbscan) printf("Minimum points: %u\n", minpts);
        print_points(points, num_points, &clusters);
    }
    free(points);
    
    if (debug_dbscan) printf("\n\n");
    for(int j=0; j < clusters.size() && debug_dbscan; j++ ){
        printf("Cluster dbscan_id=%d\n", clusters.at(j).dbscan_id);
        
        for(int jj=0; jj < clusters.at(j).points.size(); jj++){
            printf("  -> (%f, %f, %f)\n", clusters.at(j).points[jj]->coord.x, clusters.at(j).points[jj]->coord.y, clusters.at(j).points[jj]->coord.z);
        }
        printf("\n");
    }
    
    //Add random color and compute the geometric center for each cluster
    for(int i=0; i < clusters.size(); i++){
        Color c;
        c.R = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        c.G = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        c.B = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        clusters[i].color = c;
        clusters[i].id = i;
        vector<Point *> nls; // ref o the points that will be used in the computation of the geometric center of the current cluster
        for(int j=0; j < clusters[i].points.size(); j++){
            nls.push_back((clusters[i].points[j]));
        }
        //Vec3 gc = Utils::Generic::computeGeometricCenter(nls);
        //Vector4 gc;
        //gc.computeGeometricCenter(nls);
        //clusters[i].gc = gc;
        //points->push_back(fClusters[i].gc);
    }
    
    
    return clusters;
}