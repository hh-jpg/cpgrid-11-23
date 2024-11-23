/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-04-30 16:32:51
 * @LastEditTime: 2024-06-19 14:53:49
 * @FilePath: /cpgrid/include/geo.h
 * @Description:
 *
 */

#ifndef GEO_H
#define GEO_H

#include <vector>
#include <cmath>
#include <petscsys.h>
#include <petscsystypes.h>
#include <petscsnes.h>
#include "point.h"
#include "grid.h"
#include "config.h"

// Macro definitions
#define MIN(i, j) ((i) < (j) ? (i) : (j))
#define MAX(i, j) ((i) > (j) ? (i) : (j))
#define LINE_INTERSECTION(a1, a2, b1, b2) (((a1 > b1) && (a2 < b2)) || ((a1 < b1) && (a2 > b2)))

// 单元描述
//  *
//  *   HEX8: 7        6
//  *         o--------z
//  *        /:       /|         z
//  *       / :      / |          ^   y (into page)
//  *    4 /  :   5 /  |          | /
//  *     o--------o   |          |/
//  *     |   o....|...o 2        o---> x
//  *     |  .3    |  /
//  *     | .      | /
//  *     |.       |/
//  *     o--------o
//  *     0        1
//  *
// edge_nodes 的map 关系
// edge_nodes_map[0][0] 将返回 0，表示边 0 的第一个节点的索引是 0。
// edge_nodes_map[0][1] 将返回 1，表示边 0 的第二个节点的索引是 1。

const unsigned int edge_nodes_map[12][2] =
    {
        {0, 1}, // Edge 0
        {0, 2}, // Edge 1
        {2, 3}, // Edge 2
        {0, 3}, // Edge 3
        {0, 4}, // Edge 4
        {1, 5}, // Edge 5
        {2, 6}, // Edge 6
        {3, 7}, // Edge 7
        {4, 5}, // Edge 8
        {5, 6}, // Edge 9
        {6, 7}, // Edge 10
        {4, 7}  // Edge 11
};
// 右手
const unsigned int side_nodes_map[6][4] =
    {
        {0, 3, 2, 1}, // Side 0 底面  mrst 5 unstructure mesh 4
        {0, 1, 5, 4}, // Side 1 前面 mrst 3 unstructure mesh 2
        {1, 2, 6, 5}, // Side 2 右面 mrst 2 unstructure mesh 1
        {2, 3, 7, 6}, // Side 3 后面 mrst 4 unstructure mesh 3
        {3, 0, 4, 7}, // Side 4 左面 mrst 1 unstructure mesh 0
        {4, 5, 6, 7}  // Side 5 顶面 mrst 6 unstructure mesh 5
};

// denote the face in the neigbor cell face number
const unsigned int neighbor_face[6] = {5, 3, 4, 1, 2, 0};
// 用数组简化结构定义，提高灵活性和可扩展性
typedef struct
{
  Point nodes[8]; // 每个单元有8个节点
} HexNodes;

typedef struct
{
  PetscScalar face[6]; // 每个单元有6个面
} HexFaces;

typedef struct
{
  Point normals[6]; // 假设有6个法向量
} HexFaceNormals;

typedef struct
{
  Point centroids[6]; // 每个单元有6个面心
} HexFaceCentroids;

typedef struct
{
  PetscScalar indices[6]; // 每个单元有6个邻居
} HexNeighborIndices;

typedef struct
{
  Point centroid;
  double volume;
} CellCentroidVol;


// Function declarations
Point cross(const Point &a, const Point &b);
int cell_index(const int i, const int j, const int k, const GridSize &gridsize);
Point polygonCellCentroid(const std::vector<Point> &points,
                          const Point &face_centroid,
                          const Point &cell_centroid);
Point polygonCellCentroid(const std::vector<Point *> &points,
                          const Point &face_centroid,
                          const Point &cell_centroid);

double inner(const std::vector<double> &a, const std::vector<double> &b);
double determinantOf(const Point *a);

double triarea(const Point *c);
double simplex_volume(const Point *a);

Point average(const std::vector<Point> &points);
Point average(const std::vector<Point *> &points);

Point polygonCentroid(const std::vector<Point> &Points, const Point &inPoint);
Point polygonCentroid(const std::vector<Point *> &points, const Point &inPoint);

Point polygonNormal(const std::vector<Point> &Points, const Point &centroid);
Point polygonNormal(const std::vector<Point *> &Points, const Point &centroid);

double polygonArea(const std::vector<Point> &Points, const Point &centroid);
double polygonArea(const std::vector<Point *> &Points, const Point &centroid);

double polygonCellVolume(const std::vector<Point> &Points,
                         const Point &face_centroid, const Point &cell_centroid);
double polygonCellVolume(const std::vector<Point*> &Points,
                         const Point &face_centroid, const Point &cell_centroid);

double distance(Point p1, Point p2);
HexFaces faceareas_setup(const HexNodes &current_cell);
HexFaceNormals facenormals_setup(const HexNodes &current_cell);
HexFaceCentroids facecentroids_setup(const HexNodes &current_cell);
CellCentroidVol calculateHexCellCV(const HexNodes &current_cell);
PetscErrorCode find_neighbor(Grid *grid);
void index_to_ijk(int index, const GridSize &gridsize, int &i, int &j, int &k);
int ijk_to_index(const int i, const int j, const int k, const GridSize &gridsize);
#endif // GEOMETRY_H
