/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-04-30 17:06:46
 * @LastEditTime: 2024-08-28 16:39:50
 * @FilePath: /cpgrid/src/geo.C
 * @Description:
 *
 */

#include "geo.h"
#include "grid.h"
#include <cassert>
#include <numeric>
#include <algorithm>
#include <cmath>

Point cross(const Point &p1, const Point &p2)
{
  return Point(
      p1.y() * p2.z() - p1.z() * p2.y(),
      p1.z() * p2.x() - p1.x() * p2.z(),
      p1.x() * p2.y() - p1.y() * p2.x());
}

double inner(const std::vector<double> &a, const std::vector<double> &b)
{
  return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
}

double determinantOf(const Point *a)
{
  return a[0].x() * (a[1].y() * a[2].z() - a[2].y() * a[1].z()) -
         a[0].y() * (a[1].x() * a[2].z() - a[2].x() * a[1].z()) +
         a[0].z() * (a[1].x() * a[2].y() - a[2].x() * a[1].y());
}

Point average(const std::vector<Point> &points)
{
  int num_points = points.size();
  assert(num_points > 0);
  Point pt = points[0];
  for (int i = 1; i < num_points; ++i)
  {
    pt += points[i];
  }
  pt /= double(num_points);
  return pt;
}

Point average(const std::vector<Point *> &points)
{
  int num_points = points.size();
  assert(num_points > 0);
  Point pt(0.0, 0.0, 0.0);
  for (int i = 0; i < num_points; ++i)
  {
    pt += *(points[i]);
  }
  pt /= double(num_points);
  return pt;
}

double triarea(const Point *c)
{
  Point d0 = c[1] - c[0];
  Point d1 = c[2] - c[0];
  Point crossprod = cross(d0, d1);
  return 0.5 * crossprod.two_norm();
}

double simplex_volume(const Point *a)
{
  Point tmp[3];
  for (int i = 0; i < 3; ++i)
  {
    tmp[i] = a[i + 1] - a[i];
  }
  return determinantOf(tmp) / (6.0);
}

Point polygonCentroid(const std::vector<Point> &Points, const Point &inPoint)
{
  double tot_area = 0.0;
  Point tot_centroid;
  int num_Points = Points.size();
  for (int i = 0; i < num_Points; ++i)
  {
    Point tri[3] = {inPoint, Points[i], Points[(i + 1) % num_Points]};
    double tri_area = triarea(tri);
    Point tri_w_mid = (tri[0] + tri[1] + tri[2]);
    tri_w_mid *= tri_area / 3.0;
    tot_area += tri_area;
    tot_centroid += tri_w_mid;
  }
  if (std::abs(tot_area) > 0.0)
  {
    tot_centroid /= tot_area;
  }
  else
  {
    tot_centroid = inPoint;
  }
  return tot_centroid;
}

Point polygonCentroid(const std::vector<Point *> &points, const Point &inPoint)
{
  double tot_area = 0.0;
  Point tot_centroid;
  int num_points = points.size();
  for (int i = 0; i < num_points; ++i)
  {
    Point tri[3] = {inPoint, *(points[i]), *(points[(i + 1) % num_points])};
    double tri_area = triarea(tri);
    Point tri_w_mid = (tri[0] + tri[1] + tri[2]);
    tri_w_mid *= tri_area / 3.0;
    tot_area += tri_area;
    tot_centroid += tri_w_mid;
  }
  if (std::abs(tot_area) > 0.0)
  {
    tot_centroid /= tot_area;
  }
  else
  {
    tot_centroid = inPoint;
  }
  return tot_centroid;
}

Point polygonNormal(const std::vector<Point *> &face, const Point &centroid)
{
  Point tot_normal(0.0, 0.0, 0.0);
  double tot_area = 0.0;
  int num_points = face.size();
  for (int i = 0; i < num_points; ++i)
  {
    Point tri[3] = {centroid, *(face[i]), *(face[(i + 1) % num_points])};
    Point d0 = tri[1] - tri[0];
    Point d1 = tri[2] - tri[0];
    Point w_normal = cross(d0, d1) / 2.0;
    tot_area += triarea(tri);
    tot_normal += w_normal;
  }
  double length = tot_normal.two_norm();
  if (length > 0.0)
  {
    tot_normal /= tot_area;
  }
  return tot_normal;
}

Point polygonNormal(const std::vector<Point> &face, const Point &centroid)
{
  Point tot_normal(0.0, 0.0, 0.0);
  double tot_area = 0.0;
  int num_points = face.size();
  for (int i = 0; i < num_points; ++i)
  {
    Point tri[3] = {centroid, face[i], face[(i + 1) % num_points]};
    Point d0 = tri[1] - tri[0];
    Point d1 = tri[2] - tri[0];
    Point w_normal = cross(d0, d1) / 2.0;
    tot_area += triarea(tri);
    tot_normal += w_normal;
  }
  double length = tot_normal.two_norm();
  if (length > 0.0)
  {
    tot_normal /= tot_area;
  }
  return tot_normal;
}

Point polygonCellCentroid(const std::vector<Point> &points,
                          const Point &face_centroid,
                          const Point &cell_centroid)
{
  Point centroid(0.0);            // 初始化质心为零向量
  double tot_volume = 0.0;        // 初始化总体积为零
  int num_points = points.size(); // 获取多边形顶点数

  for (int i = 0; i < num_points; ++i)
  {
    // 定义一个四面体，顶点分别是 cell_centroid, face_centroid, points[i], points[(i+1)%num_points]
    Point tet[4] = {cell_centroid, face_centroid, points[i], points[(i + 1) % num_points]};

    // 计算该四面体的体积
    double small_volume = std::fabs(simplex_volume(tet));
    assert(small_volume >= 0);

    // 计算该四面体的质心
    Point small_centroid = (tet[0] + tet[1] + tet[2] + tet[3]) / 4.0;

    // 乘以体积进行加权
    small_centroid *= small_volume;

    // 累加加权质心
    centroid += small_centroid;

    // 累加体积
    tot_volume += small_volume;
  }

  // 归一化质心
  centroid /= tot_volume;
  assert(tot_volume >= 0);

  return centroid; // 返回计算得到的多面体质心
}

Point polygonCellCentroid(const std::vector<Point *> &points,
                          const Point &face_centroid,
                          const Point &cell_centroid)
{
  Point centroid(0.0);            // 初始化质心为零向量
  double tot_volume = 0.0;        // 初始化总体积为零
  int num_points = points.size(); // 获取多边形顶点数

  for (int i = 0; i < num_points; ++i)
  {
    // 定义一个四面体，顶点分别是 cell_centroid, face_centroid, points[i], points[(i+1)%num_points]
    Point tet[4] = {cell_centroid, face_centroid, *(points[i]), *(points[(i + 1) % num_points])};

    // 计算该四面体的体积
    double small_volume = std::fabs(simplex_volume(tet));
    assert(small_volume >= 0);

    // 计算该四面体的质心
    Point small_centroid = (tet[0] + tet[1] + tet[2] + tet[3]) / 4.0;

    // 乘以体积进行加权
    small_centroid *= small_volume;

    // 累加加权质心
    centroid += small_centroid;

    // 累加体积
    tot_volume += small_volume;
  }

  // 归一化质心
  centroid /= tot_volume;
  assert(tot_volume >= 0);

  return centroid; // 返回计算得到的多面体质心
}

double polygonArea(const std::vector<Point> &Points, const Point &centroid)
{
  double tot_area = 0.0;
  int num_Points = Points.size();
  for (int i = 0; i < num_Points; ++i)
  {
    Point tri[3] = {centroid, Points[i], Points[(i + 1) % num_Points]};
    tot_area += triarea(tri);
  }
  return tot_area;
}

double polygonArea(const std::vector<Point *> &Points, const Point &centroid)
{
  double tot_area = 0.0;
  int num_Points = Points.size();
  for (int i = 0; i < num_Points; ++i)
  {
    Point tri[3] = {centroid, *(Points[i]), *(Points[(i + 1) % num_Points])};
    tot_area += triarea(tri);
  }
  return tot_area;
}

double polygonCellVolume(const std::vector<Point> &Points,
                         const Point &face_centroid, const Point &cell_centroid)
{
  double tot_volume = 0.0;
  int num_Points = Points.size();
  for (int i = 0; i < num_Points; ++i)
  {
    Point tet[4] = {cell_centroid, face_centroid, Points[i], Points[(i + 1) % num_Points]};
    double small_volume = std::fabs(simplex_volume(tet));
    assert(small_volume >= 0);
    tot_volume += small_volume;
  }
  assert(tot_volume >= 0);
  return tot_volume;
}
double polygonCellVolume(const std::vector<Point *> &Points,
                         const Point &face_centroid, const Point &cell_centroid)
{
  double tot_volume = 0.0;
  int num_Points = Points.size();
  for (int i = 0; i < num_Points; ++i)
  {
    Point tet[4] = {cell_centroid, face_centroid, *(Points[i]), *(Points[(i + 1) % num_Points])};
    double small_volume = std::fabs(simplex_volume(tet));
    assert(small_volume >= 0);
    tot_volume += small_volume;
  }
  assert(tot_volume >= 0);
  return tot_volume;
}

int cell_index(const int i, const int j, const int k, const GridSize &gridsize)
{
  int nx = gridsize.nx;
  int ny = gridsize.ny;
  int nz = gridsize.nz;
  if (i < 0 || i >= nx || j < 0 || j >= ny || k < 0 || k >= nz)
  {
    return -1;
  }
  // 计算索引并返回
  return i + j * nx + k * nx * ny;
}

void index_to_ijk(const int index, const GridSize &gridsize, int &i, int &j, int &k)
{
  int nx = gridsize.nx;
  int ny = gridsize.ny;
  // 计算i, j, k
  k = index / (nx * ny);       // 确定第k层
  int rem = index % (nx * ny); // 第k层内的一维索引
  j = rem / nx;                // 确定第j行
  i = rem % nx;                // 确定第i列
}

int ijk_to_index(const int i, const int j, const int k, const GridSize &gridsize)
{
  int nx = gridsize.nx;
  int ny = gridsize.ny;
  return k * (nx * ny) + j * nx + i;
}

double distance(Point p1, Point p2)
{
  // 计算两个点之间的差向量
  Point diff = p1 - p2;
  // 计算差向量的模，即两点之间的距离
  return diff.two_norm();
}

HexFaces faceareas_setup(const HexNodes &current_cell)
{
  HexFaces hex_face_areas;
  for (int i = 0; i < 6; i++)
  {
    std::vector<Point> current_face;
    for (int j = 0; j < 4; j++)
    {
      current_face.push_back(current_cell.nodes[side_nodes_map[i][j]]);
    }
    Point face_avg = average(current_face);
    hex_face_areas.face[i] = polygonArea(current_face, face_avg);
  }
  return hex_face_areas;
}

HexFaceNormals facenormals_setup(const HexNodes &current_cell)
{
  HexFaceNormals hex_face_normals;
  for (int i = 0; i < 6; i++)
  {
    std::vector<Point> current_face;
    for (int j = 0; j < 4; j++)
    {
      current_face.push_back(current_cell.nodes[side_nodes_map[i][j]]);
    }
    Point avg = average(current_face);
    // Point centroid = polygonCentroid(current_face, avg);
    hex_face_normals.normals[i] = polygonNormal(current_face, avg);
  }
  return hex_face_normals;
}

HexFaceCentroids facecentroids_setup(const HexNodes &current_cell)
{
  HexFaceCentroids hex_face_centroids;
  for (int i = 0; i < 6; i++)
  {
    std::vector<Point> current_face;
    for (int j = 0; j < 4; j++)
    {
      current_face.push_back(current_cell.nodes[side_nodes_map[i][j]]);
    }
    Point avg = average(current_face);
    hex_face_centroids.centroids[i] = polygonCentroid(current_face, avg);
  }
  return hex_face_centroids;
}

Point calculateMean(const std::vector<Point> &vec)
{
  if (vec.empty())
  {
    throw std::invalid_argument("Vector is empty.");
  }
  Point sum(0.0, 0.0, 0.0);
  for (unsigned int i = 0; i < vec.size(); i++)
  {
    sum += vec[i];
  }
  Point avg = sum / vec.size();
  return avg;
}

CellCentroidVol calculateHexCellCV(const HexNodes &current_cell)
{
  double tot_cell_vol = 0.0;
  Point cell_centroid(0.0, 0.0, 0.0);
  std::vector<Point> current_cell_vector;
  for (unsigned int i = 0; i < 8; i++)
  {
    current_cell_vector.push_back(current_cell.nodes[i]);
  }
  Point cell_avg = average(current_cell_vector);
  for (int i = 0; i < 6; i++)
  {
    std::vector<Point> current_face;
    for (int j = 0; j < 4; j++)
    {
      current_face.push_back(current_cell.nodes[side_nodes_map[i][j]]);
    }
    Point avg = average(current_face);
    Point face_centroids = polygonCentroid(current_face, avg);
    double small_vol = polygonCellVolume(current_face, avg, cell_avg);
    tot_cell_vol += small_vol;
    Point face_contrib = polygonCellCentroid(current_face, avg, cell_avg);
    face_contrib *= small_vol;
    cell_centroid += face_contrib;
  }
  if (tot_cell_vol > 0.)
  {
    cell_centroid /= tot_cell_vol;
  }
  return {cell_centroid, tot_cell_vol};
}

PetscErrorCode find_neighbor(Grid *grid)
{
  PetscFunctionBegin;
  HexNeighborIndices ***neighbor_indices_array;
  int i, j, k, xl, yl, zl, nxl, nyl, nzl;
  Vec global_neighbor_indices;
  PetscCall(DMCreateGlobalVector(grid->dmdas.da6, &global_neighbor_indices));
  PetscCall(DMDAVecGetArray(grid->dmdas.da6, global_neighbor_indices, &neighbor_indices_array));

  PetscCall(DMDAGetCorners(grid->dmdas.da6, &xl, &yl, &zl, &nxl, &nyl, &nzl));

  for (k = 0; k < nzl; k++)
  {
    for (j = 0; j < nyl; j++)
    {
      for (i = 0; i < nxl; i++)
      {
        neighbor_indices_array[k + zl][j + yl][i + xl].indices[0] = cell_index(i + xl, j + yl, k + zl - 1, grid->size);
        neighbor_indices_array[k + zl][j + yl][i + xl].indices[1] = cell_index(i + xl, j + yl - 1, k + zl, grid->size);
        neighbor_indices_array[k + zl][j + yl][i + xl].indices[2] = cell_index(i + xl + 1, j + yl, k + zl, grid->size);
        neighbor_indices_array[k + zl][j + yl][i + xl].indices[3] = cell_index(i + xl, j + yl + 1, k + zl, grid->size);
        neighbor_indices_array[k + zl][j + yl][i + xl].indices[4] = cell_index(i + xl - 1, j + yl, k + zl, grid->size);
        neighbor_indices_array[k + zl][j + yl][i + xl].indices[5] = cell_index(i + xl, j + yl, k + zl + 1, grid->size);
      }
    }
  }
  PetscCall(DMDAVecRestoreArray(grid->dmdas.da6, global_neighbor_indices, &neighbor_indices_array));
  PetscCall(DMGlobalToLocalBegin(grid->dmdas.da6, global_neighbor_indices, INSERT_VALUES, grid->vectors.neighbor_indices));
  PetscCall(DMGlobalToLocalEnd(grid->dmdas.da6, global_neighbor_indices, INSERT_VALUES, grid->vectors.neighbor_indices));
  PetscCall(VecDestroy(&global_neighbor_indices));
  PetscFunctionReturn(PETSC_SUCCESS);
}
