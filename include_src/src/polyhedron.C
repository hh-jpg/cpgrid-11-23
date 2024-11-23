/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-12 15:59:52
 * @LastEditTime: 2024-08-28 17:27:07
 * @FilePath: /cpgrid/src/polyhedron.C
 * @Description:
 *
 */

#include "polyhedron.h"
#include "geo.h"
#include "utils.h"
// 参数构造函数
Polyhedron::Polyhedron(const std::vector<Face> &faces, unsigned int id)
    : _faces(faces), _id(id) {}

std::vector<Face> &Polyhedron::get_faces()
{
  return _faces;
}

const std::vector<Face> &Polyhedron::get_faces() const
{
  return _faces;
}
unsigned int Polyhedron::id() const
{
  return _id;
}

double Polyhedron::perm() const
{
  return _perm;
}

unsigned int Polyhedron::n_faces() const
{
  return _faces.size();
}

Point Polyhedron::compute_vertex_avg() const
{
  Point cell_avg(0.0, 0.0, 0.0);
  unsigned int fsize = _faces.size();
  for (unsigned int i = 0; i < fsize; i++)
  {
    cell_avg += _faces[i].compute_face_avg() / fsize;
  }
  return cell_avg;
}

Point Polyhedron::compute_cell_centroid() const
{
#if MESH_DIM == 2
  double tot_area = 0.0;
  Point tot_centroid;
  Point avg = compute_vertex_avg();
  unsigned int fsize = _faces.size();
  for (unsigned int i = 0; i < fsize; i++)
  {
    Face current_face = _faces[i];
    Point *p1 = current_face.get_points()[0];
    Point *p2 = current_face.get_points()[1];
    Point tri[3] = {avg, *p1, *p2};
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
    tot_centroid = avg;
  }
  return tot_centroid;
#else
  double tot_cell_vol = 0.0;
  Point cell_centroid(0.0, 0.0, 0.0);
  Point cell_avg = compute_vertex_avg();
  unsigned int fsize = _faces.size();
  for (unsigned int i = 0; i < fsize; i++)
  {
    Face current_face = _faces[i];
    Point fc = current_face.compute_face_centroid();
    Point favg = current_face.compute_face_avg();
    std::vector<Point *> face_points = current_face.get_points();
    double small_vol = polygonCellVolume(face_points, favg, cell_avg);
    tot_cell_vol += small_vol;
    Point face_contrib = polygonCellCentroid(face_points, favg, cell_avg);
    face_contrib *= small_vol;
    cell_centroid += face_contrib;
  }
  if (tot_cell_vol > 0.)
  {
    cell_centroid /= tot_cell_vol;
  }
  return cell_centroid;
#endif
}

double Polyhedron::compute_cell_volume() const
{
#if MESH_DIM == 2
  double tot_area = 0.0;
  Point cell_centroid(0.0, 0.0, 0.0);
  Point cell_avg = compute_vertex_avg();
  unsigned int fsize = _faces.size();
  for (unsigned int i = 0; i < fsize; i++)
  {
    Face current_face = _faces[i];
    Point *p1 = current_face.get_points()[0];
    Point *p2 = current_face.get_points()[1];
    Point tri[3] = {cell_avg, *p1, *p2};
    tot_area += triarea(tri);
  }
  return tot_area;
#else
  double tot_cell_vol = 0.0;
  Point cell_centroid(0.0, 0.0, 0.0);
  Point cell_avg = compute_vertex_avg();
  unsigned int fsize = _faces.size();
  for (unsigned int i = 0; i < fsize; i++)
  {
    Face current_face = _faces[i];
    Point fc = current_face.compute_face_centroid();
    Point favg = current_face.compute_face_avg();
    std::vector<Point *> face_points = current_face.get_points();
    double small_vol = polygonCellVolume(face_points, favg, cell_avg);
    tot_cell_vol += small_vol;
  }
  return tot_cell_vol;

#endif
}

void Polyhedron::set_elem_perm(double perm)
{
  _perm = perm;
}

// 将 Polyhedron 对象输出到流
std::ostream &operator<<(std::ostream &os, const Polyhedron &poly)
{

  os << "Polyhedron ID: " << poly._id << "\nFaces:\n";

  for (const auto &face : poly._faces)
  {
    os << face << "\n";
  }
  return os;
}
