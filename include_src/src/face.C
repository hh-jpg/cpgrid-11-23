/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-12 16:04:21
 * @LastEditTime: 2024-08-28 21:59:29
 * @FilePath: /cpgrid/src/face.C
 * @Description:
 *
 */
#include "face.h"
#include "polyhedron.h"
#include "geo.h"
#include "utils.h"
// 构造函数
Face::Face(std::vector<Point *> points, unsigned int pos, unsigned int id, Polyhedron *elem, Polyhedron *neighbor)
    : _points(points), _pos(pos), _id(id), _elem(elem), _neighbor(neighbor) {}

Point Face::compute_face_centroid() const
{
#if MESH_DIM == 2
  if (_points.size() != 2)
  {
    throw std::runtime_error("In 2D mesh, the face must be a line segment with exactly 2 points.");
  }
  // 计算线段的中点
  double x = 0.5 * ((*_points[0]).x() + (*_points[1]).x());
  double y = 0.5 * ((*_points[0]).y() + (*_points[1]).y());
  return Point(x, y);
#else
  Point avg = average(_points);
  return polygonCentroid(_points, avg);
#endif
}

double Face::compute_face_area() const
{
#if MESH_DIM == 2
  if (_points.size() != 2)
  {
    throw std::runtime_error("In 2D mesh, the face must be a line segment with exactly 2 points.");
  }
  // 计算线段的长度
  double dx = (*_points[1]).x() - (*_points[0]).x();
  double dy = (*_points[1]).y() - (*_points[0]).y();
  return std::sqrt(dx * dx + dy * dy);
#else
  Point avg = average(_points);
  return polygonArea(_points, avg);
#endif
}

Point Face::compute_face_avg() const
{
  return average(_points);
}

Point Face::compute_face_normal() const
{
// do not guarantee that the normal vector is outward
//  we next compute trans set cc2fc.dot(fn) to positive.
#if MESH_DIM == 2
  if (_points.size() != 2)
  {
    throw std::runtime_error("In 2D mesh, the face must be a line segment with exactly 2 points.");
  }
  double dx = (*_points[1]).x() - (*_points[0]).x();
  double dy = (*_points[1]).y() - (*_points[0]).y();
  double length = std::sqrt(dx * dx + dy * dy);
  if (length == 0)
  {
    throw std::runtime_error("Zero-length line segment is not valid.");
  }
  // 计算单位法向量
  Point normal(-dy / length, dx / length);
  return normal;
#else
  if (_points.size() < 3)
  {
    throw std::runtime_error("In 3D mesh, the face must be a polygon with at least 3 points.");
  }
  Point avg = average(_points);
  return polygonNormal(_points, avg);
#endif
}

Point Face::compute_add_cell(const Polyhedron &elem) const // cc is the cellcentroid;
{
  Point cc = elem.compute_cell_centroid();
  Point fc = compute_face_centroid();
  Point ac = fc + fc - cc;
  // PetscPrintf(PETSC_COMM_WORLD,"cc=(%.2f , %.2f), fc=(%.2f , %.2f), ac=(%.2f , %.2f)\n",cc.x(),cc.y(),fc.x(),fc.y(),ac.x(),ac.y());
  return ac;
}

double Face::compute_half_trans(const Polyhedron &elem) const // cc is the cellcentroid;
{
  Point cc = elem.compute_cell_centroid();
  Point fc = compute_face_centroid();
  Point cc2fc = fc - cc;
  double farea = compute_face_area();
  Point fn = compute_face_normal();
  double perm = elem.perm();
  double hT = farea * perm * cc2fc.dot(fn) / (cc2fc.two_norm() * cc2fc.two_norm());
  return std::abs(hT);
}

double Face::compute_trans() const
{
  double cT = compute_half_trans(*_elem); // Compute half transmissibility of current element

  if (_neighbor == nullptr)
  {
    return cT; // No neighbor, return current element transmissibility
  }

  // Compute transmissibility using the harmonic mean
  double nT = compute_half_trans(*_neighbor);
  return half_harmonicMean(cT, nT);
}

double Face::compute_efftrans(const PetscScalar mob_up) const
{
  double T = compute_trans();
  return mob_up * T;
}

std::vector<Point *> Face::get_points() const
{
  return _points;
}

unsigned short Face::pos() const
{
  return _pos;
}

unsigned int Face::id() const
{
  return _id;
}

Polyhedron *Face::neighbor() const
{
  return _neighbor;
}

Polyhedron *Face::elem() const
{
  return _elem;
}
unsigned int Face::n_nodes() const
{
  return _points.size();
}

unsigned short Face::convertpos()
{
  unsigned int convert[6] = {4, 2, 1, 3, 0, 5};
  return convert[_pos];
}
void Face::set_neighbor(Polyhedron *n)
{
  _neighbor = n;
}
void Face::set_elem(Polyhedron *c)
{
  _elem = c;
}
// 将 Face 对象输出到流
std::ostream &operator<<(std::ostream &os, const Face &face)
{
  os << "Face ID: " << face._id << "\n";
  os << "Position: " << face._pos << "\n";
  os << " ";
  os << "Neighbor ID: ";

  if (face._neighbor == nullptr)
  {
    os << "this face is boundary face \n";
  }
  else
  {
    os << face._neighbor->id() << "\n";
  }

  os << "Points: ";
  for (const auto &point : face._points)
  {
    os << *point << " ";
  }
  os << "\n";

  return os;
}
