/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-12 15:55:30
 * @LastEditTime: 2024-08-28 15:16:28
 * @FilePath: /cpgrid/include/face.h
 * @Description:
 *
 */
#ifndef FACE_H
#define FACE_H

#include <vector>
#include <iostream>
#include "point.h"
#include <petscsys.h>
#include "config.h"

class Polyhedron;
class Face
{
public:
  unsigned short convertpos(); // convert to the right hand coord

  // 优化的构造函数
  Face(std::vector<Point *> points, unsigned int pos,
       unsigned int id, Polyhedron *elem, Polyhedron *neighbor);

  Face() = default; // 默认构造函数
  unsigned int n_nodes() const;

  // 获取私有变量的公共成员函数
  std::vector<Point *> get_points() const;
  unsigned short pos() const;//we now use the mrst position - 1;
  unsigned int id() const;
  Point compute_face_centroid() const;
  Point compute_face_avg() const;
  Point compute_face_normal() const;
  double compute_face_area() const;
  Point compute_add_cell(const Polyhedron &elem) const;
  double compute_half_trans(const Polyhedron &elem) const;
  double compute_trans() const;
  double compute_efftrans(const PetscScalar mob_up) const;
  // face的两侧是elem a 和elem b，在a中返回的是b，在b中返回的是a
  Polyhedron *neighbor() const;
  Polyhedron *elem() const;
  void set_neighbor(Polyhedron *n);
  void set_elem(Polyhedron *c);
  friend std::ostream &operator<<(std::ostream &os, const Face &face);
  
private:
  std::vector<Point *> _points;
  unsigned short _pos = 0; // this is the mrst pos - 1
  unsigned int _id = 0;
  int _dim = 2;
  Polyhedron *_neighbor = nullptr;
  Polyhedron *_elem = nullptr;
};

#endif // FACE_H
