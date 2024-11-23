/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-12 15:59:52
 * @LastEditTime: 2024-08-28 10:18:07
 * @FilePath: /cpgrid/include/polyhedron.h
 * @Description:
 *
 */
#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <vector>
#include <iostream>
#include "face.h"
#include "config.h"

class Face;
class Polyhedron
{
public:
  // 参数构造函数
  Polyhedron(const std::vector<Face> &faces, unsigned int id);

  Polyhedron() = default; // 默认构造函数

  // 获取私有变量的公共成员函数
  std::vector<Face> &get_faces(); // need to modify the faces, not const
  const std::vector<Face> &get_faces() const;
  unsigned int n_faces() const;
  unsigned int id() const;
  double perm() const;
  // 计算相关函数
  Point compute_cell_centroid() const;
  Point compute_vertex_avg() const;
  double compute_cell_volume() const;
  // 设置相关函数
  void set_elem_perm(double perm);
  inline int &processor_id();
  inline const int &processor_id() const;
  inline void set_id(unsigned int id);
  // 重载输出运算符
  friend std::ostream &operator<<(std::ostream &os, const Polyhedron &poly);

private:
  unsigned int _id = 0;
  std::vector<Face> _faces;
  int _processor_id = 0;
  double _perm = 0.0; // 初始化成员变量
};

inline void Polyhedron::set_id(unsigned int id)
{
  _id = id;
}
inline int &Polyhedron::processor_id()
{
  return _processor_id;
}
inline const int &Polyhedron::processor_id() const
{
  return _processor_id;
}
#endif // POLYHEDRON_H
