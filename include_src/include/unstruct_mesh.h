/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-11 11:01:14
 * @LastEditTime: 2024-08-28 10:09:04
 * @FilePath: /cpgrid/include/unstruct_mesh.h
 * @Description:  arbitrary polygonal elements
 *
 */
#ifndef UNSTRUCT_MESH_H
#define UNSTRUCT_MESH_H
#include <vector>
#include <memory>
#include "point.h"
#include "face.h"
#include "polyhedron.h"
//#include "metis_partitioner.h"
#include <array>
#include <petscsys.h>
#include <iostream>
#include <string>
#include "utils.h"
#include "elem_iterator.h"
#include "config.h"

// 假设 Point 和 Polyhedron 的定义在其他地方
class Mesh
{
public:
  Mesh() = default; // 默认构造函数

  explicit Mesh(const std::string &fpath); // 带文件路径参数的构造函数

  // 读取文件函数
  void read(const std::string &fpath);

  // 输出型函数
  void printInfo() const;                                   // 打印网格信息
  void writeNodesToFile(const std::string &filename) const; // 将节点数据写入文件
  void writeElemsToFile(const std::string &filename) const; // 将单元数据写入文件
  // 其他有用的函数
  size_t getNumberOfNodes() const;          // 获取节点数
  size_t getNumberOfPolyhedrons() const;    // 获取多面体数
  size_t getNumberOfFaces() const;          // 内部面只算一次,from mrst
  const Point &getNode(size_t index) const; // 根据索引获取节点
  Polyhedron &getPolyhedron(size_t index);
  void init();                               // this just set the mesh topology
  
  void prepare_for_use();                    // partition and reset elem id
  void partition();
  // 设置渗透率的函数
  void read_set_perms();
  void setFaceNeighbors();                   // 设置邻居关系
  const communicator &comm() const
  {
    return PETSC_COMM_WORLD;
  }
  int n_processors() const;
  int processor_id() const;
  unsigned int &set_n_partitions()
  {
    return _n_parts;
  }
  auto elem_range() const
  {
    return FilterRange(_polyhedras.begin(), _polyhedras.end(), is_global_range, 0);
  }

  auto elem_processor_range(int processor_id)
  {
    return FilterRange(_polyhedras.begin(), _polyhedras.end(), is_processor_range, processor_id);
  }
  auto elem_local_range() const
  {
    if (!_partitioned)
    {
      throw std::runtime_error("Mesh must be partitioned before calling elem_local_range.");
    }
    int local_processor_id = processor_id();
    return FilterRange(_polyhedras.begin(), _polyhedras.end(), is_processor_range, local_processor_id);
  }
  int &n_local_elem()
  {
    return _n_local_elem;
  }
  bool partitioned();
  void reset_elem_id(); // after the partition,we need reset the elem_id, need to distribute for the dofs;
private:
  void read_nodes(); // 读取节点信息
  void read_elems();
  void addPoint(const Point &p, int id = -1); // 添加节点
  bool _partitioned = false;                  // 读取单元信息
  communicator _comm;
  int _n_local_elem = 0; // set use the n_local_elem() in the partition;
  unsigned int _n_parts = 0;
  std::string filepath;                // 文件路径
  std::vector<Point> _nodes;           // 节点集合
  std::vector<Polyhedron> _polyhedras; // 多面体集合
  static bool is_global_range(const Polyhedron &e, int processor_id)
  {
    return true;
  }

  static bool is_processor_range(const Polyhedron &e, int processor_id)
  {
    return e.processor_id() == processor_id;
  }
};

// 内联函数的定义应放在头文件中
inline size_t Mesh::getNumberOfNodes() const
{
  return _nodes.size();
}
inline bool Mesh::partitioned()
{
  return this->_partitioned;
}

inline int Mesh::n_processors() const
{
  int size;
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  return size;
}
inline int Mesh::processor_id() const
{
  int rank;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  return rank;
}
inline size_t Mesh::getNumberOfPolyhedrons() const
{
  return _polyhedras.size();
}

const inline Point &Mesh::getNode(size_t index) const
{
  return _nodes.at(index);
}

inline Polyhedron &Mesh::getPolyhedron(size_t index)
{
  if (index < _polyhedras.size())
  {
    return _polyhedras[index];
  }
  throw std::out_of_range("Polyhedron index out of range");
}

template <typename T>
std::vector<std::vector<T>> readCSV(const std::string &filename);

#endif