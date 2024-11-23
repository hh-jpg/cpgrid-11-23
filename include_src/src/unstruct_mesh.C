/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-11 15:01:13
 * @LastEditTime: 2024-08-28 15:30:47
 * @FilePath: /cpgrid/src/unstruct_mesh.C
 * @Description:
 *
 */
#include "unstruct_mesh.h"
#include "point.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <petscsys.h>
#include "metis_csr_graph.h"
#include "metis_partitioner.h"
Mesh::Mesh(const std::string &fpath) : filepath(fpath)
{
  read(fpath); // 在构造函数中读取文件
}

void Mesh::read(const std::string &fpath)
{
  // 确保 fpath 以 '/' 结尾
  if (!fpath.empty() && fpath.back() != '/')
  {
    this->filepath = fpath + '/';
  }
  else
  {
    this->filepath = fpath;
  }
  this->filepath = fpath;
  this->read_nodes();

  this->read_elems();
}

void Mesh::read_nodes()
{
  std::string node_coords = filepath + "/G/nodes/coords/data.csv";
  std::vector<std::vector<double>> mrst_point;
  mrst_point = readCSV<double>(node_coords);
  for (unsigned int i = 0; i < mrst_point.size(); i++)
  {
    const unsigned int node_id = i;
    double x = 0, y = 0, z = 0;
    x = mrst_point[i][0];
    y = mrst_point[i][1];
#if MESH_DIM == 3
    z = mrst_point[i][2];
#endif
    this->addPoint(Point(x, y, z), node_id);
  }
}

void Mesh::read_elems()
{
  std::string mrst_cells_index_path = filepath + "/G/cells/indexMap/data.csv";
  std::vector<std::vector<int>> mrst_elem_ids = readCSV<int>(mrst_cells_index_path);

  std::string mrst_cells_faces_path = filepath + "/G/cells/faces/data.csv";
  std::vector<std::vector<int>> mrst_cell_faces = readCSV<int>(mrst_cells_faces_path);

  std::string mrst_cells_facepos_path = filepath + "/G/cells/facePos/data.csv";
  std::vector<std::vector<int>> mrst_cell_facepos = readCSV<int>(mrst_cells_facepos_path);

  std::string mrst_faces_nodes_path = filepath + "/G/faces/nodes/data.csv";
  std::vector<std::vector<int>> mrst_faces_nodes = readCSV<int>(mrst_faces_nodes_path);

  std::string mrst_faces_nodepos_path = filepath + "/G/faces/nodePos/data.csv";
  std::vector<std::vector<int>> mrst_faces_nodepos = readCSV<int>(mrst_faces_nodepos_path);
  // std::string mrst_faces_neighbor_path = filepath + "/G/faces/neighbor/data.csv";
  // std::vector<std::vector<unsigned int>> mrst_faces_neighbor = readCSV<int>(mrst_faces_neighbor_path);

  _polyhedras.resize(mrst_elem_ids.size());
  for (unsigned int i = 0; i < mrst_elem_ids.size(); i++)
  {
    unsigned int n_faces = mrst_cell_facepos[i + 1][0] - mrst_cell_facepos[i][0];

    std::vector<Face> faces(n_faces);
    for (unsigned int j = 0; j < n_faces; j++)
    {
      unsigned int face_index = mrst_cell_faces[mrst_cell_facepos[i][0] + j][0];
      unsigned int face_pos = mrst_cell_faces[mrst_cell_facepos[i][0] + j][1];
      unsigned int n_nodes = mrst_faces_nodepos[face_index + 1][0] - mrst_faces_nodepos[face_index][0];

      std::vector<Point *> points(n_nodes);
      for (unsigned int k = 0; k < n_nodes; k++)
      {
        unsigned int node_index = mrst_faces_nodes[mrst_faces_nodepos[face_index][0] + k][0];
        points[k] = &_nodes[node_index];
      }
      Face face(points, face_pos, face_index, nullptr, nullptr); // the neighbor will be set after
      faces[j] = face;
    }
    Polyhedron cell(faces, i);
    _polyhedras[i] = cell;
  }
}

size_t Mesh::getNumberOfFaces() const
{
  std::string mrst_faces_num_path = filepath + "/G/faces/num/data.csv";
  std::vector<std::vector<int>> mrst_faces_num = readCSV<int>(mrst_faces_num_path);
  return mrst_faces_num[0][0] + 1; // readCSV时自动会减一；
}

void Mesh::setFaceNeighbors()
{
  // 读取邻居关系的CSV文件
  std::string mrst_faces_neighbor_path = filepath + "/G/faces/neighbors/data.csv";
  std::vector<std::vector<int>> mrst_faces_neighbor = readCSV<int>(mrst_faces_neighbor_path);

  // 遍历所有多面体
  for (Polyhedron &current_elem : _polyhedras)
  {
    unsigned int elem_id = current_elem.id();
    std::vector<Face> &faces = current_elem.get_faces();

    // 遍历当前多面体的所有面
    for (Face &face : faces)
    {
      face.set_elem(&current_elem);
      unsigned int face_id = face.id();
      const auto &neighbors = mrst_faces_neighbor[face_id];
      if (neighbors[0] == elem_id)
      {
        unsigned int neighbor_id = neighbors[1];
        if (neighbor_id == -1)
          face.set_neighbor(nullptr);
        else
          face.set_neighbor(&_polyhedras[neighbor_id]);
      }
      else if (neighbors[1] == elem_id)
      {
        unsigned int neighbor_id = neighbors[0];
        if (neighbor_id == -1)
          face.set_neighbor(nullptr);
        else
          face.set_neighbor(&_polyhedras[neighbor_id]);
      }
      else
      {
        std::cerr << "Error: Wrong neighbor set for face ID " << face_id << " in element ID " << elem_id << std::endl;
      }
    }
  }
}

void Mesh::addPoint(const Point &p, int id)
{
  if (id != -1)
  {
    if (id < _nodes.size())
      _nodes[id] = p;
    else
    {
      _nodes.resize(id + 1);
      _nodes[id] = p;
    }
  }
  else
    _nodes.push_back(p);
}

// 输出网格信息的函数实现
void Mesh::printInfo() const
{
  PetscPrintf(PETSC_COMM_WORLD, "Mesh Information:\n");
  PetscPrintf(PETSC_COMM_WORLD, "Filepath: %s\n", filepath.c_str());
  PetscPrintf(PETSC_COMM_WORLD, "Number of nodes: %zu\n", _nodes.size());
  PetscPrintf(PETSC_COMM_WORLD, "Number of polyhedras: %zu\n", _polyhedras.size());
  PetscPrintf(PETSC_COMM_WORLD, "The mesh dimension is %d\n", MESH_DIM);
}

// 将节点数据写入文件的函数实现
void Mesh::writeNodesToFile(const std::string &filename) const
{
  std::ofstream file(filename);
  if (!file.is_open())
  {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }
  for (const auto &node : _nodes)
  {
    // 假设 Point 类有适当的输出运算符
    file << node << std::endl;
  }
  file.close();
}

// 将单元数据写入文件的函数实现
void Mesh::writeElemsToFile(const std::string &filename) const
{
  std::ofstream file(filename);

  if (!file.is_open())
  {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }
  for (const auto &polyhedron : _polyhedras)
  {
    // 假设 Polyhedron 类有适当的输出运算符
    file << polyhedron << std::endl;
  }
  file.close();
}

void Mesh::init()
{
  read_set_perms();
  setFaceNeighbors();
}
// 读取matlab数据，其中整型-1,因为matlab 从1开始，c++从0开始
template <typename T>
std::vector<std::vector<T>> readCSV(const std::string &filename)
{
  std::ifstream file(filename);
  std::vector<std::vector<T>> data;
  if (!file.is_open())
  {
    std::cerr << "无法打开文件: " << filename << std::endl;
    return data; // 返回一个空的向量
  }
  std::string line;
  while (std::getline(file, line))
  {
    std::vector<T> row;
    std::stringstream lineStream(line);
    std::string value;

    while (std::getline(lineStream, value, ','))
    {
      try
      {
        if constexpr (std::is_same_v<T, int>)
        {
          row.push_back(static_cast<int>(std::stod(value) - 1)); // MATLAB 从 1 开始,我们遵守c++从0开始的规则
        }
        else if constexpr (std::is_same_v<T, double>)
        {
          row.push_back(static_cast<double>(std::stod(value)));
        }
        else
        {
          throw std::invalid_argument("Unsupported data type for conversion");
        }
      }
      catch (const std::invalid_argument &e)
      {
        std::cerr << "无效的数据值: " << value << std::endl;
      }
    }
    data.push_back(row);
  }
  file.close();
  return data;
}

void Mesh::reset_elem_id() // 也许能多线程进行优化
{
  if (!_partitioned)
  {
    throw std::runtime_error("Mesh must be partitioned before calling reset_elem_id.");
  }
  int current_id = 0;
  for (int p_id = 0; p_id < n_processors(); p_id++)
  {
    for (Polyhedron &elem : elem_processor_range(p_id))
    {
      elem.set_id(current_id);
      current_id++;
    }
  }
}
void Mesh::partition()
{
  METISPartitioner partitioner;
  unsigned int n_pieces = n_processors();
  METIS_CSR_Graph csr_graph = CSRGraphFactory::createCSRGraph(*this);
  partitioner.partition(*this, csr_graph, n_pieces);
  _partitioned = true;
}
void Mesh::prepare_for_use()
{
  partition();
  reset_elem_id();
}
void Mesh::read_set_perms()
{
  const std::string perms_path = filepath + "/rock/perm/data.csv";
  std::vector<std::vector<double>> mrst_perms;
  mrst_perms = readCSV<double>(perms_path);
  for (int i = 0; i < _polyhedras.size(); i++)
  {
    double perm = mrst_perms[i][0];
    _polyhedras[i].set_elem_perm(perm);
  }
}