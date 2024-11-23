/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-01 09:41:05
 * @LastEditTime: 2024-08-16 16:31:53
 * @FilePath: /cpgrid/src/dof_map.C
 * @Description:
 *
 */
#include "dof_map.h"
#include <algorithm>
#include <memory>
#include <set>
#include <sstream>
#include <unordered_map>

DofMap::DofMap(Mesh &mesh)
    : _mesh(mesh), _variables()
{
}
Mesh & DofMap::get_mesh() const
{
  return _mesh;
}
void DofMap::init()
{
  distribute_dofs();
  prepare_send_list();
  build_ghosted_to_local_map();
}
const std::vector<unsigned int> &DofMap::get_send_list() const
{
  return _send_list;
}

void DofMap::build_ghosted_to_local_map()
{
  int ghosted_size = _send_list.size();
  int proc_id = processor_id();
  for (int i = 0; i < ghosted_size; i++)
  {
    _ghosted_to_local_map[_send_list[i]] = i + (_end_dof[proc_id] - _first_dof[proc_id]);
  }
}

int DofMap::processor_id() const
{
  return _mesh.processor_id();
}

communicator DofMap::comm() const
{
  return _mesh.comm();
}

int DofMap::n_local_dofs() const
{
  return n_vars() * (_mesh.n_local_elem());
}

int DofMap::n_dofs() const
{
  return n_vars() * _mesh.getNumberOfPolyhedrons();
}

int DofMap::n_processors() const
{
  return _mesh.n_processors();
}
void DofMap::set_dofarrange_type(DofArrangementType type)
{
  _dof_arrange_type = type;
}
void DofMap::add_variable(const std::string_view &var_name)
{
  std::string varname(var_name);
  _variables.add_variable(varname);
}
void DofMap::clear_send_list()
{
  _send_list.clear();
}

DofArrangementType DofMap::dof_arrange_type() const
{
  return _dof_arrange_type;
}

Variables DofMap::variables() const
{
  return _variables;
}

std::string DofMap::variable_name(const int vars_num)
{
  return _variables.variable_name(vars_num);
}

int DofMap::variable_num(std::string_view var) const
{
  return _variables.variable_num(var);
}

int DofMap::n_vars() const
{
  return _variables.size();
}
bool DofMap::distributed()
{
  return _distributed;
}
int DofMap::first_dof()
{
  int processor_id = _mesh.processor_id();
  return _first_dof[processor_id];
}

int DofMap::end_dof()
{
  int processor_id = _mesh.processor_id();
  return _end_dof[processor_id];
}

void DofMap::prepare_send_list() // var dof send_list
{
  clear_send_list();
  int n_var = _variables.size();
  // 如果只有一个处理器，直接返回
  if (this->n_processors() == 1)
    return;
  if (!_mesh.partitioned())
    std::cout << "prepare_send_list need the mesh partitioned " << std::endl;
  for (const auto &elem : _mesh.elem_local_range())
  {
    const auto &faces = elem.get_faces();
    for (const Face &f : faces)
    {
      Polyhedron *neighbor = f.neighbor();
      if (neighbor == nullptr)
        continue;
      int processor_id = _mesh.processor_id();
      // 如果这个单元不属于当前处理器
      for (int i = 0; i < n_var; i++)
      {
        int neighbor_dof = dof_indices(*neighbor, i);

        if (neighbor_dof < _first_dof[processor_id] || neighbor_dof >= _end_dof[processor_id])
        {
          _send_list.push_back(neighbor_dof);
        }
      }
    }
  }

  std::sort(_send_list.begin(), _send_list.end());
  // Now use std::unique to remove duplicate entries
  std::vector<unsigned int>::iterator new_end =
      std::unique(_send_list.begin(), _send_list.end());

  std::vector<unsigned int>(_send_list.begin(), new_end).swap(_send_list);
}

// this function return global dof index of the vn-th var in elem e
int DofMap::dof_indices(const Polyhedron &elem, const int vn)
{
  if (!distributed())
    throw std::invalid_argument("DOF must be distributed before dof_indices");

  int elem_id = elem.id();
  int proc_id = elem.processor_id();
  int first_dof = _first_dof[proc_id];
  int end_dof = _end_dof[proc_id];
  int num_proc_elem = (end_dof - first_dof) / n_vars();
  int first_elem_proc = first_dof / n_vars();

  switch (_dof_arrange_type)
  {
  case PVE:
    return first_dof + vn * num_proc_elem + elem_id - first_elem_proc;
  case PEV:
    return elem_id * n_vars() + vn;
  default:
    throw std::runtime_error("Undefined DOF arrangement type.");
  }
}

int DofMap::dof_local_indices(const Polyhedron &elem, const int vn)
{
  // 获取全局索引
  int global_index = dof_indices(elem, vn);
  int proc_id = _mesh.processor_id();

  // 检查全局索引是否在有效范围内
  if (global_index >= _first_dof[proc_id] && global_index < _end_dof[proc_id])
  {
    // 计算并返回局部索引
    return global_index - _first_dof[proc_id];
  }

  // 查找全局索引是否在 ghosted to local map 中
  auto it = _ghosted_to_local_map.find(global_index);
  if (it != _ghosted_to_local_map.end())
  {
    // 返回找到的局部索引
    return it->second;
  }
  else
  {
    // 处理未找到的情况，抛出异常或返回错误码
    throw std::runtime_error("dof index not found in local ");
    // 或者返回一个特定的错误码，例如 -1
    // return -1;
  }
}

void DofMap::distribute_dofs()
{
  int n_elems = _mesh.getNumberOfPolyhedrons();
  int n_local_elems = _mesh.n_local_elem();
  int n_processors = _mesh.n_processors();
  std::vector<int> n_dofs_on_proc(n_processors, 0);
  int n_vars = this->n_vars();
  int n_local_dofs = n_vars * n_local_elems;
  MPI_Allgather(&n_local_dofs, 1, MPI_INT, n_dofs_on_proc.data(), 1, MPI_INT, PETSC_COMM_WORLD);
  _first_dof.resize(n_processors);
  _end_dof.resize(n_processors);
  _first_dof[0] = 0;
  for (int i = 1; i < n_processors; ++i)
    _first_dof[i] = _end_dof[i - 1] = _first_dof[i - 1] + n_dofs_on_proc[i - 1];
  _end_dof[n_processors - 1] = _first_dof[n_processors - 1] + n_dofs_on_proc[n_processors - 1];
  _distributed = true;
}

std::unordered_map<int, int> &DofMap::ghosted_to_local_map()
{
  return _ghosted_to_local_map;
}

bool DofMap::is_in_local_range(PetscInt index)
{
  int proc_id = processor_id();
  return _first_dof[proc_id] <= index && index < _end_dof[proc_id];
}