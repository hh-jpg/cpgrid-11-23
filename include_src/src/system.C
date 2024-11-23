/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-03 15:10:01
 * @LastEditTime: 2024-08-29 11:11:22
 * @FilePath: /cpgrid/src/system.C
 * @Description:
 *
 */

#include "system.h"
#include "sparsity.h"
#include "petscdmcp.h"
#include "utils.h"
#include "parameters.h"
System::System(Mesh &mesh) : _mesh(mesh),
                             _dof_map(std::make_unique<DofMap>(mesh)), _para(){}

void System::init()
{
  _dof_map->init();
  _initialized = true;
  create_local_vec(&_old_solution);
}

Mesh &System::get_mesh()
{
  return _mesh;
}

int System::n_local_dofs()
{
  return _mesh.n_local_elem();
}

int System::n_dofs()
{
  return _mesh.getNumberOfPolyhedrons();
}
const std::string &System::name()
{
  return _name;
}
void System::set_name(const std::string_view &name)
{
  std::string sysname(name);
  _name = sysname;
}
int System::n_processors()
{
  return _mesh.n_processors();
}

void System::add_variable(const std::string_view &varname)
{
  _dof_map->add_variable(varname);
}
int System::variable_number(std::string_view var) const
{

  return _dof_map->variable_num(var);
}
int System::n_vars()
{
  return _dof_map->n_vars();
}

DofMap &System::dof_map() const
{
  return *_dof_map;
}

const communicator System::comm() const
{
  return _mesh.comm();
}

PetscErrorCode System::create_global_vec(Vec *globalvec)
{
  if (!_initialized) // 假设 _mesh 有 is_initialized 方法
  {
    throw std::runtime_error("System is not initialized before creat local vector.");
  }

  PetscInt petsc_n_local = static_cast<PetscInt>(_dof_map->n_local_dofs());

  PetscInt petsc_n = static_cast<PetscInt>(_dof_map->n_dofs());
  PetscCall(VecCreate(PETSC_COMM_WORLD, globalvec));
  PetscCall(VecSetSizes(*globalvec, petsc_n_local, petsc_n));
  PetscCall(VecSetType(*globalvec, VECMPI));
  // PetscFunctionReturn(PETSC_SUCCESS);
  return PETSC_SUCCESS; // 更正返回语句
}

PetscErrorCode System::create_local_vec(Vec *localvec) // 更正拼写
{
  if (!_initialized) // 假设 _mesh 有 is_initialized 方法
  {
    throw std::runtime_error("System is not initialized before create local vector."); // 更正拼写
  }

  if (_mesh.n_processors() == 1)
  {
    create_global_vec(localvec); // 更正拼写
  }

  std::vector<unsigned int> send_list = _dof_map->get_send_list();
  PetscInt petsc_n_local = static_cast<PetscInt>(_dof_map->n_local_dofs());
  PetscInt petsc_n = static_cast<PetscInt>(_dof_map->n_dofs());
  PetscInt petsc_n_ghost = static_cast<PetscInt>(send_list.size());
  PetscInt *petsc_ghost = send_list.empty() ? nullptr : const_cast<PetscInt *>(reinterpret_cast<const PetscInt *>(send_list.data()));

  PetscCall(VecCreateGhost(PETSC_COMM_WORLD, petsc_n_local, petsc_n,
                           petsc_n_ghost, petsc_ghost, localvec));
  PetscCall(VecSetFromOptions(*localvec));
  PetscCall(VecSetType(*localvec, VECMPI));

  return PETSC_SUCCESS; // 更正返回语句
}

void transform_preallocation_arrays(const PetscInt blocksize,
                                    const std::vector<int> &n_nz,
                                    const std::vector<int> &n_oz,
                                    std::vector<int> &b_n_nz,
                                    std::vector<int> &b_n_oz)
{
  b_n_nz.clear(); /**/
  b_n_nz.reserve(n_nz.size() / blocksize);
  b_n_oz.clear(); /**/
  b_n_oz.reserve(n_oz.size() / blocksize);

  for (std::size_t nn = 0, nnzs = n_nz.size(); nn < nnzs; nn += blocksize)
  {
    b_n_nz.push_back(n_nz[nn] / blocksize);
    b_n_oz.push_back(n_oz[nn] / blocksize);
  }
}

PetscErrorCode System::create_mat(Mat &mat)
{
  Sparsity sp(*_dof_map);
  sp.creatmat(mat);
}

Vec& System::get_old_solution()
{
  return _old_solution;
}

Vec& System::get_older_solution()
{
  return _older_solution;
}

PetscErrorCode System::set_old_solution(const Vec &old_solution)
{
  PetscErrorCode ierr;
  if (_old_solution == nullptr)
  {
    // 如果 _old_solution 为空，则创建它
    ierr = VecDuplicate(old_solution, &_old_solution);
    if (ierr)
      return ierr;
  }
  // 将 old_solution 复制到 _old_solution
  ierr = VecCopy(old_solution, _old_solution);

  return ierr;
}

PetscErrorCode System::set_older_solution(const Vec &older_solution)
{
  PetscErrorCode ierr;
  if (_older_solution == nullptr)
  {
    // 如果 _old_solution 为空，则创建它
    ierr = VecDuplicate(older_solution, &_older_solution);
    if (ierr)
      return ierr;
  }
  // 将 old_solution 复制到 _old_solution
  ierr = VecCopy(older_solution, _older_solution);
  return ierr;
}

Para &System::para()
{
  return _para;
}