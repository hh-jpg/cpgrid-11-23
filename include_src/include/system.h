/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-28 14:39:46
 * @LastEditTime: 2024-08-29 11:11:12
 * @FilePath: /cpgrid/include/system.h
 * @Description:
 *
 */

#ifndef SYSTEM_H
#define SYSTEM_H

#include "unstruct_mesh.h"
#include <vector>
#include <string>
#include <map>
#include <memory>
#include "petscsys.h"
#include "petsc.h"
#include "utils.h"
#include "dof_map.h"
#include "parameters.h"
#include "config.h"
class System {
public:
  System(Mesh &mesh);
  System(const System &)            = delete; // 禁用拷贝构造函数
  System &operator=(const System &) = delete; // 禁用拷贝赋值操作符
  // 添加一个移动构造函数和移动赋值操作符，如果需要
  System(System &&) noexcept            = default;
  System &operator=(System &&) noexcept = default;
  ~System()
  {
    if (_old_solution) {
      VecDestroy(&_old_solution); // Destroy the old solution vector
    }

    if (_older_solution) {
      VecDestroy(&_older_solution); // Destroy the older solution vector
    }
  }
  void  set_name(const std::string_view &name);
  void  init();
  Mesh &get_mesh();
  int   n_local_dofs();
  int   n_dofs();
  int   n_processors();
  void  add_variable(const std::string_view &varname);
  int   variable_number(std::string_view var) const;

  int                n_vars();
  DofMap            &dof_map() const;
  const communicator comm() const;
  PetscErrorCode     create_global_vec(Vec *globalvec);
  PetscErrorCode     create_local_vec(Vec *localvec);

  PetscErrorCode     create_mat(Mat &mat);
  const std::string &name();
  Vec               &get_old_solution();
  Vec               &get_older_solution();
  PetscErrorCode     set_old_solution(const Vec &old_solution);
  PetscErrorCode     set_older_solution(const Vec &older_solution);
  Para              &para();
  double             time = 0.0;

private:
  Mesh                   &_mesh;
  std::unique_ptr<DofMap> _dof_map;
  std::string             _name        = "cp_system";
  bool                    _initialized = false;
  Vec                     _old_solution;
  Vec                     _older_solution;
  Para                    _para;
};

#endif // SYSTEM_H
