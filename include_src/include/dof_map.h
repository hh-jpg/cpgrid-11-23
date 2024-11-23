/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-06-28 15:09:39
 * @LastEditTime: 2024-08-28 10:16:12
 * @FilePath: /cpgrid/include/dof_map.h
 * @Description:
 *
 */
#ifndef DOFMAP_H
#define DOFMAP_H

#include <vector>
#include <string>
#include <unordered_map>
#include <petscsys.h>
#include <petsc.h>
#include "unstruct_mesh.h"
#include "variables.h"
#include "utils.h"
#include "config.h"
class DofMap
{
public:
  DofMap(Mesh &mesh);
  void init();
  const std::vector<unsigned int> &get_send_list() const;

  communicator comm() const;

  int n_local_dofs() const;
  int n_dofs() const;
  int n_processors() const;
  void distribute_dofs();
  void add_variable(const std::string_view &varname);
  void prepare_send_list(); // var dof send_list
  int dof_indices(const Polyhedron &elem, const int vn);
  int dof_local_indices(const Polyhedron &elem, const int vn);
  void clear_send_list();
  void build_ghosted_to_local_map();
  int processor_id() const;
  std::string variable_name(const int vars_num);
  int variable_num (std::string_view var) const;
  DofArrangementType dof_arrange_type() const;
  Variables variables() const;
  int n_vars() const;
  int first_dof();
  int end_dof();
  bool distributed();
  void set_dofarrange_type(DofArrangementType type);// this need to before init();
  std::unordered_map<int, int> &  ghosted_to_local_map();
  bool is_in_local_range(PetscInt index);
  Mesh & get_mesh() const;
private:
  Mesh &_mesh;
  DofArrangementType _dof_arrange_type = PVE;
  std::vector<unsigned int> _send_list;
  std::unordered_map<int, int> _ghosted_to_local_map;
   /**
   * First DOF index on processor \p p.
   */
  std::vector<unsigned int> _first_dof;
  std::vector<unsigned int> _end_dof; //[_first_position, _end_position);
  Variables _variables;
  bool _distributed = false;
};

#endif
