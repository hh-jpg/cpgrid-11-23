/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-23 14:22:28
 * @LastEditTime: 2024-09-01 16:53:22
 * @FilePath: /cpgrid/example/twophase/convergence_test_1/init.C
 * @Description:
 *
 */
#include <petscvec.h>
#include <petscdm.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscdmtypes.h>
#include "utils.h"
#include "petscdmcp.h"
#include "vec_func.h"
#include "system.h"
#include "dof_map.h"
#include "ghosted_array_wrapper.h"
#include "fluid.h"
#include "mat_func.h"
#include "case_fluid.h"

PetscErrorCode set_para(Para &para)
{
  PetscFunctionBeginUser;
  // Adding simulation parameters with default values the final value is on the option file
  para.add_parameter<PetscScalar>("-dt", 86400.0);
  para.add_parameter<PetscInt>("-n_steps", 10);
  para.add_parameter<PetscScalar>("-phi", 0.2);
  para.add_parameter<PetscScalar>("-beta", 2.0);
  para.add_parameter<PetscScalar>("-mu_w", 0.002);
  para.add_parameter<PetscScalar>("-mu_o", 0.003);
  para.add_parameter<PetscScalar>("-srw", 0.0);
  para.add_parameter<PetscScalar>("-srn", 0.2);

  para.add_parameter<PetscScalar>("-rho_w", 1000.);
  para.add_parameter<PetscScalar>("-cf_w", 0.0);
  para.add_parameter<PetscScalar>("-rho_o", 0.);
  para.add_parameter<PetscScalar>("-cf_o", 0.);
  para.add_parameter<PetscScalar>("-ref_pw", 0.);
  para.add_parameter<PetscScalar>("-ref_po", 0.);
  para.add_parameter<PetscScalar>("-Bc", 0.);
  para.add_parameter<PetscInt>("-time_order", 1);
  para.print_parameters();
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormInitialGuess(SNES snes, Vec x, void *ctx)
{
  PetscFunctionBeginUser;
  System *sys = (System *)ctx;
  Mesh &mesh = sys->get_mesh();
  DofMap &dof_map = sys->dof_map();
  int p_var = dof_map.variable_num("p");
  int s_var = dof_map.variable_num("s");
  PetscScalar *x_array;
  PetscCall(VecGetArray(x, &x_array));
  double t = sys->time;
  Para &para = sys->para();
  CaseFluid fluid(para);
  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int s_dof = dof_map.dof_local_indices(e, s_var);
    Point a = e.compute_cell_centroid();
    PetscReal s_exact = fluid.Rsw(t, a.x(), a.y());
    PetscReal p_exact = fluid.Rpw(t, a.x(), a.y());
    x_array[p_dof] = p_exact;
    x_array[s_dof] = s_exact;
  }
  PetscCall(VecRestoreArray(x, &x_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}
