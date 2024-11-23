/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-01 15:43:34
 * @LastEditTime: 2024-08-25 21:16:51
 * @FilePath: /cpgrid/example/twophase/buckley_leverett/main.C
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

PetscErrorCode FormFunction(SNES, Vec, Vec, void *);
PetscErrorCode FormInitialGuess(SNES, Vec, void *);
PetscErrorCode FormJacobian(SNES, Vec, Mat, Mat, void *);
double wellindex(double K, double d1, double d2, double d3);

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  // Mesh settings
  Mesh mesh;
  std::string case_name = "CFL";
  std::string filepath = "./grid/" + case_name;

  mesh.read(filepath);
  mesh.init();
  mesh.printInfo();
  mesh.partition();
  mesh.reset_elem_id();

  // System settings
  System sys(mesh);
  sys.add_variable("p");
  sys.add_variable("s");
  sys.init();

  sys.set_name(case_name);
  Para &para = sys.para();

  // 添加模拟参数 值为默认值，实际的设置要在option file 里设置
  para.add_parameter<PetscScalar>("-dt", 86400.0);
  para.add_parameter<PetscInt>("-n_steps", 10);
  para.add_parameter<PetscScalar>("-phi", 0.2);
  para.add_parameter<PetscScalar>("-well_dx", 1.0);
  para.add_parameter<PetscScalar>("-well_dy", 1.0);
  para.add_parameter<PetscScalar>("-well_dz", 1.0);
  para.add_parameter<PetscScalar>("-beta", 2.0);
  para.add_parameter<PetscScalar>("-mu_w", 0.002);
  para.add_parameter<PetscScalar>("-mu_o", 0.003);
  para.add_parameter<PetscScalar>("-srw", 0.0);
  para.add_parameter<PetscScalar>("-srn", 0.2);
  para.add_parameter<PetscScalar>("-p_out", 1.0e5);
  para.add_parameter<PetscScalar>("-f_in", 3.47222e-7);
  para.add_parameter<PetscScalar>("-t_end", 8640000.0);
  para.print_parameters();

  // DM (Distributed Mesh) settings
  DM dm_cp;
  PetscCall(DMRegister(DMCP, DMCreate_CP));
  PetscCall(DMCreate(PETSC_COMM_WORLD, &dm_cp));
  PetscCall(DMSetType(dm_cp, DMCP));

  // 将 System 与 DM 关联
  PetscErrorCode (*f)(DM, System &) = nullptr;
  PetscCall(PetscObjectQueryFunction((PetscObject)dm_cp, "DMCPSetSystem_C", &f));
  PetscCall((*f)(dm_cp, sys));
  PetscCall(DMSetFromOptions(dm_cp));
  PetscCall(DMSetUp(dm_cp));

  // SNES (Nonlinear solver) settings
  SNES snes;
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes));
  PetscCall(SNESSetDM(snes, dm_cp));

  // 创建全局向量
  Vec x, res;
  PetscCall(DMCreateGlobalVector(dm_cp, &x));
  PetscCall(VecDuplicate(x, &res));

  // 设置初始猜测和函数
  // PetscCall(SNESSetComputeInitialGuess(snes, FormInitialGuess, &sys));
  PetscCall(FormInitialGuess(snes, x, &sys));

  PetscCall(SNESSetFunction(snes, res, FormFunction, &sys));

  Mat j;
  sys.create_mat(j);
  PetscCall(SNESSetJacobian(snes, j, j, FormJacobian, NULL));

  // 从命令行选项设置 SNES
  PetscCall(SNESSetFromOptions(snes));
  PetscCall(SNESSetUp(snes));
  // 时间步进

  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscInt n_timesteps = para.get_parameter<PetscInt>("-n_steps");
  PetscScalar t_end = dt * n_timesteps;
  PetscInt t_step = 1;
  char filename[256];

  while (sys.time < t_end)
  {
    PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Solving time step %d, time = %g, dt = %g ***\n", t_step, sys.time, dt);
    sys.time += dt;

    // 保存旧的解并进行求解
    sys.set_old_solution(x);
    PetscCall(SNESSolve(snes, NULL, x));

    // 获取并保存解
    PetscCall(SNESGetSolution(snes, &x));

    // if (t_step % 10 == 0)
    // {
    //   snprintf(filename, sizeof(filename), "./output/cfl_bl_solution_%d_%.2f.m", t_step,dt);
    //   PetscCall(SaveVecToMatlab(x, filename, "p_s"));
    // }
    if ((t_step * dt) == t_end)
    {
      snprintf(filename, sizeof(filename), "./output/cfl_bl_solution_%d_%.0e.m", t_step,dt);
      PetscCall(SaveVecToMatlab(x, filename, "p_s"));
    }
    t_step++;
  }
  // 清理资源
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&res));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm_cp));
  PetscCall(PetscFinalize());
  return 0;
}

PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecZeroEntries(f));
  DM dm_cp;
  System *sys;
  Vec x_local, x_old_local;
  PetscScalar *f_array, *x_array, *x_old_array;
  PetscErrorCode ierr;
  
  // 获取DM和System
  PetscCall(SNESGetDM(snes, &dm_cp));
  PetscCall(DMCPGetSystem_CP(dm_cp, sys));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = sys->dof_map();
  Mesh &mesh = sys->get_mesh();
  Para &para = sys->para();
  Fluid fluid(para);

  // 获取参数
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscScalar phi = para.get_parameter<PetscScalar>("-phi");
  PetscScalar mu_w = para.get_parameter<PetscScalar>("-mu_w");
  PetscScalar mu_o = para.get_parameter<PetscScalar>("-mu_o");
  PetscScalar f_in = para.get_parameter<PetscScalar>("-f_in");
  PetscScalar p_out = para.get_parameter<PetscScalar>("-p_out");

  // 获取变量编号
  int p_var = dof_map.variable_num("p");
  int s_var = dof_map.variable_num("s");

  // 将全局向量转换为局部向量
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  sys->create_local_vec(&x_local);
  VecDuplicate(x_local, &x_old_local);
  PetscCall(globalvec_to_local(x, ghosted_index, &x_local));
  PetscCall(globalvec_to_local(sys->get_old_solution(), ghosted_index, &x_old_local));

  // 获取向量数组
  PetscCall(VecGetArray(x_local, &x_array));
  PetscCall(VecGetArray(x_old_local, &x_old_array));
  PetscCall(VecGetArray(f, &f_array));

  // 遍历网格中的单元
  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int p_global = dof_map.dof_indices(e, p_var);
    int s_dof = dof_map.dof_local_indices(e, s_var);
    PetscScalar e_volume = e.compute_cell_volume();
    PetscScalar p = x_array[p_dof];
    PetscScalar s = x_array[s_dof];
   /*
   if (s<0.0)
    {
      s=0.0;
    }else if (s>1.0)
    {
      s=1.0-1.e-24;
    }
    */

    PetscScalar s_old = x_old_array[s_dof];

    // 计算残差
    f_array[p_dof] += phi * (s_old - s) * e_volume;
    f_array[s_dof] += phi * (s - s_old) * e_volume;

    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 0)
        {
          f_array[s_dof] -= dt * f_in * face.compute_face_area();
        }
        else if (face.pos() == 1)
        {
          PetscScalar mob_oc = fluid.func_lambda_o(s);
          PetscScalar mob_wc = fluid.func_lambda_w(s);
          // f_array[p_dof] += dt * face.compute_efftrans(mob_oc) * (p - p_out);
          // f_array[s_dof] += dt * face.compute_efftrans(mob_wc) * (p - p_out);
          f_array[p_dof] -= - dt * f_in * face.compute_face_area();
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);
        int s_n_dof = dof_map.dof_local_indices(*neighbor, s_var);

        PetscScalar p_n = x_array[p_n_dof];
        PetscScalar s_n = x_array[s_n_dof];
        PetscScalar mob_oc = fluid.func_lambda_o(s);
        PetscScalar mob_on = fluid.func_lambda_o(s_n);
        PetscScalar mob_wc = fluid.func_lambda_w(s);
        PetscScalar mob_wn = fluid.func_lambda_w(s_n);
        PetscScalar mob_oup, mob_wup;
        if (p > p_n)
        {
          mob_oup = mob_oc;
          mob_wup = mob_wc;
        }
        else
        {
          mob_oup = mob_on;
          mob_wup = mob_wn;
        }
        f_array[p_dof] += dt * face.compute_efftrans(mob_oup) * (p - p_n);
        f_array[s_dof] += dt * face.compute_efftrans(mob_wup) * (p - p_n);
      }
    }
  }

  // 恢复向量数组
  PetscCall(VecRestoreArray(x_local, &x_array));
  PetscCall(VecRestoreArray(x_old_local, &x_old_array));
  PetscCall(VecRestoreArray(f, &f_array));
  PetscCall(VecDestroy(&x_local));
  PetscCall(VecDestroy(&x_old_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat jac, Mat B, void *ctx)
{
  PetscFunctionBeginUser;
  DM dm_cp;
  PetscCall(MatZeroEntries(jac));
  System *sys;
  Vec x_local, x_old_local;
  PetscScalar *x_array;
  PetscErrorCode ierr;

  // 获取DM和System
  PetscCall(SNESGetDM(snes, &dm_cp));
  PetscCall(DMCPGetSystem_CP(dm_cp, sys));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = sys->dof_map();
  Mesh &mesh = sys->get_mesh();
  Para &para = sys->para();
  Fluid fluid(para);
  // 获取参数
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscScalar phi = para.get_parameter<PetscScalar>("-phi");
  PetscScalar mu_w = para.get_parameter<PetscScalar>("-mu_w");
  PetscScalar mu_o = para.get_parameter<PetscScalar>("-mu_o");
  PetscScalar f_in = para.get_parameter<PetscScalar>("-f_in");
  PetscScalar p_out = para.get_parameter<PetscScalar>("-p_out");
  // 获取变量编号
  int p_var = dof_map.variable_num("p");
  int s_var = dof_map.variable_num("s");

  // 将全局向量转换为局部向量
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  sys->create_local_vec(&x_local);
  VecDuplicate(x_local, &x_old_local);
  PetscCall(globalvec_to_local(x, ghosted_index, &x_local));
  // 获取向量数组
  PetscCall(VecGetArray(x_local, &x_array));

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    PetscScalar Kpp = 0.0, Kps = 0.0, Ksp = 0.0, Kss = 0.0;
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int s_dof = dof_map.dof_local_indices(e, s_var);
    int p_g_dof = dof_map.dof_indices(e, p_var);
    int s_g_dof = dof_map.dof_indices(e, s_var);
    PetscScalar e_volume = e.compute_cell_volume();

    PetscScalar p = x_array[p_dof];
    PetscScalar s = x_array[s_dof];
    Kps += -phi * e_volume;
    Kss += phi * e_volume;
    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 1)
        {
          PetscScalar mob_oc = fluid.func_lambda_o(s);
          PetscScalar mob_wc = fluid.func_lambda_w(s);
          Kpp += dt * face.compute_efftrans(mob_oc);
          Kps += dt * face.compute_trans() * fluid.func_dlambda_o_ds(s) * (p - p_out);
          Ksp += dt * face.compute_efftrans(mob_wc);
          Kss += dt * face.compute_trans() * fluid.func_dlambda_w_ds(s) * (p - p_out);
          // double a=1;
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);
        int s_n_dof = dof_map.dof_local_indices(*neighbor, s_var);
        int p_n_g_dof = dof_map.dof_indices(*neighbor, p_var);
        int s_n_g_dof = dof_map.dof_indices(*neighbor, s_var);
        PetscScalar Kpnp = 0.0, Kpns = 0.0, Knpp = 0.0, knps = 0.0, Ksnp = 0.0, Ksns = 0.0;
        PetscScalar p_n = x_array[p_n_dof];
        PetscScalar s_n = x_array[s_n_dof];
        PetscScalar mob_oc = fluid.func_lambda_o(s);
        PetscScalar mob_on = fluid.func_lambda_o(s_n);
        PetscScalar mob_wc = fluid.func_lambda_w(s);
        PetscScalar mob_wn = fluid.func_lambda_w(s_n);
        PetscScalar mob_oup, mob_wup;
        if (p > p_n)
        {
          Kpp += dt * face.compute_efftrans(mob_oc);
          Kpnp += -dt * face.compute_efftrans(mob_oc);
          Kps += dt * face.compute_trans() * fluid.func_dlambda_o_ds(s) * (p - p_n);
          Ksp += dt * face.compute_efftrans(mob_wc);
          Ksnp += -dt * face.compute_efftrans(mob_wc);
          Kss += dt * face.compute_trans() * fluid.func_dlambda_w_ds(s) * (p - p_n);
        }
        else
        {
          Kpp += dt * face.compute_efftrans(mob_on);
          Kpnp += -dt * face.compute_efftrans(mob_on);
          Kpns += dt * face.compute_trans() * fluid.func_dlambda_o_ds(s_n) * (p - p_n);
          Ksp += dt * face.compute_efftrans(mob_wn);
          Ksnp += -dt * face.compute_efftrans(mob_wn);
          Ksns += dt * face.compute_trans() * fluid.func_dlambda_w_ds(s_n) * (p - p_n);
        }
        MatSetValue(jac, p_g_dof, s_n_g_dof, Kpns, ADD_VALUES);
        MatSetValue(jac, p_g_dof, p_n_g_dof, Kpnp, ADD_VALUES);
        MatSetValue(jac, s_g_dof, p_n_g_dof, Ksnp, ADD_VALUES);
        MatSetValue(jac, s_g_dof, s_n_g_dof, Ksns, ADD_VALUES);
      }
    }
    MatSetValue(jac, p_g_dof, p_g_dof, Kpp, ADD_VALUES);
    MatSetValue(jac, p_g_dof, s_g_dof, Kps, ADD_VALUES);
    MatSetValue(jac, s_g_dof, p_g_dof, Ksp, ADD_VALUES);
    MatSetValue(jac, s_g_dof, s_g_dof, Kss, ADD_VALUES);
  }
  PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(VecRestoreArray(x_local, &x_array));
  PetscCall(VecDestroy(&x_local));
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

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int s_dof = dof_map.dof_local_indices(e, s_var);
    x_array[p_dof] = 1.e5;
    // x_array[p_dof] = 0.0;
    x_array[s_dof] = 0.0;
  }
  PetscCall(VecRestoreArray(x, &x_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}