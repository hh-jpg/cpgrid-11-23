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
PetscErrorCode FormFunction_p(SNES, Vec, Vec, void *);
PetscErrorCode FormFunction_s(SNES, Vec, Vec, void *);
PetscErrorCode FormInitialGuess(SNES, Vec, void *);
PetscErrorCode FormInitialGuess_p(SNES, Vec, void *);
PetscErrorCode FormInitialGuess_s(SNES, Vec, void *);
PetscErrorCode FormJacobian_p(SNES, Vec, Mat, Mat, void *);
PetscErrorCode FormJacobian_s(SNES, Vec, Mat, Mat, void *);

#define pow2(a) a*a
#define pow8(a) pow2(a)*pow2(a)*pow2(a)*pow2(a)

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  // Mesh settings
  Mesh mesh;
  std::string case_name = "bl";
  std::string filepath = "./grid/" + case_name;

  mesh.read(filepath);
  mesh.init();
  mesh.printInfo();
  mesh.partition();
  mesh.reset_elem_id();

  // System settings
  System sysp(mesh),syss(mesh);
  sysp.add_variable("p");
  sysp.init();
  syss.add_variable("s");
  syss.init();

  sysp.set_name("pressure");
  syss.set_name("saturation");

  sysp.para().add_parameter<PetscScalar>("-dt", 86400.0);
  sysp.para().add_parameter<PetscInt>("-n_steps", 10);
  sysp.para().add_parameter<PetscScalar>("-phi", 0.2);
  sysp.para().add_parameter<PetscScalar>("-beta", 2.0);
  sysp.para().add_parameter<PetscScalar>("-mu_w", 0.002);
  sysp.para().add_parameter<PetscScalar>("-mu_o", 0.003);
  sysp.para().add_parameter<PetscScalar>("-srw", 0.0);
  sysp.para().add_parameter<PetscScalar>("-srn", 0.2);
  sysp.para().add_parameter<PetscScalar>("-p_out", 1.0e5);
  sysp.para().add_parameter<PetscScalar>("-f_in", 3.47222e-7);
  sysp.para().add_parameter<PetscScalar>("-t_end", 8640000.0);
  sysp.para().print_parameters();

  // DM (Distributed Mesh) settings
  DM dm_p,dm_s;
  PetscCall(DMRegister(DMCP, DMCreate_CP));
  PetscCall(DMCreate(PETSC_COMM_WORLD, &dm_p));
  PetscCall(DMCreate(PETSC_COMM_WORLD, &dm_s));
  PetscCall(DMSetType(dm_p, DMCP));
  PetscCall(DMSetType(dm_s, DMCP));

  // 将 System 与 DM 关联
  PetscErrorCode (*fp)(DM, System &) = nullptr;
  PetscCall(PetscObjectQueryFunction((PetscObject)dm_p, "DMCPSetSystem_C", &fp));
  PetscCall((*fp)(dm_p, sysp));
  PetscCall(DMSetFromOptions(dm_p));
  PetscCall(DMSetUp(dm_p));
  PetscErrorCode (*fs)(DM, System &) = nullptr;
  PetscCall(PetscObjectQueryFunction((PetscObject)dm_s, "DMCPSetSystem_C", &fs));
  PetscCall((*fs)(dm_s, syss));
  PetscCall(DMSetFromOptions(dm_s));
  PetscCall(DMSetUp(dm_s));

  // SNES (Nonlinear solver) settings
  SNES snes_p,snes_s;
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes_p));
  PetscCall(SNESCreate(PETSC_COMM_WORLD, &snes_s));
  PetscCall(SNESSetDM(snes_p, dm_p));
  PetscCall(SNESSetDM(snes_s, dm_s));

  // 创建全局向量
  Vec xp, xs, resp, ress;
  PetscCall(DMCreateGlobalVector(dm_p, &xp));
  PetscCall(DMCreateGlobalVector(dm_s, &xs));
  PetscCall(VecDuplicate(xp, &resp));
  PetscCall(VecDuplicate(xs, &ress));

  // 设置初始猜测和函数
  PetscCall(FormInitialGuess_p(snes_p, xp, &sysp));
  PetscCall(FormInitialGuess_s(snes_s, xs, &syss));

  PetscCall(SNESSetFunction(snes_p, resp, FormFunction_p, &syss));
  PetscCall(SNESSetFunction(snes_s, ress, FormFunction_s, &sysp));
  // 计算F(x)
  // PetscCall(SNESComputeFunction(snes_p, xp, resp));
  // PetscCall(SNESComputeFunction(snes_s, xs, ress));

  // 设置Jacobi函数
  Mat jp;
  sysp.create_mat(jp);
  PetscCall(SNESSetJacobian(snes_p, jp, jp, FormJacobian_p, &syss));
  Mat js;
  sysp.create_mat(js);
  PetscCall(SNESSetJacobian(snes_s, js, js, FormJacobian_s, &sysp));

  // 从命令行选项设置 SNES
  SNESSetOptionsPrefix(snes_p, "p_");
  SNESSetOptionsPrefix(snes_s, "s_");

  // 获取 KSP 对象
  KSP ksp_p,ksp_s;
  SNESGetKSP(snes_p, &ksp_p);
  KSPSetOptionsPrefix(ksp_p, "p_");
  SNESGetKSP(snes_s, &ksp_s);
  KSPSetOptionsPrefix(ksp_s, "s_");

  // 获取 PC 对象并设置预处理器类型
  // PC pc_p,pc_s;
  // KSPGetPC(ksp_p, &pc_p);
  // PCSetType(pc_p, PCLU);  // 设置为 LU 预处理器类型
  // KSPGetPC(ksp_s, &pc_s);
  // PCSetType(pc_s, PCLU);  // 设置为 LU 预处理器类型

  PetscCall(KSPSetFromOptions(ksp_p));
  PetscCall(KSPSetFromOptions(ksp_s));
  PetscCall(SNESSetFromOptions(snes_p));
  PetscCall(SNESSetFromOptions(snes_s));
  PetscCall(SNESSetUp(snes_p));
  PetscCall(SNESSetUp(snes_s));
  // 时间步进

  PetscScalar dt = sysp.para().get_parameter<PetscScalar>("-dt");
  PetscInt n_timesteps = sysp.para().get_parameter<PetscInt>("-n_steps");
  PetscScalar t_end = dt * n_timesteps;
  PetscInt t_step = 1;
  char filename[256];
  
  while (sysp.time < t_end)
  {
    PetscPrintf(PETSC_COMM_WORLD, "\n\n*** Solving time step %d, time = %g, dt = %g ***\n", t_step, sysp.time, dt);
    sysp.time += dt;

    // snprintf(filename, sizeof(filename), "./Test_output/xsini_%d.m", t_step);
    // PetscCall(SaveVecToMatlab(xs, filename, "s"));

    // 保存旧的解并进行求解
    PetscPrintf(PETSC_COMM_WORLD, "----------------- snes_p -----------------\n");
    syss.set_old_solution(xs);

  // ierr = SNESComputeFunction(snes_p,xp,resp);CHKERRQ(ierr);
  // SaveVecToMatlab(resp, "./Test_output/resp.m", "res");

    PetscCall(SNESSolve(snes_p, NULL, xp));
    PetscCall(SNESGetSolution(snes_p, &xp));// 获取并保存解, x---p^{n+1} 、xold---S_w^{n}
    // if (t_step % 5 == 0 || t_step <=10)
    // {
    // snprintf(filename, sizeof(filename), "./output/xp_%d.m", t_step);
    // PetscCall(SaveVecToMatlab(xp, filename, "pressure"));
    // }
    if (t_step % 10 == 0 )
    {
    snprintf(filename, sizeof(filename), "./Test_output/xp_%d.m", t_step);
    PetscCall(SaveVecToMatlab(xp, filename, "p"));
    }
// return 0;

    // 保存旧的解并进行求解
    PetscPrintf(PETSC_COMM_WORLD, "----------------- snes_s -----------------\n");
    sysp.set_old_solution(xp);
    PetscCall(SNESSolve(snes_s, NULL, xs));
    PetscCall(SNESGetSolution(snes_s, &xs));// 获取并保存解, xold---p^{n+1} 、x---S_w^{n+1}
    if (t_step % 10 == 0 )
    {
    snprintf(filename, sizeof(filename), "./Test_output/xs_%d.m", t_step);
    PetscCall(SaveVecToMatlab(xs, filename, "s"));
    }
// return 0;
    t_step++;
  }
  
  // 清理资源  
  PetscCall(VecDestroy(&xp));
  PetscCall(VecDestroy(&xs));
  PetscCall(VecDestroy(&resp));
  PetscCall(VecDestroy(&ress));
  PetscCall(SNESDestroy(&snes_p));
  PetscCall(SNESDestroy(&snes_s));
  PetscCall(DMDestroy(&dm_p));
  PetscCall(DMDestroy(&dm_s));
  PetscCall(PetscFinalize());
  return 0;
}
PetscErrorCode FormFunction_p(SNES snes, Vec p, Vec f, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecZeroEntries(f));
  System *syss = (System*) ctx, *sysp;
  Vec p_local, s_old_local;
  PetscScalar *f_array, *p_array, *s_old_array;
  DM dm_p;

  // 获取DM和System
  // 根据 snes 和 dm_p 获取 sysp 的信息
  PetscCall(SNESGetDM(snes, &dm_p));
  PetscCall(DMCPGetSystem_CP(dm_p, sysp));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = sysp->dof_map();
  Mesh &mesh = sysp->get_mesh();
  Para &para = sysp->para(); 
  Fluid fluid(para);

  // 获取参数
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt"); 
  PetscScalar phi = para.get_parameter<PetscScalar>("-phi");
  PetscScalar mu_w = para.get_parameter<PetscScalar>("-mu_w");
  PetscScalar mu_o = para.get_parameter<PetscScalar>("-mu_o");
  PetscScalar f_in = para.get_parameter<PetscScalar>("-f_in");
  PetscScalar p_out = para.get_parameter<PetscScalar>("-p_out");
  // std::cout<<"Right pressure = "<<p_out<<std::endl;

  // 获取变量编号
  int p_var = dof_map.variable_num("p");

  // 将全局向量转换为局部向量
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  sysp->create_local_vec(&p_local);
  VecDuplicate(p_local, &s_old_local);
  PetscCall(globalvec_to_local(p, ghosted_index, &p_local));
  PetscCall(globalvec_to_local(syss->get_old_solution(), ghosted_index, &s_old_local));

  // PetscCall(SaveVecToMatlab(syss->get_old_solution(), "./Test_output/xsgetini.m", "s"));
  // PetscCall(SaveVecToMatlab(s_old_local, "./Test_output/xsgetinilocal.m", "s"));

  // 获取向量数组
  PetscCall(VecGetArray(p_local, &p_array));
  PetscCall(VecGetArray(s_old_local, &s_old_array));
  PetscCall(VecGetArray(f, &f_array));

  // 遍历网格中的单元
  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int p_dof = dof_map.dof_local_indices(e, p_var);
    PetscScalar e_volume = e.compute_cell_volume();
    PetscScalar p = p_array[p_dof];
    PetscScalar s = s_old_array[p_dof];
    // PetscScalar s = 0.0;

    // PetscPrintf(PETSC_COMM_WORLD,"sw_%d=%.5e\n",p_dof,(s));
    // PetscPrintf(PETSC_COMM_WORLD,"p_%d=%.5e\n",p_dof,(p));
    // 计算残差
    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 0)
        {
          f_array[p_dof] -= f_in * face.compute_face_area() * pow8(10);
          // PetscPrintf(PETSC_COMM_WORLD,"p_dof=%d\n",p_dof);
          // PetscPrintf(PETSC_COMM_WORLD,"Left q =%.5e\n",f_in * face.compute_face_area());   // 1.157407407407407e-07
        }
        else if (face.pos() == 1)
        {
          PetscScalar mob_oc = fluid.func_lambda_o(s);
          PetscScalar mob_wc = fluid.func_lambda_w(s);
          f_array[p_dof] += face.compute_efftrans(mob_oc + mob_wc) * (p - p_out) * pow8(10);
          // PetscPrintf(PETSC_COMM_WORLD,"Right q =%.5e\n",face.compute_efftrans(mob_oc + mob_wc) * (p_out));  // 6.579488444773420e-08
          // PetscPrintf(PETSC_COMM_WORLD,"A_%d(1)=%.5e\n",p_dof,face.compute_efftrans(mob_oc + mob_wc));
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);

      // PetscPrintf(PETSC_COMM_WORLD,"p_n_dof=%d\n",p_n_dof);
        PetscScalar p_n = p_array[p_n_dof];
        PetscScalar s_n = s_old_array[p_n_dof];
        PetscScalar mob_oc = fluid.func_lambda_o(s);
        PetscScalar mob_on = fluid.func_lambda_o(s_n);
        PetscScalar mob_wc = fluid.func_lambda_w(s);
        PetscScalar mob_wn = fluid.func_lambda_w(s_n);
        PetscScalar mob_up;
        if (p > p_n)
        {
          mob_up = mob_oc + mob_wc;
        }
        else
        {
          mob_up = mob_on + mob_wn;
        }
        f_array[p_dof] += face.compute_efftrans(mob_up) * (p - p_n) * pow8(10);
       }
      //  std::cout<<"f_array["<<p_dof<<"]="<<f_array[p_dof]<<std::endl;
    }

    // PetscPrintf(PETSC_COMM_WORLD,"pppp_%d=%.5e\n",p_dof,(p));
  }


  // 恢复向量数组
  PetscCall(VecRestoreArray(p_local, &p_array));
  PetscCall(VecRestoreArray(s_old_local, &s_old_array));
  PetscCall(VecRestoreArray(f, &f_array));
  PetscCall(VecDestroy(&p_local));
  PetscCall(VecDestroy(&s_old_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormFunction_s(SNES snes, Vec s, Vec f, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecZeroEntries(f));
  Vec s_local, p_old_local, s_old_local;
  PetscScalar *f_array, *s_array, *p_old_array, *s_old_array;
  System *sysp = (System*) ctx,*syss;
  DM dm_s;
  // 获取DM和System
  PetscCall(SNESGetDM(snes, &dm_s));
  PetscCall(DMCPGetSystem_CP(dm_s, syss));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = syss->dof_map();
  Mesh &mesh = syss->get_mesh();
  Para &para = sysp->para();
  Fluid fluid(para);

  // 获取参数
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscScalar phi = para.get_parameter<PetscScalar>("-phi");
  PetscScalar mu_w = para.get_parameter<PetscScalar>("-mu_w");
  PetscScalar mu_o = para.get_parameter<PetscScalar>("-mu_o");
  PetscScalar f_in = para.get_parameter<PetscScalar>("-f_in");
  PetscScalar p_out = para.get_parameter<PetscScalar>("-p_out");

  // 获取变量编号
  int s_var = dof_map.variable_num("s");

  // 将全局向量转换为局部向量
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  syss->create_local_vec(&s_local);
  VecDuplicate(s_local, &p_old_local);
  VecDuplicate(s_local, &s_old_local);
  PetscCall(globalvec_to_local(s, ghosted_index, &s_local));
  PetscCall(globalvec_to_local(sysp->get_old_solution(), ghosted_index, &p_old_local));
  PetscCall(globalvec_to_local(syss->get_old_solution(), ghosted_index, &s_old_local));

  // 获取向量数组
  PetscCall(VecGetArray(s_local, &s_array));
  PetscCall(VecGetArray(p_old_local, &p_old_array));
  PetscCall(VecGetArray(s_old_local, &s_old_array));
  PetscCall(VecGetArray(f, &f_array));

  // 遍历网格中的单元
  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int s_dof = dof_map.dof_local_indices(e, s_var);
    PetscScalar e_volume = e.compute_cell_volume();
    PetscScalar p = p_old_array[s_dof];
    PetscScalar s = s_array[s_dof];
    PetscScalar s_old = s_old_array[s_dof];
    /*
    if (s<0.0)
    {
      s=0.0;
    }else if (s>1.0)
    {
      s=1.0-1.e-24;
    }
    */
    // 计算残差
    f_array[s_dof] += phi * (s - s_old) * e_volume * pow8(10);

    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 0)
        {
          f_array[s_dof] -= dt * f_in * face.compute_face_area()* pow8(10);
          // PetscPrintf(PETSC_COMM_WORLD,"q=%.6e\n",-f_in);
        }
        else if (face.pos() == 1)
        {
          PetscScalar mob_wc = fluid.func_lambda_w(s);
          PetscScalar mob_oc = fluid.func_lambda_o(s);
          // PetscPrintf(PETSC_COMM_WORLD,"mob_wc = %.4g, mob_oc = %.4g\n",mob_wc,mob_oc);
          PetscScalar fw = mob_wc / (mob_wc + mob_oc);
          PetscScalar q = face.compute_efftrans(mob_wc) * (p_out - p);
          // PetscPrintf(PETSC_COMM_WORLD,"q=%.5g\n",q);
          f_array[s_dof] -= dt * ( (q>0)?q:0 + ((q<0)?q:0) * fw)* pow8(10);
          // f_array[s_dof] += dt * face.compute_efftrans(mob_wc) * (p - p_out);
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int s_n_dof = dof_map.dof_local_indices(*neighbor, s_var);

        PetscScalar s_n = s_array[s_n_dof];
        PetscScalar p_n = p_old_array[s_n_dof];
        PetscScalar mob_wc = fluid.func_lambda_w(s);
        PetscScalar mob_wn = fluid.func_lambda_w(s_n);
        PetscScalar mob_up = (p > p_n) ? mob_wc : mob_wn;
        f_array[s_dof] += dt * face.compute_efftrans(mob_up) * (p - p_n)* pow8(10);
      }
    }
  }

  // 恢复向量数组
  PetscCall(VecRestoreArray(s_local, &s_array));
  PetscCall(VecRestoreArray(s_old_local, &s_old_array));
  PetscCall(VecRestoreArray(p_old_local, &p_old_array));
  PetscCall(VecRestoreArray(f, &f_array));
  PetscCall(VecDestroy(&s_local));
  PetscCall(VecDestroy(&s_old_local));
  PetscCall(VecDestroy(&p_old_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

// PetscErrorCode FormJacobian_p(SNES snes, Vec x, Mat jac, Mat B, void *ctx){
//   PetscFunctionReturn(PETSC_SUCCESS);}
// PetscErrorCode FormJacobian_s(SNES snes, Vec x, Mat jac, Mat B, void *ctx){
//   PetscFunctionReturn(PETSC_SUCCESS);}

PetscErrorCode FormJacobian_p(SNES snes, Vec x, Mat jac, Mat B, void *ctx)
{
  PetscFunctionBeginUser;
  DM dm_p;
  PetscCall(MatZeroEntries(jac));
  System *syss = (System*) ctx,*sysp;
  Vec p_local, s_old_local;
  PetscScalar *p_array,*s_old_array;
  PetscErrorCode ierr;

  // 获取DM和System
  // 根据 snes 和 dm_p 获取 sysp 的信息
  PetscCall(SNESGetDM(snes, &dm_p));
  PetscCall(DMCPGetSystem_CP(dm_p, sysp));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = sysp->dof_map();
  Mesh &mesh = sysp->get_mesh();
  Para &para = sysp->para();
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

  // 将全局向量转换为局部向量
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  sysp->create_local_vec(&p_local);
  VecDuplicate(p_local, &s_old_local);
  PetscCall(globalvec_to_local(x, ghosted_index, &p_local));
  PetscCall(globalvec_to_local(syss->get_old_solution(), ghosted_index, &s_old_local));
  // 获取向量数组
  PetscCall(VecGetArray(p_local, &p_array));
  PetscCall(VecGetArray(s_old_local, &s_old_array));

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    PetscScalar Kpp = 0.0, Kps = 0.0, Ksp = 0.0, Kss = 0.0;
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int p_g_dof = dof_map.dof_indices(e, p_var);
    PetscScalar e_volume = e.compute_cell_volume();

    PetscScalar p = p_array[p_dof];
    PetscScalar s = s_old_array[p_dof];
    
    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 1)
        {
          PetscScalar mob_oc = fluid.func_lambda_o(s);
          PetscScalar mob_wc = fluid.func_lambda_w(s);
          PetscScalar mob_c = mob_oc + mob_wc;
          Kpp += face.compute_efftrans(mob_c);
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);
        int p_n_g_dof = dof_map.dof_indices(*neighbor, p_var);
        PetscScalar Kpnp = 0.0, Kpns = 0.0, knps = 0.0, Ksnp = 0.0, Ksns = 0.0;
        PetscScalar p_n = p_array[p_n_dof];
        PetscScalar s_n = s_old_array[p_n_dof];
        PetscScalar mob_oc = fluid.func_lambda_o(s);
        PetscScalar mob_on = fluid.func_lambda_o(s_n);
        PetscScalar mob_wc = fluid.func_lambda_w(s);
        PetscScalar mob_wn = fluid.func_lambda_w(s_n);
        PetscScalar mob_c = mob_oc + mob_wc, mob_n = mob_on + mob_wn;
        if (p > p_n)
        {
          Kpp += face.compute_efftrans(mob_c);
          Kpnp += -face.compute_efftrans(mob_c);
        }
        else
        {
          Kpp += face.compute_efftrans(mob_n);
          Kpnp += -face.compute_efftrans(mob_n);
        }
        MatSetValue(jac, p_g_dof, p_n_g_dof, Kpnp, ADD_VALUES);
      }
    }
    MatSetValue(jac, p_g_dof, p_g_dof, Kpp, ADD_VALUES);
  }
  PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(VecRestoreArray(p_local, &p_array));
  PetscCall(VecRestoreArray(s_old_local, &s_old_array));
  PetscCall(VecDestroy(&p_local));
  PetscCall(VecDestroy(&s_old_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode FormJacobian_s(SNES snes, Vec x, Mat jac, Mat B, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(MatZeroEntries(jac));
  // Vec s_local, p_old_local, s_old_local;
  // PetscScalar *s_array, *p_old_array, *s_old_array;
  Vec s_local, p_old_local;
  PetscScalar *s_array, *p_old_array;
  System *sysp = (System*) ctx,*syss;
  DM dm_s;
  PetscErrorCode ierr;

  // 获取DM和System
  PetscCall(SNESGetDM(snes, &dm_s));
  PetscCall(DMCPGetSystem_CP(dm_s, syss));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = syss->dof_map();
  Mesh &mesh = syss->get_mesh();
  Para &para = sysp->para();
  Fluid fluid(para);
  // 获取参数
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscScalar phi = para.get_parameter<PetscScalar>("-phi");
  PetscScalar mu_w = para.get_parameter<PetscScalar>("-mu_w");
  PetscScalar mu_o = para.get_parameter<PetscScalar>("-mu_o");
  PetscScalar f_in = para.get_parameter<PetscScalar>("-f_in");
  PetscScalar p_out = para.get_parameter<PetscScalar>("-p_out");
  // 获取变量编号
  int s_var = dof_map.variable_num("s");

  // 将全局向量转换为局部向量
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  syss->create_local_vec(&s_local);
  // VecDuplicate(s_local, &s_old_local);
  VecDuplicate(s_local, &p_old_local);
  PetscCall(globalvec_to_local(x, ghosted_index, &s_local));
  PetscCall(globalvec_to_local(sysp->get_old_solution(), ghosted_index, &p_old_local));
  // PetscCall(globalvec_to_local(syss->get_old_solution(), ghosted_index, &s_old_local));

  // 获取向量数组
  PetscCall(VecGetArray(s_local, &s_array));
  // PetscCall(VecGetArray(s_old_local, &s_old_array));
  PetscCall(VecGetArray(p_old_local, &p_old_array));

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    PetscScalar Kpp = 0.0, Kps = 0.0, Ksp = 0.0, Kss = 0.0;
    int s_dof = dof_map.dof_local_indices(e, s_var);
    int s_g_dof = dof_map.dof_indices(e, s_var);
    PetscScalar e_volume = e.compute_cell_volume();

    PetscScalar p = p_old_array[s_dof];
    PetscScalar s = s_array[s_dof];
    Kss += phi * e_volume;
    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 1)
        {
          PetscScalar mob_wc = fluid.func_lambda_w(s);
          Kss += dt * face.compute_trans() * fluid.func_dlambda_w_ds(s) * (p - p_out);
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int s_n_dof = dof_map.dof_local_indices(*neighbor, s_var);
        int s_n_g_dof = dof_map.dof_indices(*neighbor, s_var);
        PetscScalar Kpnp = 0.0, Kpns = 0.0, Knpp = 0.0, knps = 0.0, Ksnp = 0.0, Ksns = 0.0;
        PetscScalar p_n = p_old_array[s_n_dof];
        PetscScalar s_n = s_array[s_n_dof];
        if (p > p_n)
        {
          Kss += dt * face.compute_trans() * fluid.func_dlambda_w_ds(s) * (p - p_n);
        }
        else
        {
          Ksns += dt * face.compute_trans() * fluid.func_dlambda_w_ds(s_n) * (p - p_n);
        }
        MatSetValue(jac, s_g_dof, s_n_g_dof, Ksns, ADD_VALUES);
      }
    }
    MatSetValue(jac, s_g_dof, s_g_dof, Kss, ADD_VALUES);
  }
  PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(VecRestoreArray(s_local, &s_array));
  PetscCall(VecRestoreArray(p_old_local, &p_old_array));
  PetscCall(VecDestroy(&s_local));
  PetscCall(VecDestroy(&p_old_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormInitialGuess_p(SNES snes, Vec x, void *ctx)
{
  PetscFunctionBeginUser;
  System *sys = (System *)ctx;
  Mesh &mesh = sys->get_mesh();
  DofMap &dof_map = sys->dof_map();
  int p_var = dof_map.variable_num("p");
  PetscScalar *x_array;
  PetscCall(VecGetArray(x, &x_array));

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int p_dof = dof_map.dof_local_indices(e, p_var);
    x_array[p_dof] = 1.e5;
  }
  PetscCall(VecRestoreArray(x, &x_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormInitialGuess_s(SNES snes, Vec x, void *ctx)
{
  PetscFunctionBeginUser;
  System *sys = (System *)ctx;
  Mesh &mesh = sys->get_mesh();
  DofMap &dof_map = sys->dof_map();
  int s_var = dof_map.variable_num("s");
  PetscScalar *x_array;
  PetscCall(VecGetArray(x, &x_array));

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int s_dof = dof_map.dof_local_indices(e, s_var);
    x_array[s_dof] = 0.0;
  }
  PetscCall(VecRestoreArray(x, &x_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}
