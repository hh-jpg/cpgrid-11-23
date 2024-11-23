/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-23 14:22:43
 * @LastEditTime: 2024-09-01 16:55:04
 * @FilePath: /cpgrid/example/twophase/convergence_test_1/res_jacobian.C
 * @Description:
 *
 */
#include <petscvec.h>
#include <petscdm.h>
#include <petscsnes.h>
#include <petscsys.h>
#include <petscdmtypes.h>
#include <iomanip>
#include "utils.h"
#include "petscdmcp.h"
#include "vec_func.h"
#include "system.h"
#include "dof_map.h"
#include "ghosted_array_wrapper.h"
#include "fluid.h"
#include "mat_func.h"
#include "case_fluid.h"
#define Ghost 1
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecZeroEntries(f));
  DM dm_cp;
  System *sys;
  Vec x_local, x_old_local;
  PetscScalar *f_array, *x_array, *x_old_array;

  // Get DM and System objects
  PetscCall(SNESGetDM(snes, &dm_cp));
  PetscCall(DMCPGetSystem_CP(dm_cp, sys));
  // Access system components
  DofMap &dof_map = sys->dof_map();
  Mesh &mesh = sys->get_mesh();
  Para &para = sys->para();
  CaseFluid fluid(para);

  // Get time and parameters
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscScalar phi = para.get_parameter<PetscScalar>("-phi");
  PetscInt t_order = para.get_parameter<PetscInt>("-time_order"); 
  double t = sys->time;

  // Get variable numbers for pressure and saturation
  int p_var = dof_map.variable_num("p");
  int s_var = dof_map.variable_num("s");

  // Convert global vectors to local vectors
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  PetscCall(sys->create_local_vec(&x_local));
  PetscCall(VecDuplicate(x_local, &x_old_local));
  PetscCall(globalvec_to_local(x, ghosted_index, &x_local));
  PetscCall(globalvec_to_local(sys->get_old_solution(), ghosted_index, &x_old_local));
  
  // Access vector arrays
  PetscCall(VecGetArray(x_local, &x_array));
  PetscCall(VecGetArray(x_old_local, &x_old_array));
  PetscCall(VecGetArray(f, &f_array));

  // Loop over local elements in the mesh
  for (const Polyhedron &e : mesh.elem_local_range())
  {
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int s_dof = dof_map.dof_local_indices(e, s_var);
    PetscScalar e_volume = e.compute_cell_volume();
    Point a = e.compute_cell_centroid();

    // Exact solutions for pressure, saturation, and fluxes
    PetscReal s_exact = fluid.Rsw(t, a.x(), a.y());
    PetscReal p_exact = fluid.Rpw(t, a.x(), a.y());
    PetscReal qw_exact = fluid.Rqw(t, a.x(), a.y());
    PetscReal qo_exact = fluid.Rqo(t, a.x(), a.y());
    PetscReal pc_exact = fluid.Rpc(t, a.x(), a.y());
    PetscReal rhow_exact = fluid.Rrhow(t, a.x(), a.y());
    PetscReal rhoo_exact = fluid.Rrho_o(t, a.x(), a.y());
    PetscScalar mobo_exact = fluid.func_lambda_o(s_exact) * rhoo_exact;
    PetscScalar mobw_exact = fluid.func_lambda_w(s_exact) * rhow_exact;

    PetscScalar p = x_array[p_dof];
    PetscScalar s = x_array[s_dof];
    PetscScalar s_old = x_old_array[s_dof];
    PetscScalar p_old = x_old_array[p_dof];
    PetscScalar so = 1.0 - s;
    PetscScalar so_old = 1.0 - s_old;

    // Mobility and pressure calculations
    PetscScalar mobo = fluid.func_lambda_o(s);
    PetscScalar mobw = fluid.func_lambda_w(s);

    PetscScalar pc = fluid.func_pc(s, s_old);
    PetscScalar po = p + pc;

    PetscScalar rhow = fluid.func_rho_w(p);
    PetscScalar rhow_old = fluid.func_rho_w(p_old);
    PetscScalar rhoo = fluid.func_rho_o(p);
    PetscScalar rhoo_old = fluid.func_rho_o(p_old);

    // Compute residuals
    // time term
   if (t_order == 1)
   {
    f_array[p_dof] += phi / dt * (rhoo * so - rhoo_old * so_old) * e_volume;
    f_array[s_dof] += phi / dt * (rhow * s - rhow_old * s_old) * e_volume;
   }   
   
    // source term
    f_array[p_dof] -= qo_exact * e_volume;
    f_array[s_dof] -= qw_exact * e_volume;

// PetscPrintf(PETSC_COMM_WORLD,"[%d].Rqn=%.5f,[%d].Rqw=%.5f\n",p_dof,qo_exact,s_dof-15,qw_exact);//No problem

    // Loop over faces of the element
    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        #if !Ghost
        // Boundary condition
        Point fc = face.compute_face_centroid();

        PetscReal p_bc = fluid.Rpw(t, fc.x(), fc.y());
        PetscReal pc_bc = fluid.Rpc(t, fc.x(), fc.y());
        f_array[p_dof] += face.compute_efftrans(mobo_exact) * (p + pc - pc_bc - p_bc);
        f_array[s_dof] += face.compute_efftrans(mobw_exact) * (p - p_bc);
        #elif Ghost
        // Boundary condition
        Point ac = face.compute_add_cell(e); 
        // PetscPrintf(PETSC_COMM_WORLD,"ac.x(%d)=%.5f, ac.y(%d)=%.5f\n",p_dof,ac.x(),p_dof,ac.y());//扩张节点没有问题
        PetscReal s_a_exact    = fluid.Rsw(t, ac.x(), ac.y());
        PetscReal p_a_exact    = fluid.Rpw(t, ac.x(), ac.y());
        PetscReal pc_a_exact   = fluid.Rpc(t, ac.x(), ac.y());
        PetscReal po_a_exact   = p_a_exact + pc_a_exact;
        PetscReal rhow_a_exact = fluid.func_rho_w(p_a_exact);
        PetscReal rhoo_a_exact = fluid.func_rho_o(p_a_exact);

        PetscScalar mobo_a_exact = fluid.func_lambda_o(s_a_exact); //mobility (lamada) = kr / mu
        PetscScalar mobw_a_exact = fluid.func_lambda_w(s_a_exact);

        PetscScalar mobo_a_aver = (mobo + mobo_a_exact) / 2;
        PetscScalar mobw_a_aver = (mobw + mobw_a_exact) / 2;
        // 1.2' 验证 mobw_a_aver 和 mobo_a_aver √
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_dof/2,mobw_a_aver);
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_dof/2,mobo_a_aver);

        //边界处也进行迎风处理，鬼点处取值为真解，边界点为数值解
        PetscScalar mobwup = (p > p_a_exact) ? rhow * mobw_a_aver : rhow_a_exact * mobw_a_aver; // 当前cell的压力大于邻居cell的压力时, 取当前cell的mob
        PetscScalar moboup = (po > po_a_exact) ? rhoo * mobo_a_aver : rhoo_a_exact * mobo_a_aver;
        // 1.3' 验证 mobwup √ 和 moboup √
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_dof/2,mobwup);
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_dof/2,moboup);

        // std::cout<<face.compute_trans()<<std::endl;
        f_array[p_dof] += face.compute_efftrans(moboup) * (po - po_a_exact) / 2; //p+pc=po
        f_array[s_dof] += face.compute_efftrans(mobwup) * (p - p_a_exact) / 2;
        // PetscScalar a=1./2 * 0.5;
        // std::cout<<a<<std::endl;
        // f_array[s_dof] += (a) * face.compute_efftrans(mobwup) * (p - p_a_exact);
        //注意，计算trans不能在前面乘1/2，结果会直接为0。
        //应该在后面 /2, 
        //或者写成1./2[但后者与lr师兄的代码出入较大，但不是数量级的]

        // 1.4' trans * 压力差
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_dof/2,face.compute_efftrans(mobwup) * (p - p_a_exact)/2);
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_dof/2,face.compute_efftrans(moboup) * (po - po_a_exact)/2);

        // 1.1' 验证 p - p_n 和 po - po_n √`
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_dof/2,p - p_a_exact);
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",(p_dof)/2,p_dof/2,po - po_a_exact);
        #endif
      }
      else
      {
        // Internal face, handle neighbor contributions
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);
        int s_n_dof = dof_map.dof_local_indices(*neighbor, s_var);

        PetscScalar p_n = x_array[p_n_dof];
        PetscScalar s_n = x_array[s_n_dof];
        PetscScalar pc_n = fluid.func_pc(s_n, s_old);
        PetscScalar po_n = p_n + pc_n;

        PetscScalar rho_on = fluid.func_rho_o(p_n);
        PetscScalar rho_wn = fluid.func_rho_w(p_n);
        PetscScalar mobon = fluid.func_lambda_o(s_n);
        PetscScalar mobwn = fluid.func_lambda_w(s_n);

        PetscScalar mobo_aver = (mobo + mobon) / 2;
        PetscScalar mobw_aver = (mobw + mobwn) / 2;
        // 1.2 验证 mobw_aver 和 mobo_aver √
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_n_dof/2,mobw_aver);
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_n_dof/2,mobo_aver);

        // Compute upstream mobilities and residual updates
        PetscScalar moboup = (po > po_n) ? rhoo * mobo_aver : rho_on * mobo_aver;
        PetscScalar mobwup = (p > p_n) ? rhow * mobw_aver : rho_wn * mobw_aver;
        // 1.3 验证 mobwup √ 和 moboup √
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_n_dof/2,mobwup);
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_n_dof/2,moboup);

        // std::cout<<face.compute_trans()<<std::endl;
        f_array[p_dof] += face.compute_efftrans(moboup) * (po - po_n);
        f_array[s_dof] += face.compute_efftrans(mobwup) * (p - p_n);

        // 1.4 trans * 压力差
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_n_dof/2,face.compute_efftrans(mobwup) * (p - p_n));
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_n_dof/2,face.compute_efftrans(moboup) * (po - po_n));
        
        // 1.1 验证 p - p_n 和 po - po_n √`
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",p_dof/2,p_n_dof/2,p - p_n);
        // PetscPrintf(PETSC_COMM_WORLD,"e[%d]-f[%d]=%.5f\n",(p_dof)/2,(p_n_dof)/2,po - po_n);
      }
    }
  }

  // Restore vector arrays
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
  PetscCall(MatZeroEntries(jac));
  DM dm_cp;
  System *sys;
  Vec x_local, x_old_local;
  PetscScalar *x_array, *x_old_array;

  // Get DM and System objects
  PetscCall(SNESGetDM(snes, &dm_cp));
  PetscCall(DMCPGetSystem_CP(dm_cp, sys));

  // Access system components
  DofMap &dof_map = sys->dof_map();
  Mesh &mesh = sys->get_mesh();
  Para &para = sys->para();
  CaseFluid fluid(para);

  // Get time and parameters
  PetscScalar dt = para.get_parameter<PetscScalar>("-dt");
  PetscScalar phi = para.get_parameter<PetscScalar>("-phi");
  double t = sys->time;

  // Get variable numbers for pressure and saturation
  int p_var = dof_map.variable_num("p");
  int s_var = dof_map.variable_num("s");

  // Convert global vectors to local vectors
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  sys->create_local_vec(&x_local);
  PetscCall(VecDuplicate(x_local, &x_old_local));
  PetscCall(globalvec_to_local(x, ghosted_index, &x_local));
  PetscCall(globalvec_to_local(sys->get_old_solution(), ghosted_index, &x_old_local));

  // Access vector arrays
  PetscCall(VecGetArray(x_local, &x_array));
  PetscCall(VecGetArray(x_old_local, &x_old_array));

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    PetscScalar Kpp = 0.0, Kps = 0.0, Ksp = 0.0, Kss = 0.0;
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int s_dof = dof_map.dof_local_indices(e, s_var);
    int p_g_dof = dof_map.dof_indices(e, p_var);
    int s_g_dof = dof_map.dof_indices(e, s_var);
    PetscScalar e_volume = e.compute_cell_volume();
    Point a = e.compute_cell_centroid();

    // Exact solutions for pressure, saturation, and fluxes
    PetscReal s_exact = fluid.Rsw(t, a.x(), a.y());
    PetscReal p_exact = fluid.Rpw(t, a.x(), a.y());
    PetscReal qw_exact = fluid.Rqw(t, a.x(), a.y());
    PetscReal qo_exact = fluid.Rqo(t, a.x(), a.y());
    PetscReal rhow_exact = fluid.Rrhow(t, a.x(), a.y());
    PetscReal rhoo_exact = fluid.Rrho_o(t, a.x(), a.y());
    PetscScalar mobo_exact = fluid.func_lambda_o(s_exact) * rhoo_exact;
    PetscScalar mobw_exact = fluid.func_lambda_w(s_exact) * rhow_exact;

    PetscScalar p = x_array[p_dof];
    PetscScalar s = x_array[s_dof];
    PetscScalar s_old = x_old_array[s_dof];
    PetscScalar p_old = x_old_array[p_dof];

    PetscScalar so = 1.0 - s;
    PetscScalar so_old = 1.0 - s_old;

    // Mobility and pressure calculations
    PetscScalar mobo = fluid.func_lambda_o(s);
    PetscScalar mobw = fluid.func_lambda_w(s);

    PetscScalar pc = fluid.func_pc(s, s_old);
    PetscScalar po = p + pc;

    PetscScalar rhow = fluid.func_rho_w(p);
    PetscScalar rhow_old = fluid.func_rho_w(p_old);
    PetscScalar rhoo = fluid.func_rho_o(p);
    PetscScalar rhoo_old = fluid.func_rho_o(p_old);
    // time term
    Kpp += phi / dt * fluid.func_drhoo_dp(p) * so * e_volume;
    Kps -= phi / dt * rhoo * e_volume;
    Ksp += phi / dt * fluid.func_drhow_dp(p) * s * e_volume;
    Kss += phi / dt * rhow * e_volume;

    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        #if !Ghost
        Kpp += face.compute_efftrans(mobo_exact);
        Kps += face.compute_efftrans(mobo_exact) * fluid.func_dpc_ds(s);
        Ksp += face.compute_efftrans(mobw_exact);
        #elif Ghost
        Point ac = face.compute_add_cell(e); 
        PetscReal s_a_exact    = fluid.Rsw(t, ac.x(), ac.y());
        PetscReal p_a_exact    = fluid.Rpw(t, ac.x(), ac.y());
        PetscReal pc_a_exact   = fluid.Rpc(t, ac.x(), ac.y());
        PetscReal po_a_exact   = p_a_exact + pc_a_exact;
        PetscReal rhow_a_exact = fluid.func_rho_w(p_a_exact);
        PetscReal rhoo_a_exact = fluid.func_rho_o(p_a_exact);

        PetscScalar mobo_a_exact = fluid.func_lambda_o(s_a_exact); //mobility (lamada) = kr / mu
        PetscScalar mobw_a_exact = fluid.func_lambda_w(s_a_exact);

        PetscScalar mobo_a_aver = (mobo + mobo_a_exact) / 2;
        PetscScalar mobw_a_aver = (mobw + mobw_a_exact) / 2;

        //边界处也进行迎风处理，鬼点处取值为真解，边界点为数值解
        PetscScalar mobwup = (p > p_a_exact) ? rhow * mobw_a_aver : rhow_a_exact * mobw_a_aver; // 当前cell的压力大于邻居cell的压力时, 取当前cell的mob
        PetscScalar moboup = (po > po_a_exact) ? rhoo * mobo_a_aver : rhoo_a_exact * mobo_a_aver;

        Kpp += face.compute_efftrans(moboup) / 2;
        Kps += face.compute_efftrans(moboup) * fluid.func_dpc_ds(s) / 2;
        Ksp += face.compute_efftrans(mobwup) / 2;
        #endif
      }
      else
      {
        PetscScalar Kpnp = 0.0, Kpns = 0.0, Knpp = 0.0, knps = 0.0, Ksnp = 0.0, Ksns = 0.0;
        // Internal face, handle neighbor contributions
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);
        int s_n_dof = dof_map.dof_local_indices(*neighbor, s_var);
        int p_n_g_dof = dof_map.dof_indices(*neighbor, p_var);
        int s_n_g_dof = dof_map.dof_indices(*neighbor, s_var);
        PetscScalar p_n = x_array[p_n_dof];
        PetscScalar s_n = x_array[s_n_dof];
        PetscScalar pc_n = fluid.func_pc(s_n);
        PetscScalar po_n = p_n + pc_n;

        PetscScalar rho_on = fluid.func_rho_o(p_n);
        PetscScalar rho_wn = fluid.func_rho_w(p_n);
        PetscScalar mobon = fluid.func_lambda_o(s_n);
        PetscScalar mobwn = fluid.func_lambda_w(s_n);

        PetscScalar mobo_aver = (mobo + mobon) / 2;
        PetscScalar mobw_aver = (mobw + mobwn) / 2;

        // Compute upstream mobilities and residual updates
        PetscScalar moboup = (po > po_n) ? rhoo * mobo_aver : rho_on * mobo_aver;
        PetscScalar mobwup = (p > p_n) ? rhow * mobw_aver : rho_wn * mobw_aver;
        
        Kpp += face.compute_efftrans(moboup);
        Kps += face.compute_efftrans(moboup) * fluid.func_dpc_ds(s);
        Kpnp -= face.compute_efftrans(moboup);
        Kpns -= face.compute_efftrans(moboup) * fluid.func_dpc_ds(s_n);
        Ksp += face.compute_efftrans(mobwup);
        Ksnp -= face.compute_efftrans(mobwup);

        if (po > po_n)
        {
          Kpp += fluid.func_drhoo_dp(p) * mobo * face.compute_trans() * (po - po_n);
          Kps += rhoo * fluid.func_dlambda_o_ds(s) * face.compute_trans() * (po - po_n);
        }
        else
        {
          Kpnp += fluid.func_drhoo_dp(p_n) * mobon * face.compute_trans() * (po - po_n);
          Kpns += rho_on * fluid.func_dlambda_o_ds(s_n) * face.compute_trans() * (po - po_n);
        }
        if (p > p_n)
        {
          Ksp += fluid.func_drhow_dp(p) * mobw * face.compute_trans() * (p - p_n);
          Kss += rhow * fluid.func_dlambda_w_ds(s) * face.compute_trans() * (p - p_n);
        }
        else
        {
          Ksnp += fluid.func_drhow_dp(p_n) * mobwn * face.compute_trans() * (p - p_n);
          Ksns += rho_wn * fluid.func_dlambda_w_ds(s_n) * face.compute_trans() * (p - p_n);
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
