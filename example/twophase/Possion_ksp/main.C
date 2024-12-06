/*
1D Incompressible Single Phase Flow Equation

1D Possion

div(grad p) = f,  0 < x,y < 1
     with
       forcing function f = -cos(m*pi*x),
       Neuman boundary conditions
        dp/dx = 0 for x = 0, x = 1.
The exact solution is \frac{1}{(m^2)\pi^2}\cos(m\pi x)
Here m = n = 1;

故左右边界处的 \int_{\Gamma_{-1,0}}\nabla p\cdot \vec n ds
            =\int_{\Gamma_{99,100}}\nabla p\cdot \vec n ds
            =0

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

#define pi 3.141596524
#define fr(x) PetscCosScalar(x*pi)
#define fe(x) PetscCosScalar(x*pi)/(pi*pi)

PetscErrorCode ComputeRHS(KSP , Vec , void *);
PetscErrorCode FormFunction(SNES, Vec, Vec, void *);
PetscErrorCode FormInitialGuess(SNES, Vec, void *);
PetscErrorCode FormJacobian(KSP, Mat, Mat, void *);
PetscErrorCode FormComputeExact(SNES nes, Vec , void *);
int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  char filename[256];
  PetscViewer viewer;

  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  // Mesh settings
  Mesh mesh;
  std::string case_name = "grid_possion/possion_128";
  std::string filepath = "./grid/" + case_name;

  mesh.read(filepath);
  mesh.init();
  mesh.printInfo();
  mesh.partition();
  mesh.reset_elem_id();

  // System settings
  System sys(mesh);
  sys.add_variable("p");
  sys.init();
  sys.set_name(case_name);

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
  PetscCall(SNESSetType(snes, SNESNEWTONLS));
  PetscCall(SNESSetDM(snes, dm_cp));

  // 创建全局向量
  Vec x, res, xexact;
  PetscCall(DMCreateGlobalVector(dm_cp, &x));
  PetscCall(VecDuplicate(x, &res));
  PetscCall(VecDuplicate(x, &xexact));

  // 设置初始猜测和函数
  PetscCall(FormInitialGuess(snes, x, &sys));
  PetscCall(SNESSetFunction(snes, res, FormFunction, &sys));

  // ierr = SNESComputeFunction(snes, x, res);
  // CHKERRQ(ierr);
  // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.m", &viewer);
  // VecView(x, viewer);
  // PetscViewerASCIIOpen(PETSC_COMM_WORLD, "resinit.m", &viewer);
  // VecView(res, viewer);

  // Mat j;
  // sys.create_mat(j);
  // PetscCall(SNESSetJacobian(snes, j, j, FormJacobian, NULL));

  // 从命令行选项设置 SNES
  PetscCall(SNESSetFromOptions(snes));
  PetscCall(SNESSetUp(snes));
  // return 0;
  // PetscCall(SNESSolve(snes, NULL, x));

  KSP ksp;
  SNESGetKSP(snes,&ksp);
  KSPSetUp(ksp);
  PetscCall(KSPSetComputeRHS(ksp, ComputeRHS, NULL));
  PetscCall(KSPSetComputeOperators(ksp,FormJacobian , NULL));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(KSPSolve(ksp, NULL, x));





  // // 获取并保存解
  // PetscCall(SNESGetSolution(snes, &x));
  // ierr = SNESComputeFunction(snes, x, res);
  // CHKERRQ(ierr);

  // ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "res.m", &viewer);
  // CHKERRQ(ierr);
  // VecView(res, viewer);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.m", &viewer);
  CHKERRQ(ierr);
  VecView(x, viewer);

  // 计算真解
  PetscCall(FormComputeExact(snes, xexact, &sys));
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD, "xexact.m", &viewer);
  CHKERRQ(ierr);
  VecView(xexact, viewer);

  PetscScalar residual_x,residual_xexact,norm2,norm1,norminf;
	ierr = VecAXPY(x,-1.0,xexact); CHKERRQ(ierr);
	ierr = VecNorm(x,NORM_2,&residual_x);CHKERRQ(ierr);
	ierr = VecNorm(xexact,NORM_2,&residual_xexact);CHKERRQ(ierr);
  norm2 = residual_x / residual_xexact;
	PetscPrintf(PETSC_COMM_WORLD,"The residual NORM_2 of solution =%.5e\n",residual_x/residual_xexact);
  
  ierr = VecNorm(x,NORM_1,&residual_x);CHKERRQ(ierr);
	ierr = VecNorm(xexact,NORM_1,&residual_xexact);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The residual NORM_1 of solution =%.5e\n",residual_x/residual_xexact);
  norm1 = residual_x / residual_xexact;

  ierr = VecNorm(x,NORM_INFINITY,&residual_x);CHKERRQ(ierr);
	ierr = VecNorm(xexact,NORM_INFINITY,&residual_xexact);CHKERRQ(ierr);
	PetscPrintf(PETSC_COMM_WORLD,"The residual NORM_inf of solution =%.5e\n",residual_x/residual_xexact);
  norminf = residual_x / residual_xexact;

  PetscPrintf(PETSC_COMM_WORLD,"& %.5e &  & %.5e &  & %.5e &  \\\\ \n",norminf,norm2,norm1);

  // 清理资源
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&xexact));
  PetscCall(VecDestroy(&res));
  PetscCall(SNESDestroy(&snes));
  PetscCall(DMDestroy(&dm_cp));
  PetscCall(PetscFinalize());
  return 0;
}

PetscErrorCode ComputeRHS(KSP ksp, Vec b, void *ctx)
{
  PetscFunctionBeginUser;
  // PetscCall(VecZeroEntries(f));
  DM dm_cp;
  System *sys;
  Vec x_local;
  PetscScalar *b_array;
  PetscErrorCode ierr;
  MatNullSpace nullspace;

  // 获取DM和System
  PetscCall(KSPGetDM(ksp, &dm_cp));
  PetscCall(DMCPGetSystem_CP(dm_cp, sys));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = sys->dof_map();
  Mesh &mesh = sys->get_mesh();

  // 获取变量编号
  int p_var = dof_map.variable_num("p");

  // 获取向量数组
  PetscCall(VecGetArray(b, &b_array));

  // 遍历网格中的单元
  for (const Polyhedron &e : mesh.elem_local_range())
  {
    // PetscPrintf(PETSC_COMM_WORLD,"-----======-----\n");
    double perm = e.perm();
    PetscScalar e_volume = e.compute_cell_volume();
    Point cent=e.compute_cell_centroid();
    // std::cout<<cent.x()<<std::endl;   // 输出当前cell的x坐标
    int p_dof = dof_map.dof_local_indices(e, p_var);

    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 1)
        {
          b_array[p_dof] =  -face.compute_efftrans(1) / perm * (1/(pi*pi));
        }
      }
    }
    b_array[p_dof] = fr(cent.x()) * e_volume;

    // PetscPrintf(PETSC_COMM_WORLD,"p_%d=%.5e\n",p_dof,(p));
  }

  // 恢复向量数组
  PetscCall(VecRestoreArray(b, &b_array));
  PetscCall(VecAssemblyBegin(b));
  PetscCall(VecAssemblyEnd(b));
  
  // 处理奇异矩阵的右端项
  PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace));
  PetscCall(MatNullSpaceRemove(nullspace, b));
  PetscCall(MatNullSpaceDestroy(&nullspace));
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecZeroEntries(f));
  DM dm_cp;
  System *sys;
  Vec x_local;
  PetscScalar *f_array, *x_array;
  PetscErrorCode ierr;

  // 获取DM和System
  PetscCall(SNESGetDM(snes, &dm_cp));
  PetscCall(DMCPGetSystem_CP(dm_cp, sys));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = sys->dof_map();
  Mesh &mesh = sys->get_mesh();

  // 获取变量编号
  int p_var = dof_map.variable_num("p");
  // std::cout<<"p_var="<<p_var<<std::endl;

  // 将全局向量转换为局部向量
  const std::vector<unsigned int> &ghosted_index = dof_map.get_send_list();
  sys->create_local_vec(&x_local);
  PetscCall(globalvec_to_local(x, ghosted_index, &x_local));

  // 获取向量数组
  PetscCall(VecGetArray(x_local, &x_array));
  PetscCall(VecGetArray(f, &f_array));

  // 遍历网格中的单元
  for (const Polyhedron &e : mesh.elem_local_range())
  {
    // PetscPrintf(PETSC_COMM_WORLD,"-----======-----\n");
    double perm = e.perm();
    PetscScalar e_volume = e.compute_cell_volume();
    Point cent=e.compute_cell_centroid();
    // std::cout<<cent.x()<<std::endl;   // 输出当前cell的x坐标
    int p_dof = dof_map.dof_local_indices(e, p_var);
    PetscScalar p = x_array[p_dof];

    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 1)
        {
          f_array[p_dof] +=  face.compute_efftrans(1) / perm * (p + 1/(pi*pi));
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);
        PetscScalar p_n = x_array[p_n_dof];
        // PetscPrintf(PETSC_COMM_WORLD,"p_n=%g\n",p_n);
        f_array[p_dof] += face.compute_efftrans(1) / perm * (p - p_n);
      }
    }
        f_array[p_dof] -= fr(cent.x()) * e_volume;

    // PetscPrintf(PETSC_COMM_WORLD,"p_%d=%.5e\n",p_dof,(p));
  }

  // 恢复向量数组
  PetscCall(VecRestoreArray(x_local, &x_array));
  PetscCall(VecRestoreArray(f, &f_array));
  PetscCall(VecDestroy(&x_local));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormJacobian(KSP ksp,  Mat J, Mat jac, void *ctx)
{
  PetscFunctionBeginUser;
  DM dm_cp;
  PetscCall(MatZeroEntries(jac));
  System *sys;
  Vec x_local, x_old_local;
  PetscScalar *x_array;
  PetscErrorCode ierr;
  MatNullSpace nullspace;

  // 获取DM和System
  PetscCall(KSPGetDM(ksp, &dm_cp));
  PetscCall(DMCPGetSystem_CP(dm_cp, sys));

  // 获取系统的DofMap和Mesh
  DofMap &dof_map = sys->dof_map();
  Mesh &mesh = sys->get_mesh();
  // 获取变量编号
  int p_var = dof_map.variable_num("p");

  for (const Polyhedron &e : mesh.elem_local_range())
  {
    PetscScalar Kpp = 0.0, Kps = 0.0, Ksp = 0.0, Kss = 0.0;
    int p_dof = dof_map.dof_local_indices(e, p_var);
    int p_g_dof = dof_map.dof_indices(e, p_var);
    PetscScalar e_volume = e.compute_cell_volume();
    double perm = e.perm();

    PetscScalar p = x_array[p_dof];
    for (const Face &face : e.get_faces())
    {
      if (face.neighbor() == nullptr)
      {
        if (face.pos() == 1)
        {
          Kpp += face.compute_efftrans(1)/perm;
        }
      }
      else
      {
        Polyhedron *neighbor = face.neighbor();
        int p_n_dof = dof_map.dof_local_indices(*neighbor, p_var);
        int p_n_g_dof = dof_map.dof_indices(*neighbor, p_var);
        PetscScalar Kpnp = 0.0, Kpns = 0.0, Knpp = 0.0, knps = 0.0, Ksnp = 0.0, Ksns = 0.0;
        PetscScalar p_n = x_array[p_n_dof];

        Kpp += face.compute_efftrans(1)/perm;
        Kpnp += - face.compute_efftrans(1)/perm;
        // std::cout<<Kpnp<<" ";
        
        MatSetValue(jac, p_g_dof, p_n_g_dof, Kpnp, ADD_VALUES);
      }
    }
    // std::cout<<Kpp<<std::endl;
    MatSetValue(jac, p_g_dof, p_g_dof, Kpp, ADD_VALUES);
  }
  PetscCall(MatAssemblyBegin(jac, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(jac, MAT_FINAL_ASSEMBLY));

  PetscCall(MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, 0, &nullspace));
  PetscCall(MatSetNullSpace(J, nullspace));
  PetscCall(MatNullSpaceDestroy(&nullspace));
  PetscFunctionReturn(PETSC_SUCCESS);
}


PetscErrorCode FormInitialGuess(SNES snes, Vec x, void *ctx)
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
    x_array[p_dof] = 1;
  }
  PetscCall(VecRestoreArray(x, &x_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode FormComputeExact(SNES snes, Vec x, void *ctx)
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
    Point cent=e.compute_cell_centroid();
    x_array[p_dof] = fe(cent.x());
  }
  PetscCall(VecRestoreArray(x, &x_array));
  PetscFunctionReturn(PETSC_SUCCESS);
}

