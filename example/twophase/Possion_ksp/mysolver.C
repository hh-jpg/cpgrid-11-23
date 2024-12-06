/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-10-22 09:30:50
 * @LastEditTime: 2024-11-13 11:29:27
 * @FilePath: /cpgrid/src/mysolver.C
 * @Description:
 *
 */
#include "mysolver.h"
#include <petscsnes.h>
#include <petscksp.h>
#include "mat_func.h"
#include "vec_func.h"
#include <petscsys.h>
#include <petscvec.h>
#include <petscksp.h>
#include <petscsnes.h>
//int a = 0;
PetscErrorCode mynewtonsolver(SNES snes, Vec x)
{
  PetscFunctionBeginUser;
  Vec F, Jdx;
  Mat J, Jpre;
  KSP ksp;
  PetscInt nit = 0, maxits = 100;
  PetscReal tol = 1e-5, normF = 1.e99;
  KSPConvergedReason reason;

  // 从 SNES 获取雅可比矩阵
  PetscCall(SNESGetJacobian(snes, &J, &Jpre, NULL, NULL));
  PetscCall(SNESGetKSP(snes, &ksp));

  // 创建向量 F 和 Jdx
  PetscCall(VecDuplicate(x, &F));
  PetscCall(VecDuplicate(x, &Jdx));

  // 迭代求解
  while (normF > tol && nit <= maxits)
  {
    // 计算雅可比矩阵和函数值
    PetscCall(SNESComputeJacobian(snes, x, J, Jpre));
    PetscCall(SNESComputeFunction(snes, x, F));
    /*if (!a)
    {
      PetscPrintf(PETSC_COMM_WORLD, "start write the jac and vec ");
      SaveVecToMatlab(x, "/Users/jltu/code/reservoir_simulation/mymatlab/twophase_tpfa/egg_nonuniform/x_1.m", "vec_1");
      SaveVecToMatlab(F, "/Users/jltu/code/reservoir_simulation/mymatlab/twophase_tpfa/egg_nonuniform/f_1.m", "res_1");
      SaveMatToMatlab(J, "/Users/jltu/code/reservoir_simulation/mymatlab/twophase_tpfa/egg_nonuniform/exact_jac_1.m", "jac_1");
    }*/
    // 求解线性系统 J * dx = -F
    PetscCall(KSPSetOperators(ksp, J, Jpre)); // 使用 J 和 Jpre 作为操作矩阵
    PetscCall(KSPSolve(ksp, F, Jdx));
    /*if (!a)
    {
      SaveVecToMatlab(Jdx, "/Users/jltu/code/reservoir_simulation/mymatlab/twophase_tpfa/egg_nonuniform/jdx_1.m", "jx_1");
      a++;
    }*/
    // 检查 KSP 的收敛状态
    PetscCall(KSPGetConvergedReason(ksp, &reason));
    if (reason < 0)
    {
      // 如果 KSP 不收敛，返回错误
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "KSP did not converge at iteration %D\n", nit));
      PetscCall(VecDestroy(&F));
      PetscCall(VecDestroy(&Jdx));
      return PETSC_ERR_NOT_CONVERGED; // 返回未收敛错误
    }

    // 更新解 x = x - dx
    PetscCall(VecAXPY(x, -1.0, Jdx));
    PetscCall(VecNorm(F, NORM_2, &normF));
    PetscPrintf(PETSC_COMM_WORLD, "Nonlinear Iteration %D: Res = %1.15e\n", nit, (double)normF);
    // 增加迭代次数
    nit++;
  }
  SNESComputeFunction(snes, x, F);
  PetscCall(VecNorm(F, NORM_2, &normF));

  if (normF <= tol)
  {
    PetscPrintf(PETSC_COMM_WORLD, "Converged at iteration %D with norm %g\n", nit, (double)normF);
  }
  else
  {
    PetscPrintf(PETSC_COMM_WORLD, "Reached maximum iterations without convergence\n");
  }
  PetscPrintf(PETSC_COMM_WORLD, "Nonlinear Iteration %D: true Res = %1.15e\n", nit, (double)normF);

  // 清理
  PetscCall(VecDestroy(&F));
  PetscCall(VecDestroy(&Jdx));
  PetscFunctionReturn(PETSC_SUCCESS);
}
