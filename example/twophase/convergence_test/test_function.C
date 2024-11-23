/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-25 10:23:11
 * @LastEditTime: 2024-08-25 10:54:25
 * @FilePath: /cpgrid/example/twophase/convergence_test_1/test_function.C
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
#include "lirui.h"
PetscErrorCode set_para(Para &para);
void test_Rqw_Rqo_equality(CaseFluid &fluid, PetscReal t, PetscReal x, PetscReal y, PetscReal tolerance = 1e-6)
{
  // 计算Rqw和Rqo的值
  PetscReal result_Rqw = fluid.Rqw(t, x, y);
  PetscReal result_Rqo = fluid.Rqo(t, x, y);
  PetscReal x_loc = x;
  PetscReal y_loc = y;

  PetscReal qw = phi * (Rrhowdt(t, x_loc, y_loc) * Rsw(t, x_loc, y_loc) + Rrhow(t, x_loc, y_loc) * Rswdt(t, x_loc, y_loc)) +
                 (-2 * Rsw(t, x_loc, y_loc) / mu_w * Rswdx(t, x_loc, y_loc) * Rrhow(t, x_loc, y_loc) + (-pow2(Rsw(t, x_loc, y_loc)) / mu_w) * Rrhowdx(t, x_loc, y_loc)) * Rpwdx(t, x_loc, y_loc) + Rrhow(t, x_loc, y_loc) * (-pow2(Rsw(t, x_loc, y_loc)) / mu_w) * Rpwdxx(t, x_loc, y_loc) +
                 (-2 * Rsw(t, x_loc, y_loc) / mu_w * Rswdy(t, x_loc, y_loc) * Rrhow(t, x_loc, y_loc) + (-pow2(Rsw(t, x_loc, y_loc)) / mu_w) * Rrhowdy(t, x_loc, y_loc)) * Rpwdy(t, x_loc, y_loc) + Rrhow(t, x_loc, y_loc) * (-pow2(Rsw(t, x_loc, y_loc)) / mu_w) * Rpwdyy(t, x_loc, y_loc);
  PetscReal qo = phi * (Rrhondt(t, x_loc, y_loc) * (1 - Rsw(t, x_loc, y_loc)) + Rrhon(t, x_loc, y_loc) * (-Rswdt(t, x_loc, y_loc))) + (2 * (1 - Rsw(t, x_loc, y_loc)) / mu_n * Rswdx(t, x_loc, y_loc) * Rrhon(t, x_loc, y_loc) + (-pow2(1 - Rsw(t, x_loc, y_loc)) / mu_n) * Rrhondx(t, x_loc, y_loc)) * Rpndx(t, x_loc, y_loc) + Rrhon(t, x_loc, y_loc) * (-pow2(1 - Rsw(t, x_loc, y_loc)) / mu_n) * Rpndxx(t, x_loc, y_loc) +
                 (2 * (1 - Rsw(t, x_loc, y_loc)) / mu_n * Rswdy(t, x_loc, y_loc) * Rrhon(t, x_loc, y_loc) + (-pow2(1 - Rsw(t, x_loc, y_loc)) / mu_n) * Rrhondy(t, x_loc, y_loc)) * Rpndy(t, x_loc, y_loc) + Rrhon(t, x_loc, y_loc) * (-pow2(1 - Rsw(t, x_loc, y_loc)) / mu_n) * Rpndyy(t, x_loc, y_loc);

  // 输出结果
  std::cout << "Rqw: " << result_Rqw << " qw " << qw << std::endl;
  std::cout << "Rqo: " << result_Rqo << " qo " << qo << std::endl;

  // 判断两者是否相等
  if (std::fabs(result_Rqw - qw) < tolerance)
  {
    std::cout << "Test passed! Rqw and qw are equal within tolerance." << std::endl;
  }
  else
  {
    std::cout << "Test failed! Rqw and qw are not equal." << std::endl;
  }
}

int main(int argc, char **argv)
{
  PetscErrorCode ierr;
  PetscFunctionBeginUser;
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  // Mesh settings
  Mesh mesh;
  std::string case_name, mesh_name;
  char casename_char[256], meshname_char[256];
  PetscCall(PetscOptionsGetString(nullptr, nullptr, "-case_name", casename_char, sizeof(casename_char), nullptr));
  case_name = std::string(casename_char);
  PetscCall(PetscOptionsGetString(nullptr, nullptr, "-mesh_name", meshname_char, sizeof(meshname_char), nullptr));
  mesh_name = std::string(meshname_char);

  std::string filepath = "./grid/" + mesh_name;
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
  PetscCall(set_para(para));
  CaseFluid fluid(para);
  // 调用测试函数并增加输出
  std::cout << "Calling test_Rqw_Rqo_equality with t=0.01, x=0.5, y=0.6" << std::endl;
  test_Rqw_Rqo_equality(fluid, 0.01, 0.5, 0.6);

  std::cout << "Calling test_Rqw_Rqo_equality with t=0.02, x=0.3, y=0.7" << std::endl;
  test_Rqw_Rqo_equality(fluid, 0.02, 0.3, 0.7);

  std::cout << "Calling test_Rqw_Rqo_equality with t=0.05, x=0.6, y=0.4" << std::endl;
  test_Rqw_Rqo_equality(fluid, 0.05, 0.6, 0.4);

  std::cout << "Calling test_Rqw_Rqo_equality with t=0.1, x=0.4, y=0.5" << std::endl;
  test_Rqw_Rqo_equality(fluid, 0.1, 0.4, 0.5);

  std::cout << "Calling test_Rqw_Rqo_equality with t=0.15, x=0.2, y=0.8" << std::endl;
  test_Rqw_Rqo_equality(fluid, 0.15, 0.2, 0.8);

  std::cout << "Calling test_Rqw_Rqo_equality with t=0.25, x=0.7, y=0.3" << std::endl;
  test_Rqw_Rqo_equality(fluid, 0.25, 0.7, 0.3);

  PetscFunctionReturn(PETSC_SUCCESS);
}