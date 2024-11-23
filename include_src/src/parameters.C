/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-05-15 09:44:49
 * @LastEditTime: 2024-08-06 16:12:46
 * @FilePath: /cpgrid/src/parameters.C
 * @Description:
 *
 */
#include "parameters.h"
#include <petscsys.h>
#include <array>
#include <cstring>


PetscErrorCode Para::print_parameters() const
{
  PetscFunctionBeginUser;
  for (const auto &param : _parameters)
  {
    PetscPrintf(PETSC_COMM_WORLD, "%s: ", param.first.c_str());
    if (param.second.type() == typeid(PetscInt))
    {
      PetscPrintf(PETSC_COMM_WORLD, "%d\n", std::any_cast<PetscInt>(param.second));
    }
    else if (param.second.type() == typeid(PetscReal) || param.second.type() == typeid(double) || param.second.type() == typeid(PetscScalar))
    {
      PetscPrintf(PETSC_COMM_WORLD, "%g\n", std::any_cast<PetscReal>(param.second));
    }
    else if (param.second.type() == typeid(std::string))
    {
      PetscPrintf(PETSC_COMM_WORLD, "%s\n", std::any_cast<std::string>(param.second).c_str());
    }
    else
    {
      PetscPrintf(PETSC_COMM_WORLD, "Unknown type\n");
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}
