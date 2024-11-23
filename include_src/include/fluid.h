/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-03 15:53:37
 * @LastEditTime: 2024-08-28 20:14:32
 * @FilePath: /cpgrid/include/fluid.h
 * @Description:
 *
 */
#ifndef FLUID_H
#define FLUID_H

#include "parameters.h"
#include <cmath>
#include "config.h"
class Fluid
{
public:
  Fluid(Para &para);

  virtual ~Fluid() = default;

  // Capillary pressure and its derivatives
  virtual PetscReal func_pc(const PetscReal s, const PetscReal s_old);
  virtual PetscReal func_dpc_ds(const PetscReal s, const PetscReal s_old);
  virtual PetscReal func_ddpc_ds(const PetscReal s, const PetscReal s_old);

  // Water phase mobility and its derivative
  virtual PetscReal func_lambda_w(const PetscReal s);
  virtual PetscReal func_dlambda_w_ds(const PetscReal s);

  // Oil phase mobility and its derivative
  virtual PetscReal func_lambda_o(const PetscReal s);
  virtual PetscReal func_dlambda_o_ds(const PetscReal s);

protected:
  Para &_para;

  // Cached parameters
  PetscReal _srw;
  PetscReal _srn;
  PetscReal _mu_w;
  PetscReal _mu_o;
  PetscReal _beta;

  // Helper methods for common expressions
  PetscReal func_se(const PetscReal s);
  PetscReal func_dse_ds() const;
};
#endif // FLUID_H
