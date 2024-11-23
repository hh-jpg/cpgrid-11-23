/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-04 16:47:31
 * @LastEditTime: 2024-08-28 20:15:21
 * @FilePath: /cpgrid/src/fluid.C
 * @Description: 
 * 
 */
#include "fluid.h"
Fluid::Fluid(Para &para)
    : _para(para)
{
  // Cache the parameters on initialization to avoid repeated lookups
  _srw = _para.get_parameter<PetscReal>("-srw");
  _srn = _para.get_parameter<PetscReal>("-srn");
  _mu_w = _para.get_parameter<PetscReal>("-mu_w");
  _mu_o = _para.get_parameter<PetscReal>("-mu_o");
  _beta = _para.get_parameter<PetscReal>("-beta");
}

PetscReal Fluid::func_se(const PetscReal s)
{
  return (s - _srw) / (1. - _srw - _srn);
}

PetscReal Fluid::func_dse_ds() const
{
  return 1. / (1. - _srw - _srn);
}

PetscReal Fluid::func_pc(const PetscReal s, const PetscReal s_old)
{
  return 50. * std::pow(s, -0.5);
}

PetscReal Fluid::func_dpc_ds(const PetscReal s, const PetscReal s_old)
{
  return -25.0 * std::pow(s, -1.5);
}

PetscReal Fluid::func_ddpc_ds(const PetscReal s, const PetscReal s_old)
{
  return 37.5 * std::pow(s, -2.5);
}

PetscReal Fluid::func_lambda_w(const PetscReal s)
{
  PetscReal se = func_se(s);
  PetscReal krw = std::pow(se, _beta);
  return krw / _mu_w;
}

PetscReal Fluid::func_dlambda_w_ds(const PetscReal s)
{
  PetscReal se = func_se(s);
  PetscReal dse_ds = func_dse_ds();
  PetscReal dkrw_dse = _beta * std::pow(se, _beta - 1);
  return dkrw_dse * dse_ds / _mu_w;
}

PetscReal Fluid::func_lambda_o(const PetscReal s)
{
  PetscReal se = func_se(s);
  PetscReal krn = std::pow((1 - se), _beta);
  return krn / _mu_o;
}

PetscReal Fluid::func_dlambda_o_ds(const PetscReal s)
{
  PetscReal se = func_se(s);
  PetscReal dse_ds = func_dse_ds();
  PetscReal dkrn_dse = _beta * std::pow(1 - se, _beta - 1) * (-1);
  return dkrn_dse * dse_ds / _mu_o;
}