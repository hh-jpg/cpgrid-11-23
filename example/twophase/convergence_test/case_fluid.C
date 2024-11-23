/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-22 11:24:25
 * @LastEditTime: 2024-08-29 09:47:18
 * @FilePath: /cpgrid/example/twophase/convergence_test_1/case_fluid.C
 * @Description:
 *
 */
#include "case_fluid.h"
#include "fluid.h"
#include <math.h>
CaseFluid::CaseFluid(Para &para)
    : Fluid(para) // 调用父类的构造函数
{
  // 初始化子类独有的成员变量
  rho_w = _para.get_parameter<PetscScalar>("-rho_w");
  cf_w = _para.get_parameter<PetscScalar>("-cf_w");
  rho_o = _para.get_parameter<PetscScalar>("-rho_o");
  cf_o = _para.get_parameter<PetscScalar>("-cf_o");
  ref_pw = _para.get_parameter<PetscScalar>("-ref_pw");
  ref_po = _para.get_parameter<PetscScalar>("-ref_po");
  Bc = _para.get_parameter<PetscScalar>("-Bc");
}

PetscReal CaseFluid::func_pc(const PetscReal s, const PetscReal s_old)
{
  PetscReal se = std::max(func_se(s), 1.e-10);
  return -Bc * log(se);
}

PetscReal CaseFluid::func_rho_w(const PetscReal p)
{
  PetscScalar p_ref = 0.0;
  // Output the values
  
  return rho_w * std::exp(cf_w * (p - p_ref));  //-cf_w = 1.e-12, -rho_w = 1000   
}
     
PetscReal CaseFluid::func_rho_o(const PetscReal p)
{
  return rho_o; //-rho_o = 1000 
}

PetscReal CaseFluid::func_dpc_ds(const PetscReal s, const PetscReal s_old)
{
  PetscReal se = func_se(s);
  PetscReal dse_ds = func_dse_ds();
  return -Bc * (1.0 / se) * dse_ds;
}

PetscReal CaseFluid::func_drhow_dp(const PetscReal p)
{
  PetscScalar p_ref = 0.0;
  return rho_w * cf_w * std::exp(cf_w * (p - p_ref));
}

PetscReal CaseFluid::func_drhoo_dp(const PetscReal p)
{
  return 0.;
}

PetscReal CaseFluid::Rsw(PetscReal t, PetscReal x, PetscReal y) const
{
  return sin(t + x - y + 1);
}

PetscReal CaseFluid::Rdswdx(PetscReal t, PetscReal x, PetscReal y) const
{
  return cos(t + x - y + 1);
}

PetscReal CaseFluid::Rdswdy(PetscReal t, PetscReal x, PetscReal y) const
{
  return -cos(t + x - y + 1);
}

PetscReal CaseFluid::Rdswdt(PetscReal t, PetscReal x, PetscReal y) const
{
  return cos(t + x - y + 1);
}

PetscReal CaseFluid::Rddswdx(PetscReal t, PetscReal x, PetscReal y) const
{
  return -sin(t + x - y + 1);
}

PetscReal CaseFluid::Rddswdy(PetscReal t, PetscReal x, PetscReal y) const
{
  return -sin(t + x - y + 1);
}

// Rdpw 系列函数
PetscReal CaseFluid::Rpw(PetscReal t, PetscReal x, PetscReal y) const
{
  return cos(t + x - y);
}

PetscReal CaseFluid::Rdpwdx(PetscReal t, PetscReal x, PetscReal y) const
{
  return -sin(t + x - y);
}

PetscReal CaseFluid::Rddpwdx(PetscReal t, PetscReal x, PetscReal y) const
{
  return -cos(t + x - y);
}

PetscReal CaseFluid::Rdpwdy(PetscReal t, PetscReal x, PetscReal y) const
{
  return sin(t + x - y);
}

PetscReal CaseFluid::Rddpwdy(PetscReal t, PetscReal x, PetscReal y) const
{
  return -cos(t + x - y);
}

PetscReal CaseFluid::Rdpwdt(PetscReal t, PetscReal x, PetscReal y) const
{
  return -sin(t + x - y);
}

// Rpc 系列函数
PetscReal CaseFluid::Rpc(PetscReal t, PetscReal x, PetscReal y) const
{
  return -Bc * log(Rsw(t, x, y));
}

PetscReal CaseFluid::Rdpcdx(PetscReal t, PetscReal x, PetscReal y) const
{
  return -Bc * Rdswdx(t, x, y) / Rsw(t, x, y);
}

PetscReal CaseFluid::Rddpcdx(PetscReal t, PetscReal x, PetscReal y) const
{
  return Bc * pow(Rdswdx(t, x, y), 2) / pow(Rsw(t, x, y), 2) - Bc * Rddswdx(t, x, y) / Rsw(t, x, y);
}

PetscReal CaseFluid::Rdpcdy(PetscReal t, PetscReal x, PetscReal y) const
{
  return -Bc * Rdswdy(t, x, y) / Rsw(t, x, y);
}

PetscReal CaseFluid::Rddpcdy(PetscReal t, PetscReal x, PetscReal y) const
{
  return Bc * pow(Rdswdy(t, x, y), 2) / pow(Rsw(t, x, y), 2) - Bc * Rddswdy(t, x, y) / Rsw(t, x, y);
}

PetscReal CaseFluid::Rdpcdt(PetscReal t, PetscReal x, PetscReal y) const
{
  return -Bc * Rdswdt(t, x, y) / Rsw(t, x, y);
}

// Rrhow 系列函数
PetscReal CaseFluid::Rrhow(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_w * std::exp(cf_w * (Rpw(t, x, y) - ref_pw));
}

PetscReal CaseFluid::Rdrhowdt(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_w * std::exp(cf_w * (Rpw(t, x, y) - ref_pw)) * cf_w * Rdpwdt(t, x, y);
}

PetscReal CaseFluid::Rdrhowdx(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_w * std::exp(cf_w * (Rpw(t, x, y) - ref_pw)) * cf_w * Rdpwdx(t, x, y);
}

PetscReal CaseFluid::Rdrhowdy(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_w * std::exp(cf_w * (Rpw(t, x, y) - ref_pw)) * cf_w * Rdpwdy(t, x, y);
}

// Rrho_o 系列函数
PetscReal CaseFluid::Rrho_o(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_o * std::exp(cf_o * (Rpo(t, x, y) - ref_po));
}

PetscReal CaseFluid::Rdrho_odt(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_o * std::exp(cf_o * (Rpo(t, x, y) - ref_po)) * cf_o * Rdpodt(t, x, y);
}

PetscReal CaseFluid::Rdrho_odx(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_o * std::exp(cf_o * (Rpo(t, x, y) - ref_po)) * cf_o * Rdpodx(t, x, y);
}

PetscReal CaseFluid::Rdrho_ody(PetscReal t, PetscReal x, PetscReal y) const
{
  return rho_o * std::exp(cf_o * (Rpo(t, x, y) - ref_po)) * cf_o * Rdpody(t, x, y);
}

// Rpo 系列函数
PetscReal CaseFluid::Rpo(PetscReal t, PetscReal x, PetscReal y) const
{
  return Rpw(t, x, y) + Rpc(t, x, y);
}

PetscReal CaseFluid::Rdpodt(PetscReal t, PetscReal x, PetscReal y) const
{
  return Rdpwdt(t, x, y) + Rdpcdt(t, x, y);
}

PetscReal CaseFluid::Rdpodx(PetscReal t, PetscReal x, PetscReal y) const
{
  return Rdpwdx(t, x, y) + Rdpcdx(t, x, y);
}

PetscReal CaseFluid::Rddpodx(PetscReal t, PetscReal x, PetscReal y) const
{
  return Rddpwdx(t, x, y) + Rddpcdx(t, x, y);
}

PetscReal CaseFluid::Rdpody(PetscReal t, PetscReal x, PetscReal y) const
{
  return Rdpwdy(t, x, y) + Rdpcdy(t, x, y);
}

PetscReal CaseFluid::Rddpody(PetscReal t, PetscReal x, PetscReal y) const
{
  return Rddpwdy(t, x, y) + Rddpcdy(t, x, y);
}

PetscReal CaseFluid::Rqw(PetscReal t, PetscReal x, PetscReal y) const
{
  PetscReal phi = 0.8;

  // 第一项
  PetscReal term1 = phi * (Rdrhowdt(t, x, y) * Rsw(t, x, y) + Rrhow(t, x, y) * Rdswdt(t, x, y));

  // 第二项
  PetscReal term2a = -2 * Rsw(t, x, y) / _mu_w * Rdswdx(t, x, y) * Rrhow(t, x, y);
  PetscReal term2b = (-pow(Rsw(t, x, y), 2) / _mu_w) * Rdrhowdx(t, x, y);
  PetscReal term2 = (term2a + term2b) * Rdpwdx(t, x, y);

  // 第三项
  PetscReal term3 = Rrhow(t, x, y) * (-pow(Rsw(t, x, y), 2) / _mu_w) * Rddpwdx(t, x, y);

  // 第四项
  PetscReal term4a = -2 * Rsw(t, x, y) / _mu_w * Rdswdy(t, x, y) * Rrhow(t, x, y);
  PetscReal term4b = (-pow(Rsw(t, x, y), 2) / _mu_w) * Rdrhowdy(t, x, y);
  PetscReal term4 = (term4a + term4b) * Rdpwdy(t, x, y);

  // 第五项
  PetscReal term5 = Rrhow(t, x, y) * (-pow(Rsw(t, x, y), 2) / _mu_w) * Rddpwdy(t, x, y);

  // 总和
  return term1 + term2 + term3 + term4 + term5;
}

PetscReal CaseFluid::Rqo(PetscReal t, PetscReal x, PetscReal y) const
{
  PetscReal phi = 0.8;

  // 第一项
  PetscReal term1 = phi * (Rdrho_odt(t, x, y) * (1 - Rsw(t, x, y)) + Rrho_o(t, x, y) * (-Rdswdt(t, x, y)));

  // 第二项
  PetscReal term2a = 2 * (1 - Rsw(t, x, y)) / _mu_o * Rdswdx(t, x, y) * Rrho_o(t, x, y);
  PetscReal term2b = (-pow(1 - Rsw(t, x, y), 2) / _mu_o) * Rdrho_odx(t, x, y);
  PetscReal term2 = (term2a + term2b) * Rdpodx(t, x, y);

  // 第三项
  PetscReal term3 = Rrho_o(t, x, y) * (-pow(1 - Rsw(t, x, y), 2) / _mu_o) * Rddpodx(t, x, y);

  // 第四项
  PetscReal term4a = 2 * (1 - Rsw(t, x, y)) / _mu_o * Rdswdy(t, x, y) * Rrho_o(t, x, y);
  PetscReal term4b = (-pow(1 - Rsw(t, x, y), 2) / _mu_o) * Rdrho_ody(t, x, y);
  PetscReal term4 = (term4a + term4b) * Rdpody(t, x, y);

  // 第五项
  PetscReal term5 = Rrho_o(t, x, y) * (-pow(1 - Rsw(t, x, y), 2) / _mu_o) * Rddpody(t, x, y);

  // 总和
  return term1 + term2 + term3 + term4 + term5;
}
