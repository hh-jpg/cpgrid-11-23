/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-22 11:20:35
 * @LastEditTime: 2024-08-28 20:53:32
 * @FilePath: /cpgrid/example/twophase/convergence_test_1/case_fluid.h
 * @Description:this class include the fluid
 *
 */
#ifndef CASEFLUID_H
#define CASEFLUID_H



#include "fluid.h"

class CaseFluid : public Fluid
{
public:
  CaseFluid(Para &para);
  // Destructor
  virtual ~CaseFluid() = default;

  PetscReal func_pc(const PetscReal s, const PetscReal s_old = 0.) override;
  PetscReal func_rho_w(const PetscReal p);
  PetscReal func_rho_o(const PetscReal p);
  PetscReal func_dpc_ds(const PetscReal s, const PetscReal s_old = 0.) override;
  PetscReal func_drhow_dp(const PetscReal p);
  PetscReal func_drhoo_dp(const PetscReal p);
  // Rdsw 系列函数
  PetscReal Rsw(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdswdx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdswdy(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdswdt(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddswdx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddswdy(PetscReal t, PetscReal x, PetscReal y) const;

  // Rpw 系列函数
  PetscReal Rpw(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpwdx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddpwdx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpwdy(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddpwdy(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpwdt(PetscReal t, PetscReal x, PetscReal y) const;

  // Rpc 系列函数
  PetscReal Rpc(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpcdx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddpcdx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpcdy(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddpcdy(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpcdt(PetscReal t, PetscReal x, PetscReal y) const;

  // Rrhow 系列函数
  PetscReal Rrhow(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdrhowdx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdrhowdy(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdrhowdt(PetscReal t, PetscReal x, PetscReal y) const;

  // Rrhon 系列函数
  PetscReal Rrho_o(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdrho_odx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdrho_ody(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdrho_odt(PetscReal t, PetscReal x, PetscReal y) const;
  // Rpo 系列函数
  PetscReal Rpo(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpodt(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpodx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddpodx(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rdpody(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rddpody(PetscReal t, PetscReal x, PetscReal y) const;
 
  // Rqw and Rqo
  PetscReal Rqw(PetscReal t, PetscReal x, PetscReal y) const;
  PetscReal Rqo(PetscReal t, PetscReal x, PetscReal y) const;
  
private:
  // 成员变量
  PetscScalar rho_w;
  PetscScalar cf_w;
  PetscScalar rho_o;
  PetscScalar cf_o;
  PetscScalar ref_pw;
  PetscScalar ref_po;
  PetscScalar Bc;
};

#endif // CASEFLUID_H
