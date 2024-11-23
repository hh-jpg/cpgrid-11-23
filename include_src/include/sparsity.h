/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-07-09 09:30:31
 * @LastEditTime: 2024-08-28 10:18:16
 * @FilePath: /cpgrid/include/sparsity.h
 * @Description:
 *
 */

#include "petscsys.h"
#include "petsc.h"
#include "utils.h"
#include "dof_map.h"
#include <iostream>
#include "config.h"
class Sparsity
{
private:
  DofMap &_dof_map;
  std::vector<PetscInt> _d_nnz;
  std::vector<PetscInt> _o_nnz;
  /* data */
public:
  Sparsity(DofMap &dofmap) : _dof_map(dofmap){prepare_sparsity();};
  void prepare_sparsity();
  PetscErrorCode creatmat(Mat&mat);
};
