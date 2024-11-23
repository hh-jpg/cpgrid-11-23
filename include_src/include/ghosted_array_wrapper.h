/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-03 10:59:13
 * @LastEditTime: 2024-08-03 11:22:38
 * @FilePath: /cpgrid/include/ghosted_array_wrapper.h
 * @Description:
 *
 */

#ifndef GHOSTEDARRAYWRAPPER_H
#define GHOSTEDARRAYWRAPPER_H

#include "petscsys.h"
#include "dof_map.h"
#include "petscvec.h"
#include <unordered_map>
#include <stdexcept>
#include <iostream>
#include "config.h"

class GhostedArrayWrapper
{
public:
  GhostedArrayWrapper(Vec x_local, DofMap &dof_map);
  ~GhostedArrayWrapper();

  const PetscScalar &operator[](PetscInt global_index) const;
  void print(int num_elements) const;

private:
  Vec x_local;
  const PetscScalar *ghosted_array;
  DofMap &dofmap;
};

#endif // GHOSTEDARRAYWRAPPER_H
