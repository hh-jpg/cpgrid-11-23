/*
 * @Author: jltu && jltu.siat.ac.cn, 1620808124@qq.com
 * @Date: 2024-08-02 14:28:20
 * @LastEditTime: 2024-08-24 16:47:57
 * @FilePath: /cpgrid/src/vec_func.C
 * @Description:
 *
 */
#include "vec_func.h"
#include <iostream>
PetscErrorCode globalvec_to_local(Vec &global_vec, const std::vector<unsigned int> &ghosted_index, Vec *local_vec)
{
  PetscFunctionBegin;
 
  // 拷贝global_vec的数据到local_vec
  PetscCall(VecCopy(global_vec, *local_vec));
  PetscCall(VecGhostUpdateBegin(*local_vec, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGhostUpdateEnd(*local_vec, INSERT_VALUES, SCATTER_FORWARD));
  PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode localvec_to_global(Vec &local_vec, Vec *global_vec)
{
  PetscErrorCode ierr;
  MPI_Comm comm = PetscObjectComm((PetscObject)*global_vec);

  // 从local_vec拷贝数据到global_vec
  const PetscScalar *local_array;
  PetscScalar *global_array;

  // 获取local_vec的数据
  ierr = VecGetArrayRead(local_vec, &local_array);
  CHKERRQ(ierr);
  // 获取global_vec的数据
  ierr = VecGetArray(*global_vec, &global_array);
  CHKERRQ(ierr);

  // 复制数据到global_vec的本地部分
  PetscInt local_start, local_end;
  ierr = VecGetOwnershipRange(*global_vec, &local_start, &local_end);
  CHKERRQ(ierr);
  for (PetscInt i = local_start; i < local_end; ++i)
  {
    global_array[i - local_start] = local_array[i - local_start];
  }

  // 释放local_vec的数据
  ierr = VecRestoreArrayRead(local_vec, &local_array);
  CHKERRQ(ierr);

  // 释放global_vec的数据
  ierr = VecRestoreArray(*global_vec, &global_array);
  CHKERRQ(ierr);

  return ierr;
}

PetscErrorCode local_index_in_vec(Vec &ghosted_vec, PetscInt global_index, PetscInt *local_index)
{
  PetscErrorCode ierr;
  PetscInt local_start, local_end;

  // 获取本地拥有的全局索引范围
  ierr = VecGetOwnershipRange(ghosted_vec, &local_start, &local_end);
  CHKERRQ(ierr);

  // 检查global_index是否在本地拥有的范围内
  if (global_index >= local_start && global_index < local_end)
  {
    *local_index = global_index - local_start;
    return 0; // 成功找到
  }

  // 获取幽灵节点数量
  PetscInt local_size, total_size;
  ierr = VecGetLocalSize(ghosted_vec, &local_size);
  CHKERRQ(ierr);
  ierr = VecGetSize(ghosted_vec, &total_size);
  CHKERRQ(ierr);
  PetscInt nghost = local_size - (local_end - local_start);

  // 获取幽灵节点的全局索引
  const PetscInt *ghosted_indices;
  ierr = VecGetOwnershipRanges(ghosted_vec, &ghosted_indices);
  CHKERRQ(ierr);

  // 检查global_index是否在幽灵节点索引范围内
  for (PetscInt i = 0; i < nghost; ++i)
  {
    if (ghosted_indices[local_end - local_start + i] == global_index)
    {
      *local_index = local_end - local_start + i;
      return 0; // 成功找到
    }
  }

  // 如果没有找到，返回一个错误
  *local_index = -1;
  return PETSC_ERR_ARG_OUTOFRANGE; // 索引未找到
}

PetscErrorCode SaveVecToMatlab(Vec &vec, const char *filename, const char *vecname)
{
  PetscViewer viewer;
  PetscFunctionBegin;
  PetscCall(PetscObjectSetName((PetscObject)vec, vecname));
  PetscCall(PetscViewerASCIIOpen(PETSC_COMM_WORLD, filename, &viewer));
  PetscCall(PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB));
  PetscCall(VecView(vec, viewer));
  PetscCall(PetscViewerPopFormat(viewer));
  PetscCall(PetscViewerDestroy(&viewer));
  PetscFunctionReturn(PETSC_SUCCESS);
}
